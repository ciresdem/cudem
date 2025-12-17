#!/usr/bin/env python3
"""
Canonical metadata table from .gpkg files (recursive) with cleaning.

Columns (default, in order):
    Source, Title, Date, Data_Type, Resolution, Hdatum, Vdatum, URL
Optional:
    --include-gpkg  → prepend 'GPKG' (filename only) as first column

Defaults:
  • Unique rows (dedupe) and drop rows where all canonical cols are empty
  • Sort by Source (empty last), tiebreak by Title
  • Write CSV and XLSX (prompts to install pandas/openpyxl if missing)

Fix-ups (heuristics + strict):
  • Swap Date <-> Source when values appear reversed
  • Swap Data_Type <-> Resolution when values appear reversed
  • STRICT: Source must be text-only (no digits) & Date numeric-ish;
            if Source looks numeric-ish and Date looks text-only → swap

Flags:
  • --include-gpkg     include GPKG column (filename only) as first column
  • --keep_all         keep duplicates (turn off deduplication)
  • --keep_empty       keep rows even if all canonical columns are empty
  • --no-xlsx          skip writing Excel
  • --out, --xlsx-out  custom output paths
"""

import os
import sys
import csv
import re
import argparse
import subprocess
from osgeo import ogr

ogr.UseExceptions()

# ---------------- Excel dependency handling (interactive) ----------------

def ensure_excel_deps_interactive():
    def _try():
        try:
            import pandas  # noqa
            has_pd = True
        except Exception:
            has_pd = False
        try:
            import openpyxl  # noqa
            has_ox = True
        except Exception:
            has_ox = False
        return has_pd, has_ox

    has_pd, has_ox = _try()
    missing = []
    if not has_pd: missing.append("pandas")
    if not has_ox: missing.append("openpyxl")
    if not missing:
        return True, True

    try:
        resp = input(f"Excel output requires {', '.join(missing)}. Install now with pip? [Y/n]: ").strip().lower()
    except EOFError:
        resp = "n"

    if resp in ("", "y", "yes"):
        pkgs = []
        if not has_pd: pkgs.append("pandas")
        if not has_ox: pkgs.append("openpyxl")
        try:
            print(f"→ Installing: {' '.join(pkgs)}")
            subprocess.run([sys.executable, "-m", "pip", "install", *pkgs], check=True)
        except Exception as e:
            print(f"⚠️  Installation failed: {e}")
        has_pd, has_ox = _try()
        if not (has_pd and has_ox):
            print("ℹ️  Excel dependencies still missing. XLSX will be skipped.")
    else:
        print("ℹ️  Skipping Excel dependency install. XLSX will be skipped.")

    return has_pd, has_ox

# ---------------- Canonicalization config ----------------

CANON_COLS = ["Source", "Title", "Date", "Data_Type", "Resolution", "Hdatum", "Vdatum", "URL"]

ALIAS = {
    "Source": [
        "source","data_source","datasource","provider","agency","organization","organisation",
        "org","originator","credit","publisher","owner","maintainer","author","contact_org",
        "contact","responsibleparty","responsible_party","producer","supplied_by"
    ],
    "Title": [
        "title","name","dataset","dataset_name","layer","layer_name","mapname","product","theme"
    ],
    "Date": [
        "date","pub_date","publication_date","publish_date","issued","issue_date","updated",
        "update_date","last_update","last_updated","revised","revision_date","created",
        "creation_date","date_created","date_published","year","yyyy"
    ],
    "Data_Type": [
        "data_type","datatype","type","layer_type","category","format","geom_type","feature_type","raster_type"
    ],
    "Resolution": [
        "resolution","spatial_resolution","nominal_resolution","cellsize","cell_size",
        "pixel_size","grid_res","x_res","y_res","ground_sample_distance","gsd"
    ],
    "Hdatum": [
        "hdatum","horizontal_datum","h_datum","horiz_datum","datum_h","geodetic_datum",
        "reference_frame","crs_horizontal","hdatum_name"
    ],
    "Vdatum": [
        "vdatum","vertical_datum","v_datum","vert_datum","vertical_reference","tidal_datum","vdatum_name"
    ],
    "URL": [
        "url","link","weblink","website","homepage","source_url","data_url","download",
        "download_url","doi","portal","landing_page"
    ],
}

EXTRA_SUBSTR_TO_CANON = {
    "agency": "Source",
    "name":   "Title",
}

URL_REGEX = re.compile(r"https?://[^\s\)\]]+", re.IGNORECASE)

# ---------------- Heuristic classifiers ----------------

_org_tokens = re.compile(r"(noaa|usgs|usace|ncei|nerr|ocm|epa|usfs|blm|esri|navy|coast guard|univ|university|dept|department|state|county|city|survey|hydrographic|natural resources|dnr|dwr|boem|nps|fhwa|dot|transport|oceans|geological|academy|institute|center|centre|authority)", re.I)
_date_tokens = re.compile(r"(^\s*\d{4}(\-\d{1,2}(\-\d{1,2})?)?\s*$)|([A-Za-z]{3,9}\s+\d{1,2},\s*\d{4})|(Q[1-4]\s*\d{4})", re.I)
_res_tokens  = re.compile(r"(\b\d+(\.\d+)?\s*(m|meter|metre|ft|feet|km|cm)\b)|(\b\d+(\.\d+)?\s*(arc-?sec(ond)?|arcmin(ute)?|deg(ree)?)\b)|(\b1/\d+\s*arc-?sec\b)", re.I)
_dtype_tokens = re.compile(r"(lidar|topo|bathym|multibeam|singlebeam|sonar|chart|raster|vector|contour|photogrammetry|imagery|dem|dtm|dsm|elevation|soundings?)", re.I)

def looks_like_date(s):       return bool(s) and (_date_tokens.search(s) is not None)
def looks_like_source(s):     return bool(s) and (_org_tokens.search(s) is not None) and not looks_like_date(s)
def looks_like_resolution(s): return bool(s) and (_res_tokens.search(s) is not None)
def looks_like_dtype(s):      return bool(s) and (_dtype_tokens.search(s) is not None)

# -------- STRICT checks: Source text-only; Date numeric-ish (improved) --------

# Accept year ranges like "2018 to 2019", "2018–2019", "2018-2019", "2018/2019"
_YEAR = r"(?:19|20)\d{2}"
_NUMISH_PATTERNS = [
    re.compile(rf"^\s*{_YEAR}\s*(?:to|–|—|-|/)\s*{_YEAR}\s*$", re.I),
    re.compile(rf"^\s*{_YEAR}\s*$", re.I),
    re.compile(r"^\s*\d{1,2}[/-]\d{1,2}[/-]\d{2,4}\s*$"),
    re.compile(r"^\s*Q[1-4]\s*(?:'?\d{{2}}|\d{{4}})\s*$", re.I),
    re.compile(r"^\s*\d[\d\s./-]*\d?\s*$"),  # general numeric-ish fallback
]

_TXT_ONLY = re.compile(r"^[^\d]*[A-Za-z][^\d]*$")  # at least one letter, no digits

def is_text_only_no_digits(s: str) -> bool:
    s = (s or "").strip()
    return bool(_TXT_ONLY.match(s))

def is_numericish_date(s: str) -> bool:
    s = (s or "").strip()
    if not s:
        return False
    # allow the word 'to' (date range) and separators
    s_norm = re.sub(r"\bto\b", "", s, flags=re.I)
    # quick win: contains a four-digit year?
    has_year = re.search(_YEAR, s_norm) is not None
    # no other letters allowed after removing 'to'
    has_letters = re.search(r"[A-Za-z]", s_norm) is not None
    # match common numeric-ish forms
    matches_numeric = any(p.match(s_norm) for p in _NUMISH_PATTERNS)
    return (has_year and not has_letters) or matches_numeric

# ---------------- field mapping helpers ----------------

def normalize_fieldname(s: str) -> str:
    s = s or ""
    s = s.strip().lower()
    return re.sub(r"[^a-z0-9_]+", "", s)

def classify_field_to_canonical(fieldname: str):
    nf = normalize_fieldname(fieldname)
    for sub, canon in EXTRA_SUBSTR_TO_CANON.items():
        if sub in nf:
            return canon
    for canon, patterns in ALIAS.items():
        if nf in patterns:
            return canon
    if nf.endswith("url") or "http" in nf or "link" in nf:
        return "URL"
    if "datum" in nf and ("v" in nf or "vertical" in nf or "tidal" in nf):
        return "Vdatum"
    if "datum" in nf and ("h" in nf or "horizontal" in nf or "geodetic" in nf):
        return "Hdatum"
    if "res" in nf or "pixel" in nf or "cell" in nf:
        return "Resolution"
    if nf.endswith("type") or "type" in nf:
        return "Data_Type"
    return None

# ---------------- I/O helpers ----------------

def iter_gpkg_paths(rootdir):
    for dirpath, _, filenames in os.walk(rootdir):
        for fn in filenames:
            if fn.lower().endswith(".gpkg") and not fn.endswith(("-wal", "-journal")):
                yield os.path.join(dirpath, fn)

def get_layer_field_names(layer):
    defn = layer.GetLayerDefn()
    return [defn.GetFieldDefn(i).GetName() for i in range(defn.GetFieldCount())]

def as_text(v):
    if v is None:
        return ""
    if isinstance(v, (list, tuple)):
        return " ".join(map(str, v))
    return str(v)

def de_dupe_preserve_order(values):
    seen = set()
    out = []
    for v in values:
        vv = v.strip()
        if not vv:
            continue
        key = vv.lower()
        if key not in seen:
            seen.add(key)
            out.append(vv)
    return out

# ---------------- Core logic ----------------

def collect_rows(rootdir):
    rows = []
    for gpkg_path in iter_gpkg_paths(rootdir):
        gpkg_name = os.path.basename(gpkg_path)
        try:
            ds = ogr.Open(gpkg_path, 0)
        except Exception:
            continue
        if ds is None:
            continue

        for i in range(ds.GetLayerCount()):
            layer = ds.GetLayerByIndex(i)
            if layer is None:
                continue

            fnames = get_layer_field_names(layer)
            mapping = {fn: classify_field_to_canonical(fn) for fn in fnames}
            layer.ResetReading()

            for feat in layer:
                buckets = {c: [] for c in CANON_COLS}

                # 1) Direct matches via mapping
                for fn in fnames:
                    canon = mapping.get(fn)
                    if canon:
                        try:
                            val = feat.GetField(fn)
                        except KeyError:
                            val = None
                        if val not in (None, ""):
                            buckets[canon].append(as_text(val))

                # 2) URL regex fallback across all attributes
                for fn in fnames:
                    try:
                        s = as_text(feat.GetField(fn))
                    except KeyError:
                        continue
                    for m in URL_REGEX.findall(s):
                        buckets["URL"].append(m)

                # 3) Consolidate into a single row dict
                row = {}
                for c in CANON_COLS:
                    row[c] = " | ".join(de_dupe_preserve_order(buckets[c]))

                # 4) Heuristic & strict fix-ups
                row = clean_swaps(row)

                # Attach provenance (filename) for optional output
                row["_GPKG"] = gpkg_name

                rows.append(row)
    return rows

def clean_swaps(row):
    """Heuristic swaps first, then strict enforcement for Source/Date."""
    src = row.get("Source", "").strip()
    dat = row.get("Date", "").strip()
    dty = row.get("Data_Type", "").strip()
    res = row.get("Resolution", "").strip()

    # --- heuristic Source/Date swaps ---
    if looks_like_date(src) and (looks_like_source(dat) or not looks_like_date(dat)):
        row["Source"], row["Date"] = dat, src
        src, dat = row["Source"], row["Date"]
    elif looks_like_source(dat) and looks_like_date(src):
        row["Source"], row["Date"] = dat, src
        src, dat = row["Source"], row["Date"]

    # --- heuristic Data_Type/Resolution swaps ---
    if (not looks_like_resolution(res)) and looks_like_resolution(dty):
        row["Data_Type"], row["Resolution"] = res, dty
        dty, res = row["Data_Type"], row["Resolution"]
    elif (not looks_like_dtype(dty)) and looks_like_resolution(res):
        pass
    elif (not looks_like_resolution(res)) and (not looks_like_resolution(dty)):
        if any(t in dty.lower() for t in [" m", "meter", "metre", "ft", "arc", "pixel", "cell", "gsd", "1/"]):
            row["Data_Type"], row["Resolution"] = res, dty
            dty, res = row["Data_Type"], row["Resolution"]

    # --- STRICT enforcement: Source must be text-only; Date numeric-ish ---
    src = row.get("Source", "").strip()
    dat = row.get("Date", "").strip()
    if is_numericish_date(src) and is_text_only_no_digits(dat):
        row["Source"], row["Date"] = dat, src

    return row

def unique_rows(rows, include_gpkg):
    key_cols = (["_GPKG"] if include_gpkg else []) + CANON_COLS
    seen = set()
    out = []
    for r in rows:
        key = tuple((r.get(c if c != "_GPKG" else "_GPKG", "").strip().lower()) for c in key_cols)
        if key not in seen:
            seen.add(key)
            out.append(r)
    return out

def drop_empty_rows(rows):
    out = []
    for r in rows:
        keep = any((r.get(c, "").strip() != "") for c in CANON_COLS)
        if keep:
            out.append(r)
    return out

def sort_rows(rows):
    def sort_key(r):
        src = r.get("Source", "")
        ttl = r.get("Title", "")
        empty = (src.strip() == "")
        return (empty, src.casefold(), ttl.casefold())
    return sorted(rows, key=sort_key)

def write_csv(rows, out_csv, include_gpkg):
    cols = (["GPKG"] if include_gpkg else []) + CANON_COLS
    os.makedirs(os.path.dirname(os.path.abspath(out_csv)) or ".", exist_ok=True)
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in rows:
            rr = {c: (r.get("_GPKG", "") if c == "GPKG" else r.get(c, "")) for c in cols}
            w.writerow(rr)

def write_xlsx(rows, out_xlsx, include_gpkg, force=False):
    has_pd, has_ox = ensure_excel_deps_interactive()
    if not (has_pd and has_ox):
        if force:
            print("⚠️  Could not satisfy Excel deps; XLSX not written.")
        return
    import pandas as pd
    cols = (["GPKG"] if include_gpkg else []) + CANON_COLS
    data = []
    for r in rows:
        data.append({c: (r.get("_GPKG", "") if c == "GPKG" else r.get(c, "")) for c in cols})
    df = pd.DataFrame(data, columns=cols)
    os.makedirs(os.path.dirname(os.path.abspath(out_xlsx)) or ".", exist_ok=True)
    with pd.ExcelWriter(out_xlsx, engine="openpyxl") as writer:
        df.to_excel(writer, index=False, sheet_name="canonical")

# ---------------- CLI ----------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("root", nargs="?", default=os.getcwd(),
                    help="Root directory to scan (default: current directory).")
    ap.add_argument("--out", default="gpkg_canonical_table.csv",
                    help="CSV output path (default: gpkg_canonical_table.csv).")
    ap.add_argument("--xlsx-out", default="gpkg_canonical_table.xlsx",
                    help="XLSX output path (default: gpkg_canonical_table.xlsx).")
    ap.add_argument("--include-gpkg", action="store_true",
                    help="Include 'GPKG' column (filename only) as first column.")
    ap.add_argument("--keep_all", action="store_true",
                    help="Keep duplicates (default removes duplicates).")
    ap.add_argument("--keep_empty", action="store_true",
                    help="Keep rows where all canonical columns are empty (default drops).")
    ap.add_argument("--no-xlsx", action="store_true",
                    help="Do not write the Excel workbook.")
    args = ap.parse_args()

    print(f"🔍 scanning: {args.root}")
    rows = collect_rows(args.root)
    print(f"📦 collected rows: {len(rows)}")

    if not args.keep_all:
        before = len(rows)
        rows = unique_rows(rows, include_gpkg=args.include_gpkg)
        after = len(rows)
        print(f"🧹 unique default: {before} → {after} rows")

    if not args.keep_empty:
        before = len(rows)
        rows = drop_empty_rows(rows)
        after = len(rows)
        print(f"🚮 drop_empty default: {before} → {after} rows")

    rows = sort_rows(rows)

    write_csv(rows, args.out, include_gpkg=args.include_gpkg)
    print(f"✅ wrote CSV: {args.out} (rows={len(rows)}, cols={(8 + (1 if args.include_gpkg else 0))})")

    if not args.no_xlsx:
        write_xlsx(rows, args.xlsx_out, include_gpkg=args.include_gpkg, force=True)

if __name__ == "__main__":
    main()
