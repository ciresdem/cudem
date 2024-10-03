import argparse
import ast
import dateparser
import geopandas
import numpy
import os
import pandas
import pyproj
import rasterio, rasterio.warp, rasterio.crs
import shapely
import sys
import typing

from cudem import dlim, regions
import utils.pickle_blosc


def get_photon_dataframe(polygon_bbox_or_dem_fname: typing.Union[shapely.geometry, list, tuple, str],
                         dem_horz_reference_frame: typing.Union[str, None, pyproj.CRS] = None,
                         dem_vert_reference_frame: typing.Union[str, None, pyproj.CRS] = None,
                         start_date: str = "a year ago midnight",
                         end_date: str = "midnight today",
                         other_columns: typing.Union[dict, None] = None,
                         classify_bathymetry: bool = True,
                         classify_buildings: bool = True,
                         classifications_to_keep: typing.Union[list, tuple, str, set, numpy.ndarray] = (1, 2, 3, 4),
                         conf_levels_to_keep: typing.Union[list, tuple, str, set, numpy.ndarray] = "4",
                         ) -> typing.Union[pandas.DataFrame, geopandas.GeoDataFrame, None]:
    """Return a dataframe of all classified photons within the polygon/bounding-box and timeframe, from NSIDC data.

    Photons will be returned in the vertical and horizontal reference frame specified by dem_horz_reference_frame and
    dem_vert_reference_frame, respectively.

    The original X,Y,Z coordinates (in WGS84 (ITRF 2014) lat/lon and EGM2008 height) will be in columns

    Parameters:
        polygon_bbox_or_dem_fname: A shapely.geometry.Polygon object, or bounding-box is a 4-value list/tuple
            of (xmin,xmax,ymin,ymax). Coordinates should be in the reference frame outlined in dem_horz_reference_frame.
            Alternately, can also be the name of a DEM file from which to derive the bounding box.
        dem_horz_reference_frame: The horizontal reference frame of the DEM. Default None.
            If None, attempt to derive the horizontal reference frame from the DEM file if a filename is provided in
            the parameter polygon_bbox_or_dem_fname. This is None and no filename is provided, throws a ValueError.
        dem_vert_reference_frame: The vertical reference frame of the DEM. Can be any value suppored by cudems.vdatum.
            Run "vdatums --list-epsg" for a full list of options. Default None.
            If None, attempt to derive the vertical reference frame from the DEM file if a filename is provided in
            the parameter polygon_bbox_or_dem_fname. If unsuccessful, it will attempt to get it from the horizontal
            reference frame, in case it is a compound reference frame. If unsuccessful, it will throw a ValueError.
        start_date: The start date of the time period to query. Can be any string parsed by dateparser.parse().
            Default 'a year ago midnight'.
        end_date: The end date of the time period to query. Can be any string parsed by dateparser.parse().
            Default 'midnight today'. If end_date is before start_date, will throw a ValueError.
        other_columns: A dictionary of other columns to include in the dataframe. It should include the ATL03 variable
            path followed by the column name you want to assign it. For example {'/heights/h_ph':'h_ph'} Default None.
        classify_bathymetry: Whether to classify the bathymetry. Default True, which will run CShelph.
        classify_buildings: Whether to classify the buildings. Default True, which will use Bing building footprints.
        classifications_to_keep: A list of the photon classifications to keep. Default (1, 2, 3, 4, 5, 7).
            See "dlim --modules icesat2" for a list of classification numbers.
        conf_levels_to_keep: A list of the photon confidence levels to keep.
            Default "4" (only keep highest-confidence photons). Can also be a set, tuple, or numpy.ndarray. If a string,
            values should be integers 0-4 seperated by forward slashes (no spaces).

    Raises:
        ValueError: if dem_horz_reference_frame is None, and no filename is provided in the
            parameter polygon_bbox_or_dem_fname. Or if the start_date is on or after the end_date.

    Returns:
        A geopandas.GeoDataFrame of all classified photons within the polygon/bounding-box and timeframe.
        None if no ICESat-2 data is returned. The x,y,z points will be in the horizontal and vertical reference frames
        specified in the parameters. If merge_dataframes is False, return a list of the resultant dataframes.
    """

    # Extract the dem_horz and def_vert reference frames from either the file or the user inputs.
    horz_proj, vert_proj = _get_dem_horz_and_vert_reference_frame(polygon_bbox_or_dem_fname,
                                                                  dem_horz_reference_frame,
                                                                  dem_vert_reference_frame)
    assert horz_proj is not None and vert_proj is not None

    # Make sure the date strings are valid and aren't on the same day or in reverse order.
    start_datestr = _process_input_date_str(start_date)
    end_datestr = _process_input_date_str(end_date)
    if start_datestr >= end_datestr:
        raise ValueError("The start date must occur one or more days before the end date.")

    # Get the list of classifications to keep as a "0/1/2/3" type of string of integers separated by forward-slashes.
    classes_str = _get_classifications_str(classifications_to_keep)
    conf_str = _get_classifications_str(conf_levels_to_keep)

    dem_wgs84_bbox = get_wgs84_bounding_box(polygon_bbox_or_dem_fname, horz_proj)
    region = regions.Region().from_list(list(dem_wgs84_bbox))

    fetches_module = f'icesat2:subset=True:time_start={start_datestr}:time_end={end_datestr}'

    # Get a short string for the srs of the DEM, from the two projections above.
    dem_srs_string = get_dem_srs_string(horz_proj, vert_proj)
    # ICESat-2, by default using dlim, comes in WGS84 lat-lon coordinates and EGM2008 vertical reference frame.
    is2_srs_string = "EPSG:4326+3855"

    # print(f"fn={fetches_module}")
    # print(f"src_region=regions.Region().from_list({dem_wgs84_bbox})")
    # print(f"classes={classes_str}")
    # print(f"confidence_levels={conf_str}")
    # print(f'columns={ {} if other_columns is None else other_columns}')
    # print(f"classify_bathymetry={classify_bathymetry}")
    # print(f"cshelph={classify_bathymetry}")
    # print(f"classify_buildings={classify_buildings}")

    # import sys
    # sys.exit(0)

    # Initialize the icesat2 query.
    ds = dlim.IceSat2Fetcher(fn=fetches_module,
                             src_region=region,
                             src_srs=is2_srs_string,
                             dst_srs=dem_srs_string,
                             classes=classes_str,
                             confidence_levels=conf_str,
                             columns=other_columns if other_columns else {},
                             classify_bathymetry=classify_bathymetry,
                             cshelph=classify_bathymetry,
                             classify_buildings=classify_buildings).initialize()

    # Run the icesat2 query over all files and lasers in that box.
    list_of_dfs = []

    for i, df in enumerate(ds.yield_points()):
        # Drop any rows that have NaN values.
        df.dropna(inplace=True)
        # Drop columns we don't need.
        df.drop(columns=["ref_elevation", "ref_azimuth", "ref_sat_alt", "ph_h_watermask"], inplace=True)
        # Rename the column "ph_h_classed" to "class_code"
        df.rename(columns={"ph_h_classed": "class_code"}, inplace=True)
        # Add a column that is the unique ID for each individual laser track in the dataframe.
        # TODO: Later put the actual granule name in here.
        df["unique_laser_id"] = i

        # If there are any points remaining after this, append the dataframe to the list.
        if len(df) > 0:
            list_of_dfs.append(df)

    if len(list_of_dfs) == 0:
        return None

    # Merge the dataframes together.
    df_out = pandas.concat(list_of_dfs, ignore_index=True)

    # Now filter the points again, this time using the dem's native bbox or polygon provided by the user.
    gdf = _filter_points_against_dem(df_out, polygon_bbox_or_dem_fname, horz_proj)

    if len(gdf) == 0:
        return None
    else:
        return gdf


def get_dem_srs_string(horz_reference: pyproj.CRS,
                       vert_reference: pyproj.CRS) -> str:
    """Get the short string for the projection from the projections.

    This assumes a valid SRS projection has been found for both the horizontal and vertical components. If it is a 3D
    reference frame they should be the same (as returned by _get_dem_horz_and_vert_reference_frame), and that the
    authority (EPSG, etc) of both reference frames should be the same.

    Raises:
        ValueError if the two datums are based on different authorities.

    Returns:
        string in the format ("AUTH:HORZ+VERT")
    """
    horz_auth = horz_reference.list_authority()[0].auth_name
    vert_auth = vert_reference.list_authority()[0].auth_name

    if horz_auth != vert_auth:
        raise ValueError("Reference authorites for the horizontal datum and vertical datum must match.")

    if horz_reference.equals(vert_reference):
        return horz_reference.srs
    else:
        return f"{horz_auth}:{horz_reference.list_authority()[0].code}+{vert_reference.list_authority()[0].code}"


def _get_dem_horz_and_vert_reference_frame(polygon_bbox_or_dem_fname: str,
                                           dem_horz_reference_frame: typing.Union[str, None],
                                           dem_vert_reference_frame: typing.Union[str, None]) -> typing.Tuple:
    """From the user inputs, get the horizontal and vertical reference frame of the DEM.

    Use the user-inputs if they're defined, otherwise attempt to get them from the DEM file (if provided).

    Parameters:
        polygon_bbox_or_dem_fname: A shapely.geometry.Polygon object, or bounding-box is a 4-value list/tuple
            of (xmin,xmax,ymin,ymax). Coordinates should be in the reference frame outlined in dem_horz_reference_frame.
            Alternately, can also be the name of a DEM file from which to derive the bounding box.
        dem_horz_reference_frame: The horizontal reference frame of the DEM.
        dem_vert_reference_frame: The vertical reference frame of the DEM. Can be any value suppored by cudems.vdatum.
            Run "vdatums --list-epsg" for a full list of options.

    Raises:
        ValueError: if dem_horz_reference_frame is None, and no filename is provided in the
            parameter polygon_bbox_or_dem_fname.

    Returns:
        A tuple of (dem_horz_reference_frame, dem_vert_reference_frame) as pyproj.CRS objects.
    """

    # First, get a lat-lon bounding box from the "polygon_bbox_or_dem_fname" parameter.
    if ((dem_horz_reference_frame is None or dem_vert_reference_frame is None)
            and not os.path.exists(polygon_bbox_or_dem_fname)):
        raise ValueError("If no filename is provided in the parameter polygon_bbox_or_dem_fname, "
                         "dem_horz_reference_frame and dem_vert_reference_frame must be specified.")

    # Get the DEM reference frames, both horizontal and vertical
    if (dem_horz_reference_frame is None) and (dem_vert_reference_frame is None):
        dem_horz_reference_frame, dem_vert_reference_frame = get_dem_reference_frame_from_file(
            polygon_bbox_or_dem_fname,
            vert_horz_or_both="both")

    # If a value is None, then attempt to get it from the file.
    # If a value is a caller-created string, then derive the reference frame from that.
    elif dem_horz_reference_frame is None:
        # If only the DEM's vertical reference frame is given here, then get the horizontal one from the file.
        dem_horz_reference_frame = get_dem_reference_frame_from_file(polygon_bbox_or_dem_fname,
                                                                     vert_horz_or_both="horz")
        # Get the vertical reference frame from the user input
        dem_vert_reference_frame = get_dem_reference_frame_from_user_input(dem_vert_reference_frame,
                                                                           vert_horz_or_both="vert")

    elif dem_vert_reference_frame is None:
        # If only the DEM's horizontal reference frame is given here, then get the vertical one from the file.
        dem_vert_reference_frame = get_dem_reference_frame_from_file(polygon_bbox_or_dem_fname,
                                                                     vert_horz_or_both="vert")
        # Get the horizontal reference frame from the user input
        dem_horz_reference_frame = get_dem_reference_frame_from_user_input(dem_horz_reference_frame,
                                                                           vert_horz_or_both="horz")
    else:
        # Get the horizontal and vertical reference frames from the user input
        dem_horz_reference_frame = get_dem_reference_frame_from_user_input(dem_horz_reference_frame,
                                                                           vert_horz_or_both="horz")
        dem_vert_reference_frame = get_dem_reference_frame_from_user_input(dem_vert_reference_frame,
                                                                           vert_horz_or_both="vert")

    # At this point, we should have a valid vertical and horizontal reference frame to be working from.
    if dem_horz_reference_frame is None or dem_vert_reference_frame is None:
        prefix_text = ""
        if dem_horz_reference_frame is None:
            prefix_text = prefix_text + "dem_horz_reference_frame is not specified or could not be derived.\n"
        if dem_vert_reference_frame is None:
            prefix_text = prefix_text + "dem_vert_reference_frame is not specified or coult not be derived.\n"
        raise ValueError(
            prefix_text + "dem_horz_reference_frame and dem_vert_reference_frame must both be specified or "
                          "able to be derived from a DEM file.")

    return dem_horz_reference_frame, dem_vert_reference_frame


def get_polygon(polygon_bbox_or_dem_fname: typing.Union[shapely.geometry, list, tuple, str]) -> shapely.geometry.Polygon:
    """From a filename, a shapely geometry, or a bounding box, return a shapely polygon object."""
    if type(polygon_bbox_or_dem_fname) in (shapely.geometry.Polygon, shapely.geometry.MultiPolygon):
        return polygon_bbox_or_dem_fname
    elif type(polygon_bbox_or_dem_fname) is list or type(polygon_bbox_or_dem_fname) is tuple:
        if len(polygon_bbox_or_dem_fname) == 4:
            bbox = polygon_bbox_or_dem_fname
            return shapely.geometry.box(bbox[0], bbox[2], bbox[1], bbox[3])
        elif len(polygon_bbox_or_dem_fname) > 4 and len(polygon_bbox_or_dem_fname) % 2 == 0:
            return shapely.geometry.Polygon(polygon_bbox_or_dem_fname)
        else:
            raise TypeError("polygon_bbox_or_dem_fname must be a 4-value list of (xmin, xmax, ymin, ymax) or a greater "
                            "length list/tuple defining a points of a polygon outline.")
    elif type(polygon_bbox_or_dem_fname) is str:
        if os.path.exists(polygon_bbox_or_dem_fname):
            bbox = rasterio.open(polygon_bbox_or_dem_fname).bounds
            return shapely.geometry.box(*bbox)
        else:
            raise FileNotFoundError(f"If polygon_bbox_or_fname is a string, it must point to a filename that exists. "
                                    f"Did not find {polygon_bbox_or_dem_fname}.")
    else:
        raise TypeError("polygon_bbox_or_dem_fname must be either a 4-item bounding box, a list of coordinates "
                        "defining a polygon, a shapely geometry Polygon, or a filename.")


def _filter_points_against_dem(df: pandas.DataFrame,
                               polygon_bbox_or_dem_fname: typing.Union[shapely.geometry, list, tuple, str],
                               dem_crs: pyproj.CRS,
                               xy_cols: typing.Union[list, tuple] = ("x", "y")) -> pandas.DataFrame:
    """Filter the dataframe of points against the DEM outline or user-defined polygon in its native coordinates."""
    print("Converting to GeoDataFrame and subsetting photons...", end=" ", flush=True)
    gdf = geopandas.GeoDataFrame(df, geometry=geopandas.points_from_xy(df[xy_cols[0]], df[xy_cols[1]]), crs=dem_crs)

    gdf_subset = gdf[gdf.intersects(get_polygon(polygon_bbox_or_dem_fname))]
    print("Done.", flush=True)

    return gdf_subset


def get_wgs84_bounding_box(polygon_bbox_or_dem_fname: typing.Union[shapely.geometry, list, tuple, str],
                           dem_horz_reference_frame: typing.Union[str, pyproj.CRS, None] = None) -> list:
    """From a filename, a shapely geometry, or a bounding box, return a 4-value list of (xmin,xmax,ymin,ymax) in
    WGS84 (ITRF 2014) coordinates appropriate for querying NSIDC.

    Parameters:
        polygon_bbox_or_dem_fname: A shapely.geometry.Polygon object, or bounding-box is a 4-value list/tuple
            of (xmin,xmax,ymin,ymax). Coordinates should be in the reference frame outlined in dem_horz_reference_frame.
            Alternately, can also be the name of a DEM file from which to derive the bounding box.
        dem_horz_reference_frame: The horizontal reference frame of the DEM. Default None.
            If None, attempt to derive the horizontal reference frame from the DEM file if a filename is provided in
            the parameter polygon_bbox_or_dem_fname. This is None and no filename is provided, throws a ValueError.

    Raises:
        ValueError: if dem_horz_reference_frame is None, and no filename is provided in the
            parameter polygon_bbox_or_dem_fname.
        TypeError: If either parameter is of an unhandled type.

    Returns:
        A 4-value list of (xmin,xmax,ymin,ymax) in WGS84 (EPSG: 4326) coordinates appropriate for querying NSIDC.
    """
    polygon = None
    if type(polygon_bbox_or_dem_fname) is shapely.geometry.Polygon:
        polygon = shapely.Polygon(polygon_bbox_or_dem_fname.exterior.coords[:])

        dem_horz_reference_frame = get_dem_reference_frame_from_user_input(dem_horz_reference_frame, "horz")

    elif type(polygon_bbox_or_dem_fname) is list or type(polygon_bbox_or_dem_fname) is tuple:
        if len(polygon_bbox_or_dem_fname) == 4:
            bbox = polygon_bbox_or_dem_fname
            # Convert the bounds from (xmin, xmax, ymin, ymax) to (xmin, ymin, xmax, ymax) work with shapely.
            bbox = [bbox[0], bbox[2], bbox[1], bbox[3]]
            polygon = shapely.geometry.box(*bbox)
        elif len(polygon_bbox_or_dem_fname) > 4 and len(polygon_bbox_or_dem_fname) % 2 == 0:
            polygon = shapely.geometry.Polygon(polygon_bbox_or_dem_fname)
        else:
            raise TypeError("polygon_bbox_or_dem_fname must be a 4-value list of (xmin, xmax, ymin, ymax) or a greater "
                            "length list/tuple defining a points of a polygon outline.")

        dem_horz_reference_frame = get_dem_reference_frame_from_user_input(dem_horz_reference_frame, "horz")

    elif type(polygon_bbox_or_dem_fname) is str:
        if os.path.exists(polygon_bbox_or_dem_fname):
            if dem_horz_reference_frame is None:
                dem_horz_reference_frame = get_dem_reference_frame_from_file(polygon_bbox_or_dem_fname, "horz")
            else:
                dem_horz_reference_frame = get_dem_reference_frame_from_user_input(dem_horz_reference_frame, "horz")

            bbox = rasterio.open(polygon_bbox_or_dem_fname).bounds
            polygon = shapely.geometry.box(*bbox)

        else:
            raise FileNotFoundError(f"If polygon_bbox_or_fname is a string, it must point to a filename that exists. "
                                    f"Did not find {polygon_bbox_or_dem_fname}.")

    if dem_horz_reference_frame is None:
        raise ValueError("dem_horz_reference_frame not defined.")

    # These should both the true if we got here.
    assert type(polygon) is shapely.geometry.Polygon
    assert type(dem_horz_reference_frame) is pyproj.CRS
    assert not dem_horz_reference_frame.is_compound

    wgs84_crs = pyproj.CRS.from_user_input("EPSG:4326")

    # This DEM is already in WGS84, just return the bounding box of the polygon in (xmin,xmax,ymin,ymax) format.
    if dem_horz_reference_frame.equals(wgs84_crs):
        return polygon.bounds[0], polygon.bounds[2], polygon.bounds[1], polygon.bounds[3]

    # Otherwise, translate the polygon to WGS84 and return the bounding box of the polygon in (xmin,xmax,ymin,ymax) format.
    transformer = pyproj.Transformer.from_crs(dem_horz_reference_frame, wgs84_crs, always_xy=True)
    polygon_wgs84 = shapely.geometry.Polygon(shell=transformer.itransform(polygon.exterior.coords[:]))

    return polygon_wgs84.bounds[0], polygon_wgs84.bounds[2], polygon_wgs84.bounds[1], polygon_wgs84.bounds[3]


def get_dem_reference_frame_from_file(dem_fname: str, vert_horz_or_both: str = "both") -> \
        typing.Union[pyproj.CRS, tuple, None]:
    """Get the reference frame of the DEM.

    Parameters:
        dem_fname: The name of the DEM file.
        vert_horz_or_both: Whether to return the vertical or horizontal reference frame. Will accept any value that
            starts in "v", "h", or "b" to mean "vertical", "horizontal", or "both", respectively. Default "both".

    Raises:
        ValueError: if vert_horz_or_both does not begin with one of "v", "h", or "b".
        FileNotFoundError: if the DEM file does not exist or has no crs attached.

    Returns:
        The vertical or horizontal reference frame of the DEM as a pyproj.CRS, or None if no reference frame could be found.
        If "vert_horz_or_both" is "both", 2-value a tuple of the vertical and horizontal reference frames is returned.
    """
    if not os.path.exists(dem_fname):
        raise FileNotFoundError(f"DEM file {dem_fname} does not exist.")
    dem_ds = rasterio.open(dem_fname)

    if dem_ds.crs is None:
        dem_ds_str = ""
    else:
        dem_ds_str = dem_ds.crs

    return get_dem_reference_frame_from_user_input(dem_ds_str, vert_horz_or_both)


def get_dem_reference_frame_from_user_input(crs: typing.Union[pyproj.CRS, rasterio.crs.CRS, str, None],
                                            vert_horz_or_both: str = "both") -> typing.Union[pyproj.CRS, tuple, None]:
    """Get the horizontal reference frame from an input string.

    Parameters:
        crs: The CRS string or object. If None, this will return None.
        vert_horz_or_both: Whether to return the vertical or horizontal reference frame. Will accept any value that
            starts in "v", "h", or "b" to mean "vertical", "horizontal", or "both", respectively. Default "both".

    Raises:
        ValueError: if vert_horz_or_both doesn't begin with one of "v", "h", or "b".

    Returns:
        The vertical or horizontal reference frame of the DEM as a pyproj.CRS, or None if no reference frame could be found.
        If "vert_horz_or_both" is "both", 2-value a tuple of the vertical and horizontal reference frames is returned.
    """
    if crs is None or crs == "":
        crs_obj = None
    elif isinstance(crs, rasterio.crs.CRS) or isinstance(crs, pyproj.CRS):
        crs_obj = pyproj.CRS(crs)
    else:
        crs_obj = pyproj.CRS.from_user_input(crs)

    if crs_obj is None:
        horz, vert = None, None
    elif crs_obj.is_compound:
        horz, vert = crs_obj.sub_crs_list
    elif len(crs_obj.axis_info) == 3:
        # If it's a 3D axis it's not "compound" but has 2 horz axes and a vert axis (3 total) which counts here.
        horz, vert = crs_obj, crs_obj
    elif crs_obj.is_vertical:
        horz, vert = None, crs_obj
    else:
        horz, vert = crs_obj, None

    # Get the lower-case first non-space letter from the vert_horz_or_both input string.
    choice_letter = vert_horz_or_both.strip().lower()[0]
    if choice_letter == "b":
        retval = horz, vert
    elif choice_letter == "h":
        retval = horz
    elif choice_letter == "v":
        retval = vert
    else:
        raise ValueError(
            f"Uknown choice '{vert_horz_or_both}' for 'vert_horz_or_both'. Should begin with 'h', 'v', or 'b'")

    return retval


def _process_input_date_str(date_str: str) -> str:
    """Process a date string and return a 'YYYY-MM-DD' string.

    Input date string can be anything that can be parsed by dateparser.parse().
    """
    return dateparser.parse(date_str).strftime("%Y-%m-%d")


def _get_classifications_str(classifications: typing.Union[str, list, tuple, numpy.ndarray]) -> str:
    """Given a list of photon classifications codes, return a string separated by forward-slashes, compatible with dlim.

    Parameters:
        classifications: A string or list of strings of photon classifications codes.

    Returns:
        A string of classifications codes separated by forward-slashes.
    """
    if type(classifications) is str:
        return classifications.replace(" ", "").replace(",", "/")

    else:
        return "/".join([str(int(c)) for c in classifications])


def export_as_vector(gdf: geopandas.GeoDataFrame,
                     outfile: str):
    """Export a GeoDataFrame as a shapefile or geopackage."""
    if os.path.splitext(outfile)[-1] in (".gpkg", ".shp"):
        gdf.to_file(outfile)
    elif os.path.splitext(outfile)[-1] in (".blosc", ".blosc2"):
        utils.pickle_blosc.write(gdf, outfile)

    print(outfile, f"written with {len(gdf)}",
          f"{'photons' if isinstance(gdf.geometry.iloc[0], shapely.geometry.Point) else 'lines'}.")


def export_lines(points_gdf_or_fname: typing.Union[geopandas.GeoDataFrame, str],
                 outfile: str,
                 tolerance=0.00001):
    """Take a points file produced from the ICESat-2 photon data and export a simplified set of lines that follows
    the path of each laser.

    Export as a vector file."""
    if isinstance(points_gdf_or_fname, str):
        print("Reading", os.path.basename(points_gdf_or_fname), end="...", flush=True)
        if os.path.splitext(points_gdf_or_fname)[-1].lower() in (".shp", ".gpkg"):
            points_gdf = geopandas.read_file(points_gdf_or_fname)
        elif os.path.splitext(points_gdf_or_fname)[-1].lower() in (".blosc", ".blosc2"):
            points_gdf = utils.pickle_blosc.read(points_gdf_or_fname)
        else:
            raise ValueError(f"Unrecognized file extension in {points_gdf_or_fname}.")

        print(" Done.", flush=True)
    else:
        assert isinstance(points_gdf_or_fname, geopandas.GeoDataFrame)
        points_gdf = points_gdf_or_fname

    # Group by the points by each unique laser ID, sort them by y-coordinate, and export as LineStrings into a new GDF.
    lines_gdf = (points_gdf.groupby("unique_laser_id")["geometry"]
                 .apply(lambda x: shapely.geometry.LineString(sorted(x.tolist(),
                                                                     key=lambda p: p.coords[0][1])).simplify(tolerance=tolerance)))

    export_as_vector(lines_gdf, outfile)


def define_and_parse_args():
    """Define and parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Query ICESat-2 data.")
    parser.add_argument("outfile", help="The name of the output file.")
    parser.add_argument("-reg", "--region", target="region", type=str, default=None,
                        help="The region to query (in output coordinates). A comma-separated string (no spaces) "
                             "of xmin,xmax,ymin,ymax, or a longer string of x1,y1,x2,y2,x3,y3,etc coordinates "
                             "defining a polygon.")
    parser.add_argument("-rast", "--raster", dest="raster", default=None,
                        help="The name of a raster file to use to define the polygon extent. "
                             "One of 'region' or 'input_raster' must be set.")
    parser.add_argument("-crs", "--crs", dest="crs", default=None,
                        help="The coordinate reference system of the polygon or raster. Default: if 'raster' is "
                             "set, get the horizontal reference frame from the raster. If used, will override any crs "
                             "defined in the raster."
                             "At least one of --raster or --crs must be set.")
    parser.add_argument("-vd", "--vdatum", dest="vdatum", default=None,
                        help="The vertical datum of the polygon or raster. Default: if --raster is set, attempt"
                             "to get the horizontal datum from the input raster. Else, if --crs is set and is a "
                             "compound or 3D coordinate system, get the vertical datum from that. "
                             "If neither of these conditions is true, --vdatum must be set.")
    parser.add_argument("--xyz_cols", dest="xyz_cols", default=None,
                        help="The names of the columns containing the transformed x, y, and z coordinates, in a comma-separated string."
                             "Default: 'x_dem,y_dem,z_dem'. Setting to x,y,z will overwrite the original coordinate columns.")
    parser.add_argument("-s", "--start_date", target="start_date", default="a year ago midnight",
                        help="The start date in a format that python.dateparser can read. Default is 'a year ago midnight'.")
    parser.add_argument("-e", "--end_date", default="midnight today",
                        help="The end date in a format that python.dateparser can read. Default is 'midnight today'.")
    parser.add_argument("-cols", "--other_columns", dest="other_columns", default=None,
                        help="A comma-separated list of extra ICESat-2 ATL03 data columns to include in the output. "
                             "Default is only necessary columns.")
    parser.add_argument("-conf", "--confidence_levels", dest="conf_levels", default="4",
                        help="A comma-separated (or /-separated) list of confidence levels to include in the output. "
                             "Values can be 1, 2, 3, 4. Default is 4 (highest-confidence photons only).")
    parser.add_argument("-cls", "--classes", dest="classes", default="1,2,3,4",
                        help="A comma-separated (or /-separated) list of classifications to include in the output. "
                             "Values can be -1, 0, 1, 2, 3, 4, 5, 6, 7. Run 'dlim --modules -214' for a full list."
                             "Default is 1,2,3,4.")
    parser.add_argument("-nobuild", "--no_buildings", dest="buildings", default=True, actions="store_false",
                        help="Do not include building classifications in the output. Default: Use Bing to classify "
                             "buildings as class_code 7.")
    parser.add_argument("-nobathy", "--no_bathymetry", dest="bathymetry", default=True, actions="store_false",
                        help="Do not include bathymetry classifications in the output. Default: "
                             "Use CShelph to classify bathymetry as class_code 4 (bathy floor) and 5 (bathy surface).")
    parser.add_argument("-l", "--lines", dest="lines", default=None,
                        help="The name of a vector file where simplified lines following the path of each laser will be exported.")

    args = parser.parse_args()
    return args


def main():
    args = define_and_parse_args()
    if args.region is None and args.input_raster is None:
        print("One of --region or --input_raster must be set.", file=sys.stderr)
        sys.exit(1)
    if args.region is not None and args.input_raster is not None:
        print("Only one of --region or --input_raster can be set.", file=sys.stderr)
        sys.exit(1)
    if not args.crs and not args.raster:
        print("At least one of --crs or --raster must be set.", file=sys.stderr)
        sys.exit(1)

    if args.region:
        args.region = [float(x) for x in args.region.split(",")]

    if args.other_columns:
        args.other_columns = ast.literal_eval(args.other_columns)

    gdf = get_photon_dataframe(args.raster if args.raster else args.region,
                               dem_horz_reference_frame=args.crs,
                               dem_vert_reference_frame=args.vdatum,
                               start_date=args.start_date,
                               end_date=args.end_date,
                               other_columns=args.other_columns,
                               classify_bathymetry=args.bathymetry,
                               classify_buildings=args.buildings,
                               classifications_to_keep=args.classes,
                               conf_levels_to_keep=args.conf_levels
                               )

    export_as_vector(gdf, args.outfile)

    if args.lines:
        export_lines(gdf, args.lines)


if __name__ == "__main__":
    main()
