### nceithredds.py
##
## Copyright (c) 2010 - 2026 Regents of the University of Colorado
##
## nceithredds.py is part of CUDEM
##
## Permission is hereby granted, free of charge, to any person obtaining a copy 
## of this software and associated documentation files (the "Software"), to deal 
## in the Software without restriction, including without limitation the rights 
## to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies 
## of the Software, and to permit persons to whom the Software is furnished to do so, 
## subject to the following conditions:
##
## The above copyright notice and this permission notice shall be included in all
## copies or substantial portions of the Software.
##
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
## INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR 
## PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE 
## FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
## ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
## SOFTWARE.
##
### Commentary:
##
## Fetch NOAA NCEI DEMs via THREDDS Catalog.
##
### Code:

import os
from typing import List, Dict, Optional, Any
from urllib.parse import urljoin, urlencode
from cudem import utils
from cudem.fetches import fetches
from cudem.fetches import FRED

## ==============================================
## Constants
## ==============================================
NGDC_BASE_URL = "https://www.ngdc.noaa.gov"
THREDDS_CATALOG_URL = "https://www.ngdc.noaa.gov/thredds/catalog/demCatalog.xml"

## ==============================================
## NCEI THREDDS Module
## ==============================================
class NCEIThreddsCatalog(fetches.FetchModule):
    """NOAA NCEI DEMs via THREDDS

    Fetch DEMs from NCEI THREDDS Catalog.
    
    Digital Elevation Models around the world at various resolutions and extents.

    https://www.ngdc.noaa.gov/thredds/demCatalog.xml

    Configuration Example:
    < ncei_thredds:where=None:want_wcs=False >
    """

    def __init__(self, where: str = '', want_wcs: bool = False, datatype: Optional[str] = None, **kwargs):
        super().__init__(name='ncei_thredds', **kwargs)
        self.where = [where] if where else []
        self.want_wcs = want_wcs
        self.datatype = datatype
        
        if self.datatype:
            self.where.append(f"ID LIKE '%{self.datatype}%'")

        ## Initialize FRED
        self.FRED = FRED.FRED(name=self.name, verbose=self.verbose)
        self.update_if_not_in_FRED()

        
    def update_if_not_in_FRED(self):
        """Update the vector reference in FRED if not present."""
        
        self.FRED._open_ds()
        try:
            self.FRED._attribute_filter([f"DataSource = '{self.name}'"])
            if len(self.FRED.layer) == 0:
                ## Close before update to allow write access
                self.FRED._close_ds()
                self.update()
        finally:
            if self.FRED.ds is not None:
                self.FRED._close_ds()

                
    def _get_full_url(self, base_url: str, relative_url: str) -> str:
        """Construct full URL handling absolute/relative paths correctly."""
        
        if relative_url.startswith("/"):
            return f"{NGDC_BASE_URL}{relative_url}"
        else:
            return f"{os.path.dirname(base_url)}/{relative_url}"

        
    def _extract_services(self, dataset_node, services_node) -> Dict[str, str]:
        """Extract service URLs (ISO, WCS, HTTP) for a dataset."""
        
        urls = {'iso': None, 'wcs': None, 'http': None}
        try:
            ds_path = dataset_node.attrib['urlPath']
        except KeyError:
            return urls

        for service in services_node:
            s_name = service.attrib.get('name')
            s_base = service.attrib.get('base')
            
            if s_name in urls and s_base:
                urls[s_name] = f"{NGDC_BASE_URL}{s_base}{ds_path}"
                
        return urls

    
    def _parse_dataset(self, catalog_url: str) -> List[Dict]:
        """Parse a specific dataset XML catalog."""
        
        surveys = []
        ntCatXml = fetches.iso_xml(catalog_url)
        
        if ntCatXml.xml_doc is None:
            return surveys

        datasets = ntCatXml.xml_doc.findall('.//th:dataset', namespaces=ntCatXml.namespaces)
        services = ntCatXml.xml_doc.findall('.//th:service', namespaces=ntCatXml.namespaces)

        if not datasets:
            return surveys

        desc = datasets[0].attrib.get('name', 'Unknown Dataset')
        
        with utils.ccp(total=len(datasets), desc=f'Scanning {desc}', leave=self.verbose) as pbar:
            for node in datasets:
                pbar.update(1)
                
                this_id = node.attrib.get('ID')
                if not this_id:
                    continue

                ## Check if already exists in FRED
                self.FRED._attribute_filter([f"ID = '{this_id}'"])
                if self.FRED.layer is not None and len(self.FRED.layer) > 0:
                    continue

                ## Handle Sub-Catalogs (recursion)
                subCatRefs = node.findall('.//th:catalogRef', namespaces=ntCatXml.namespaces)
                if subCatRefs:
                    self._parse_catalog(catalog_url) 
                    break

                ## Extract Service URLs
                service_urls = self._extract_services(node, services)
                iso_url = service_urls['iso']

                if not iso_url:
                    continue

                ## Parse ISO Metadata
                try:
                    this_xml = fetches.iso_xml(iso_url)
                    if this_xml.xml_doc is None:
                        continue
                        
                    title = this_xml.title()
                    h_epsg, v_epsg = this_xml.reference_system()
                    geom = this_xml.bounds(geom=True)
                    
                    ## Determine Z Variable
                    zvar = 'z'
                    zv_nodes = this_xml.xml_doc.findall(
                        './/gmd:dimension/gmd:MD_Band/gmd:sequenceIdentifier/gco:MemberName/gco:aName/gco:CharacterString',
                        namespaces=this_xml.namespaces
                    )
                    if zv_nodes:
                        for zvs in zv_nodes:
                            if zvs.text in ['bathy', 'Band1', 'z']:
                                zvar = zvs.text
                                break

                    if geom is not None:
                        surveys.append({
                            'Name': title,
                            'ID': this_id,
                            'Agency': 'NOAA',
                            'Date': this_xml.date(),
                            'MetadataLink': this_xml.url,
                            'MetadataDate': this_xml.xml_date(),
                            'DataLink': service_urls['http'],
                            'IndexLink': service_urls['wcs'],
                            'Link': THREDDS_CATALOG_URL,
                            'DataType': 'raster',
                            'DataSource': 'ncei_thredds',
                            'HorizontalDatum': h_epsg,
                            'VerticalDatum': v_epsg,
                            'Etcetra': zvar,
                            'Info': this_xml.abstract(),
                            'geom': geom
                        })
                except Exception as e:
                    if self.verbose:
                        utils.echo_warning_msg(f"Error parsing metadata for {this_id}: {e}")
        return surveys

    def _parse_catalog(self, catalog_url: str):
        """Recursively parse THREDDS catalogs."""
        
        ntCatalog = fetches.iso_xml(catalog_url)
        if ntCatalog.xml_doc is None:
            return

        ## Find Catalog References
        ntCatRefs = ntCatalog.xml_doc.findall('.//th:catalogRef', namespaces=ntCatalog.namespaces)
        
        for ntCatRef in ntCatRefs:
            href = ntCatRef.attrib.get('{http://www.w3.org/1999/xlink}href')
            if not href:
                continue

            full_cat_url = self._get_full_url(catalog_url, href)
            
            ## Parse Datasets in this catalog
            surveys = self._parse_dataset(full_cat_url)
            if surveys:
                self.FRED._add_surveys(surveys)

                
    def update(self):
        """Scan the THREDDS Catalog and fill the reference vector in FRED."""
        
        self.FRED._open_ds(1)
        try:
            self._parse_catalog(THREDDS_CATALOG_URL)
        except Exception as e:
            utils.echo_error_msg(f"Error updating NCEI THREDDS FRED: {e}")
        finally:
            self.FRED._close_ds()

            
    def run(self):
        """Search for data in the reference vector file."""
        
        _results = FRED._filter_FRED(self)
        
        for surv in _results:
            if self.want_wcs and surv.get('IndexLink'):
                ## WCS Request Construction
                wcs_params = {
                    'request': 'GetCoverage',
                    'version': '1.0.0',
                    'service': 'WCS',
                    'coverage': surv.get('Etcetra', 'z'),
                    'bbox': self.region.format('bbox'),
                    'format': 'geotiff_float'
                }
                
                query_string = urlencode(wcs_params)
                wcs_url = f"{surv['IndexLink']}?{query_string}"
                
                ## Infer filename from DataLink (NetCDF) -> GeoTIFF
                out_fn = surv['DataLink'].split(',')[0].split('/')[-1].replace('.nc', '.tif')
                
                self.add_entry_to_results(wcs_url, out_fn, surv['DataType'])
            
            elif surv.get('DataLink'):
                ## Direct HTTP Download
                for link in surv['DataLink'].split(','):
                    if link.strip():
                        self.add_entry_to_results(
                            link, 
                            link.split('/')[-1], 
                            surv['DataType']
                        )

        return self

    
### End
