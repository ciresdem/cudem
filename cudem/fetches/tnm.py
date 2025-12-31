### tnm.py
##
## Copyright (c) 2010 - 2025 Regents of the University of Colorado
##
## tnm.py is part of CUDEM
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
## Fetch elevation data from The National Map (TNM) API.
##
### Code:

import datetime
from typing import Optional, List, Any
from cudem import utils
from cudem.fetches import fetches

## ==============================================
## Constants
## ==============================================
TNM_API_PRODUCTS_URL = 'http://tnmaccess.nationalmap.gov/api/v1/products?'

DATASET_CODES = [
    "National Boundary Dataset (NBD)",
    "National Elevation Dataset (NED) 1 arc-second",
    "Digital Elevation Model (DEM) 1 meter",
    "National Elevation Dataset (NED) 1/3 arc-second",
    "National Elevation Dataset (NED) 1/9 arc-second",
    "National Elevation Dataset (NED) Alaska 2 arc-second",
    "Alaska IFSAR 5 meter DEM",
    "National Elevation Dataset (NED) 1/3 arc-second - Contours",
    "Original Product Resolution (OPR) Digital Elevation Model (DEM)",
    "Ifsar Digital Surface Model (DSM)",
    "Ifsar Orthorectified Radar Image (ORI)",
    "Lidar Point Cloud (LPC)",
    "Historical Topographic Maps",
    "National Hydrography Dataset Plus High Resolution (NHDPlus HR)",
    "National Hydrography Dataset (NHD) Best Resolution",
    "National Watershed Boundary Dataset (WBD)",
    "Map Indices",
    "National Geographic Names Information System (GNIS)",
    "Small-scale Datasets - Boundaries",
    "Small-scale Datasets - Contours",
    "Small-scale Datasets - Hydrography",
    "Small-scale Datasets - Transportation",
    "National Structures Dataset (NSD)",
    "Combined Vector",
    "National Transportation Dataset (NTD)",
    "US Topo Current",
    "US Topo Historical",
    "Land Cover - Woodland",
    "3D Hydrography Program (3DHP)",
]

FORMAT_KEYWORDS = [
    "ArcExport", "ArcGrid", "BIL", "FileGDB", "FileGDB 10.1", "FileGDB 10.2",
    "GeoPDF", "GeoTIFF", "GridFlow", "IMG", "JPEG2000", "LAS,LAZ", "NLAPS",
    "PDF", "SDE Export", "Shapefile", "Text", "TIFF", "TXT (pipes)",
]

DATASET_KEYWORDS = [
    "7.5", "Airports", "Bend", "Bridge", "Building", "Canal", "Cape", "Cave",
    "Cemetery", "Census", "Channel", "Chamber", "Church", "Civil", "Cliff",
    "Coast", "Coastline", "Conduit", "Contour", "Crossing", "Dam", "Dams",
    "Ditch", "Elevation", "Falls", "Flat", "Flume", "Forest", "Gaging", "Gaging Station",
    "Gap", "Gate", "Glacier", "Gut", "HUC", "Harbor", "Hospital", "Hydrography",
    "Hydrologic", "Images", "Intake", "Island", "Isthmus", "Lake", "Lakes", "Lava",
    "Levee", "Lidar", "Lock", "Map", "Maps", "Military", "Mine", "Oilfield", "Outflow",
    "Park", "Pillar", "Pipeline", "Plain", "Populated Place", "Post Office", "Range",
    "Rapids", "Reef", "Reserve", "Reservoir", "Ridge", "Rise", "River", "Rock",
    "School", "Sea", "Seep", "Shore", "Sink", "Slope", "Spring", "Stream", "Summit",
    "Swamp", "Tower", "Trail", "Tunnel", "Valley", "Water", "Waterfall", "Watershed",
    "Weir", "Well", "Woods",
]

DATASET_EXTENTS = [
    "10000 x 10000 meter", "1500 x 1500 meter", "15 x 15 minute", "1 x 1 degree",
    "1 x 2 degree", "1 x 3 degree", "1 x 4 degree", "2 x 1 degree", "30 x 30 minute",
    "30 x 60 minute", "3.75 minute x 3.75 minute", "3 x 3 degree", "7.5 x 15 minute",
    "7.5 x 7.5 minute", "Contiguous US", "HU-2 Region", "HU-4 Subregion", "HU-8 Subbasin",
    "National", "North America", "State", "Varies",
]

DATE_TYPES = ["dateCreated", "lastUpdated", "Publication"]

## ==============================================
## The National Map Module
## ==============================================
class TheNationalMap(fetches.FetchModule):
    """The National Map:

    Fetch elevation data from The National Map.
    
    http://tnmaccess.nationalmap.gov/

    < tnm:datasets=None:formats=None:extents=None:q=None:date_type=None:date_start=None:date_end=None >
    """
    
    ## Generate dynamic docstring for help menus
    __doc__ = f'''{__doc__}
    Dataset codes (datasets):
    {utils.list_str(DATASET_CODES)}

    Format keywords (formats):
    {FORMAT_KEYWORDS}

    Dataset keywords (q):
    {DATASET_KEYWORDS}

    Dataset extents (extents):
    {DATASET_EXTENTS}

    Date types (date_types):
    {DATE_TYPES}
    '''

    def __init__(
            self, datasets: Optional[str] = None, formats: Optional[str] = None, 
            extents: Optional[str] = None, q: Optional[str] = None,
            date_type: Optional[str] = None, date_start: Optional[str] = None, 
            date_end: Optional[str] = None, **kwargs
    ):
        super().__init__(name='tnm', **kwargs)
        self.q = q
        self.formats = formats
        self.extents = extents
        self.datasets = datasets
        self.date_type = date_type
        self.date_start = date_start
        self.date_end = date_end

        
    def run(self):
        """Run the TNM fetching module."""
        
        if self.region is None:
            return []

        offset = 0
        total = 0
        
        while True:
            # Base Parameters
            params = {
                'bbox': self.region.format('bbox'),
                'max': 100,
                'offset': offset
            }

            ## Handle Datasets
            ## User passes slash-separated indices (e.g., "1/3") mapping to DATASET_CODES
            if self.datasets is not None:
                try:
                    ## Parse indices
                    ds_indices = [int(x) for x in self.datasets.split('/')]
                    ## Map to names
                    selected_datasets = [DATASET_CODES[i] for i in ds_indices if 0 <= i < len(DATASET_CODES)]
                    params['datasets'] = ','.join(selected_datasets)
                except (ValueError, IndexError) as e:
                    utils.echo_warning_msg(f"Could not parse datasets '{self.datasets}': {e}")
                    ## Fallback default
                    params['datasets'] = "National Elevation Dataset (NED) 1 arc-second"

            ## Handle Other Filters
            if self.q is not None:
                params['q'] = str(self.q)
            
            if self.formats is not None:
                ## API expects comma-separated
                params['prodFormats'] = ','.join(self.formats.split('/'))
            
            if self.extents is not None:
                params['prodExtents'] = ','.join(self.extents.split('/'))

            ## Handle Dates
            if self.date_start is not None:
                params['start'] = self.date_start
                params['end'] = self.date_end if self.date_end else datetime.datetime.now().strftime('%Y-%m-%d')
                params['dateType'] = self.date_type if self.date_type else 'dateCreated'

            ## Execute Request
            req = fetches.Fetch(
                TNM_API_PRODUCTS_URL, 
                verbose=self.verbose
            ).fetch_req(params=params, timeout=60, read_timeout=60)
            
            if req is not None and req.status_code == 200:
                utils.echo_debug_msg(req.url)
                
                ## Check for API error message in body
                if req.text.startswith("{errorMessage"):
                    utils.echo_error_msg(f"TNM API Error: {req.text}")
                    break
                
                try:
                    features = req.json()
                    total = features.get('total', 0)
                    items = features.get('items', [])

                    ## Double-Check the bounding-box in case TNM
                    ## returned data outside of our region.
                    if 'boundingBox' in feature:
                        bb = feature['boundingBox']
                        feature_region = regions.Region().from_list([
                            bb.get('minX'), bb.get('maxX'), 
                            bb.get('minY'), bb.get('maxY')
                        ])
                        
                        if not regions.regions_intersect_p(self.region, feature_region):
                            continue
                    
                    for feature in items:
                        self.add_entry_to_results(
                            feature['downloadURL'],
                            feature['downloadURL'].split('/')[-1],
                            feature['format']
                        )
                except Exception as e:
                     utils.echo_error_msg(f"Error parsing TNM JSON: {e}")
                     break

            ## Pagination
            offset += 100
            if offset >= total:
                break
        
        return self

    
## ==============================================
## Shortcuts
## ==============================================
class NED(TheNationalMap):
    """National Elevation Dataset (NED) 1 & 1/3 arc-second via TNM."""

    def __init__(self, **kwargs):
        ## Indices 1 and 3 in DATASET_CODES correspond to NED 1 and NED 1/3
        super().__init__(datasets='1/3', **kwargs)
        self.data_format = 200


class NED1(TheNationalMap):
    """National Elevation Dataset (NED) 1 meter via TNM."""
    
    def __init__(self, **kwargs):
        ## Index 2 in DATASET_CODES corresponds to DEM 1 meter
        super().__init__(datasets='2', **kwargs)
        self.data_format = 200


class TNM_LAZ(TheNationalMap):
    """Lidar (LAZ) via TNM."""

    def __init__(self, **kwargs):
        super().__init__(formats="LAZ", **kwargs)
        self.data_format = 300

### End
