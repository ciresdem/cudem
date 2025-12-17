### emodnet.py
##
## Copyright (c) 2010 - 2025 Regents of the University of Colorado
##
## fetches.py is part of CUDEM
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
###############################################################################
### Commentary:
##
##
### Code:

import lxml.etree
from cudem.fetches import fetches

## EMODNet - EU data
class EMODNet(fetches.FetchModule):
    """EU elevation data extracts from EMOD DTM.

    Fetch raster data from the EMODNET DTM
    
    https://portal.emodnet-bathymetry.eu/
    https://erddap.emodnet.eu

    erddap formats (default is nc):
    https://erddap.emodnet.eu/erddap/griddap/documentation.html#fileType

    Data
    fileTypes	Description
    .asc	View OPeNDAP-style ISO-8859-1 comma-separated text.
    .csv	Download a ISO-8859-1 comma-separated text table (line 1: names; line 2: units; ISO 8601 times).
    .csvp	Download a ISO-8859-1 .csv file with line 1: name (units). Times are ISO 8601 strings.
    .csv0	Download a ISO-8859-1 .csv file without column names or units. Times are ISO 8601 strings.
    .das	View the dataset's metadata via an ISO-8859-1 OPeNDAP Dataset Attribute Structure (DAS).
    .dds	View the dataset's structure via an ISO-8859-1 OPeNDAP Dataset Descriptor Structure (DDS).
    .dods	OPeNDAP clients use this to download the data in the DODS binary format.
    .esriAscii	Download an ISO-8859-1 ESRI ASCII file (latitude longitude data only; longitude must be all below or all above 180).
    .fgdc	View the dataset's UTF-8 FGDC .xml metadata.
    .graph	View a Make A Graph web page.
    .help	View a web page with a description of griddap.
    .html	View an OPeNDAP-style HTML Data Access Form.
    .htmlTable	View a UTF-8 .html web page with the data in a table. Times are ISO 8601 strings.
    .iso19115	View the dataset's ISO 19115-2/19139 UTF-8 .xml metadata.
    .itx	Download an ISO-8859-1 Igor Text File. Each axis variable and each data variable becomes a wave.
    .json	View a table-like UTF-8 JSON file (missing value = 'null'; times are ISO 8601 strings).
    .jsonlCSV1	View a UTF-8 JSON Lines CSV file with column names on line 1 (mv = 'null'; times are ISO 8601 strings).
    .jsonlCSV	View a UTF-8 JSON Lines CSV file without column names (mv = 'null'; times are ISO 8601 strings).
    .jsonlKVP	View a UTF-8 JSON Lines file with Key:Value pairs (missing value = 'null'; times are ISO 8601 strings).
    .mat	Download a MATLAB binary file.
    .nc	Download a NetCDF-3 binary file with COARDS/CF/ACDD metadata.
    .ncHeader	View the UTF-8 header (the metadata) for the NetCDF-3 .nc file.
    .ncml	View the dataset's structure and metadata as a UTF-8 NCML .xml file.
    .nccsv	Download a NetCDF-3-like 7-bit ASCII NCCSV .csv file with COARDS/CF/ACDD metadata.
    .nccsvMetadata	View the dataset's metadata as the top half of a 7-bit ASCII NCCSV .csv file.
    .ncoJson	Download a UTF-8 NCO lvl=2 JSON file with COARDS/CF/ACDD metadata.
    .odvTxt	Download time,latitude,longitude,otherVariables as an ODV Generic Spreadsheet File (.txt).
    .timeGaps	View a UTF-8 list of gaps in the time values which are larger than the median gap.
    .tsv	Download a ISO-8859-1 tab-separated text table (line 1: names; line 2: units; ISO 8601 times).
    .tsvp	Download a ISO-8859-1 .tsv file with line 1: name (units). Times are ISO 8601 strings.
    .tsv0	Download a ISO-8859-1 .tsv file without column names or units. Times are ISO 8601 strings.
    .wav	Download a .wav audio file. All columns must be numeric and of the same type.
    .xhtml	View a UTF-8 XHTML (XML) file with the data in a table. Times are ISO 8601 strings.

    < emodnet:want_erddap=True:erddap_format=nc >
    """

    def __init__(self, want_erddap = False, erddap_format = 'nc', **kwargs):
        super().__init__(name='emodnet', **kwargs)
        self.want_erddap = want_erddap
        self.erddap_format = erddap_format

        ## Emodnet URLs
        self._emodnet_grid_url = 'https://ows.emodnet-bathymetry.eu/wcs?'
        self._emodnet_grid_url_erddap = ('https://erddap.emodnet.eu/erddap/griddap/'
                                         f'dtm_2020_v2_e0bf_e7e4_5b8f.{erddap_format}?')

        ## for dlim
        self.src_srs = 'epsg:4326'

        
    def run(self):
        """Run the EMODNET fetching module"""
        
        if self.region is None: return([])

        if self.want_erddap:
            suff = 'elevation%5B({}):1:({})%5D%5B({}):1:({})%5D'.format(
                self.region.ymin, self.region.ymax, self.region.xmin, self.region.xmax
            )
            erddap_url = self._emodnet_grid_url_erddap + suff
            outf = 'emodnet_{}.{}'.format(self.region.format('fn'), self.erddap_format)
            self.add_entry_to_results(erddap_url, outf, self.erddap_format)
            
        else:
            _data = {
                'request': 'DescribeCoverage',
                'version': '2.0.1',
                'CoverageID': 'emodnet:mean',
                'service': 'WCS'
            }
            _req = fetches.Fetch(self._emodnet_grid_url).fetch_req(params=_data)
            _results = lxml.etree.fromstring(_req.text.encode('utf-8'))
            g_env = _results.findall(
                './/{http://www.opengis.net/gml/3.2}GridEnvelope',
                namespaces=fetches.namespaces
            )[0]
            hl = [float(x) for x in g_env.find(
                '{http://www.opengis.net/gml/3.2}high'
            ).text.split()]
            g_bbox = _results.findall(
                './/{http://www.opengis.net/gml/3.2}Envelope'
            )[0]
            lc = [float(x) for x in  g_bbox.find(
                '{http://www.opengis.net/gml/3.2}lowerCorner'
            ).text.split()]
            uc = [float(x) for x in g_bbox.find(
                '{http://www.opengis.net/gml/3.2}upperCorner'
            ).text.split()]
            ds_region = regions.Region().from_list(
                [lc[1], uc[1], lc[0], uc[0]]
            )
            resx = (uc[1] - lc[1]) / hl[0]
            resy = (uc[0] - lc[0]) / hl[1]
            if regions.regions_intersect_ogr_p(self.region, ds_region):
                _wcs_data = {
                    'service': 'WCS',
                    'request': 'GetCoverage',
                    'version': '1.0.0',
                    'Identifier': 'emodnet:mean',
                    'coverage': 'emodnet:mean',
                    'format': 'GeoTIFF',
                    'bbox': self.region.format('bbox'),
                    'resx': resx,
                    'resy': resy,
                    'crs': 'EPSG:4326',
                }
                emodnet_wcs = '{}{}'.format(self._emodnet_grid_url, urlencode(_wcs_data))
                outf = 'emodnet_{}.tif'.format(self.region.format('fn'))
                self.add_entry_to_results(emodnet_wcs, outf, 'emodnet')
            
        return(self)

### End
