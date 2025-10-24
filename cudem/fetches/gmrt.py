### gmrt.py
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

from cudem.fetches import fetches

###############################################################################
## GMRT
###############################################################################
def gmrt_fetch_point(latitude=None, longitude=None):
    gmrt_point_url = "https://www.gmrt.org:443/services/PointServer?"
    headers = {
        'User-Agent': ('Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:89.0) '
                       'Gecko/20100101 Firefox/89.0')
    }
    data = {'longitude':longitude, 'latitude':latitude}
    req = fetches.Fetch(gmrt_point_url).fetch_req(
        params=data, tries=10, timeout=2
    )
    if req is not None:
        return(req.text)
    else:
        return(None)

    
class GMRT(fetches.FetchModule):
    """The Global Multi-Resolution Topography synthesis.
    
    The Global Multi-Resolution Topography (GMRT) synthesis is a multi-resolutional 
    compilation of edited multibeam sonar data collected by scientists and 
    institutions worldwide, that is reviewed, processed and gridded by the GMRT 
    Team and merged into a single continuously updated compilation of global elevation 
    data. The synthesis began in 1992 as the Ridge Multibeam Synthesis (RMBS), 
    was expanded to include multibeam bathymetry data from the Southern Ocean, 
    and now includes bathymetry from throughout the global and coastal oceans.

    Data Formats
    GMT v3 Compatible NetCDF (GMT id=cf)
    COARDS/CF1.6 Compliant NetCDF (GMT id=nd)
    ESRI ArcASCII
    GeoTIFF

    Metadata Formats
    XML (metadata)
    JSON (metadata)
    Plain text (metadata)
    
    layers: 'topo' or 'topo-mask'
    fmt: 'geotiff', 'netcdf'
    
    Data is assumed instantaneous MSL (5773?)
    
    https://www.gmrt.org

    < gmrt:res=max:fmt=geotiff:layer=topo >
    """
    
    def __init__(
            self,
            res='default',
            fmt='netcdf',
            layer='topo',
            want_swath=False,
            **kwargs
    ):
        super().__init__(name='gmrt', **kwargs)
        # GMRT resolution
        self.res = res 
        # GMRT format
        self.fmt = fmt 
        # fetch the swath vector along with the data, used to clip non-swath data
        self.want_swath = want_swath
        # GMRT layer
        self.layer = 'topo' \
            if (layer != 'topo' and layer != 'topo-mask') \
               else layer 

        ## The various urls to use for GMRT
        self._gmrt_url = 'https://www.gmrt.org'
        self._gmrt_grid_url = "https://www.gmrt.org:443/services/GridServer?"
        self._gmrt_grid_urls_url = "https://www.gmrt.org:443/services/GridServer/urls?"
        self._gmrt_grid_metadata_url = "https://www.gmrt.org/services/GridServer/metadata?"        
        self._gmrt_swath_poly_url = "https://www.gmrt.org/shapefiles/gmrt_swath_polygons.zip"

        ## buffer the input region and correct to wgs extremes
        self.gmrt_region = self.region.copy()
        self.gmrt_region.buffer(pct=2.33,x_inc=.0088,y_inc=.0088)
        self.gmrt_region._wgs_extremes(just_below=True)

        ## dlim variables, parse with GDAL and set to WGS84/MSL
        self.data_format = 200
        self.src_srs = 'epsg:4326+3855'
        self.title = 'GMRT'
        self.source = 'GMRT'
        self.date = None
        self.data_type = 'Raster'
        self.resolution = None
        self.hdatum = 'wgs84'
        self.vdatum = 'msl'
        self.url = self._gmrt_url

        ## Firefox on windows for this one.
        self.headers = {
            'User-Agent': ('Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:89.0) '
                           'Gecko/20100101 Firefox/89.0')
        }

        
    def run(self):
        '''Run the GMRT fetching module'''

        if self.region is None:
            return([])

        self.data = {
            'north': self.gmrt_region.ymax,
            'west': self.gmrt_region.xmin,
            'south': self.gmrt_region.ymin,
            'east': self.gmrt_region.xmax,
            'mformat': 'json',
            'resolution': self.res,
            'format': self.fmt,
            'layer': self.layer
        }
        req = fetches.Fetch(
            self._gmrt_grid_url
        ).fetch_req(
            params=self.data, tries=10, timeout=2
        )

        if req is not None:
            outf = 'gmrt_{}_{}_{}.{}'.format(
                self.layer,
                self.res,
                self.region.format('fn_full'),
                'tif' if self.fmt == 'geotiff' else 'grd'
            )
            self.add_entry_to_results(req.url, outf, 'gmrt')
        # else:
        #     ## we got multiple URLs, so lets loop through those
        #     ## and fetch them individually
        #     gmrt_urls = req.json()
        #     for url in gmrt_urls:
        #         if self.layer == 'topo-mask':
        #             url = url.replace('topo', 'topo-mask')

        #         opts = {}
        #         for url_opt in url.split('?')[1].split('&'):
        #             opt_kp = url_opt.split('=')
        #             opts[opt_kp[0]] = opt_kp[1]

        #         url_region = regions.Region().from_list([
        #             float(opts['west']),
        #             float(opts['east']),
        #             float(opts['south']),
        #             float(opts['north'])
        #         ])
        #         outf = 'gmrt_{}_{}.{}'.format(
        #             opts['layer'],
        #             url_region.format('fn'),
        #             'tif' if self.fmt == 'geotiff' else 'grd'
        #         )
        #         self.add_entry_to_results(url, outf, 'gmrt')

        #         ## if want_swath is True, we will download the swath
        #         ## polygons so that we can clip the data to that in
        #         ## dlim or elsewhere.
        #         if self.want_swath:
        #             self.add_entry_to_results(
        #                 self._gmrt_swath_poly_url,
        #                 'gmrt_swath_polygons.zip',
        #                 'gmrt'
        #             )
                
        return(self)

### End
