import json
from cudem import utils
from cudem.fetches import fetches

## NOAA DEMs
## doesn't really work well, use ncei_thredds or digital_coast instead...
class DEMMosaic(fetches.FetchModule):
    """
    Fields:

    OBJECTID ( type: esriFieldTypeOID, alias: OBJECTID )
    Shape ( type: esriFieldTypeGeometry, alias: Shape )
    Name ( type: esriFieldTypeString, alias: Name, length: 200 )
    MinPS ( type: esriFieldTypeDouble, alias: MinPS )
    MaxPS ( type: esriFieldTypeDouble, alias: MaxPS )
    LowPS ( type: esriFieldTypeDouble, alias: LowPS )
    HighPS ( type: esriFieldTypeDouble, alias: HighPS )
    Category ( type: esriFieldTypeInteger, alias: Category , 
      Coded Values: [0: Unknown] , [1: Primary] , [2: Overview] , ...6 more... )
    Tag ( type: esriFieldTypeString, alias: Tag, length: 100 )
    GroupName ( type: esriFieldTypeString, alias: GroupName, length: 100 )
    ProductName ( type: esriFieldTypeString, alias: ProductName, length: 100 )
    CenterX ( type: esriFieldTypeDouble, alias: CenterX )
    CenterY ( type: esriFieldTypeDouble, alias: CenterY )
    ZOrder ( type: esriFieldTypeInteger, alias: ZOrder )
    Shape_Length ( type: esriFieldTypeDouble, alias: Shape_Length )
    Shape_Area ( type: esriFieldTypeDouble, alias: Shape_Area )
    ZOrder_1 ( type: esriFieldTypeInteger, alias: ZOrder_1 )
    DEM_ID ( type: esriFieldTypeSmallInteger, alias: DEM_ID )
    DateCompleted ( type: esriFieldTypeDate, alias: DateCompleted, length: 8 )
    CellsizeArcseconds ( type: esriFieldTypeSingle, alias: CellsizeArcseconds )
    DemName ( type: esriFieldTypeString, alias: DemName, length: 100 )
    MetadataURL ( type: esriFieldTypeString, alias: MetadataURL, length: 250 )
    VerticalDatum ( type: esriFieldTypeString, alias: VerticalDatum, length: 50 )

    https://gis.ngdc.noaa.gov/arcgis/rest/services/DEM_mosaics/DEM_global_mosaic/ImageServer

    ** doesn't really work well, use ncei_thredds or digital_coast instead...
    """
    
    def __init__(self, where = '1=1', layer = 1, index = False, **kwargs):
        super().__init__(name='demmosaic', **kwargs)
        self.where = where
        self.index = index
        
        ## The various DEMMosaic URLs
        self._dem_mosaic_url = ('https://gis.ngdc.noaa.gov/arcgis/rest/services/'
                                'DEM_mosaics/DEM_all/ImageServer')
        #self._dem_mosaic_query_url = '{0}/{1}/query?'.format(self._dem_mosaic_url, layer)
        self._dem_mosaic_query_url = '{0}/query?'.format(self._dem_mosaic_url)
        
        ## for dlim
        self.data_format = None # bag/xyz data are different, reset later
        self.src_srs = None # dems vary, set later

        
    def run(self):
        """Run the DEMMosaic fetching module"""
        
        if self.region is None:
            return([])

        _data = {
            'where': self.where,
            'outFields': '*',
            'geometry': self.region.format('bbox'),
            'inSR':4326,
            'outSR':4326,
            'f':'pjson',
            'returnGeometry':'False',
            'geometryType':'esriGeometryEnvelope',
            'spatialRel':'esriSpatialRelIntersects'
        }
        _req = fetches.Fetch(
            self._dem_mosaic_query_url,
            verbose=self.verbose
        ).fetch_req(params=_data)
        if _req is not None:
            features = _req.json()
            for feature in features['features']:
                if self.index:
                    print(json.dumps(feature['attributes'], indent=4))
                else:
                    #utils.echo_msg(feature['attributes']['MetadataURL'])
                    print(feature)
                    Name = feature['attributes']['Name']
                    ID = feature['attributes']['DEM_ID']
                    link = feature['attributes']['MetadataURL']
                    if link is not None:
                        utils.echo_msg(Name)
                        utils.echo_msg(ID)
                        utils.echo_msg(link)                        
                        page = fetches.Fetch(link).fetch_xml()
                        print(page)
                        sys.exit()
