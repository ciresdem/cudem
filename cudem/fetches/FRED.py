### FRED.py
##
## Copyright (c) 2010 - 2021 CIRES Coastal DEM Team
##
## FRED.py is part of CUDEM
##
## Permission is hereby granted, free of charge, to any person obtaining a copy 
## of this software and associated documentation files (the "Software"), to deal 
## in the Software without restriction, including without limitation the rights 
## to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies 
## of the Software, and to permit persons to whom the Software is furnished to do so, 
## subject to the following conditions:
##
## The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
##
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
## INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR 
## PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE 
## FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
## ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
##
### Commentary:
### Code:

import os
import json
import lxml.etree
from osgeo import ogr
from cudem import utils
from cudem import regions
import cudem.fetches.utils as f_utils

## =============================================================================
##
## XML Metadata parsing
##
## =============================================================================
def xml2py(node):
    '''parse an xml file into a python dictionary'''
    texts = {}
    if node is None: return(None)
    for child in list(node):
        child_key = lxml.etree.QName(child).localname
        if 'name' in child.attrib.keys(): child_key = child.attrib['name']
        if '{http://www.w3.org/1999/xlink}href' in child.attrib.keys():
            href = child.attrib['{http://www.w3.org/1999/xlink}href']
        else: href = None
        if child.text is None or child.text.strip() == '':
            if href is not None:
                if child_key in texts.keys():
                    texts[child_key].append(href)
                else: texts[child_key] = [href]
            else:
                if child_key in texts.keys():
                    ck = xml2py(child)
                    texts[child_key][list(ck.keys())[0]].update(ck[list(ck.keys())[0]])
                else: texts[child_key] = xml2py(child)
        else:
            if child_key in texts.keys():
                texts[child_key].append(child.text)
            else: texts[child_key] = [child.text]
    return(texts)

class iso_xml:
    def __init__(self, xml_url, timeout = 2, read_timeout = 10):
        self.url = xml_url
        self.xml_doc = self._fetch(timeout = timeout, read_timeout = read_timeout)
        self.namespaces = {
            'gmd': 'http://www.isotc211.org/2005/gmd', 
            'gmi': 'http://www.isotc211.org/2005/gmi', 
            'gco': 'http://www.isotc211.org/2005/gco',
            'gml': 'http://www.isotc211.org/2005/gml',
            'th': 'http://www.unidata.ucar.edu/namespaces/thredds/InvCatalog/v1.0',
            'wms': 'http://www.opengis.net/wms',
        }
        
    def _fetch(self, timeout = 2, read_timeout = 10):
        
        return(f_utils.Fetch(self.url).fetch_xml(timeout=timeout, read_timeout=read_timeout))

    def title(self):
        t = self.xml_doc.find('.//gmd:MD_DataIdentification/gmd:citation/gmd:CI_Citation/gmd:title/gco:CharacterString', namespaces = self.namespaces)
        return(t.text if t is not None else 'Unknown')
        
    def bounds(self, geom = True):
        wl = self.xml_doc.find('.//gmd:westBoundLongitude/gco:Decimal', namespaces = self.namespaces)
        el = self.xml_doc.find('.//gmd:eastBoundLongitude/gco:Decimal', namespaces = self.namespaces)
        sl = self.xml_doc.find('.//gmd:southBoundLatitude/gco:Decimal', namespaces = self.namespaces)                            
        nl = self.xml_doc.find('.//gmd:northBoundLatitude/gco:Decimal', namespaces = self.namespaces)                           
        if wl is not None and el is not None and sl is not None and nl is not None:
            region = [float(wl.text), float(el.text), float(sl.text), float(nl.text)]
            if geom: return(regions.Region().from_list([float(wl.text), float(el.text), float(sl.text), float(nl.text)]).export_as_geom())
            else: return(region)
        else: return(None)

    def polygon(self, geom = True):
        opoly = []
        polygon = self.xml_doc.find('.//{*}Polygon', namespaces = self.namespaces)
        if polygon is not None:
            nodes = polygon.findall('.//{*}pos', namespaces = self.namespaces)
            [opoly.append([float(x) for x in node.text.split()]) for node in nodes]
            if geom: return(utils.wkt2geom(regions.create_wkt_polygon(opoly)))
            else: return(opoly)
        else: return(None)
        
    def date(self):
        dt = self.xml_doc.find('.//gmd:date/gco:Date', namespaces = self.namespaces)
        if dt is None:
            dt = self.xml_doc.find('.//gmd:MD_DataIdentification/gmd:citation/gmd:CI_Citation/gmd:date/gmd:CI_Date/gmd:date/gco:Date', namespaces = self.namespaces)
        return(dt.text[:4] if dt is not None else '0000')

    def xml_date(self):
        mddate = self.xml_doc.find('.//gmd:dateStamp/gco:DateTime', namespaces = self.namespaces)
        return(utils.this_date() if mddate is None else mddate.text)
        
    def reference_system(self):
        ref_s = self.xml_doc.findall('.//gmd:MD_ReferenceSystem', namespaces = self.namespaces)
        if ref_s is None or len(ref_s) == 0: return(None, None)
        h_epsg = ref_s[0].find('.//gmd:code/gco:CharacterString', namespaces = self.namespaces)
        if h_epsg is not None: h_epsg = h_epsg.text.split(':')[-1]
        if len(ref_s) > 1:
            v_epsg = ref_s[1].find('.//gmd:code/gco:CharacterString', namespaces = self.namespaces)
            if v_epsg is not None: v_epsg = v_epsg.text.split(':')[-1]
        else: v_epsg = None
            
        return(h_epsg, v_epsg)

    def abstract(self):
        try:
            abstract = self.xml_doc.find('.//gmd:abstract/gco:CharacterString', namespaces = self.namespaces)
            abstract = '' if abstract is None else abstract.text
        except: abstract = ''
        return(abstract)

    def linkages(self):
        linkage = self.xml_doc.find('.//{*}linkage/{*}URL', namespaces = self.namespaces)
        if linkage is not None: linkage = linkage.text
        return(linkage)
    
    def data_links(self):
        dd = {}        
        dfs = self.xml_doc.findall('.//gmd:MD_Format/gmd:name/gco:CharacterString', namespaces = self.namespaces)
        dus = self.xml_doc.findall('.//gmd:onLine/gmd:CI_OnlineResource/gmd:linkage/gmd:URL', namespaces =  self.namespaces)

        if dfs is not None:
            for i,j in enumerate(dfs):
                if j.text in dd.keys():
                    dd[j.text].append(dus[i].text)
                else:
                    dd[j.text] = [dus[i].text]
        return(dd)

class WCS:
    def __init__(self, url):
        self.url = url
        self.namespaces = {
            'wms': 'http://www.opengis.net/wms', 'wcs': 'http://www.opengis.net/wcs/2.0',
            'ows': 'http://www.opengis.net/ows/2.0', 'gml': 'http://www.opengis.net/gml/3.2',
            'gmlcov': 'http://www.opengis.net/gmlcov/1.0'}
        self._get_capabilities()
        self._s_version = self._si()['ServiceTypeVersion'][0]

    def _get_capabilities(self):
        _data = {'request': 'GetCapabilities', 'service': 'WCS'}
        c = f_utils.Fetch(self.url).fetch_req(params=_data)
        cx = lxml.etree.fromstring(c.text.encode('utf-8'))
        self.service_provider = cx.find('.//ows:ServiceProvider', namespaces = self.namespaces)
        self.service_identification = cx.find('.//ows:ServiceIdentification', namespaces = self.namespaces)
        self.operations_metadata = cx.find('.//ows:OperationsMetadata', namespaces = self.namespaces)
        self.service_metadata = cx.find('.//wcs:ServiceMetadata', namespaces = self.namespaces)
        self.contents = cx.find('.//wcs:Contents', namespaces = self.namespaces)

    def _contents(self):
        c = []
        for coverage in self.contents.xpath('//wcs:CoverageSummary', namespaces = self.namespaces):
            c.append(xml2py(coverage))
        return(c)

    def _om(self):
        return(xml2py(self.operations_metadata))

    def _sp(self):
        return(xml2py(self.service_provider))
    
    def _si(self):
        return(xml2py(self.service_identification))
    
    def fix_coverage_id(self, coverage):
        return(':'.join(coverage.split('__')))

    def unfix_coverage_id(self, coverage):
        return('__'.join(coverage.split(':')))

    def _describe_coverage(self, coverage):
        c_d = {}
        valid = False
        c = self._contents()
        for cc in c:
            if coverage == cc['CoverageId']:
                valid = True
                c_d = cc
                break

        om = self._om()
        url = om['DescribeCoverage']['DCP']['HTTP']['Get'][0]
        _data = {'request': 'DescribeCoverage', 'service': 'WCS',
            'version': self._s_version, 'CoverageID': self.unfix_coverage_id(coverage)}
        d = f_utils.Fetch(url).fetch_req(params=_data)
        d_r = lxml.etree.fromstring(d.text.encode('utf-8'))
        cd = d_r.find('.//wcs:CoverageDescription', namespaces = self.namespaces)
        return(xml2py(d_r.find('.//wcs:CoverageDescription', namespaces = self.namespaces)))

    def _get_coverage_region(self, cov_desc):
        uc = [float(x) for x in cov_desc["boundedBy"]["Envelope"]["upperCorner"][0].split()]
        lc = [float(x) for x in cov_desc["boundedBy"]["Envelope"]["lowerCorner"][0].split()]
        return(regions.Region().from_list([lc[1], uc[1], lc[0], uc[0]]))
    
    def _get_coverage_url(self, coverage, region = None):
        dl_coverage = self.fix_coverage_id(coverage)
        cov_desc = self._describe_coverage(coverage)
        fmt = cov_desc["ServiceParameters"]["nativeFormat"][0]        
        hl = [float(x) for x in cov_desc["domainSet"]["RectifiedGrid"]["limits"]["GridEnvelope"]['high'][0].split()]
        uc = [float(x) for x in cov_desc["boundedBy"]["Envelope"]["upperCorner"][0].split()]
        lc = [float(x) for x in cov_desc["boundedBy"]["Envelope"]["lowerCorner"][0].split()]
        ds_region = regions.Region().from_list([lc[1], uc[1], lc[0], uc[0]])
        resx = (uc[1] - lc[1]) / hl[0]
        resy = (uc[0] - lc[0]) / hl[1]
        data = {'request': 'GetCoverage', 'version': '1.0.0', 'service': 'WCS',
                'resx': resx, 'resy': resy, 'crs': 'EPSG:4326', 'format': fmt,
                'coverage': coverage, 'Identifier': coverage}
        if region is not None: data['bbox'] = region.format('bbox')
        enc_data = f_utils.urlencode(data)
        #try:
        #    enc_data = urllib.urlencode(data)
        #except: enc_data = urllib.parse.urlencode(data)
        return('{}{}'.format(self.url, enc_data))
    
    def fetch_coverage(coverage, region = None):
        c_url = self._get_coverage_url(coverage, region)
        return(f_utils.Fetch(c_url, verbose=True).fetch_file('{}_{}.tif'.format(coverage, region.format('fn')), params=data))
    
## =============================================================================
##
## Fetches Remote Elevation Datalist (FRED)
##
## the reference vector location and related functions
##
## =============================================================================
this_dir, this_filename = os.path.split(__file__)
fetchdata = os.path.join(this_dir, 'data')

class FRED:
    def __init__(self, name='FRED', verbose=False, local=False):
        self._verbose = verbose
        self.fetchdata = os.path.join(this_dir, 'data')
        self.driver = ogr.GetDriverByName('GeoJSON')
        self.fetch_v = '{}.geojson'.format(name)
        if local:
            self.FREDloc = self.fetch_v
        elif os.path.exists(self.fetch_v):
            self.FREDloc = self.fetch_v
        elif os.path.exists(os.path.join(self.fetchdata, self.fetch_v)):
            self.FREDloc = os.path.join(self.fetchdata, self.fetch_v)
        else: self.FREDloc = self.fetch_v
        if self._verbose: utils.echo_msg('using {}'.format(self.FREDloc))
        self.ds = None
        self.layer = None
        self.open_p = False
        
        self._fields = ['Name', 'ID', 'Date', 'Agency', 'MetadataLink',
                        'MetadataDate', 'DataLink', 'IndexLink', 'Link',
                        'DataType', 'DataSource', 'Resolution', 'HorizontalDatum',
                        'VerticalDatum', 'LastUpdate', 'Etcetra', 'Info']
        
    def _create_ds(self):
        utils.remove_glob(self.FREDloc)
        self.ds = self.driver.CreateDataSource(self.FREDloc)        
        self.layer = self.ds.CreateLayer('FRED', None, ogr.wkbMultiPolygon)
        ldfn = self.layer.GetLayerDefn()
        self.layer.CreateField(ogr.FieldDefn('Name', ogr.OFTString))
        self.layer.CreateField(ogr.FieldDefn('ID', ogr.OFTString))
        self.layer.CreateField(ogr.FieldDefn('Date', ogr.OFTInteger))
        self.layer.CreateField(ogr.FieldDefn('Agency', ogr.OFTString))
        self.layer.CreateField(ogr.FieldDefn('MetadataLink', ogr.OFTString))
        self.layer.CreateField(ogr.FieldDefn('MetadataDate', ogr.OFTString))
        self.layer.CreateField(ogr.FieldDefn('DataLink', ogr.OFTString))
        self.layer.CreateField(ogr.FieldDefn('IndexLink', ogr.OFTString))
        self.layer.CreateField(ogr.FieldDefn('Link', ogr.OFTString))
        self.layer.CreateField(ogr.FieldDefn('DataType', ogr.OFTString))
        self.layer.CreateField(ogr.FieldDefn('DataSource', ogr.OFTString))
        self.layer.CreateField(ogr.FieldDefn('Resolution', ogr.OFTString))
        self.layer.CreateField(ogr.FieldDefn('HorizontalDatum', ogr.OFTString))
        self.layer.CreateField(ogr.FieldDefn('VerticalDatum', ogr.OFTString))
        self.layer.CreateField(ogr.FieldDefn('LastUpdate', ogr.OFTString))
        self.layer.CreateField(ogr.FieldDefn('Etcetra', ogr.OFTString))
        self.layer.CreateField(ogr.FieldDefn('Info', ogr.OFTString))

    def _add_survey(self, Name = 'fetches', ID = None, Date = None, Agency = None, MetadataLink = None,
                    MetadataDate = None, DataLink = None, IndexLink = None, Link = None, DataType = None,
                    DataSource = None, Resolution = None, HorizontalDatum = None, VerticalDatum = None,
                    Lastupdate = None, Etcetra = None, Info = None, geom = None):
        if not self.open_p: return(None)
        if geom is not None and self.layer is not None:
            self._add_feature([
                {'Name': Name, 'ID': ID, 'Agency': Agency,
                 'Date': Date, 'MetadataLink': MetadataLink,
                 'MetadataDate': MetadataDate, 'DataLink': DataLink,
                 'IndexLink': IndexLink, 'Link': Link, 'DataType': DataType,
                 'DataSource': DataSource, 'Resolution': Resolution,
                 'HorizontalDatum': HorizontalDatum, 'VerticalDatum': VerticalDatum,
                 'LastUpdate': utils.this_date(), 'Etcetra': Etcetra, 'Info': Info},
                geom.ExportToJson()
            ])
            [self.layer.SetFeature(ff) for ff in self.layer]
        else: return(None)
        
    def _open_ds(self, mode = 0):
        if not self.open_p:
            try:
                self.ds = self.driver.Open(self.FREDloc, mode)
            except: self.ds = None

            if self.ds is None or len(self.ds.GetLayer()) == 0:
                self.ds = None
                self._create_ds()
                
            self.layer = self.ds.GetLayer()
            self.open_p = True
            return(0)
        else: return(-1)
        
    def _close_ds(self):
        self.layer = self.ds = None
        self.open_p = False
        return(0)
        
    def _get_fields(self):
        if self.open_p:
            schema = []
            ldefn = self.layer.GetLayerDefn()
            for n in range(ldefn.GetFieldCount()):
                fdefn = ldefn.GetFieldDefn(n)
                schema.append(fdefn.name)
            return(schema)
        else: return(-1)
        
    def _add_feature(self, survey):
        '''add a survey to the reference vector layer'''
        if self.open_p:
            layer_defn = self.layer.GetLayerDefn()
            feat = ogr.Feature(layer_defn)
            geom = ogr.CreateGeometryFromJson(survey[1])
            feat.SetGeometry(geom)
            for field in self._fields:
                try:
                    feat.SetField(field, survey[0][field])
                except: feat.SetField(field, -1)
            self.layer.CreateFeature(feat)
            feat = None
            return(0)
        else: return(-1)

    def _edit_feature(self, feature, survey):
        if self.open_p:
            geom = ogr.CreateGeometryFromJson(survey[1])
            feature.SetGeometry(geom)
            for field in self._fields:
                try:
                    feature.SetField(field, survey[0][field])
                except: feature.SetField(field, -1)
            self.layer.SetFeature(feature)
            return(0)
        else: return(-1)

    def _add_surveys(self, surveys):
        '''update or create a reference vector using a list of surveys'''
        if self.open_p:
            if self.layer is not None:
                for survey in surveys:
                    self._add_survey(**survey)
            return(0)
        else: return(-1)

    def _get_region(self, where = [], layers = []):
        out_regions = []
        if self._verbose: _prog = _progress('gathering regions from {}...'.format(self.FREDloc))
        for i, layer in enumerate(layers):
            if self._verbose: _prog.update_perc((i, len(layers)))
            this_layer = self.layer
            this_layer.SetAttributeFilter("DataSource = '{}'".format(layer))
            [this_layer.SetAttributeFilter('{}'.format(filt)) for filt in where]
            for feat in this_layer:
                geom = feat.GetGeometryRef()
                wkt = geom.ExportToWkt()
                this_region = ogr.CreateGeometryFromWkt(wkt).GetEnvelope()
                if len(out_regions) > 0:
                    out_regions = regions_merge(out_regions, this_region)
                else: out_regions = this_region
        return(out_regions)

    def _attribute_filter(self, where = []):
        if self.open_p:
            attr_str = []
            [attr_str.append(f) for f in where]
            wf = ' AND '.join(attr_str)
            self.layer.SetAttributeFilter(wf)
            return(0)
        else: return(-1)
        
    def _filter(self, region = None, where = [], layers = []):
        '''Search for data in the reference vector file'''
        _results = []
        if region is not None:
            _boundsGeom = region.export_as_geom()#region2geom(region)
        else: _boundsGeom = None

        if self._verbose: _prog = utils.CliProgress('filtering {}...'.format(self.FREDloc))
        if not self.open_p:
            self._open_ds()
            close_p = True
        else: close_p = False

        for i, layer in enumerate(layers):
            if self._verbose: _prog.update_perc((i, len(layers)))
            #this_layer = self.layer
            where.append("DataSource = '{}'".format(layer))
            if self._verbose: utils.echo_msg('FRED filter: {}'.format(where))
            self._attribute_filter(where = where)
            for feat in self.layer:
                if _boundsGeom is not None:
                    geom = feat.GetGeometryRef()
                    if geom is not None:
                        if _boundsGeom.Intersects(geom):
                            _results.append({})
                            f_j = json.loads(feat.ExportToJson())
                            for key in f_j['properties'].keys():
                                _results[-1][key] = feat.GetField(key)
                else:
                    _results.append({})
                    f_j = json.loads(feat.ExportToJson())
                    for key in f_j['properties'].keys():
                        _results[-1][key] = feat.GetField(key)
                        
            #this_layer = None
        if close_p: self._close_ds()
        if self._verbose: _prog.end(0, 'filtered \033[1m{}\033[m data records from FRED'.format(len(_results)))
        return(_results)

## ==============================================
## lambdas for the FRED using the module object `mod`
## ==============================================
_filter_FRED = lambda mod: mod.FRED._filter(mod.region, mod.where, [mod.name])
_update_FRED = lambda mod, s: mod.FRED._add_surveys(s)
_filter_FRED_index = lambda mod: [utils.echo_msg(json.dumps(f, indent = 2)) for f in _filter_FRED(mod)]

### End
