#!/usr/bin/env python
#
#  Description:  Return the number of rows and columns in a given gdal-supported grid file.
#

import sys
from osgeo import gdal
#from xml.dom import minidom
import xml.etree.ElementTree as ET
import datetime 

gm_version = 0.1

def Usage():
    print('Usage: gdal_metadata.py grdfile')
    print('')
    print('gdal_metadata v.%s' %(gm_version))

def grdinfo(ingrd):
    # Process the grid file
    ds = gdal.Open(ingrd)
    comp_geot = ds.GetGeoTransform()
    
    xrow = ds.RasterXSize
    yrow = ds.RasterYSize

    cellsize = [float(comp_geot[1]), float(comp_geot[5])]
    xll = float(comp_geot[0])
    yur = float(comp_geot[3])
    xur = xll + (cellsize[0] * xrow)
    yll = yur - (cellsize[1] * yrow)

    return cellsize,xll,xur,yll,yur,xrow,yrow

if __name__ == "__main__":

    # Process the command-line
    infile = None
    verbose = False
    mtitle = ""
    mgrdnum = 0
    vdatum = "MHW"
    sxml = None
    dxml = None

    i = 1
    argv = sys.argv
    while i < len(argv):
        arg = argv[i]

        if arg == '-help' or arg == '-h':
            Usage()
            sys.exit(0)

        elif arg == '-title':
            mtitle = argv[i+1]
            i = i + 1

        elif arg == '-groovy_number':
            grdnum = argv[i+1]
            i = i + 1

        elif arg == '-vdatum':
            vdatum = argv[i+1]
            i = i + 1

        elif arg == '-sxml':
            sxml = argv[i+1]
            i = i + 1

        elif arg == '-dxml':
            dxml = argv[i+1]
            i = i + 1

        elif arg == '-verbose':
            verbose = True

        elif infile is None:
            infile = arg

        elif arg[0] == '-':
            Usage()

        #else:
         #   Usage()

        i = i + 1    
    print(infile)
    if infile is None:
        print(Usage())
        sys.exit(0)

    try:
        gdal_infos = grdinfo(infile)
    except:
        print(Usage())
        sys.exit(0)

    if sxml is None or dxml is None:
        print(Usage())
        print("you must supply an input and output xml file")
        sys.exit(0)

    ## Metadata Namespace
    ET._namespace_map["http://www.isotc211.org/2005/gmi"] = "gmi"
    ET._namespace_map["http://www.isotc211.org/2005/gco"] = "gco"
    ET._namespace_map["http://www.isotc211.org/2005/gmd"] = "gmd"
    ET._namespace_map["http://www.isotc211.org/2005/gmx"] = "gmx"
    ET._namespace_map["http://www.isotc211.org/2005/gsr"] = "gsr"
    ET._namespace_map["http://www.isotc211.org/2005/gss"] = "gss"
    ET._namespace_map["http://www.isotc211.org/2005/gts"] = "gts"
    ET._namespace_map["http://www.opengis.net/gml/3.2"] = "gml"
    ET._namespace_map["http://www.isotc211.org/2005/gmi"] = "xsi"
    ET._namespace_map["http://www.isotc211.org/2005/srv"] = "srv"
    ET._namespace_map["http://www.w3.org/1999/xlink"] = "xlink"
   

    gmi="{http://www.isotc211.org/2005/gmi}" 
    gco="{http://www.isotc211.org/2005/gco}" 
    gmd="{http://www.isotc211.org/2005/gmd}" 
    gmx="{http://www.isotc211.org/2005/gmx}" 
    gsr="{http://www.isotc211.org/2005/gsr}" 
    gss="{http://www.isotc211.org/2005/gss}" 
    gts="{http://www.isotc211.org/2005/gts}" 
    gml="{http://www.opengis.net/gml/3.2}" 
    xsi="{http://www.isotc211.org/2005/gmi http://www.ngdc.noaa.gov/metadata/published/xsd/schema.xsd}"
    srv="{http://www.isotc211.org/2005/srv}"
    xlink="{http://www.w3.org/1999/xlink}"

    print(sxml)
    tree = ET.parse(sxml)
    doc = tree.getroot()

    # Date Stamp
    ds_elem = doc.findall("{0}dateStamp".format(gmd))
    print(ds_elem[0].find("{0}Date".format(gco)).text)
    today = datetime.date.today()
    ds_elem[0].find("{0}Date".format(gco)).text = str(today)
    
    # Spatial Representation
    sp_elem = doc.findall("{0}spatialRepresentationInfo/{0}MD_GridSpatialRepresentation/{0}axisDimensionProperties".format(gmd))

    ## Resolution
    for i in sp_elem:
        se = i.findall("{0}MD_Dimension/{0}resolution".format(gmd))
        for j in se:
            print(j.find("{0}Measure".format(gco)).text)
            j.find("{0}Measure".format(gco)).text = str(gdal_infos[0][0])

    ## Row/Col
    for i in sp_elem:
        se = i.findall("{0}MD_Dimension".format(gmd))
        for j in se:
            if j.find("{0}dimensionName/{0}MD_DimensionNameTypeCode".format(gmd)).text == "row":
                j.find("{0}dimensionSize".format(gmd)).find("{0}Integer".format(gco)).text = str(gdal_infos[5])
            if j.find("{0}dimensionName/{0}MD_DimensionNameTypeCode".format(gmd)).text == "column":
                j.find("{0}dimensionSize".format(gmd)).find("{0}Integer".format(gco)).text = str(gdal_infos[6])


    # Data Identification Info
    id_elem = doc.findall("{0}identificationInfo/{0}MD_DataIdentification/".format(gmd))
    
    ## Citation
    ### Title:
    id_elem[0].find("{0}citation/{0}CI_Citation/{0}title/{1}CharacterString".format(gmd,gco)).text = mtitle
    ### Date:
    id_elem[0].find("{0}citation/{0}CI_Citation/{0}date/{0}CI_Date/{0}date/{1}Date".format(gmd,gco)).text = str(today)
    ### Responsible Party
    rp = id_elem[0].findall("{0}citation/{0}CI_Citation/{0}citedResponsibleParty".format(gmd))
    #### online resources
    for i in rp:
        if i.find("{0}CI_ResponsibleParty/{0}positionName/{1}CharacterString".format(gmd,gco)) is not None:
            j = i.find("{0}CI_ResponsibleParty/{0}contactInfo/{0}CI_Contact/{0}onlineResource/{0}CI_OnlineResource".format(gmd))
            # url
            j.find("{0}linkage/{0}URL".format(gmd)).text = "http://www.ngdc.noaa.gov/dem/squareCellGrid/getReport/0".format(mgrdnum)
            # title
            j.find("{0}name/{1}CharacterString".format(gmd,gco)).text = mtitle
    ### series
    id_elem[0].find("{0}citation/{0}CI_Citation/{0}series/{0}CI_Series/{0}issueIdentification/{1}CharacterString".format(gmd,gco)).text = mtitle

    ## Graphic Overview
    ### URL
    id_elem[0].find("{0}graphicOverview/{0}MD_BrowseGraphic/{0}fileName/{1}CharacterString".format(gmd,gco)).text = "http://www.ngdc.noaa.gov/dem/squareCellGrid/getGraphic/{0}".format(mgrdnum)
    
    ##
    id_elem[0].find("{0}aggregationInfo/{0}MD_AggregateInformation/{0}aggregateDataSetName/{0}CI_Citation/{0}title/{1}CharacterString".format(gmd,gco)).text = mtitle
    id_elem[0].find("{0}aggregationInfo/{0}MD_AggregateInformation/{0}aggregateDataSetName/{0}CI_Citation/{0}date/{0}CI_Date/{0}date/{1}Date".format(gmd,gco)).text = str(today)
    

    ## Report Info

    ## Extent
    gext = id_elem[0].findall("{0}extent/{0}EX_Extent/{0}geographicElement/{0}EX_GeographicBoundingBox".format(gmd))
    gext[0].find("{0}westBoundLongitude/{1}Decimal".format(gmd,gco)).text = str(gdal_infos[1])
    gext[0].find("{0}eastBoundLongitude/{1}Decimal".format(gmd,gco)).text = str(gdal_infos[2])
    gext[0].find("{0}southBoundLatitude/{1}Decimal".format(gmd,gco)).text = str(gdal_infos[3])
    gext[0].find("{0}northBoundLatitude/{1}Decimal".format(gmd,gco)).text = str(gdal_infos[4])

    # Write tree to file
    tree.write(dxml)

# END
