#!/usr/bin/env python

#import urllib2
import os
import sys
import glob
from xml.dom import minidom

di_version = '0.0.1'

dem_xml=sys.argv[1]

def dem_info(dem_xml):
    f = open(dem_xml, 'r')
    xmldoc = minidom.parseString(f.read())

    # Title
    try:
        xml_ident = xmldoc.getElementsByTagName("citation")
        for t in xml_ident:
            xml_citation = t.getElementsByTagName("citeinfo")
            title = xml_citation[0].getElementsByTagName("title")[0].firstChild.nodeValue.strip().replace(",", " ").replace("/", "_").replace("(", "_").replace(")", "_")
            date = xml_citation[0].getElementsByTagName("pubdate")[0].firstChild.nodeValue.strip()
    except:
        title = "Unknown"
        date = "Unknown"

    xml_spdom = xmldoc.getElementsByTagName("spdom")
    xml_bounding = xml_spdom[0].getElementsByTagName("bounding")
    wb = xml_bounding[0].getElementsByTagName("westbc")[0].firstChild.nodeValue.strip()
    eb = xml_bounding[0].getElementsByTagName("eastbc")[0].firstChild.nodeValue.strip()
    nb = xml_bounding[0].getElementsByTagName("northbc")[0].firstChild.nodeValue.strip()
    sb = xml_bounding[0].getElementsByTagName("southbc")[0].firstChild.nodeValue.strip()

    try:
        xml_spref = xmldoc.getElementsByTagName("spref")
        xml_horiz = xml_spref[0].getElementsByTagName("horizsys")[0].getElementsByTagName("geograph")
        latres = xml_horiz[0].getElementsByTagName("latres")[0].firstChild.nodeValue.strip()
        lonres = xml_horiz[0].getElementsByTagName("longres")[0].firstChild.nodeValue.strip()
    except:
        latres = "NaN"
        lonres = "NaN"
    
    return title,date,wb,eb,sb,nb,latres

dem_infos = dem_info(dem_xml)

#print dem_infos[0]
xy_region_string = "xy-region.scm -R" + dem_infos[2] + "/" + dem_infos[3] + "/" + dem_infos[4] + "/" + dem_infos[5] + " -A title:" + dem_infos[0].replace(" ", "_") + ":string" + ",date:" + dem_infos[1] + ":string" + ",resolution:" + dem_infos[6] + ":double -a >> dems.gmt"

os.system(xy_region_string)
