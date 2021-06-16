### lasfun.py
##
## Copyright (c) 2010 - 2021 CIRES Coastal DEM Team
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
import sys
import json
import struct
import time
import re
import numpy as np
from scipy import spatial
import cudem
from cudem import utils
from cudem import regions
from cudem import xyzs

## ==============================================
## las-file processing (datalists fmt:300)
## ==============================================
def las_inf(src_las):

    pts = []
    lasi = {}
    lasi['name'] = src_las
    lasi['numpts'] = 0
    lasi['minmax'] = [0, 0, 0, 0, 0, 0]
    utils.echo_msg('generating inf file for {}'.format(src_las))
    
    for i, l in enumerate(las_yield_entry([src_las])):
        if i == 0:
            lasi['minmax'] = [l[0], l[0], l[1], l[1], l[2], l[2]]
        else:
            try:
                if l[0] < lasi['minmax'][0]: lasi['minmax'][0] = l[0]
                elif l[0] > lasi['minmax'][1]: lasi['minmax'][1] = l[0]
                if l[1] < lasi['minmax'][2]: lasi['minmax'][2] = l[1]
                elif l[1] > lasi['minmax'][3]: lasi['minmax'][3] = l[1]
                if l[2] < lasi['minmax'][4]: lasi['minmax'][4] = l[2]
                elif l[2] > lasi['minmax'][5]: lasi['minmax'][5] = l[2]
            except: pass
        pts.append(l)
        lasi['numpts'] = i

    try:
        out_hull = [pts[i] for i in spatial.ConvexHull(pts, qhull_options='Qt').vertices]
        out_hull.append(out_hull[0])
        lasi['wkt'] = utils.create_wkt_polygon(out_hull, xpos = 0, ypos = 1)
    except: lasi['wkt'] = regions.region2wkt(lasi['minmax'])

    if lasi['numpts'] > 0:
        with open('{}.inf'.format(src_las), 'w') as inf:
            inf.write(json.dumps(lasi))
        
    return(lasi)

def las_inf2(src_las):
    '''scan an xyz file and find it's min/max values and
    write an associated inf file for the src_xyz file.

    returns region [xmin, xmax, ymin, ymax, zmin, zmax] of the src_xyz file.'''

    minmax = []
    out, status = utils.run_cmd('lasinfo -nc -nv -stdout -i {}'.format(src_las), verbose = False)
    for i in out.split('\n'):
        if 'min x y z' in i:
            xyz_min = [float(y) for y in [x.strip() for x in i.split(':')][1].split()]
        if 'max x y z' in i:
            xyz_max = [float(y) for y in [x.strip() for x in i.split(':')][1].split()]

    minmax = [xyz_min[0], xyz_max[0], xyz_min[1], xyz_max[1], xyz_min[2], xyz_max[2]]

    with open('{}.inf'.format(src_las), 'w') as inf:
        utils.echo_msg('generating inf file for {}'.format(src_las))
        inf.write('{}\n'.format(' '.join([str(x) for x in minmax])))
    
    return(minmax)
    
def las_inf_entry(entry):
    '''find the region of the xyz datalist entry
    
    returns the region [xmin, xmax, ymin, ymax, zmin, zmax] of the xyz entry'''

    return(las_inf(entry[0]))
    
def las_yield_entry(entry, region = None, verbose = False, z_region = None):
    '''yield the xyz data from the xyz datalist entry

    yields [x, y, z, <w, ...>]'''
    
    ln = 0
    if z_region is not None:
        min_z = None if z_region[0] is None else z_region[0]
        max_z = None if z_region[1] is None else z_region[1]
    else: min_z = max_z = None
    for line in utils.yield_cmd('las2txt -parse xyz -stdout -keep_class 2 29 -i {} {} {} {}\
    '.format(entry[0], '' if region is None else '-keep_xy {}'.format(regions.region_format(region, 'te')),\
             '' if min_z is None else '-drop_z_below {}'.format(min_z),\
             '' if max_z is None else '-drop_z_above {}'.format(max_z)), data_fun = None, verbose = False):
        ln += 1
        xyz = [float(x) for x in line.strip().split()]
        yield(xyz + [entry[2]] if entry[2] is not None else xyz)

    if verbose: utils.echo_msg('read {} data points from {}'.format(ln, entry[0]))

def las_dump_entry(entry, dst_port = sys.stdout, region = None, verbose = False, z_region = None):
    '''dump the las data from the las datalist entry to dst_port'''
    
    for xyz in las_yield_entry(entry, region, verbose, z_region):
        xyzfun.xyz_line(xyz, dst_port, True, None)        

def las_file_p(inlas):
    '''Check quickly whether file is an LAS
    file.  The first four(4) bytes of an
    LAS file should have the characters "LASF"'''

    return(open(inlas, 'rb').read(4) == b'LASF')

## ==============================================
## las-file header
## ==============================================
class LasHeader:
    
    def __init__(self, inlas=None):
        self.inlas = inlas
        try:
            if os.path.exists(self.inlas):
                self.nostream = False
            else: self.nostream = True
        except: self.nostream = True

    header_struct = (
        ('filesig',4,'c',4),
        ('filesourceid',2,'H',1),
        ('reserved',2,'H',1),
        ('guid1',4,'L',1),
        ('guid2',2,'H',1),
        ('guid3',2,'H',1),
        ('guid4',8,'B',8),
        ('vermajor',1,'B',1),
        ('verminor',1,'B',1),
        ('sysid',32,'c',32),
        ('gensoftware',32,'c',32),
        ('fileday',2,'H',1),
        ('fileyear',2,'H',1),
        ('headersize',2,'H',1),
        ('offset',4,'L',1),
        ('numvlrecords',4,'L',1),
        ('pointformat',1,'B',1),
        ('pointreclen',2,'H',1),
        ('numptrecords',4,'L',1),
        ('numptbyreturn',20,'L',5),
        ('xscale',8,'d',1),
        ('yscale',8,'d',1),
        ('zscale',8,'d',1),
        ('xoffset',8,'d',1),
        ('yoffset',8,'d',1),
        ('zoffset',8,'d',1),
        ('xmax',8,'d',1),
        ('xmin',8,'d',1),
        ('ymax',8,'d',1),
        ('ymin',8,'d',1),
        ('zmax',8,'d',1),
        ('zmin',8,'d',1),
    )

    vlheader_struct = (
        ('reserved',2,'H',1),
        ('userid',16,'c',16),
        ('recordid',2,'H',1),
        ('reclen',2,'H',1),
        ('desc',32,'c',32),
    )
    
    def parse_header(self):
        """Parse the LAS Header information into a
        python dictionary.  See the header_struct
        defenition above for the names to use to
        call certain values."""

        if self.nostream:
            return("You cannot parse what you do not have, try createHeader first...")

        ds = open(self.inlas, 'rb')
        lheader = {}
        for i in self.header_struct:
            try:
                if i[2] == 'c':
                    value = ds.read(i[1])
                elif i[3] > 1:
                    value = struct.unpack("=" + str(i[3]) + i[2], ds.read(i[1]))
                else:
                    value = struct.unpack("=" + i[2], ds.read(i[1]))[0]
            except:
                value = ds.read(i[1])
            lheader[i[0]] = value
        ds.close()
        return(lheader)

    def dump_header(self):
        """Dump the LAS file public header to the screen"""
        
        lash = self.parse_header()
        for i in self.header_struct:
            print(i[0] + ":" + "\t", lash[i[0]])

    def get_seek(self):
        """Retreive and return the proper seek
        position, which is just the value of the
        size of the header in question. - initial
        seek position, that is."""
		
        return(self.parse_header()['headersize'])

    def get_VLR_num(self):
        """Return the number of variable length
        records present in LAS file."""
		
        return(self.parse_header()['numvlrecords'])

    def parse_VLR(self, seeknum):
        """Parse a variable length header from an
        LAS file.  The existence of any variable
        length headers should be mentioned in the
        main header of the LAS file."""

        if self.nostream:
            return("You cannot parse what you do not have, try createHeader first...")
        
        ds = open(self.inlas,'rb')
        ds.seek(seeknum)
        vlhead = {}
        for i in self.vlheader_struct:
            if i[2] == 'c':
                value = ds.read(i[1])
            elif i[2] == 'H':
                value = struct.unpack("=" + i[2] , ds.read(i[1]) )[0]
            vlhead[i[0]] = value
        ds.close()
        return(vlhead)
    
    def get_VLR_recs(self):
        """Parse variable length record from an
        LAS file, based on the information in
        the LAS variable length header."""
		
        seeknum = self.get_seek()
        vlnum = self.get_VlR_num()
        vlrecs = []
        for i in range(0, vlnum):
            vlr = self.parse_VLR(seeknum)
            rid = vlr['recordid']
            if len(str(rid)) < 4:
                rid = vlr['reserved']
            seeknum = seeknum + 54 + vlr['reclen']
            recs = (rid, vlr['reclen'], seeknum)
            vlrecs.append(recs)
        return(vlrecs)

    def read_VLR_rec(self, seeknum, recid, reclen):
        ds = open(self.inlas, 'rb')
        # Las Spec Records
        if int(recid) == 0:
            '''classification lookup'''
        elif int(recid) == 2:
            '''histogram'''
        elif int(recid) == 3:
            ds.seek(seeknum+54)
            return ds.read(reclen)
        # Projection Info
        elif int(recid) == 34735:
            '''do something'''
        elif int(recid) == 34736:
            '''do something'''
        elif int(recid) == 34737:
            ds.seek(seeknum+54)
            return(ds.read(reclen))
            
    def get_scale(self):
        """Retreive and return the x, y and z scale
        information from the LAS header."""

        h = self.parse_header()
        return(h['xscale'], h['yscale'], h['zscale'])

    def get_offsets(self):
        """Return the x,y and z offsets from the header"""
        
        h = self.parse_header()
        return(h['xoffset'], h['yoffset'], h['zoffset'])

    def get_return_grp(self, instream):
        """Return the return group information.
        This is used in header creation."""

        firstret = 0
        secondret = 0
        thirdret = 0
        fourthret = 0
        fifthret = 0
        for row in instream:
            return_grp = row['r']
            return_num = return_grp & 7
            if return_num == 1 or return_num == 0:
                firstret += 1
            elif return_num == 2:
                secondret +=1
            elif return_num == 3:
                thirdret += 1
            elif return_num == 4:
                fourthret += 1
            elif return_num == 5:
                fifthret += 1
				
        return(firstret, secondret, thirdret, fourthret, fifthret)

    def get_pnt_stats(self, returns=False):
        """Return some statistics about the point-
        cloud, for header creation."""

        if self.nostream:
            instream = self.inlas
        else: instream = las_file(self.inlas).readAsArray()
        maxx = instream['x'].max()
        minx = instream['x'].min()
        maxy = instream['y'].max()
        miny = instream['y'].min()
        maxz = instream['z'].max()
        minz = instream['z'].min()
        tot = len(instream)
        if returns:
            firstret,secondret,thirdret,fourthret,fifthret = self.get_return_grp(instream)
        else:
            firstret,secondret,thirdret,fourthret,fifthret = 0,0,0,0,0
        return maxx,maxy,maxz,minx,miny,minz,tot,\
               firstret,secondret,thirdret,fourthret,fifthret

    def get_pnt_stats_as_list(self, returns=False, reclen=10):
        """Return some statistics about the point-
        cloud, for header creation."""

        if self.nostream:
            instream = self.inlas
        else:
            instream = las_file(inlas).readAsList()
        maxx = 0
        maxy = 0
        minx = 0
        miny = 0
        maxz = 0
        minz = 0
        tot = 0

        for i in range(0, len(instream)):
            for j in range(0, len(instream[i]), reclen):
                if instream[i][j] > maxx:
                    maxx = instream[i][j]
                elif instream[i][j] < minx:
                    minx = instream[i][j]
                if instream[i][j+1] > maxy:
                    maxy = instream[i][j+1]
                elif instream[i][j+1] < miny:
                    miny = instream[i][j+1]
                if instream[i][j+2] > maxz:
                    maxz = instream[i][j+2]
                elif instream[i][j+2] < minz:
                    minz = instream[i][j+2]
                tot += 1

        if returns:
            firstret,secondret,thirdret,fourthret,fifthret = self.getReturnGrp(instream)
        else:
            firstret,secondret,thirdret,fourthret,fifthret = 0,0,0,0,0
        return(maxx,maxy,maxz,minx,miny,minz,tot,\
               firstret,secondret,thirdret,fourthret,fifthret)

    def create_header(self, returns=False, scale=[0.0001,0.0001,0.0001], offsets=[0,0,0], reclen=10):
        """Create the las header, which will output
        as version 1.1 whether the input LAS file
        is 1.0 or 1.1, other versions are not yet
        supported (v1.2 & v1.3dev).  If returns equals
        False, then the return numbers wont be calculated
        for the header.
        """

        lheader = {}
        fday = time.localtime()[7]
        fyear = time.localtime()[0]
        lstats = self.get_pnt_stats_as_list(returns)
        h_size = 227
        h_offset = 229

        if reclen == 10:
            ptformat = 1
            ptlen = 28
        else:
            ptformat = 0
            ptlen = 28

        values = ["LASF", 0, 0, 0, 0, 0, (0, 0, 0, 0, 0, 0, 0, 0), 1, 1, \
                  "Python " + str(sys.version[0:5]), "GEOMODS " + _version, \
                  fday, fyear, h_size, h_offset, 0, ptformat, ptlen, lstats[6], \
                  (lstats[7], lstats[8], lstats[9], lstats[10], lstats[11]), \
                  scale[0], scale[1], scale[2], offsets[0], offsets[1], offsets[2], \
                  lstats[0] * scale[0], lstats[3] * scale[0], lstats[1] * scale[1], \
                  lstats[4] * scale[1], lstats[2] * scale[2], lstats[5] * scale[2]]

        j = 0
        for i in values:
            lheader[self.header_struct[j][0]] = i
            j += 1

        return(lheader)

    def create_np_header(self, returns=False, scale=[0.0001,0.0001,0.0001], offsets=[0,0,0]):
        """Create an LAS File header as a numpy array"""

        dt = "a4, u2, u2, u4, u2, u2, (8,)u1, u1, u1, a32, a32, u2, \
        u2, u2, u4, u4, u1, u2, u4, (5,)u4, f8, f8, f8, f8, f8, f8, \
        f8, f8, f8, f8, f8, f8"
		
        fday = time.localtime()[7]
        fyear = time.localtime()[0]
        lstats = self.getPtStats()
        h_size = 227
        h_offset = 229
        ptformat = 1
        ptlen = 28
        
        if self.nostream is False:
            offsets = self.get_offsets()
            scale = self.get_scale()
        

        values = ("LASF", 0, 0, 0, 0, 0, (0, 0, 0, 0, 0, 0, 0, 0), 1, 1, \
                  "Python " + str(sys.version[0:5]), "ML-LAS " + _version, \
                  fday, fyear, h_size, h_offset, 0, ptformat, ptlen, lstats[6], \
                  (lstats[7], lstats[8], lstats[9], lstats[10], lstats[11]), \
                  scale[0], scale[1], scale[2], offsets[0], offsets[1], offsets[2], \
                  lstats[0] * scale[0], lstats[3] * scale[0], lstats[1] * scale[1], \
                  lstats[4] * scale[1], lstats[2] * scale[2], lstats[5] * scale[2])
        
        return(np.array(values, dt))

    def create_VLR_header(self, type="proj", proj_string="NIL"):
        """Create a variable-length header for
        an LAS file in version 1.1 of the standard."""
		
        vheader = {}
        j = 0
        if type == "proj":
            rec_id = 34735
            user_id = "LASF_Projection"
            reclen = 0
            desc = proj_string
        elif type == "class":
            rec_id = 0
            user_id = "LASF_Spec"
            desc = "LASF Classifications"
            reclen = 0
            values = [0, user_id, rec_id, reclen, desc]
        for i in values:
            vheader[vlheader_struct[j][0]] = i
            j += 1
        return(vheader)

    def create_VLR_record(self):
        """do something"""

        pass

class LasFile:
    def __init__(self, inlas):
        if not las_file_p(inlas):
            sys.exit("This does not appear to be a proper LAS file. Please check your source file.")
        
        self.inlas = inlas
        self.hc = LasHeader(inlas)
        self.h = self.hc.parse_header()
        
        self.bsize = 1000
        self.bnames = ('x','y','z','i','r','c','s','u','p','g')
        self.btypes = ('i4','i4','i4','u2','u1','u1','u1','u1','u2','f8')

        if self.h['pointformat'] == 1:
            self.reclen = 10
        elif self.h['pointformat'] == 0:
            self.reclen = 9
            self.bnames = self.bnames[:-1]
            self.btypes = self.btypes[:-1]
        else:
            sys.exit("Support for point-format %s is not yet supported...") % (str(h['pointformat']))

    # =============================================================================
    def valid_point(self, pt_rec, use_recs="xyzi", clf=None):
        """Validate a point record in terms of user
        preference, for use in dumpPoints function.
        use_recs is a string containing any of the
        following letters: "xyziecaupg". clf is the
        classifcation field to use, if set, otherwise
        will validate any class records"""

        point = {}
        outpoint = []

        scale = self.hc.get_scale()
        offsets = self.hc.get_offsets()
        pt_form = tuple(re.findall(r'(\D)', use_recs))

        for form in pt_form:
            if form == 'x':
                outpoint.append(pt_rec[0] * scale[0] + offsets[0])
            elif form == 'y':
                outpoint.append(pt_rec[1] * scale[1] + offsets[1])
            elif form == 'z':
                outpoint.append(pt_rec[2] * scale[2] + offsets[2])
            else:
                outpoint.append(pt_rec[list(self.bnames).index(form)])

        return(outpoint)

    def dump_points(self, use_recs="xyz", clf=None, delim=","):

        if clf is not None:
            clfs = clf.split(",")
        
        las = open(self.inlas, 'rb')
        las.seek(self.h['offset'])

        BUFF = self.h['pointreclen']

        BUFFERSIZE = BUFF * self.bsize
        
        if BUFF == 28:
            UNPACKFORMAT = 'lllHbBBBHd'
        if BUFF == 20:
            UNPACKFORMAT = 'lllHbBBBH'

        unpackstring = "=" + UNPACKFORMAT
        
        #print >> sys.stderr, BUFF
        #while True:
        for i in range(self.h['numptrecords']):
            thisP = struct.unpack(unpackstring, las.read(BUFF))

            if clf is not None and str(thisP[5]) in clfs:
                print(str(delim).join(map(str, self.valid_point(thisP, use_recs, clf))))
            elif clf is None:
                print(str(delim).join(map(str, self.valid_point(thisP, use_recs, clf))))

    def readAsArray(self, israw=False):
        """Read a binary LAS-files point records into a
        NumPy Rec-Array.
        """

        lasf = open(self.inlas, 'rb')
        lasf.seek(self.h['offset'])
        scale = self.hc.get_scale()
        offsets = self.hc.get_offsets()

        u = np.fromfile(lasf, {'names': self.bnames, 'formats': self.btypes})

        if not israw:
            self.rawbtypes = list(self.btypes)
            self.rawbtypes[0] = 'f8'
            self.rawbtypes[1] = 'f8'
            self.rawbtypes[2] = 'f8'
            self.rawbtypes = tuple(self.rawbtypes)
            uf = np.zeros(len(u), {'names': self.bnames, 'formats': self.rawbtypes})
            for i in self.bnames:
                if i == 'x':
                    uf[i] = u[i] * scale[0] + offsets[0]
                elif i == 'y':
                    uf[i] = u[i] * scale[1] + offsets[1]
                elif i == 'z':
                    uf[i] = u[i] * scale[2] + offsets[2]
                elif i == 'r':
                    uf[i] = u[i] & 7
                else:
                    uf[i] = u[i]
            u = uf
            
        lasf.close()
        return u     

    def readAsList(self):
        """Load an LAS file into a python list.  This is
        the old default function, now depreciated, though
        it`s abilities are still in working order.
        """
		
        las = open(self.inlas, 'rb')
        u = []
        las.seek(self.h['offset'])
        if self.h['pointformat'] == 1:
            BUFFERSIZE = 28 * self.bsize
            UNPACKFORMAT = 'lllHbBBBHd'
        elif self.h['pointformat'] == 0:
            BUFFERSIZE = 20 * self.bsize
            UNPACKFORMAT = 'lllHbBBBH'
        pointlen = self.h['pointreclen']
        while True:
            try:
                data = las.read(BUFFERSIZE)
            except:
                break
            if not data:
                break
            unpackstring = "=" + UNPACKFORMAT * (len(data) / pointlen)
            u.append(struct.unpack(unpackstring, data))
        las.close()
        return(u)

    def readAsSubSet(self, col="c", clf=2):
        """'Return a subset of the las point array, based
        on the field (col) and the field-value (clf).
        """

        clfs = clf.split(",")

        lasp = self.readAsArray()

        if len(clfs) > 1:
            classindex = np.empty(0, 'i4')
            for clf in clfs:
                classindex_tmp = np.where(lasp[col] == int(clf))
                classindex = np.append(classindex, classindex_tmp)
            classindex.sort()
        else:
            classindex = np.where(lasp[col] == int(clf))[0]
        
        newArray = np.zeros(len(classindex), {'names': self.bnames, 'formats': self.rawbtypes})
        for num,item in enumerate(classindex):
            newArray[num] = lasp[item]
        return(newArray)

    def write_file(self, inarray, outfile):
        '''Write an LAS file to disk'''
        
        if isnp:
            LasWrite(inarray, outfile).write_las_file()
        else:
            LasWrite(inarray, outfile).write_las_file_from_list()

class LasWrite:
    def __init__(self, inarray, outfile):

        self.inarray = inarray
        try:
            self.inarray.shape
            self.isnp = True
        except:
            self.isnp = False
        if self.isnp:
            self.inarrayRowlen = len(self.inarray[0])
        else:
            self.inarrayRowlen = 10
        self.outlas = open(outfile, 'wb')
        self.lh = LasHeader(self.inarray)

    # =============================================================================
    def write_header(self, xyzscale=[0.01,0.01,0.01], xyzoffsets=[0,0,0], reclen=10):
        h = las_header(self.inarray).create_header(False, xyzscale, xyzoffsets, reclen)
        values = []
        for i in self.lh.header_struct:
            try:
                charfinder = h[i[0]][0] + 1
                for j in h[i[0]]:
                    values.append(j)
            except:
                values.append(h[i[0]])
        FORMATSTRING = '= 4s H H L H H 8B B B 32s \
        32s H H H L L B H L 5L d d d d d d d d d d d d'
        
        s = struct.Struct(FORMATSTRING)
        packed_data = s.pack(*values)
        self.outlas.write(packed_data)
        self.outlas.write(struct.pack('2x'))

    def write_point_records(self, nstream):
        PACKFORMAT = 'lllHbBBBHd'
        if self.isnp:
            nstream = self.inarray.tolist()
        for i in nstream:
            FORMATSTRING =  "=" + (PACKFORMAT * (len(i) / self.inarrayRowlen))
            s = struct.Struct(FORMATSTRING)
            packed_records = s.pack(*i)
            self.outlas.write(packed_records)

    def write_las_file(self, xyzscale=[0.01,0.01,0.01], xyzoffsets=[0,0,0]):
        '''Write an LAS binary file based on the
        NumPy Rec-Array of the LAS point cloud '''

        harray = self.lh.create_np_header(False,xyzscale,xyzoffsets)
        harray.tofile(self.outlas)
        self.outlas.write(struct.pack('2x'))
        self.inarray.tofile(self.outlas)
        self.outlas.close()
		
    def write_las_file_from_list(self):
        
        self.write_header()
        self.write_point_records(self.inarray)
        self.outlas.close()

    def lhead2nphead(self, lhead):
        '''Convert a header dictionary to a numpy array. '''

        dt = "a4, u2, u2, u4, u2, u2, (8,)u1, u1, u1, a32, a32, u2, \
        u2, u2, u4, u4, u1, u2, u4, (5,)u4, f8, f8, f8, f8, f8, f8, \
        f8, f8, f8, f8, f8, f8"
        nh = []
        for i in header_struct:
            nh.append(lhead[i[0]])
        return(np.array(tuple(nh), dt))

## ==============================================
##
## lasfun cli
##
## ==============================================
lasfun_usage = '''{cmd} ({version}): lasfun; Process and generate las files

usage: {cmd} [ -hqPRW [ args ] ]...

Options:

  --quiet\t\tLower the verbosity to a quiet

  --help\t\tPrint the usage text
  --version\t\tPrint the version information

Examples:
  % {cmd}

CIRES DEM home page: <http://ciresgroups.colorado.edu/coastalDEM>\
'''.format(cmd =  os.path.basename(sys.argv[0]), 
           version = cudem.__version__)

def lasfun_cli(argv = sys.argv):
    """run lasfun from command-line

    See `lasfun_usage` for full cli options.
    """

    want_verbose = True
    inlasen = []
    
    ## ==============================================
    ## parse command line arguments.
    ## ==============================================
    i = 1
    while i < len(argv):
        arg = argv[i]
        if arg == '--quiet' or arg == '-q':
            want_verbose = False
        elif arg == '--help' or arg == '-h':
            print(lasfun_usage)
            sys.exit(1)
        elif arg == '--version' or arg == '-v':
            print('{}, version {}'.format(os.path.basename(sys.argv[0]), cudem.__version__))
            sys.exit(1)
        elif arg[0] == '-':
            print(lasfun_usage)
            sys.exit(0)
        else: inlasen.append(arg)
        i = i + 1

    if len(inlasen) == 0:
        print(lasfun_usage)
        utils.echo_error_msg('you must specify at least one las file')
        sys.exit(-1)
    
    for rn, this_las in enumerate(inlasen):

        LasHeader(this_las).dump_header()
        LasFile(this_las).dump_points("xyzircsupg")
        #printHeader(infile)

    
### End
