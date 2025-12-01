### mbsfile.py - DataLists IMproved
##
## Copyright (c) 2010 - 2025 Regents of the University of Colorado
##
## mbsfile.py is part of CUDEM
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
### Examples:
##
### TODO:
##
### Code:

import os
import math
import numpy as np

from osgeo import gdal
from osgeo import ogr

from cudem import utils
from cudem import regions
from cudem import gdalfun
from cudem.datalists.dlim import ElevationDataset

class MBSParser(ElevationDataset):
    """providing an mbsystem parser

    Process MB-System supported multibeam data files.
    See MB-System for more information regarding supported file 
    formats, etc.
    
    generate_inf - generate an inf file for the MBS data
    yield_xyz - yield the MBS data as xyz
    yield_array - yield the MBS data as an array
    
    -----------
    Parameters:
    
    mb_fmt=[]
    mb_exclude=[]
    """

    def __init__(self,
                 mb_fmt=None,
                 mb_exclude='A',
                 want_mbgrid=False,
                 want_binned=False,
                 min_year=None,
                 auto_weight=True,
                 auto_uncertainty=True,
                 **kwargs):
        super().__init__(**kwargs)
        self.mb_fmt = mb_fmt
        self.mb_exclude = mb_exclude
        self.want_mbgrid = want_mbgrid
        self.want_binned = want_binned
        self.min_year = min_year
        self.auto_weight = auto_weight
        self.auto_uncertainty = auto_uncertainty
        if self.src_srs is None:
            #self.src_srs = 'epsg:4326+3855'
            self.src_srs = 'epsg:4326'
            
            
    def inf_parse(self):
        self.infos.minmax = [0,0,0,0,0,0]
        this_row = 0
        xinc = 0
        yinc = 0
        dims = []
        inf_fn = '{}.inf'.format(self.fn) \
            if self.fn.split('.')[-1] != 'inf' \
               else self.fn
        
        with open(inf_fn) as iob:
            for il in iob:
                til = il.split()
                if len(til) > 1:

                    if til[0] == 'Swath':
                        if til[2] == 'File:':
                            self.infos.name = til[3]

                    if ' '.join(til[:-1]) == 'MBIO Data Format ID:':
                        self.mb_fmt = til[-1]
                            
                    if til[0] == 'Number':
                        if til[2] == 'Records:':
                            self.infos.numpts = utils.int_or(til[3])
                            
                    if til[0] == 'Minimum':
                        if til[1] == 'Longitude:':
                            self.infos.minmax[0] = utils.float_or(til[2])
                            self.infos.minmax[1] = utils.float_or(til[5])
                        elif til[1] == 'Latitude:':
                            self.infos.minmax[2] = utils.float_or(til[2])
                            self.infos.minmax[3] = utils.float_or(til[5])
                        elif til[1] == 'Depth:':
                            self.infos.minmax[4] = utils.float_or(til[5]) * -1
                            self.infos.minmax[5] = utils.float_or(til[2]) * -1
                            
                    if til[0] == 'CM':
                        if til[1] == 'dimensions:':
                            dims = [utils.int_or(til[2]), utils.int_or(til[3])]
                            cm_array = np.zeros((dims[0], dims[1]))
                            
                    if til[0] == 'CM:':
                        for j in range(0, dims[0]):
                            cm_array[this_row][j] = utils.int_or(til[j+1])
                        this_row += 1

        mbs_region = regions.Region().from_list(self.infos.minmax)
        if len(dims) > 0:
            xinc = (mbs_region.xmax - mbs_region.xmin) / dims[0]
            yinc = (mbs_region.ymin - mbs_region.ymax) / dims[1]
        
        if abs(xinc) > 0 and abs(yinc) > 0:
            xcount, ycount, dst_gt = mbs_region.geo_transform(
                x_inc=xinc, y_inc=yinc
            )
            ds_config = {'nx': dims[0], 'ny': dims[1], 'nb': dims[1] * dims[0],
                         'geoT': dst_gt, 'proj': gdalfun.osr_wkt(self.src_srs),
                         'dt': gdal.GDT_Float32, 'ndv': 0, 'fmt': 'GTiff'}
            driver = gdal.GetDriverByName('MEM')
            ds = driver.Create(
                'tmp', ds_config['nx'], ds_config['ny'], 1, ds_config['dt']
            )
            ds.SetGeoTransform(ds_config['geoT'])
            if ds_config['proj'] is not None:
                ds.SetProjection(ds_config['proj'])
                
            ds_band = ds.GetRasterBand(1)
            ds_band.SetNoDataValue(ds_config['ndv'])
            ds_band.WriteArray(cm_array)
            tmp_ds = ogr.GetDriverByName('Memory').CreateDataSource('tmp_poly')
            tmp_layer = tmp_ds.CreateLayer('tmp_poly', None, ogr.wkbMultiPolygon)
            tmp_layer.CreateField(ogr.FieldDefn('DN', ogr.OFTInteger))
            gdal.Polygonize(ds_band, ds_band, tmp_layer, 0)
            multi = ogr.Geometry(ogr.wkbMultiPolygon)
            for feat in tmp_layer:
                feat.geometry().CloseRings()
                wkt = feat.geometry().ExportToWkt()
                multi.AddGeometryDirectly(ogr.CreateGeometryFromWkt(wkt))
                
            wkt = multi.ExportToWkt()
            tmp_ds = ds = None
        else:
            wkt = mbs_region.export_as_wkt()

        self.infos.wkt = wkt

        if self.src_srs is None:
            self.infos.src_srs = 'epsg:4326'
        else:
            self.infos.src_srs = self.src_srs
        
        return(self)

    
    def parse_(self):
        """use mblist to convert data to xyz then set the dataset as xyz and 
        use that to process...
        """
        
        mb_fn = os.path.join(self.fn)
        if self.region is not None:
            mb_region = self.region.copy()
            mb_region.buffer(pct=25)
        else:
            mb_region = None

        out_mb = utils.make_temp_fn(
            '{}_.xyz'.format(utils.fn_basename2(mb_fn)),
            self.cache_dir
        )
        out, status = utils.run_cmd('mblist -M{}{} -OXYZ -I{} > {}'.format(
            self.mb_exclude, ' {}'.format(
                mb_region.format('gmt') if mb_region is not None else ''
            ), mb_fn, out_mb
        ), verbose=True)

        if status == 0:
            data_set = DatasetFactory(
                **self._set_params(
                    mod=out_mb, data_foramt=168
                )
            )._acquire_module.initialize()
            
            # fill self.data_entries with each dataset for use outside the yield.            
            for ds in data_set.parse(): 
                self.data_entries.append(ds) 
                yield(ds)

        utils.remove_glob('{}*'.format(out_mb))

        
    def mb_inf_format(self, src_inf):
        """extract the format from the mbsystem inf file."""

        with open(src_inf, errors='ignore') as iob:
            for il in iob:
                til = il.split()
                if len(til) > 1:
                    #utils.echo_msg(til[0].strip())
                    if til[0].strip() == 'MBIO':
                        return(til[4])
        return(None)

    
    def mb_inf_data_date(self, src_inf):
        """extract the date from the mbsystem inf file."""

        if os.path.exists(src_inf):
            with open(src_inf, errors='ignore') as iob:
                for il in iob:
                    til = il.split()
                    if len(til) > 1:
                        if til[0] == 'Time:':
                            return(til[3])
        return(None)

    
    def mb_inf_perc_good(self, src_inf):
        """extract the data format from the mbsystem inf file."""

        if os.path.exists(src_inf):
            with open(src_inf, errors='ignore') as iob:
                for il in iob:
                    til = il.split(':')
                    if len(til) > 1:
                        if til[0].strip() == 'Number of Good Beams':
                            return(til[1].split()[-1].split('%')[0])
        return(None)

    
    def yield_mbgrid_ds(self):
        """process the data through mbgrid and use GDALFile to further 
        process the gridded data
        """
        
        mb_fn = os.path.join(self.fn)
        with open('_mb_grid_tmp.datalist', 'w') as tmp_dl:
            tmp_dl.write(
                '{} {} {}\n'.format(
                    self.fn,
                    self.mb_fmt if self.mb_fmt is not None else '',
                    self.weight if self.mb_fmt is not None else ''
                )
            )

        ofn = '_'.join(os.path.basename(mb_fn).split('.')[:-1])
        try:
            utils.run_cmd(
                (f'mbgrid -I_mb_grid_tmp.datalist {self.data_region.format("gmt")} ' +
                 f'-E{self.x_inc}/{self.y_inc}/degrees! -O{ofn} -A2 -F1 -C2/1 -S0 -T35 -W.1'),
                verbose=True
            )

            gdalfun.gdal2gdal('{}.grd'.format(ofn))
            utils.remove_glob(
                '_mb_grid_tmp.datalist',
                '{}.cmd'.format(ofn),
                '{}.mb-1'.format(ofn),
                '{}.grd*'.format(ofn)
            )
            
            grits_filter = grits.grits.GritsFactory(
                mod='outliers:multipass=2',
                src_dem='{}.tif'.format(ofn),
                cache_dir=self.cache_dir,
            )._acquire_module()
            
            if grits_filter is not None:
                try:
                    grits_filter()
                    os.replace(grits_filter.dst_dem, '{}'.format(ofn))
                except:
                    pass

            mbs_ds = DatasetFactory(
                **self._set_params(
                    mod='{}.tif'.format(ofn),
                    data_format=200
                )
            )._acquire_module()
            mbs_ds.initialize()
            for gdal_ds in mbs_ds.parse():
                for pts in gdal_ds.transform_and_yield_points():
                    yield(pts)
            
            #yield(mbs_ds)
            utils.remove_glob('{}*.tif*'.format(ofn))
        except Exception as e:
            #pass
            raise(e)

        
    def yield_mblist_ds(self):
        """use mblist to process the multibeam data"""

        if not self.fn.endswith('.inf'):
            mb_fn = os.path.join(self.fn)
            #xs = []
            #ys = []
            #zs = []
            #ws = []
            #us = []
            mb_points = []
            if self.region is not None:
                mb_region = self.region.copy()
                mb_region.buffer(pct=25)
            else:
                mb_region = None

            src_inf = '{}.inf'.format(self.fn)
            try:
                mb_format = self.mb_inf_format(src_inf)
                if mb_fn.split('.')[-1] == 'fbt':
                    mb_format = None
            except:
                mb_format = None

            mb_date = self.mb_inf_data_date(src_inf)
            mb_perc = self.mb_inf_perc_good(src_inf)
            this_year = int(utils.this_year())
            if self.auto_weight:
                if mb_date is not None:
                    #this_weight = min(0.9, .5*(int(mb_date)-1900)/(2020-1900))
                    #this_weight = max(0.01, 1 - ((this_year - int(mb_date)) / (this_year - 1990)))
                    this_weight = min(0.99, max(0.01, 1 - ((2024 - int(mb_date))) / (2024 - 1980)))
                    #this_weight = 1
                else:
                    this_weight = 1

                if self.weight is not None:
                    self.weight *= this_weight
                
            #'mblist -M{}{} -OXYZDSc -I{}{}'.format(
            #                                       'mblist -M{}{} -OXYZ -I{}{}'.format(
            #                    'mblist -M{}{} -OXYZDAGgFPpRrSCc -I{}{}'.format(
            for line in utils.yield_cmd(
                    'mblist -M{}{} -OXYZDS -I{}{}'.format(
                        self.mb_exclude, ' {}'.format(
                            mb_region.format('gmt') \
                            if mb_region is not None \
                            else ''
                        ), mb_fn,
                        ' -F{}'.format(mb_format) \
                        if mb_format is not None \
                        else ''
                    ),
                    verbose=False,
            ):
                this_line = [float(x) for x in line.strip().split('\t')]
                x = this_line[0]
                y = this_line[1]
                z = this_line[2]
                crosstrack_distance = this_line[3]
                speed = this_line[4]
                if speed == 0:
                    continue
                #roll = this_line[5]
                #utils.echo_msg(roll)
                #z = roll
                #u_depth = ((2+(0.02*(z*-1)))*0.51)
                #u = math.sqrt(1 + ((.023 * abs(crosstrack_distance))**2))
                if self.auto_weight:
                    #w = max(.5, 1/(1 + .005*abs(crosstrack_distance)))
                    u_depth = ((.25+(0.02*(z*-1)))*0.51)
                    u_cd =  .005*abs(crosstrack_distance)
                    #u_s = math.sqrt(1 + ((.51 * abs(speed))**2))
                    tmp_speed = min(14, abs(speed - 14))
                    u_s = 1 * (tmp_speed*.51)
                    #u_s = 2 + (.5 * abs(tmp_speed)) if tmp_speed <=0 else 1 + (.005 * abs(tmp_speed))
                    u = math.sqrt(u_depth**2 + u_cd**2 + u_s**2)
                    w = 1/u
                else:
                    w = 1
                    u = 0

                #u = u_cd
                #u = speed
                #utils.echo_msg_bold()
                out_line = [x,y,z,w,u]
                mb_points.append(out_line)
                # if self.auto_weight or self.auto_uncertainty:
                #     x = this_line[0]
                #     y = this_line[1]
                #     z = this_line[2]
                #     crosstrack_distance = this_line[3]
                #     crosstrack_slope = this_line[4]
                #     flat_bottom_grazing_angle = this_line[5]
                #     seafloor_grazing_angle = this_line[6]
                #     beamflag = this_line[7]
                #     pitch = this_line[8]
                #     draft = this_line[9]
                #     roll = this_line[10]
                #     heave = this_line[11]
                #     speed = this_line[12]
                #     sonar_alt = this_line[13]
                #     sonar_depth = this_line[14]
                #     if int(beamflag) == 0:
                #         xs.append(x)
                #         ys.append(y)
                #         zs.append(z)
                #         ws.append(1)
                #         ## uncertainty
                #         u_depth = ((2+(0.02*(z*-1)))*0.51)
                #         u_s_depth = ((2+(0.02*(sonar_depth*-1)))*0.51)
                #         u_cd = math.sqrt(1 + ((.023 * abs(crosstrack_distance))**2))
                #         if speed >= 25:
                #             u_s = math.sqrt(1 + ((.51 * abs(speed))**2))
                #             u = math.sqrt(u_depth**2 + u_cd**2 + u_s**2)
                #         else:
                #             u = math.sqrt(u_depth**2 + u_cd**2)

                #         us.append(u)
                #         #if self.auto_weight:
                #         ## weight
                #         #w = math.sqrt((1/u)) * this_weight
                #         # w = this_weight
                #         # w *= self.weight if self.weight is not None else 1
                #         # ws.append(w)
                # else:
                # xs.append(this_line[0])
                # ys.append(this_line[1])
                # zs.append(this_line[2])
                # ws.append(1)
                # us.append(0)

            if len(mb_points) > 0:
                #mb_points = np.column_stack((xs, ys, zs, ws, us))
                #xs = ys = zs = ws = us = None
                ##_ = [x.extend([1,0]) for x in mb_points]
                # for x in mb_points:
                #     x.extend([1,0])

                mb_points = np.rec.fromrecords(
                    mb_points, names='x, y, z, w, u'
                )
                if mb_points is not None:
                    yield(mb_points)

                mb_points = None

                
    def yield_mblist2_ds(self):
        """use mblist to process the multibeam data"""
        
        mb_fn = os.path.join(self.fn)
        if self.region is None and self.data_region is None:
            self.want_mbgrid = False

        if self.region is not None:
            mb_region = self.region.copy()
            mb_region.buffer(pct=25)
        else:
            mb_region = None

        #, '-F{}'.format(self.mb_fmt) if self.mb_fmt is not None else ''
        mb_points = [[float(x) for x in l.strip().split('\t')] \
                     for l in utils.yield_cmd(
                             'mblist -M{}{} -OXYZ -I{}'.format(
                                 self.mb_exclude, ' {}'.format(
                                     mb_region.format('gmt') \
                                     if mb_region is not None \
                                     else ''
                                 ), mb_fn
                             ),
                             verbose=False,
                     )]

        if len(mb_points) > 0:
            mb_points = np.rec.fromrecords(mb_points, names='x, y, z')
        else:
            mb_points = None

        if self.want_binned:
            #mb_points = self.bin_z_points(mb_points)
            point_filter = PointFilterFactory(
                mod='outlierz', points=mb_points
            )._acquire_module()
            if point_filter is not None:
                mb_points = point_filter()

        if mb_points is not None:
            yield(mb_points)

            
    def yield_points(self):        
        mb_fn = os.path.join(self.fn)
        if self.region is None or self.data_region is None:
            self.want_mbgrid = False

        ## update want_mbgrid to output points!
        if self.want_mbgrid \
           and (self.x_inc is not None \
                and self.y_inc is not None):
            for pts in self.yield_mbgrid_ds():
                yield(pts)
        else:
            for pts in self.yield_mblist_ds():
                yield(pts)


### End
