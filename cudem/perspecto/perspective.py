### perspective.py 
##
## Copyright (c) 2023 - 2025 Regents of the University of Colorado
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

from cudem.perspecto import perspecto
from cudem.perspecto import povray
from cudem import utils

class perspective(povray.POVRay):
    """Generate a perspective image

< perspective:cam_azimuth=-130:cam_elevation=30:cam_distance=265:cam_view_angle=40:light_elevation=20:light_distance=10000:vertical_exaggeration=2 >
    """
    
    def __init__(
            self,
            cam_azimuth=-130,
            cam_elevation=30,
            cam_distance=265,
            cam_view_angle=40,
            light_elevation=20,
            light_distance=10000,
            vertical_exaggeration=2,
            **kwargs
    ):
        super().__init__(**kwargs)
        self.cam_azimuth = cam_azimuth
        self.cam_elevation = cam_elevation
        self.cam_distance = cam_distance
        self.cam_view_angle = cam_view_angle
        self.light_elevation = light_elevation
        self.light_distance = light_distance
        self.vertical_exaggeration = vertical_exaggeration
        self.dem_lll = perspecto.lll(self.dem_region.ymin)
        self.output_pov = '{}_perspective.pov'.format(
            utils.fn_basename2(self.src_dem)
        )
        
        self.template = """
// DEM

global_settings {{ assumed_gamma 1 }}
#include \"colors.inc\" 

//
// Custom parameters start here
//
#declare rgb_image = \"rgb.png\"
#declare dem_image = \"dem_16bit.png\"

#declare xdim = {xdim};  //number of pixels in X-direction
#declare ydim = {ydim};  //number of pixels in y-direction
#declare max_y = {ymax}; //maximum latitude extent
#declare min_y = {ymin}; //minimum latitude extent
#declare min_z = {zmin}; //minimum elevation
#declare max_z = {zmax}; //maximum elevation

// Obtained from http://www.csgnetwork.com/degreelenllavcalc.html  
#declare deg_lat_len = {deg_lat_len}; //length of a degree of latitude in meters  
#declare deg_lon_len = {deg_lon_len}; //length of a degree of longitude in meters

// Position of camera             
#declare cam_azimuth = {cam_azimuth};
#declare cam_elevation = {cam_elevation};
#declare cam_distance = {cam_distance}; 
#declare cam_view_angle = {cam_view_angle};
                     
// Position of the \"sun\"  
#declare light_azimuth = cam_azimuth+90;
#declare light_elevation = {light_elevation};
#declare light_distance = {light_distance};             

#declare vertical_exaggeration = {vertical_exaggeration};
//
// Custom parameters end here
//
             
#declare lon_scale = deg_lon_len / deg_lat_len;
#declare z_scale = (100 * (max_z - min_z)) / (deg_lat_len * (max_y - min_y));

#declare cam_x = cam_distance * cos(radians(cam_elevation)) * sin(radians(cam_azimuth));
#declare cam_y = cam_distance * sin(radians(cam_elevation));
#declare cam_z = cam_distance * cos(radians(cam_elevation)) * cos(radians(cam_azimuth));
#declare light_x = light_distance * cos(radians(light_elevation)) * sin(radians(light_azimuth));
#declare light_y = light_distance * sin(radians(light_elevation));
#declare light_z = light_distance * cos(radians(light_elevation)) * cos(radians(light_azimuth));

#declare Texture0 = // Color elevation image (24-bit RGB PNG)
texture {{
  pigment{{
    image_map {{ 
      png "{rgb_image}" map_type 0 once interpolate 2
    }} 
  }} 
  finish {{ ambient 0.4 diffuse 0.8 }}
  rotate x*90  
}}
                  
height_field {{ // Unsigned 16-bit PNG DEM
   png "{dem_image}"
   smooth
   clipped_by {{box {{ <0, 0, 0>, <0.999, 1, 0.999> }} }}
   texture {{ Texture0 }}
   translate <-0.5, 0, -0.5>
   scale <100*lon_scale*{xdim}/{ydim},
          vertical_exaggeration*z_scale,  //Vertical exaggeration
          100>
}} 

camera {{
   angle cam_view_angle
   location <cam_x, cam_y, cam_z>
   look_at <0, 0, 0> 
}}
                                    
light_source {{ <light_x, light_y, light_z> color White shadowless parallel }}
        
background {{ White }}""".format(
    rgb_image=self.rgb_image,
    dem_image=self.dem_image,
    xdim=self.dem_infos['nx'],
    ydim=self.dem_infos['ny'],
    ymax = self.dem_region.ymax,
    ymin = self.dem_region.ymin,
    zmax = self.dem_infos['zr'][1],
    zmin = self.dem_infos['zr'][0],
    deg_lat_len = self.dem_lll[0],
    deg_lon_len = self.dem_lll[1],
    cam_azimuth = self.cam_azimuth,
    cam_elevation = self.cam_elevation,
    cam_distance = self.cam_distance,
    cam_view_angle = self.cam_view_angle,
    light_elevation = self.light_elevation,
    light_distance = self.light_distance,
    vertical_exaggeration = self.vertical_exaggeration,
)

    def run(self):
        with open(self.output_pov, 'w') as pov_temp:
            pov_temp.write(self.template)
        
        self.run_povray(
            self.output_pov, self.dem_infos['nx'], self.dem_infos['ny']
        )

### End
