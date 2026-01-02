### perspective.py 
##
## Copyright (c) 2023 - 2026 Regents of the University of Colorado
##
## perspective.py is part of CUDEM
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
## Generate a 3D perspective scene using POV-Ray.
##
### Code:

from cudem import utils
from cudem.perspecto import povray


class Perspective(povray.PovRay):
    """Generate a perspective image using POV-Ray.

    Configuration Example:
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
        
        ## Calculate degrees length at the minimum latitude
        self.dem_lll = utils.lll(self.dem_region.ymin)
        
        basename = utils.fn_basename2(self.src_dem)
        self.output_pov = f'{basename}_perspective.pov'
        self.output_png = f'{basename}_perspective.png'
        
        ## Extract variables for cleaner f-string usage
        nx = self.dem_infos['nx']
        ny = self.dem_infos['ny']
        ymax = self.dem_region.ymax
        ymin = self.dem_region.ymin
        zmax = self.dem_infos['zr'][1]
        zmin = self.dem_infos['zr'][0]
        deg_lat_len = self.dem_lll[0]
        deg_lon_len = self.dem_lll[1]

        ## POV-Ray Template
        ## Note: Double braces {{ }} are used for POV-Ray syntax. 
        ## Single braces { } are used for Python f-string variables.
        self.template = f"""
// DEM

global_settings {{ assumed_gamma 1 }}
#include "colors.inc" 

//
// Custom parameters start here
//
#declare rgb_image = "{self.rgb_image}"
#declare dem_image = "{self.dem_image}"

#declare xdim = {nx};  //number of pixels in X-direction
#declare ydim = {ny};  //number of pixels in y-direction
#declare max_y = {ymax}; //maximum latitude extent
#declare min_y = {ymin}; //minimum latitude extent
#declare min_z = {zmin}; //minimum elevation
#declare max_z = {zmax}; //maximum elevation

// Obtained from http://www.csgnetwork.com/degreelenllavcalc.html  
#declare deg_lat_len = {deg_lat_len}; //length of a degree of latitude in meters  
#declare deg_lon_len = {deg_lon_len}; //length of a degree of longitude in meters

// Position of camera             
#declare cam_azimuth = {self.cam_azimuth};
#declare cam_elevation = {self.cam_elevation};
#declare cam_distance = {self.cam_distance}; 
#declare cam_view_angle = {self.cam_view_angle};
                     
// Position of the "sun"  
#declare light_azimuth = cam_azimuth+90;
#declare light_elevation = {self.light_elevation};
#declare light_distance = {self.light_distance};             

#declare vertical_exaggeration = {self.vertical_exaggeration};
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
      png "{self.rgb_image}" map_type 0 once interpolate 2
    }} 
  }} 
  finish {{ ambient 0.4 diffuse 0.8 }}
  rotate x*90  
}}
                  
height_field {{ // Unsigned 16-bit PNG DEM
   png "{self.dem_image}"
   smooth
   clipped_by {{box {{ <0, 0, 0>, <0.999, 1, 0.999> }} }}
   texture {{ Texture0 }}
   translate <-0.5, 0, -0.5>
   scale <100*lon_scale*{nx}/{ny},
          vertical_exaggeration*z_scale,  //Vertical exaggeration
          100>
}} 

camera {{
   angle cam_view_angle
   location <cam_x, cam_y, cam_z>
   look_at <0, 0, 0> 
}}
                                    
light_source {{ <light_x, light_y, light_z> color White shadowless parallel }}
        
background {{ White }}"""

    def run(self):
        """Write the POV file and execute POV-Ray."""
        
        try:
            with open(self.output_pov, 'w') as pov_temp:
                pov_temp.write(self.template)
            
            self.run_povray(
                self.output_pov, 
                pov_width=self.dem_infos['nx'], 
                pov_height=self.dem_infos['ny']
            )
        except OSError as e:
            utils.echo_error_msg(f"Failed to run POV-Ray: {e}")

        return self.output_png

### End
