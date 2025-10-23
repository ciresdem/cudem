### sphere.py 
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

from cudem.perspecto import povray
from cudem import utils

class sphere(povray.POVRay):
    """Generate a sphere image

< sphere:cam_azimuth=310:cam_elevation=27:cam_distance=8:cam_view_angle=33:center_lat=None:center_long=None >
    """
    
    def __init__(
            self,
            cam_azimuth=310,
            cam_elevation=27,
            cam_distance=8,
            cam_view_angle=33,
            center_lat=None,
            center_long=None,
            **kwargs
    ):
        super().__init__(**kwargs)
        self.cam_azimuth = cam_azimuth
        self.cam_elevation = cam_elevation
        self.cam_distance = cam_distance
        self.cam_view_angle = cam_view_angle

        if center_lat is None or center_long is None:
            region_center = self.dem_region.center()
            self.center_long = region_center[0]
            self.center_lat = region_center[1]
        else:
            self.center_long = center_long
            self.center_lat = center_lat

        self.output_pov = '{}_sphere.pov'.format(
            utils.fn_basename2(self.src_dem)
        )            
        self.template = """
//Generates a hillshade color image on a sphere 
//CUDEM 2023MAR08

#include "colors.inc"  

#declare cam_azimuth = {cam_azimuth};
#declare cam_elevation = {cam_elevation};
#declare cam_distance = {cam_distance}; 
#declare cam_view_angle = {cam_view_angle};
        
#declare cam_x = cam_distance * cos(radians(cam_elevation)) * sin(radians(cam_azimuth));
#declare cam_y = cam_distance * sin(radians(cam_elevation));
#declare cam_z = cam_distance * cos(radians(cam_elevation)) * cos(radians(cam_azimuth));
        
// Colormap image  
#declare colorImage = 
pigment {{
  image_map {{
    png "{rgb_image}"     // Specify color map image
    map_type 1        // 0=planar, 1=spherical, 2=cylindrical, 5=torus
    interpolate 2     // 0=none, 1=linear, 2=bilinear, 4=normalized distance
  }} // image_map
}} 

// ETOPO1 topography used as a bump map
#declare hillShade = 
normal {{
  bump_map {{             // uses image color or index as bumpiness
    png "{dem_image}"  // the file to read (tiff/tga/gif/png/jpeg/tiff/sys)
    map_type 1           // 0=planar, 1=spherical, 2=cylindrical, 5=torus
    interpolate 2        // 0=none, 1=linear, 2=bilinear, 4=normalized distance
    bump_size 750        // Adjust the vertical exaggeration (typical values between 500 and 1000)
  }} // bump_map
}}

// Apply pigment to the texture 
#declare colorShade =
texture {{
  pigment {{
    colorImage
  }}
  finish {{ambient 0.3 diffuse 1}}
  normal {{ hillShade }}
}}

//Adjust values to desired cooridinates
#declare longitude = {center_long};
#declare latitude = {center_lat};

sphere {{
  <0,0,0>, 1
  texture {{ colorShade }}
  //rotate y*90 
  no_shadow
  //rotate y*longitude
  //rotate x*latitude
  //rotate x*20
}}

camera {{
   angle cam_view_angle
   location <cam_x, cam_y, cam_z>
   look_at <0, 0, 0> 
}}""".format(
    rgb_image=self.rgb_image,
    dem_image=self.dem_image,
    center_long=self.center_long,
    center_lat=self.center_lat,
    cam_azimuth = self.cam_azimuth,
    cam_elevation = self.cam_elevation,
    cam_distance = self.cam_distance,
    cam_view_angle = self.cam_view_angle,
)

    def run(self):

        with open(self.output_pov, 'w') as pov_temp:
            pov_temp.write(self.template)
        
        self.run_povray(
            self.output_pov, self.dem_infos['nx'], self.dem_infos['ny']
        )

### End
