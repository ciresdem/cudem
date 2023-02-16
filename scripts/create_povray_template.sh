#!/bin/sh

# $1 input grid

# GDAL

minmax=$(gdal_minmax.py $1)
echo $(minmax)

xdim=$(echo $minmax | awk '{print $7}')
ydim=$(echo $minmax | awk '{print $8}')
ymin=$(echo $minmax | awk '{print $4}')
ymax=$(echo $minmax | awk '{print $3}')
zmin=$(echo $minmax | awk '{print $5}')
zmax=$(echo $minmax | awk '{print $6}')

deg_lat_len=$(lll.py $ymin | awk '{print $2}')
deg_lon_len=$(lll.py $ymin | awk '{print $1}')

# GMT
# xdim=$(grdinfo $1 | grep x_min | awk '{print $11}')
# ydim=$(grdinfo $1 | grep y_min | awk '{print $11}')
# ymin=$(grdinfo $1 | grep y_min | awk '{print $3}')
# ymax=$(grdinfo $1 | grep y_min | awk '{print $5}')
# zmin=$(grdinfo $1 | grep z_min | awk '{print $3}')
# zmax=$(grdinfo $1 | grep z_min | awk '{print $5}')
# deg_lat_len=$(lonlatlen.scm $ymin | awk '{print $2}')
# deg_lon_len=$(lonlatlen.scm $ymin | awk '{print $1}')

echo "
// DEM

global_settings { assumed_gamma 1 }
#include \"colors.inc\" 
#declare Bi = 2;

//
// Custom parameters start here
//
#declare rgb_image = \"rgb.png\"
#declare dem_image = \"dem_16bit.png\"

#declare xdim = $xdim;  //number of pixels in X-direction
#declare ydim = $ydim;  //number of pixels in y-direction
#declare max_y = $ymax; //maximum latitude extent
#declare min_y = $ymin; //minimum latitude extent
#declare min_z = $zmin; //minimum elevation
#declare max_z = $zmax; //maximum elevation

// Obtained from http://www.csgnetwork.com/degreelenllavcalc.html  
#declare deg_lat_len = $deg_lat_len; //length of a degree of latitude in meters  
#declare deg_lon_len = $deg_lon_len; //length of a degree of longitude in meters

// Position of camera             
#declare cam_azimuth = 310;
#declare cam_elevation = 27;
#declare cam_distance = 235; 
#declare cam_view_angle = 33;
                     
// Position of the \"sun\"  
#declare light_azimuth = cam_azimuth+90;
#declare light_elevation = 20;
#declare light_distance = 10000;             

#declare vertical_exaggeration = 2;
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
texture {
  pigment{
    image_map { 
      png rgb_image map_type 0 once interpolate Bi 
    } 
  } 
  finish { ambient 0.4 diffuse 0.8 } 
  rotate x*90  
}

                  
height_field { // Unsigned 16-bit PNG DEM
   png dem_image 
   smooth
   clipped_by {box { <0, 0, 0>, <0.999, 1, 0.999> } }
   texture { Texture0 }
   translate <-0.5, 0, -0.5>
   scale <100*lon_scale*xdim/ydim,
          vertical_exaggeration*z_scale,  //Vertical exaggeration
          100>
} 


camera {
   angle cam_view_angle
   location <cam_x, cam_y, cam_z>
   look_at <0, 0, 0> 
}
                                    
light_source { <light_x, light_y, light_z> color White shadowless parallel }
        
background { White } 
"
