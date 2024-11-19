# DEM Uncertainty

While generating a DEM from independent scattered data sources we optionally calculate the uncertainty for each cell of the resulting DEM.
An auxiliary Uncertainty Raster can be generated along with a DEM by using the `--uncertainty` switch in the [waffles](/docs/waffles.md) command. Additionally, the waffles module 'uncertainty' can be used to only generate an uncertainty raster without generating a DEM.

## Types of uncertainty:

- [Source Data Uncertainty](#source-data-uncertainty)
  - [Source dataset uncertainty](#source-dataset-uncertainty) (dataset-wide)
  - [Source value uncertainty](#source-value-uncertainty) (per data value)
  - [Bathymetric depth uncertainty](#bathymetric-depth-uncertainty) (IHO function of depth)
- [Gridding Uncertainty](#gridding-uncertainty)
  - [Sub-pixel uncertainty](#sub-pixel-variance-and-standard-error) (variance)
  - [Interpolation uncertainty](#interpolation-uncertainty) (split-sample)
  - [Vertical datum transformation uncertainty](#vertical-datum-transformation-uncertainty)

The various uncertainty types are combined to report a Total Value Uncertainty (TVU) as supplemental raster product.

Uncertainty values are combined using the Root Sum Squared (RSS).

### Source Data Uncertainty
#### Source dataset Uncertainty

The dataset-wide source uncertainty is a single value that will be applied to an entire dataset. This is specified by the user in the 4th column of the dataset entry. Typically this value is provided by the data collector or processor in the datasets metadata. This can sometimes be reported as an RMSE or data Accuracy. If the value is reported as being in the 95th percentile confidence level, first divide that value by 1.96 to obtain an uncertainty value suitable for combining with other uncertainty values.

Specify the source data uncertainty as the 4th column in the datalist data-entry. See [dlim](/docs/dlim.md) for more information about datalist and data-entry formatting.

#### Source value Uncertainty

Some datasets, such as NOS BAG or BlueTopo, specify an uncertainty value for each data value in the dataset. Other times a user may independently calculate the uncertainty for each data value in their dataset. These uncertainty data can be used to inform the final TVU by specifying the uncertainty data of the dataset as either a seperate product or integrated into the dataset, such as with raster data or xyz data, respsectively.

For XYZ data, we can include per-point uncertainty as the 5th column in the xyz dataset, where the 4th column is the weight.

For Raster data, we can either specify a band number or a separate raster file that contains the per-pixel uncertainty data.

This is applied automatically when using fetches modules as a datalist dataset-entry.

#### Bathymetric Depth Uncertainty

For bathymetric data, the [IHO standards](https://iho.int/uploads/user/pubs/standards/s-44/S-44_Edition_6.1.0.pdf) can be used to calculate the uncertainty of each data value as a function of it's depth, where

`TVU(d) = sqrt(a**2 + (b * d)**2)`

| **Order 2** | **Order 1b** | **Order 1a** | **Special Order** | **Exclusive Order** |
|-------------|--------------|--------------|-------------------|---------------------|
| Areas where a general description of the sea floor is considered adequate. | Areas where underkeel clearance is not considered to be an issue for the type of surface shipping expected to transit the area. | Areas where underkeel clearance is considered not to be critical but features of concern to surface shipping may exist. | Areas where underkeel clearance is critical | Areas where there is strict minimum underkeel clearance and manoeuvrability criteria |
| a = 1.0m, b = 0.023 | a = 0.5m, b = 0.013 | a = 0.5m, b = 0.013 | a = 0.25m, b = 0.0075 | a = 0.15m, b = 0.0075 |

This is applied automatically when using certain fetches modules as a datalist dataset-entry, such as `hydronos` and `multibeam`.

### Gridding Uncertainty
#### Sub-pixel Variance and Standard Error
The sub-pixel variance is calculated by default when combining various datasets together into a DEM. DEM values typically represent the (optionally weighted) mean of the data measurements within an individual DEM pixel. The sub-pixel variance is the spread of the measurement values within the pixel from the mean value and the standard deviation is the square root of the variance. The uncertainty of the gridded mean value can be expressed by the standard deviation of the mean, which is also commonly known as the standard error of the mean, or simply the standard error, and is equal to the standard deviation divided by the square root of the number of measurements within the DEM pixel. The number of measurements within each DEM pixel is determined by the spatial resolution of the DEM.

The sub-pixel variance is combined with the Source Data Uncertainty of the individual measurements described in the previous sections to calculate the pooled variance. The pixel-level pooled variance is equal to the square of the weighted mean of the Source Data Uncertainty for all measurements plus the weighted variance of all measurements around the weighted mean elevation. 

The pixel-level standard error is then calculated from the pooled variance and the number of measurements within the DEM pixel. This method used to calculate the pooled variance and standard error is described in:

Amante, C. J. (2018). Estimating coastal digital elevation model
uncertainty. Journal of Coastal Research, 34(6), 1382-1397. https://doi.org/10.2112/JCOASTRES-D-17-00211.1

#### Interpolation Uncertainty

The interpolation uncertainty is calculated for all interpolated cells in a resulting DEM using a split-sample method.

The method used to calculate interpolation uncertainty is described in:

Amante, C. J. (2018). Estimating coastal digital elevation model
uncertainty. Journal of Coastal Research, 34(6), 1382-1397. https://doi.org/10.2112/JCOASTRES-D-17-00211.1

#### Vertical Datum Transformation Uncertainty

Whenever data is vertically transformed while processing a DEM, the uncertainty of that transformation is accumulated into the final TVU. The uncertainty values for these transformations are obtained from NOAA's [VDatum](https://vdatum.noaa.gov/)

[Estimation of Vertical Uncertainties in VDatum](https://vdatum.noaa.gov/docs/est_uncertainties.html)

## Examples

- [Androscoggin](#androscoggin)
- [Point Reyes](#point-reyes)

### Androscoggin

In this example we will generate tiled uncertainty rasters of the Androscoggin region using input lidar elevation data from the USGS.

#### Setup directories to hold data and dems

```bash
mkdir androscoggin
cd androscoggin
mkdir data software
cd software
```

#### Define regions

Generate a tile-set vector for the Androscoggin test region with tiles at 5000x5000 meters

```
regions -R -R379500.00/389500.00/4875000.00/4890000.00 -T 5000
```
outputs 'regions_tile_set.shp'

#### Obtain the USGS elevation data

Either gather the relevant datasets from USGS or use the ```fetches``` command to fetch them

#### gather data with fetches

Transform the region to WGS84 and buffer it for fetching

```bash
$ fetches $(regions -R379500.00/389500.00/4875000.00/4890000.00 -J epsg:6348 -P epsg:4326 -e -b 0.01) tnm:q=LPC
```

At the time of this example, we end up with 5 lidar datasets, which are fetched to the 'tnm' directory:

- ME_Maine_LiDAR_NRCS_14
- ME_SouthCoastal_2020_A20
- ME_Western_2016
- NH_Connecticut_River_2015
- NH_Umbagog

#### Move them all to '../data' for processing and storage.

```bash
mdkir ../data/ME_Maine_LiDAR_NRCS_14
mdkir ../data/ME_SouthCoastal_2020_A20
mdkir ../data/ME_Western_2016
mdkir ../data/NH_Connecticut_River_2015
mdkir ../data/NH_Umbagog

mv tnm/*ME_Maine_LiDAR_NRCS_14* ../data/ME_Maine_LiDAR_NRCS_14
mv tnm/*ME_SouthCoastal_2020_A20* ../data/ME_SouthCoastal_2020_A20
mv tnm/*ME_Western_2016* ../data/ME_Western_2016
mv tnm/*NH_Connecticut_River_2015* ../data/NH_Connecticut_River_2015
mv tnm/*NH_Umbagog* ../data/NH_Umbagog
```

#### Create Datalists

See [dlim](/docs/dlim.md) for more information of datlist formatting.

##### Create a datalist for each of the lidar datasets and generate their associated inf and geojson files

```bash
cd ../data
cd ME_Maine_LiDAR_NRCS_14
dlim -g > ME_Maine_LiDAR_NRCS_14.datalist
dlim -i ME_Maine_LiDAR_NRCS_14.datalist

cd ../ME_SouthCoastal_2020_A20
dlim -g > ME_SouthCoastal_2020_A20.datalist
dlim -i ME_SouthCoastal_2020_A20.datalist

cd ../ME_Western_2016
dlim -g > ME_Western_2016.datalist
dlim -i ME_Western_2016.datalist

cd ../NH_Connecticut_River_2015
dlim -g > NH_Connecticut_River_2015.datalist
dlim -i NH_Connecticut_River_2015.datalist

cd ../NH_Umbagog
dlim -g > NH_Umbagog.datalist
dlim -i NH_Umbagog.datalist
cd ../../software
```

##### Create the main datalist which points to each of the lidar dataset datalists, assigning weights and source uncertainty to each

```bash
nano Androscoggin.datalist
```

Add the following lines to the datalist
```
../data/ME_Maine_LiDAR_NRCS_14/ME_Maine_LiDAR_NRCS_14.datalist -1 1 0
../data/ME_SouthCoastal_2020_A20/ME_SouthCoastal_2020_A20.datalist -1 2 0
../data/ME_Western_2016/ME_Western_2016.datalist -1:mask="tmp/USGS_Maine_LiDAR_Processing_Boundary_Buffered_100m.shp":invert_mask=True 1.5 0
../data/NH_Connecticut_River_2015/NH_Connecticut_River_2015.datalist -1 1.25 0
../data/NH_Umbagog/NH_Umbagog.datalist -1 1.5 0
```

Generate the inf and geojson auxiliary files
```bash
dlim -i Androscoggin.datalist
```

#### Generate Uncertainty Rasters

See [waffles](/docs/waffles.md) for more information on the syntax of the waffles command for generating DEMs.

```waffles -R regions_tile_set.shp -E 1 -M uncertainty:waffles_module=linear:accumulate=True -O androscoggin -p res=1m -P epsg:6348+5703+geoid:g2018 Androscoggin.datalist -X 0:25 -m -k```

In this command we use the previously defined regional tile set 'regions_tile_set.shp' to generate a 1 meter Uncertainty DEM. We use the waffles 'uncertainty' module specifying the interpolation method as 'linear'. The output will be an accumulated uncertainty raster with 1 meter cell-spacing with a horizontal projection of UTM Zone 19N and a vertical projection of NAVD88 (geoid2018). With the `-m` switch we also generate an auxiliary data mask raster.

![](/media/androscoggin_unc.png)

**Figure 1.** Total estimated uncertainty (Androscoggin) 

### Point Reyes

In this example we generate 2 1/9 arc-second DEMs and associated Uncertainty rasters (Figure 2 \ref{fig2}) off the coast of Point Reyes, California using only fetches dataset modules. Instead of creating a datalist, we can just list the fetches modules in the command-line. To create 2 1/4 degree tiles, we input 2 region parameters into the waffles command-line.

#### Generate the DEM and Uncertainty

see [waffles](/docs/waffles.md) for more information on the syntax of the waffles command.

```bash
waffles \
-R -123.5/-123.25/38/38.25 -R -123.5/-123.25/37.75/38 \
-E .111111111s \
-O nocal_test \
-P epsg:4326+5703 \
-X 0:15 \
-T outliers:stacks=True:multipass=4:accumulate=True:mode=unscaled:max_weight=1.5 \
-M cudem:pre_mode=mbgrid:landmask=True:polygonize=5:pre_count=3:pre_upper_limit=-.1:min_weight=.6 \
-p -w -k -m -c -u\
mar_grav,-106:bathy_only=True,.001,.85 \
charts,-200,.01,.75 \
hydronos:datatype=xyz,-202,.1,.1 \
multibeam:exclude_survey_id=CNTL14RR/RB1604/NA085,-201,.6,.1 \
ned,-215:coast_buffer=0.0001,.61,.65 \
ned1,-215:coast_buffer=0.0001,5,.35 \
ehydro,-203,.65,0 \
hydronos:datatype=bag,-202,2,0 \
CUDEM,-210,15,0 \
CoNED,-211,14,0
```

![\label{fig2}](/media/reyes_uncertainty.png)

**Figure 2.** Uncertainty products (Point Reyes) 