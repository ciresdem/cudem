# DEM and Uncertainty

While generating a DEM from independent scattered data sources we optionally calculate the uncertainty for each cell of the resulting DEM.
An auxiliary Uncertainty Raster can be generated along with a DEM by using the `--uncertainty` switch in the [waffles](/docs/waffles.md) command. Additionally, the waffles module 'uncertainty' can be used to only generate an uncertainty raster without generating a DEM.

## Types of uncertainty:

- [Source Data Uncertainty](#source-data-uncertainty)
  - [Source dataset uncertainty](#source-dataset-uncertainty) (dataset-wide)
  - [Source value uncertainty](#source-value-uncertainty) (per data value)
  - [Bathymetric depth uncertainty](#bathymetric-depth-uncertainty) (IHO function of depth)
- [Gridding Uncertainty](#gridding-uncertainty)
  - [Sub-pixel uncertainty](sub-pixel-variance) (variance)
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
#### Sub-pixel variance

The Sub-pixel uncertainty is calculated by default when combining the various datasets together into a DEM. Whenever there is more that one data value contributing to a resulting DEM data cell, the (optionally weighted) variance of the input data is calculated.

Uncertainty is measured with a variance or its square root, which is a standard deviation. The standard deviation of a statistic is also (and more commonly) called a standard error.

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