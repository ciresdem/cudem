# DEM and Uncertainty

While generating a DEM from independent scattered data sources we optionally calculate the uncertainty for each cell of the resulting DEM.
And auxiliary Uncertainty Raster can be generated along with a DEM by using the --uncertainty switch in the [waffles](/docs/waffles.md) command. The waffles module 'uncertainty' can also be used to only generate an uncertainty raster without generating a DEM.

Types of uncertainty:

- Source uncertainty (dataset-wide)
- Source uncertainty (per data value)
- Bathymetric depth uncertainty (IHO function of depth)
- Sub-pixel uncertainty (variance)
- Interpolation uncertainty (split-sample)
- vertical datum transformation uncertainty

The various uncertainty types are combined to report a Total Value Uncertainty (TVU) as supplemental raster product.

Uncertainty values are combined using the Root Sum Squared (RSS).

1. Source Uncertainty (dataset-wide)

The dataset-wide source uncertainty is a single value that will be applied to an entire dataset. This is specified by the user in the 4th column of the dataset entry. Typically this value is provided by the data collector or processor in the datasets metadata. This can sometimes be reported as an RMSE or data Accuracy. If the value is reported as being in the 95th percentile confidence level, first divide that value by 1.96 to obtain an uncertainty value suitable for combining with other uncertainty values.

2. Source Uncertainty (per data value)

Some datasets, such as NOS BAG or BlueTopo, specify an uncertainty value for each data value in the dataset. Other times a user may independently calculate the uncertainty for each data value in their dataset. These uncertainty data can be used to inform the final TVU by specifying the uncertainty data of the dataset as either a seperate product or integrated into the dataset, such as with raster data or xyz data, respsectively.

3. Bathymetric Depth Uncertainty

For bathymetric data, the IHO standard can be used to calculate the uncertainty of each data value as a function of it's depth.

4. Sub-pixel Uncertainty

The Sub-pixel uncertainty, or sub-pixel variance, is calculated by default when combing the various datasets together into a DEM. Whenever there is more that one data value contributing a resulting DEM data cell, the (optionally weighted) variance of the input data is calculated.

5. Interpolation Uncertainty

The interpolation uncertainty is calculated for all interpolated cells in a resulting DEM using a split-sample method.

6. Vertical Datum Transformation Uncertainty

Whenever data is vertically transformed while processing a DEM, the uncertainty of that transformation is accumulated into the final TVU.

## Examples

### Androscoggin

In this example we will generate tiled uncertainty rasters of the Androscoggin region using input elevation data from the USGS.

- Define regions

- Obtain the USGS elevation data

- Create a Datalist

See [dlim](/docs/dlim.md) for more information of datlist formatting.

```
"../data/ME_Maine_LiDAR_NRCS_14/ME_Maine_LiDAR_NRCS_14.datalist" -1 1 0
"../data/ME_SouthCoastal_2020_A20/ME_SouthCoastal_2020_A20.datalist" -1 2 0
"../data/ME_Western_2016/ME_Western_2016.datalist" -1:mask="tmp/USGS_Maine_LiDAR_Processing_Boundary_Buffered_100m.shp":invert_mask=True 1.5 0
"../data/NH_Connecticut_River_2015/NH_Connecticut_River_2015.datalist" -1 1.25 0
"../data/NH_Umbagog/NH_Umbagog.datalist" -1 1.5 0
```

- Generate Uncertainty Rasters

See [waffles](/docs/waffles.md) for more information on the syntax of the waffles command for generating DEMs.

```waffles -R regions_tile_set.shp -E 1 -M uncertainty:waffles_module=linear:accumulate=True -O androscoggin -p res=1m -P epsg:6348+5703+geoid:g2018 Androscoggin.datalist -X 0:25 -m -k```

In this command we use the previously defined regional tile set 'regions_tile_set.shp' to generate a 1 meter Uncertainty DEM. We use the waffles 'uncertainty' module specifying the interpolation method as 'linear'. The output will be an accumulated uncertainty raster with 1 meter cell-spacing with a horizontal projection of UTM Zone XX and a vertical projection of NAVD88 (geoid2018). With the `-m` switch we also generate an auxiliary data mask raster.

![](/media/androscoggin_unc.png)