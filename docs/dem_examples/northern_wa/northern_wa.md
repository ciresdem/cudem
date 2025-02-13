# Overview

DEM generation example of Northern Washington.

## Specifications

| Region | Tile-Size | Cell-size | Horz Projection | Vert Projection |
|---|---|---|---|---|
| -R-125/-124/47/48.5 | .25 degrees | 1/9 Arc-Second (~3m) | NAD83 | NAVD88 |


## Generate the region vectors

### Generate a full-region buffered vector

This is for data fetching, etc.

Add a slight buffer (.1 degree) to ensure coverage of fetched data.

```bash
regions -R -125/-124/47/48.5 -B .1
```

will output [region_n48x60_w125x10.shp](region_n48x60_w125x10.geojson)

### Generate the .25 degree tiles

these will be the extents of each DEM generated.

```regions -R -121/-114/31.5/35 -T .25```

will output [regions_tile_set.shp](regions_tile_set.geojson) with 24 .25 degree tiles.

### Edit the tile set to only include desired tiles (using a GIS)

[tiles_1_9.shp](tiles_1_9.geojson)

## Generate a Coastline vector (optional)

```bash
fetch_osm_coastline.py -R region_n48x60_w125x10.geojson wa_coast.shp
```

will output [wa_coast.shp](wa_coast.geojson)

## Fetch common datasets

### Bathymetry
#### HydroNOS

```bash
fetches -R tiles_1_9.shp hydronos
```

```
hydronos/
├── bag
│   ├── D00165_MB_VR_MSL_1of1.bag
│   ├── H12219_MB_1m_MLLW_1of2.bag
│   ├── H12219_MB_2m_MLLW_2of2.bag
│   ├── H12220_MB_2m_MLLW_1of3.bag
│   ├── H12220_MB_4m_MLLW_2of3.bag
│   ├── H12220_MB_8m_MLLW_3of3.bag
│   ├── H12221_MB_1m_MLLW_1of3.bag
│   ├── H12221_MB_2m_MLLW_2of3.bag
│   ├── H12221_MB_4m_MLLW_3of3.bag
│   ├── H12222_MB_2m_MLLW_1of3.bag
│   ├── H12222_MB_4m_MLLW_2of3.bag
│   ├── H12222_MB_8m_MLLW_3of3.bag
│   ├── H12223_MB_VR_MLLW_1of1.bag
│   ├── H13412_MB_VR_Ellipsoid.bag
│   ├── H13412_MB_VR_MLLW.bag
│   ├── W00262_MB_VR_MLLW_1of1.bag
│   ├── W00442_MB_128m_MLLW_5of5.bag
│   ├── W00442_MB_128m_MLLW_Combined.bag
│   ├── W00442_MB_16m_MLLW_2of5.bag
│   ├── W00442_MB_32m_MLLW_3of5.bag
│   ├── W00442_MB_64m_MLLW_4of5.bag
│   ├── W00442_MB_8m_MLLW_1of5.bag
│   ├── W00445_MB_128m_MLLW_4of4.bag
│   ├── W00445_MB_128m_MLLW_Combined.bag
│   ├── W00445_MB_16m_MLLW_1of4.bag
│   ├── W00445_MB_32m_MLLW_2of4.bag
│   └── W00445_MB_64m_MLLW_3of4.bag
└── geodas
    ├── H02096.xyz.gz
    ├── H02170.xyz.gz
    ├── H02869.xyz.gz
    ├── H04715.xyz.gz
    ├── H04716.xyz.gz
    ├── H04728.xyz.gz
    ├── H04729.xyz.gz
    ├── H04735.xyz.gz
    ├── H05068.xyz.gz
    ├── H05069.xyz.gz
    ├── H05070.xyz.gz
    ├── H05107.xyz.gz
    ├── H05108.xyz.gz
    ├── H05109.xyz.gz
    ├── H05110.xyz.gz
    ├── H05111.xyz.gz
    ├── H05114.xyz.gz
    ├── H05146.xyz.gz
    ├── H05147.xyz.gz
    ├── H05148.xyz.gz
    ├── H05155.xyz.gz
    ├── H05157.xyz.gz
    ├── H07036.xyz.gz
    ├── H07037.xyz.gz
    ├── H08241.xyz.gz
    ├── H08242.xyz.gz
    ├── H09413.xyz.gz
    ├── H09415.xyz.gz
    ├── H09416.xyz.gz
    ├── H09418.xyz.gz
    ├── H11083.xyz.gz
    ├── H11086.xyz.gz
    ├── H12219.xyz.gz
    ├── H12220.xyz.gz
    ├── H12221.xyz.gz
    ├── H12222.xyz.gz
    └── H12223.xyz.gz
```

#### Nautical Charts
```bash
fetches -R tiles_1_9.shp charts
```

```
charts
├── US1WC01M.zip
├── US2WC03M.zip
├── US3WA01M.zip
├── US3WA03M.zip
├── US4WA36M.zip
├── US5WA01M.zip
├── US5WA04M.zip
├── US5WA20M.zip
└── US5WA60M.zip
```

#### Multibeam
```bash
fetches -R tiles_1_9.shp multibeam
```

```
multibeam/
├── 2009_Amundsen
│   ├── 0003_20090704_231701_Amundsen.all.mb58.fbt
│   ├── 0004_20090705_001701_Amundsen.all.mb58.fbt
│   └── 0005_20090705_011700_Amundsen.all.mb58.fbt
├── AII8L22
│   └── cb87oct16.rc01.mb15
├── AT11L30
│   ├── at11l30.sb.20050720.raw.mb41.fbt
│   ├── at11l30.sb.20050721.raw.mb41.fbt
│   ├── at11l30.sb.20050724.raw.mb41.fbt
│   ├── at11l30.sb.20050725.raw.mb41.fbt
│   └── at11l30.sb.20050726.raw.mb41.fbt
├── AT11L31
│   └── at11l31.sb.20050903.raw.mb41.fbt
├── AT11L32
│   └── at11l32.sb.20050919.raw.mb41.fbt
├── EX0801
│   ├── 0038_20080923_020558_EX.mb58.fbt
│   ├── 0039_20080923_024242_EX.mb58.fbt
│   └── 0040_20080923_032201_EX.mb58.fbt
├── EX1503L2
│   ├── 0294_20150608_113938_EX1503L2_MB.gsf.mb121.fbt
│   ├── 0295_20150608_123939_EX1503L2_MB.gsf.mb121.fbt
│   ├── 0296_20150608_133459_EX1503L2_MB.gsf.mb121.fbt
│   ├── 0297_20150608_134046_EX1503L2_MB.gsf.mb121.fbt
│   ├── 0298_20150608_140137_EX1503L2_MB.gsf.mb121.fbt
│   ├── 0299_20150608_142335_EX1503L2_MB.gsf.mb121.fbt
│   ├── 0300_20150608_143016_EX1503L2_MB.gsf.mb121.fbt
│   ├── 0301_20150608_150922_EX1503L2_MB.gsf.mb121.fbt
│   ├── 0302_20150608_151731_EX1503L2_MB.gsf.mb121.fbt
│   ├── 0303_20150608_152235_EX1503L2_MB.gsf.mb121.fbt
│   ├── 0304_20150608_152939_EX1503L2_MB.gsf.mb121.fbt
│   ├── 0305_20150608_155700_EX1503L2_MB.gsf.mb121.fbt
│   ├── 0306_20150608_160628_EX1503L2_MB.gsf.mb121.fbt
│   ├── 0307_20150608_161104_EX1503L2_MB.gsf.mb121.fbt
│   ├── 0308_20150608_162112_EX1503L2_MB.gsf.mb121.fbt
│   ├── 0309_20150608_164552_EX1503L2_MB.gsf.mb121.fbt
│   ├── 0310_20150608_165730_EX1503L2_MB.gsf.mb121.fbt
│   ├── 0311_20150608_175730_EX1503L2_MB.gsf.mb121.fbt
│   └── 0312_20150608_185729_EX1503L2_MB.gsf.mb121.fbt
├── FK009A
│   ├── 0045_20130823_183858_FK009A_EM710.all.mb58.fbt
│   ├── 0046_20130823_193859_FK009A_EM710.all.mb58.fbt
│   ├── 0047_20130823_203859_FK009A_EM710.all.mb58.fbt
│   ├── 0048_20130823_213859_FK009A_EM710.all.mb58.fbt
│   ├── 0049_20130823_223859_FK009A_EM710.all.mb58.fbt
│   ├── 0199_20130830_174907_FK009A_EM710.all.mb58.fbt
│   ├── 0200_20130830_184907_FK009A_EM710.all.mb58.fbt
│   ├── 0201_20130830_194907_FK009A_EM710.all.mb58.fbt
│   ├── 0202_20130830_204907_FK009A_EM710.all.mb58.fbt
│   └── 0203_20130830_214907_FK009A_EM710.all.mb58.fbt
├── FK009B
│   ├── 0028_20130909_043709_FK009B_EM710.all.mb58.fbt
│   ├── 0029_20130909_053708_FK009B_EM710.all.mb58.fbt
│   ├── 0186_20130918_055629_FK009B_EM710.all.mb58.fbt
│   ├── 0187_20130918_065628_FK009B_EM710.all.mb58.fbt
│   ├── 0188_20130918_075630_FK009B_EM710.all.mb58.fbt
│   └── 0189_20130918_085629_FK009B_EM710.all.mb58.fbt
├── FK190529
│   ├── 0001_20190604_051632_FK190529_EM710.all.mb58.fbt
│   ├── 0002_20190604_055937_FK190529_EM710.all.mb58.fbt
│   ├── 0003_20190604_065933_FK190529_EM710.all.mb58.fbt
│   ├── 0004_20190604_070537_FK190529_EM710.all.mb58.fbt
│   ├── 0005_20190604_073207_FK190529_EM710.all.mb58.fbt
│   ├── 0006_20190604_074253_FK190529_EM710.all.mb58.fbt
│   ├── 0008_20190604_175824_FK190529_EM710.all.mb58.fbt
│   ├── 0009_20190604_185825_FK190529_EM710.all.mb58.fbt
│   ├── 0010_20190604_192000_FK190529_EM710.all.mb58.fbt
│   ├── 0011_20190604_192114_FK190529_EM710.all.mb58.fbt
│   ├── 0012_20190604_195414_FK190529_EM710.all.mb58.fbt
│   ├── 0013_20190604_200055_FK190529_EM710.all.mb58.fbt
│   ├── 0014_20190604_202658_FK190529_EM710.all.mb58.fbt
│   ├── 0015_20190604_204933_FK190529_EM710.all.mb58.fbt
│   ├── 0016_20190604_211110_FK190529_EM710.all.mb58.fbt
│   ├── 0017_20190604_211218_FK190529_EM710.all.mb58.fbt
│   ├── 0018_20190604_212145_FK190529_EM710.all.mb58.fbt
│   ├── 0019_20190604_215324_FK190529_EM710.all.mb58.fbt
│   ├── 0021_20190604_225721_FK190529_EM710.all.mb58.fbt
│   ├── 0022_20190604_235602_FK190529_EM710.all.mb58.fbt
│   ├── 0023_20190605_000041_FK190529_EM710.all.mb58.fbt
│   ├── 0024_20190605_000648_FK190529_EM710.all.mb58.fbt
│   ├── 0025_20190605_000649_FK190529_EM710.all.mb58.fbt
│   ├── 0027_20190605_010902_FK190529_EM710.all.mb58.fbt
│   ├── 0028_20190605_014608_FK190529_EM710.all.mb58.fbt
│   ├── 0032_20190606_182732_FK190529_EM710.all.mb58.fbt
│   └── 0033_20190607_040857_FK190529_EM710.all.mb58.fbt
├── HEALY02
│   ├── sb200104240000.mb41.fbt
│   ├── sb200104240100.mb41.fbt
│   └── sb200104240200.mb41.fbt
├── HLY0101
│   ├── sb20011140000.mb41.fbt
│   ├── sb20011140100.mb41.fbt
│   └── sb20011140200.mb41.fbt
├── HLY0401
│   ├── sb20041220100.mb41.fbt
│   └── sb20041220200.mb41.fbt
├── HLY05TC
│   ├── sb20051530400.mb41.fbt
│   └── sb20051530500.mb41.fbt
├── HLY06TA
│   ├── sb20051231000.mb41.fbt
│   ├── sb20051231100.mb41.fbt
│   ├── sb20051231200.mb41.fbt
│   ├── sb20051231300.mb41.fbt
│   ├── sb20051231400.mb41.fbt
│   ├── sb20051231500.mb41.fbt
│   ├── sb20051231600.mb41.fbt
│   ├── sb20051240600.mb41.fbt
│   ├── sb20051240700.mb41.fbt
│   ├── sb20051240800.mb41.fbt
│   ├── sb20051240900.mb41.fbt
│   ├── sb20051241000.mb41.fbt
│   ├── sb20051241100.mb41.fbt
│   ├── sb20051250400.mb41.fbt
│   ├── sb20051250500.mb41.fbt
│   ├── sb20051250600.mb41.fbt
│   ├── sb20051250700.mb41.fbt
│   ├── sb20051250800.mb41.fbt
│   ├── sb20051270600.mb41.fbt
│   ├── sb20051270700.mb41.fbt
│   ├── sb20051270800.mb41.fbt
│   ├── sb20051270900.mb41.fbt
│   └── sb20051271000.mb41.fbt
├── HLY06TB
│   ├── sb20060930300.mb41.fbt
│   ├── sb20060930400.mb41.fbt
│   ├── sb20060940800.mb41.fbt
│   ├── sb20060940900.mb41.fbt
│   ├── sb20060941000.mb41.fbt
│   ├── sb20060941800.mb41.fbt
│   ├── sb20060941900.mb41.fbt
│   ├── sb20060950800.mb41.fbt
│   ├── sb20060950900.mb41.fbt
│   ├── sb20060951000.mb41.fbt
│   ├── sb20060961200.mb41.fbt
│   ├── sb20060961300.mb41.fbt
│   ├── sb20060961400.mb41.fbt
│   ├── sb20060961600.mb41.fbt
│   ├── sb20060961700.mb41.fbt
│   ├── sb20060961800.mb41.fbt
│   └── sb20060961900.mb41.fbt
├── HLY06TD
│   ├── sb20061210600.mb41.fbt
│   ├── sb20061210700.mb41.fbt
│   └── sb20061210800.mb41.fbt
├── HLY06TG
│   ├── sb20061892100.mb41.fbt
│   ├── sb20061892200.mb41.fbt
│   └── sb20061892300.mb41.fbt
├── HLY06TI
│   ├── sb20062460500.mb41.fbt
│   ├── sb20062460600.mb41.fbt
│   └── sb20062460700.mb41.fbt
├── HLY06TJ
│   ├── sb20063000800.mb41.fbt
│   ├── sb20063000900.mb41.fbt
│   ├── sb20063001000.mb41.fbt
│   └── sb20063091600.mb41.fbt
├── HLY07TA
│   ├── sb20070600200.mb41.fbt
│   ├── sb20070600300.mb41.fbt
│   ├── sb20070600400.mb41.fbt
│   ├── sb20070600500.mb41.fbt
│   ├── sb20070601000.mb41.fbt
│   ├── sb20070601100.mb41.fbt
│   └── sb20070601200.mb41.fbt
├── HLY07TC
│   ├── sb20070940500.mb41.fbt
│   └── sb20070940600.mb41.fbt
├── HLY07TD
│   ├── sb20071760500.mb41.fbt
│   ├── sb20071760600.mb41.fbt
│   └── sb20071760700.mb41.fbt
├── HLY07TG
│   ├── sb20072190700.mb41.fbt
│   ├── sb20072190800.mb41.fbt
│   └── sb20072190900.mb41.fbt
├── HLY08TA
│   ├── sb20080352100.mb41.fbt
│   ├── sb20080352200.mb41.fbt
│   ├── sb20080352300.mb41.fbt
│   ├── sb20080360000.mb41.fbt
│   ├── sb20080361400.mb41.fbt
│   ├── sb20080361500.mb41.fbt
│   ├── sb20080361600.mb41.fbt
│   ├── sb20080381000.mb41.fbt
│   ├── sb20080381100.mb41.fbt
│   ├── sb20080381200.mb41.fbt
│   ├── sb20080381300.mb41.fbt
│   ├── sb20080381400.mb41.fbt
│   ├── sb20080381500.mb41.fbt
│   └── sb20080381600.mb41.fbt
├── HLY08TG
│   ├── sb20081781100.mb41.fbt
│   ├── sb20081781200.mb41.fbt
│   └── sb20081781300.mb41.fbt
├── HLY09TA
│   ├── sb20090291000.mb41.fbt
│   ├── sb20090291100.mb41.fbt
│   ├── sb20090300600.mb41.fbt
│   ├── sb20090300700.mb41.fbt
│   ├── sb20090300800.mb41.fbt
│   ├── sb20090311700.mb41.fbt
│   ├── sb20090311800.mb41.fbt
│   ├── sb20090320800.mb41.fbt
│   ├── sb20090320900.mb41.fbt
│   ├── sb20090321000.mb41.fbt
│   ├── sb20090321900.mb41.fbt
│   ├── sb20090322000.mb41.fbt
│   ├── sb20090322100.mb41.fbt
│   ├── sb20090322200.mb41.fbt
│   ├── sb20090341000.mb41.fbt
│   ├── sb20090341100.mb41.fbt
│   ├── sb20090370800.mb41.fbt
│   ├── sb20090370900.mb41.fbt
│   ├── sb20090371000.mb41.fbt
│   ├── sb20090371100.mb41.fbt
│   └── sb20090371200.mb41.fbt
├── HLY10TA
│   ├── 0071_20100410_110640_Healy.all.mb58.fbt
│   ├── 0072_20100410_113640_Healy.all.mb58.fbt
│   ├── 0073_20100410_120641_Healy.all.mb58.fbt
│   ├── 0074_20100410_123641_Healy.all.mb58.fbt
│   ├── 0075_20100410_130640_Healy.all.mb58.fbt
│   ├── 0076_20100410_133640_Healy.all.mb58.fbt
│   ├── 0101_20100413_061923_Healy.all.mb58.fbt
│   ├── 0102_20100413_064923_Healy.all.mb58.fbt
│   ├── 0103_20100413_071923_Healy.all.mb58.fbt
│   └── 0104_20100413_074923_Healy.all.mb58.fbt
├── HLY10TB
│   ├── 0025_20100517_105617_ShipName.all.mb58.fbt
│   ├── 0026_20100517_112618_ShipName.all.mb58.fbt
│   ├── 0027_20100517_115617_ShipName.all.mb58.fbt
│   ├── 0028_20100517_122617_ShipName.all.mb58.fbt
│   └── 0029_20100517_125617_ShipName.all.mb58.fbt
├── HLY11TA
│   ├── 0015_20110426_034345_ShipName.all.mb58.fbt
│   ├── 0016_20110426_041346_ShipName.all.mb58.fbt
│   ├── 0017_20110426_044346_ShipName.all.mb58.fbt
│   ├── 0018_20110426_051345_ShipName.all.mb58.fbt
│   ├── 0019_20110426_054345_ShipName.all.mb58.fbt
│   ├── 0057_20110506_123628_ShipName.all.mb58.fbt
│   ├── 0058_20110506_130628_ShipName.all.mb58.fbt
│   ├── 0059_20110506_133629_ShipName.all.mb58.fbt
│   ├── 0060_20110506_140628_ShipName.all.mb58.fbt
│   ├── 0061_20110506_143629_ShipName.all.mb58.fbt
│   └── 0062_20110506_150628_ShipName.all.mb58.fbt
├── HLY11TB
│   ├── 0019_20110528_071409_Healy.all.mb58.fbt
│   ├── 0020_20110528_074409_Healy.all.mb58.fbt
│   ├── 0021_20110528_081410_Healy.all.mb58.fbt
│   ├── 0022_20110528_084409_Healy.all.mb58.fbt
│   └── 0023_20110528_091409_Healy.all.mb58.fbt
├── HLY12TA
│   ├── 0014_20120605_121450_Healy.all.mb58.fbt
│   ├── 0015_20120605_131450_Healy.all.mb58.fbt
│   ├── 0016_20120605_141449_Healy.all.mb58.fbt
│   ├── 0030_20120606_100052_Healy.all.mb58.fbt
│   ├── 0031_20120606_110051_Healy.all.mb58.fbt
│   ├── 0032_20120606_120051_Healy.all.mb58.fbt
│   └── 0033_20120606_130051_Healy.all.mb58.fbt
├── HLY12TB
│   ├── 0011_20120731_075618_Healy.all.mb58.fbt
│   ├── 0012_20120731_085618_Healy.all.mb58.fbt
│   └── 0013_20120731_095618_Healy.all.mb58.fbt
├── HLY14TF
│   ├── 0049_20140916_173214_HEALY.all.mb58.fbt
│   ├── 0050_20140916_180214_HEALY.all.mb58.fbt
│   ├── 0051_20140916_183215_HEALY.all.mb58.fbt
│   ├── 0052_20140916_190214_HEALY.all.mb58.fbt
│   ├── 0192_20140919_170214_HEALY.all.mb58.fbt
│   ├── 0193_20140919_173214_HEALY.all.mb58.fbt
│   ├── 0194_20140919_180214_HEALY.all.mb58.fbt
│   └── 0195_20140919_183214_HEALY.all.mb58.fbt
├── HLY15TD
│   ├── 0305_20151028_221231_HEALY.all.mb58.fbt
│   ├── 0306_20151028_224231_HEALY.all.mb58.fbt
│   ├── 0307_20151028_231231_HEALY.all.mb58.fbt
│   ├── 0308_20151028_234231_HEALY.all.mb58.fbt
│   └── 0309_20151029_001231_HEALY.all.mb58.fbt
├── HLY16TB
│   ├── 0006_20160608_035927_HEALY.all.mb58.fbt
│   ├── 0007_20160608_042926_HEALY.all.mb58.fbt
│   ├── 0008_20160608_045926_HEALY.all.mb58.fbt
│   └── 0009_20160608_052926_HEALY.all.mb58.fbt
├── HLY17TB
│   ├── 0011_20170628_020223_Healy.all.mb58.fbt
│   ├── 0012_20170628_030223_Healy.all.mb58.fbt
│   └── 0013_20170628_040224_Healy.all.mb58.fbt
├── HLY18TB
│   ├── 0013_20180619_070553_Healy.all.mb58.fbt
│   ├── 0014_20180619_080553_Healy.all.mb58.fbt
│   ├── 0015_20180619_090552_Healy.all.mb58.fbt
│   ├── 0229_20180630_132647_Healy.all.mb58.fbt
│   ├── 0230_20180630_142647_Healy.all.mb58.fbt
│   ├── 0231_20180630_152647_Healy.all.mb58.fbt
│   ├── 0232_20180630_162647_Healy.all.mb58.fbt
│   └── 0233_20180630_172647_Healy.all.mb58.fbt
├── HLY18TC
│   ├── 0006_20180725_053513_Healy.all.mb58.fbt
│   ├── 0007_20180725_063514_Healy.all.mb58.fbt
│   └── 0008_20180725_073513_Healy.all.mb58.fbt
├── HLY20TB
│   ├── 0012_20200706_122438_Healy.all.mb58.fbt
│   ├── 0013_20200706_125438_Healy.all.mb58.fbt
│   ├── 0014_20200706_132438_Healy.all.mb58.fbt
│   ├── 0015_20200706_135438_Healy.all.mb58.fbt
│   ├── 0016_20200706_142438_Healy.all.mb58.fbt
│   ├── 0212_20200710_162449_Healy.all.mb58.fbt
│   ├── 0213_20200710_165450_Healy.all.mb58.fbt
│   └── 9999.all.mb58.fbt
├── HLY20TC
│   ├── 0011_20200712_043652_Healy.all.mb58.fbt
│   ├── 0012_20200712_050652_Healy.all.mb58.fbt
│   ├── 0013_20200712_053652_Healy.all.mb58.fbt
│   ├── 0014_20200712_060652_Healy.all.mb58.fbt
│   └── 0015_20200712_063652_Healy.all.mb58.fbt
├── HLY21TA
│   ├── 0032_20210526_070357_Healy.all.mb58.fbt
│   ├── 0033_20210526_073357_Healy.all.mb58.fbt
│   ├── 0034_20210526_080358_Healy.all.mb58.fbt
│   ├── 0035_20210526_083358_Healy.all.mb58.fbt
│   ├── 0036_20210526_090358_Healy.all.mb58.fbt
│   ├── 0037_20210526_093357_Healy.all.mb58.fbt
│   └── 0038_20210526_100357_Healy.all.mb58.fbt
├── HLY21TG
│   ├── 0066_20211120_041220_Healy.all.mb58.fbt
│   ├── 0067_20211120_051220_Healy.all.mb58.fbt
│   └── 0068_20211120_061220_Healy.all.mb58.fbt
├── HMPR-106-2001-02
│   ├── line059.mb94.fbt
│   └── line060.mb94.fbt
├── KM0311
│   ├── em1002-174-025050-0001.mb56.fbt
│   ├── em1002-174-025427-0001.mb56.fbt
│   ├── em1002-174-030927-0002.mb56.fbt
│   ├── em1002-174-032427-0003.mb56.fbt
│   ├── em1002-174-033927-0004.mb56.fbt
│   └── em1002-174-035427-0005.mb56.fbt
├── KM1314
│   ├── 0009_20130808_000954_KM710.all.mb58.fbt
│   ├── 0010_20130808_010955_KM710.all.mb58.fbt
│   └── 0011_20130808_020955_KM710.all.mb58.fbt
├── MGL1212
│   └── 0119_20120721_143016_Langseth.all.mb58.fbt
├── MGL2102
│   ├── 0020_20210512_140203_Langseth.all.mb58.fbt
│   ├── 0021_20210512_143203_Langseth.all.mb58.fbt
│   ├── 0022_20210512_150204_Langseth.all.mb58.fbt
│   ├── 0023_20210512_153203_Langseth.all.mb58.fbt
│   ├── 0024_20210512_160202_Langseth.all.mb58.fbt
│   └── 0025_20210512_163203_Langseth.all.mb58.fbt
├── MGL2106
│   ├── 1812_20210903_044858_Langseth.all.mb58.fbt
│   ├── 1813_20210903_051858_Langseth.all.mb58.fbt
│   ├── 1814_20210903_054858_Langseth.all.mb58.fbt
│   ├── 1815_20210903_061858_Langseth.all.mb58.fbt
│   ├── 1816_20210903_064858_Langseth.all.mb58.fbt
│   └── 1817_20210903_071858_Langseth.all.mb58.fbt
├── MGL2107
│   ├── 0009_20211015_070643_Langseth.all.mb58.fbt
│   ├── 0010_20211015_073644_Langseth.all.mb58.fbt
│   ├── 0011_20211015_080643_Langseth.all.mb58.fbt
│   ├── 0012_20211015_083644_Langseth.all.mb58.fbt
│   ├── 0013_20211015_090644_Langseth.all.mb58.fbt
│   └── 0014_20211015_093643_Langseth.all.mb58.fbt
├── MGL2209
│   ├── 0143_20220908_022433_Langseth.all.mb58.fbt
│   ├── 0144_20220908_025433_Langseth.all.mb58.fbt
│   ├── 0145_20220908_032432_Langseth.all.mb58.fbt
│   ├── 0146_20220908_035432_Langseth.all.mb58.fbt
│   ├── 0147_20220908_042433_Langseth.all.mb58.fbt
│   ├── 0148_20220908_045433_Langseth.all.mb58.fbt
│   ├── 0149_20220908_052433_Langseth.all.mb58.fbt
│   └── 0150_20220908_055433_Langseth.all.mb58.fbt
├── MV1310
│   ├── 0010_20130805_020657_melville.all.mb58.fbt
│   ├── 0011_20130805_023657_melville.all.mb58.fbt
│   ├── 0012_20130805_030658_melville.all.mb58.fbt
│   └── 0013_20130805_033657_melville.all.mb58.fbt
├── MV1404
│   ├── 0001_20140611_022718_melville.all.mb58.fbt
│   ├── 0002_20140611_025718_melville.all.mb58.fbt
│   ├── 0003_20140611_032718_melville.all.mb58.fbt
│   ├── 0004_20140611_035717_melville.all.mb58.fbt
│   └── 0005_20140611_042717_melville.all.mb58.fbt
├── NA070
│   ├── 0001_20160414_015249_Nautilus.all.mb58.fbt
│   ├── 0002_20160414_025249_Nautilus.all.mb58.fbt
│   ├── 0003_20160414_035249_Nautilus.all.mb58.fbt
│   ├── 0004_20160414_045249_Nautilus.all.mb58.fbt
│   ├── 0005_20160411_053909_Nautilus.all.mb58.fbt
│   ├── 0005_20160412_195355_Nautilus.all.mb58.fbt
│   ├── 0006_20160411_063909_Nautilus.all.mb58.fbt
│   ├── 0006_20160412_205354_Nautilus.all.mb58.fbt
│   ├── 0007_20160411_073909_Nautilus.all.mb58.fbt
│   ├── 0007_20160412_215355_Nautilus.all.mb58.fbt
│   ├── 0008_20160412_225355_Nautilus.all.mb58.fbt
│   ├── 0009_20160415_083802_Nautilus.all.mb58.fbt
│   ├── 0010_20160415_093802_Nautilus.all.mb58.fbt
│   ├── 0011_20160415_103802_Nautilus.all.mb58.fbt
│   ├── 0012_20160415_113802_Nautilus.all.mb58.fbt
│   └── 0013_20160415_123802_Nautilus.all.mb58.fbt
├── NA072
│   ├── 0002_20160601_200321_Nautilus_EM302.gsf.mb121.fbt
│   ├── 0003_20160601_210321_Nautilus_EM302.gsf.mb121.fbt
│   └── 0004_20160601_220321_Nautilus_EM302.gsf.mb121.fbt
├── NA096
│   ├── 0000_20180702_233952_Nautilus.all.mb58.fbt
│   ├── 0001_20180702_041426_Nautilus.all.mb58.fbt
│   ├── 0001_20180703_003952_Nautilus.all.mb58.fbt
│   ├── 0002_20180703_013952_Nautilus.all.mb58.fbt
│   ├── 0004_20180702_061819_Nautilus.all.mb58.fbt
│   ├── 0005_20180702_063827_Nautilus.all.mb58.fbt
│   ├── 0006_20180702_064457_Nautilus.all.mb58.fbt
│   ├── 0008_20180702_070221_Nautilus.all.mb58.fbt
│   ├── 0009_20180702_071618_Nautilus.all.mb58.fbt
│   ├── 0010_20180702_072154_Nautilus.all.mb58.fbt
│   ├── 0012_20180702_073958_Nautilus.all.mb58.fbt
│   ├── 0013_20180702_075636_Nautilus.all.mb58.fbt
│   ├── 0014_20180702_080157_Nautilus.all.mb58.fbt
│   ├── 0016_20180702_082124_Nautilus.all.mb58.fbt
│   ├── 0017_20180702_083807_Nautilus.all.mb58.fbt
│   ├── 0018_20180702_084417_Nautilus.all.mb58.fbt
│   ├── 0020_20180702_090304_Nautilus.all.mb58.fbt
│   ├── 0021_20180702_092004_Nautilus.all.mb58.fbt
│   ├── 0022_20180702_092848_Nautilus.all.mb58.fbt
│   ├── 0024_20180702_094744_Nautilus.all.mb58.fbt
│   ├── 0025_20180702_100807_Nautilus.all.mb58.fbt
│   ├── 0026_20180702_101309_Nautilus.all.mb58.fbt
│   ├── 0028_20180702_103414_Nautilus.all.mb58.fbt
│   ├── 0029_20180702_105544_Nautilus.all.mb58.fbt
│   ├── 0030_20180702_110125_Nautilus.all.mb58.fbt
│   ├── 0031_20180702_111527_Nautilus.all.mb58.fbt
│   ├── 0032_20180702_112145_Nautilus.all.mb58.fbt
│   ├── 0033_20180702_114406_Nautilus.all.mb58.fbt
│   ├── 0034_20180702_115012_Nautilus.all.mb58.fbt
│   ├── 0035_20180702_120355_Nautilus.all.mb58.fbt
│   ├── 0036_20180702_121034_Nautilus.all.mb58.fbt
│   ├── 0037_20180702_123252_Nautilus.all.mb58.fbt
│   ├── 0038_20180702_123643_Nautilus.all.mb58.fbt
│   ├── 0040_20180702_125518_Nautilus.all.mb58.fbt
│   ├── 0041_20180702_131652_Nautilus.all.mb58.fbt
│   ├── 0042_20180702_132041_Nautilus.all.mb58.fbt
│   ├── 0045_20180702_135223_Nautilus.all.mb58.fbt
│   ├── 0046_20180702_142000_Nautilus.all.mb58.fbt
│   └── 9999.all.mb58.fbt
├── NA121
│   ├── 0005_20200923_081037_Nautilus.all.mb58.fbt
│   └── 0006_20200923_091038_Nautilus.all.mb58.fbt
├── REM-01MV
│   └── SBfixavg.93aug30
├── RR1303
│   ├── 0007_20130402_204700_revelle.all.mb58.fbt
│   ├── 0008_20130402_214700_revelle.all.mb58.fbt
│   ├── 0009_20130402_224701_revelle.all.mb58.fbt
│   └── 0010_20130402_234700_revelle.all.mb58.fbt
├── SKQ201607S
│   └── 0064_20160610_004937_Sikuliaq.all.mb58.fbt
├── SKQ201608S
│   ├── 0003_20160613_180640_Sikuliaq.all.mb58.fbt
│   ├── 0004_20160613_200640_Sikuliaq.all.mb58.fbt
│   └── 0144_20160625_090046_Sikuliaq.all.mb58.fbt
├── SKQ201609S
│   ├── 0003_20160629_214403_Sikuliaq.all.mb58.fbt
│   ├── 0004_20160629_234404_Sikuliaq.all.mb58.fbt
│   └── 0005_20160630_014404_Sikuliaq.all.mb58.fbt
├── SKQ201610S
│   ├── 0004_20160711_102905_Sikuliaq.all.mb58.fbt
│   ├── 0005_20160711_122905_Sikuliaq.all.mb58.fbt
│   ├── 0299_20160812_014113_Sikuliaq.all.mb58.fbt
│   ├── 0300_20160812_034113_Sikuliaq.all.mb58.fbt
│   └── 0301_20160812_054114_Sikuliaq.all.mb58.fbt
├── SKQ201704S
│   ├── 0000_20170411_141422_Sikuliaq.all.mb58.fbt
│   ├── 0001_20170411_161422_Sikuliaq.all.mb58.fbt
│   ├── 0009_20170413_154253_Sikuliaq.all.mb58.fbt
│   ├── 0010_20170413_174254_Sikuliaq.all.mb58.fbt
│   └── 0011_20170413_194254_Sikuliaq.all.mb58.fbt
├── SKQ201705S
│   ├── 0037_20170429_221204_Sikuliaq.all.mb58.fbt
│   ├── 0038_20170430_001203_Sikuliaq.all.mb58.fbt
│   └── 0039_20170430_023424_Sikuliaq.all.mb58.fbt
├── SKQ201715S
│   ├── 0000_20171006_032413_Sikuliaq.all.mb58.fbt
│   ├── 0001_20171006_034822_Sikuliaq.all.mb58.fbt
│   ├── 0008_20171003_101646_Sikuliaq.all.mb58.fbt
│   ├── 0030_20171005_063935_Sikuliaq.all.mb58.fbt
│   ├── 0031_20171005_083935_Sikuliaq.all.mb58.fbt
│   ├── 0032_20171005_103935_Sikuliaq.all.mb58.fbt
│   ├── 0033_20171005_123935_Sikuliaq.all.mb58.fbt
│   ├── 0034_20171005_143935_Sikuliaq.all.mb58.fbt
│   ├── 0035_20171005_163935_Sikuliaq.all.mb58.fbt
│   ├── 0036_20171005_183935_Sikuliaq.all.mb58.fbt
│   ├── 0037_20171005_203935_Sikuliaq.all.mb58.fbt
│   ├── 0038_20171005_223935_Sikuliaq.all.mb58.fbt
│   ├── 0039_20171006_003935_Sikuliaq.all.mb58.fbt
│   ├── 0040_20171006_023935_Sikuliaq.all.mb58.fbt
│   ├── 0046_20171006_143936_Sikuliaq.all.mb58.fbt
│   ├── 0047_20171006_163936_Sikuliaq.all.mb58.fbt
│   ├── 0048_20171006_183936_Sikuliaq.all.mb58.fbt
│   └── 0049_20171006_203936_Sikuliaq.all.mb58.fbt
├── SKQ201808S
│   ├── 0012_20180325_160323_Sikuliaq.all.mb58.fbt
│   ├── 0013_20180325_170323_Sikuliaq.all.mb58.fbt
│   ├── 0034_20180327_192835_Sikuliaq.all.mb58.fbt
│   ├── 0035_20180327_202836_Sikuliaq.all.mb58.fbt
│   ├── 0036_20180327_212835_Sikuliaq.all.mb58.fbt
│   ├── 0037_20180327_222835_Sikuliaq.all.mb58.fbt
│   ├── 0039_20180330_031342_Sikuliaq.all.mb58.fbt
│   ├── 0040_20180330_041342_Sikuliaq.all.mb58.fbt
│   ├── 0041_20180330_051343_Sikuliaq.all.mb58.fbt
│   ├── 0042_20180330_061342_Sikuliaq.all.mb58.fbt
│   ├── 0043_20180330_071343_Sikuliaq.all.mb58.fbt
│   ├── 0044_20180330_081342_Sikuliaq.all.mb58.fbt
│   ├── 0045_20180330_091342_Sikuliaq.all.mb58.fbt
│   ├── 0046_20180330_101343_Sikuliaq.all.mb58.fbt
│   ├── 0049_20180330_131343_Sikuliaq.all.mb58.fbt
│   └── 0050_20180330_141343_Sikuliaq.all.mb58.fbt
├── SKQ201822S
│   ├── 0163_20181205_140656_Sikuliaq.all.mb58.fbt
│   ├── 0164_20181205_150656_Sikuliaq.all.mb58.fbt
│   ├── 0165_20181205_160656_Sikuliaq.all.mb58.fbt
│   └── 0166_20181205_170656_Sikuliaq.all.mb58.fbt
├── SKQ201910S
│   ├── 0036_20190424_130238_Sikuliaq.all.mb58.fbt
│   ├── 0037_20190424_140238_Sikuliaq.all.mb58.fbt
│   ├── 0052_20190425_093537_Sikuliaq.all.mb58.fbt
│   ├── 0053_20190425_103537_Sikuliaq.all.mb58.fbt
│   ├── 0054_20190425_113537_Sikuliaq.all.mb58.fbt
│   ├── 0055_20190425_123537_Sikuliaq.all.mb58.fbt
│   ├── 0056_20190425_133537_Sikuliaq.all.mb58.fbt
│   ├── 0057_20190425_143537_Sikuliaq.all.mb58.fbt
│   ├── 0058_20190425_233717_Sikuliaq.all.mb58.fbt
│   └── 0059_20190426_003718_Sikuliaq.all.mb58.fbt
├── SKQ201921S
│   ├── 0028_20191012_035727_Sikuliaq.all.mb58.fbt
│   ├── 0029_20191012_045728_Sikuliaq.all.mb58.fbt
│   ├── 0032_20191012_075728_Sikuliaq.all.mb58.fbt
│   ├── 0033_20191012_085728_Sikuliaq.all.mb58.fbt
│   ├── 0034_20191012_095728_Sikuliaq.all.mb58.fbt
│   ├── 0035_20191012_105728_Sikuliaq.all.mb58.fbt
│   ├── 0036_20191012_115728_Sikuliaq.all.mb58.fbt
│   ├── 0037_20191012_125728_Sikuliaq.all.mb58.fbt
│   ├── 0038_20191012_215744_Sikuliaq.all.mb58.fbt
│   └── 0039_20191012_225744_Sikuliaq.all.mb58.fbt
├── SKQ201924S
│   ├── 0410_20191222_215322_Sikuliaq.all.mb58.fbt
│   ├── 0411_20191222_225322_Sikuliaq.all.mb58.fbt
│   ├── 0412_20191222_235322_Sikuliaq.all.mb58.fbt
│   ├── 0413_20191223_005323_Sikuliaq.all.mb58.fbt
│   └── 0414_20191223_015322_Sikuliaq.all.mb58.fbt
├── SKQ202103S
│   ├── 0008_20210313_022214_Sikuliaq.all.mb58.fbt
│   ├── 0009_20210313_032214_Sikuliaq.all.mb58.fbt
│   └── 0010_20210313_042214_Sikuliaq.all.mb58.fbt
├── SKQ202104S
│   ├── 0030_20210325_225740_Sikuliaq.all.mb58.fbt
│   ├── 0031_20210325_235740_Sikuliaq.all.mb58.fbt
│   ├── 0032_20210326_005740_Sikuliaq.all.mb58.fbt
│   ├── 0033_20210326_015740_Sikuliaq.all.mb58.fbt
│   ├── 0122_20210401_141259_Sikuliaq.all.mb58.fbt
│   ├── 0123_20210401_151251_Sikuliaq.all.mb58.fbt
│   ├── 0124_20210401_161251_Sikuliaq.all.mb58.fbt
│   ├── 0125_20210401_171251_Sikuliaq.all.mb58.fbt
│   ├── 0126_20210401_181251_Sikuliaq.all.mb58.fbt
│   ├── 0127_20210401_191251_Sikuliaq.all.mb58.fbt
│   ├── 0128_20210401_201251_Sikuliaq.all.mb58.fbt
│   ├── 0129_20210401_211251_Sikuliaq.all.mb58.fbt
│   ├── 0130_20210401_221251_Sikuliaq.all.mb58.fbt
│   ├── 0131_20210401_231251_Sikuliaq.all.mb58.fbt
│   ├── 0132_20210402_001252_Sikuliaq.all.mb58.fbt
│   ├── 0133_20210402_011252_Sikuliaq.all.mb58.fbt
│   ├── 0134_20210402_021252_Sikuliaq.all.mb58.fbt
│   ├── 0135_20210402_031252_Sikuliaq.all.mb58.fbt
│   └── 0136_20210402_041252_Sikuliaq.all.mb58.fbt
├── SKQ202205S
│   ├── 0017_20220322_164711_Sikuliaq.all.mb58.fbt
│   ├── 0018_20220322_174711_Sikuliaq.all.mb58.fbt
│   ├── 0019_20220322_184711_Sikuliaq.all.mb58.fbt
│   ├── 0020_20220322_194711_Sikuliaq.all.mb58.fbt
│   ├── 0021_20220322_204711_Sikuliaq.all.mb58.fbt
│   ├── 0055_20220324_110041_Sikuliaq.all.mb58.fbt
│   ├── 0056_20220324_120050_Sikuliaq.all.mb58.fbt
│   ├── 0057_20220324_130041_Sikuliaq.all.mb58.fbt
│   ├── 0060_20220324_160041_Sikuliaq.all.mb58.fbt
│   ├── 0061_20220324_170041_Sikuliaq.all.mb58.fbt
│   ├── 0062_20220324_180042_Sikuliaq.all.mb58.fbt
│   ├── 0063_20220324_190042_Sikuliaq.all.mb58.fbt
│   ├── 0064_20220324_200042_Sikuliaq.all.mb58.fbt
│   ├── 0065_20220324_210042_Sikuliaq.all.mb58.fbt
│   ├── 0066_20220324_220042_Sikuliaq.all.mb58.fbt
│   ├── 0067_20220324_230042_Sikuliaq.all.mb58.fbt
│   ├── 0068_20220325_000042_Sikuliaq.all.mb58.fbt
│   ├── 0069_20220325_010042_Sikuliaq.all.mb58.fbt
│   ├── 0070_20220325_020042_Sikuliaq.all.mb58.fbt
│   ├── 0071_20220325_030042_Sikuliaq.all.mb58.fbt
│   ├── 0072_20220325_040042_Sikuliaq.all.mb58.fbt
│   ├── 0073_20220325_050042_Sikuliaq.all.mb58.fbt
│   ├── 0074_20220325_060042_Sikuliaq.all.mb58.fbt
│   ├── 0075_20220325_081108_Sikuliaq.all.mb58.fbt
│   ├── 0076_20220325_091109_Sikuliaq.all.mb58.fbt
│   ├── 0077_20220325_101109_Sikuliaq.all.mb58.fbt
│   ├── 0078_20220325_111109_Sikuliaq.all.mb58.fbt
│   ├── 0079_20220325_121109_Sikuliaq.all.mb58.fbt
│   ├── 0080_20220325_131109_Sikuliaq.all.mb58.fbt
│   ├── 0081_20220325_141109_Sikuliaq.all.mb58.fbt
│   ├── 0082_20220325_151109_Sikuliaq.all.mb58.fbt
│   ├── 0083_20220325_161109_Sikuliaq.all.mb58.fbt
│   ├── 0084_20220325_171109_Sikuliaq.all.mb58.fbt
│   ├── 0092_20220326_011109_Sikuliaq.all.mb58.fbt
│   ├── 0093_20220326_021109_Sikuliaq.all.mb58.fbt
│   ├── 0094_20220326_031109_Sikuliaq.all.mb58.fbt
│   ├── 0095_20220326_041109_Sikuliaq.all.mb58.fbt
│   ├── 0136_20220327_210332_Sikuliaq.all.mb58.fbt
│   └── 0137_20220327_220332_Sikuliaq.all.mb58.fbt
├── SKQ202208S
│   ├── 0198_20220529_145312_Sikuliaq.all.mb58.fbt
│   ├── 0199_20220529_155312_Sikuliaq.all.mb58.fbt
│   ├── 0200_20220529_165312_Sikuliaq.all.mb58.fbt
│   └── 0201_20220529_175312_Sikuliaq.all.mb58.fbt
├── SO108
│   ├── O000.mb21.fbt
│   ├── P001.mb21.fbt
│   ├── P002.mb21.fbt
│   ├── P003.mb21.fbt
│   ├── P004.mb21.fbt
│   ├── R000.mb21.fbt
│   ├── R001.mb21.fbt
│   ├── S003.mb21.fbt
│   └── S004.mb21.fbt
├── SR1813
│   ├── 0019_20180918_102324_SallyRide.all.mb58.fbt
│   ├── 0020_20180918_105324_SallyRide.all.mb58.fbt
│   ├── 0021_20180918_112324_SallyRide.all.mb58.fbt
│   ├── 0022_20180918_115323_SallyRide.all.mb58.fbt
│   ├── 0023_20180918_122323_SallyRide.all.mb58.fbt
│   ├── 0024_20180918_125323_SallyRide.all.mb58.fbt
│   ├── 0025_20180918_132323_SallyRide.all.mb58.fbt
│   ├── 0033_20180918_172346_SallyRide.all.mb58.fbt
│   ├── 0034_20180918_194112_SallyRide.all.mb58.fbt
│   ├── 0074_20180920_021753_SallyRide.all.mb58.fbt
│   ├── 0075_20180920_024753_SallyRide.all.mb58.fbt
│   ├── 0077_20180920_034752_SallyRide.all.mb58.fbt
│   ├── 0078_20180920_041752_SallyRide.all.mb58.fbt
│   ├── 0079_20180920_044752_SallyRide.all.mb58.fbt
│   ├── 0080_20180920_051752_SallyRide.all.mb58.fbt
│   ├── 0081_20180920_054752_SallyRide.all.mb58.fbt
│   ├── 0082_20180920_061752_SallyRide.all.mb58.fbt
│   ├── 0083_20180920_064751_SallyRide.all.mb58.fbt
│   ├── 0084_20180920_071751_SallyRide.all.mb58.fbt
│   ├── 0085_20180920_074751_SallyRide.all.mb58.fbt
│   ├── 0086_20180920_081751_SallyRide.all.mb58.fbt
│   ├── 0087_20180920_084751_SallyRide.all.mb58.fbt
│   ├── 0088_20180920_091751_SallyRide.all.mb58.fbt
│   ├── 0089_20180920_094751_SallyRide.all.mb58.fbt
│   ├── 0090_20180920_101751_SallyRide.all.mb58.fbt
│   ├── 0091_20180920_104750_SallyRide.all.mb58.fbt
│   ├── 0092_20180920_111750_SallyRide.all.mb58.fbt
│   ├── 0093_20180920_114750_SallyRide.all.mb58.fbt
│   ├── 0094_20180920_121750_SallyRide.all.mb58.fbt
│   ├── 0095_20180920_124750_SallyRide.all.mb58.fbt
│   ├── 0096_20180920_131750_SallyRide.all.mb58.fbt
│   ├── 0097_20180920_134749_SallyRide.all.mb58.fbt
│   ├── 0098_20180920_141749_SallyRide.all.mb58.fbt
│   ├── 0099_20180920_182908_SallyRide.all.mb58.fbt
│   ├── 0100_20180920_183637_SallyRide.all.mb58.fbt
│   ├── 0101_20180920_190637_SallyRide.all.mb58.fbt
│   ├── 0102_20180920_193637_SallyRide.all.mb58.fbt
│   ├── 0111_20180921_000635_SallyRide.all.mb58.fbt
│   ├── 0112_20180921_003635_SallyRide.all.mb58.fbt
│   ├── 0193_20180924_193436_SallyRide.all.mb58.fbt
│   ├── 0194_20180924_200436_SallyRide.all.mb58.fbt
│   ├── 0195_20180924_203436_SallyRide.all.mb58.fbt
│   ├── 0196_20180924_210436_SallyRide.all.mb58.fbt
│   ├── 0197_20180924_213436_SallyRide.all.mb58.fbt
│   └── 0198_20180924_220435_SallyRide.all.mb58.fbt
├── TN144
│   ├── 0001_20020425_181335_raw.all.mb56.fbt
│   ├── 0001_20020425_202049_raw.all.mb56.fbt
│   ├── 0002_20020425_185902_raw.all.mb56.fbt
│   ├── 0003_20020425_190538_raw.all.mb56.fbt
│   ├── 0004_20020425_190638_raw.all.mb56.fbt
│   ├── 0005_20020425_190813_raw.all.mb56.fbt
│   ├── 0006_20020425_191027_raw.all.mb56.fbt
│   ├── 0007_20020425_191214_raw.all.mb56.fbt
│   ├── 0008_20020425_191437_raw.all.mb56.fbt
│   ├── 0009_20020425_191721_raw.all.mb56.fbt
│   ├── 0010_20020425_192153_raw.all.mb56.fbt
│   ├── 0011_20020425_193111_raw.all.mb56.fbt
│   ├── 0012_20020425_193323_raw.all.mb56.fbt
│   ├── 0013_20020425_193506_raw.all.mb56.fbt
│   ├── 0014_20020425_193705_raw.all.mb56.fbt
│   ├── 0015_20020425_193847_raw.all.mb56.fbt
│   ├── 0016_20020425_194054_raw.all.mb56.fbt
│   ├── 0017_20020425_194335_raw.all.mb56.fbt
│   ├── 0018_20020425_194516_raw.all.mb56.fbt
│   ├── 0019_20020425_194635_raw.all.mb56.fbt
│   ├── 0020_20020425_194758_raw.all.mb56.fbt
│   ├── 0021_20020425_194936_raw.all.mb56.fbt
│   ├── 0022_20020425_195119_raw.all.mb56.fbt
│   ├── 0023_20020425_195250_raw.all.mb56.fbt
│   ├── 0024_20020425_195422_raw.all.mb56.fbt
│   ├── 0025_20020425_195554_raw.all.mb56.fbt
│   └── 0026_20020425_195704_raw.all.mb56.fbt
├── TN146
│   ├── 0004_20020521_034644pp.mb57.fbt
│   ├── 0005_20020521_044644pp.mb57.fbt
│   └── 0006_20020521_054206pp.mb57.fbt
├── TN172
│   ├── 0001_20040902_090804_raw.all.mb56.fbt
│   └── 0002_20040902_130804_raw.all.mb56.fbt
├── TN183
│   ├── 0004_20051004_055442p.mb56.fbt
│   └── 0005_20051004_065442p.mb56.fbt
├── TN240
│   ├── 0003_20090929_053359_raw.all.mb56.fbt
│   ├── 0130_20091014_104455_raw.all.mb56.fbt
│   └── 0131_20091014_135723_raw.all.mb56.fbt
├── TN247
│   └── 0011_20100426_135319_raw.all.mb56.fbt
├── TN252
│   ├── 0006_20100727_183615_raw.all.mb56.fbt
│   ├── 0007_20100727_235217_raw.all.mb56.fbt
│   └── 0008_20100728_002748_raw.all.mb56.fbt
├── TN254
│   ├── 0003_20100909_193018_raw.all.mb56.fbt
│   ├── 0004_20100909_203019_raw.all.mb56.fbt
│   ├── 0005_20100912_164134_raw.all.mb56.fbt
│   ├── 0006_20100912_174134_raw.all.mb56.fbt
│   ├── 0007_20100912_184134_raw.all.mb56.fbt
│   ├── 0012_20100910_043019_raw.all.mb56.fbt
│   ├── 0013_20100910_053019_raw.all.mb56.fbt
│   ├── 0014_20100910_063019_raw.all.mb56.fbt
│   └── 0015_20100910_073019_raw.all.mb56.fbt
├── TN255
│   ├── 0013_20101017_005029.all.mb58.fbt
│   ├── 0014_20101017_012029.all.mb58.fbt
│   ├── 0015_20101017_013932.all.mb58.fbt
│   ├── 0016_20101017_014511.all.mb58.fbt
│   ├── 0017_20101017_014637.all.mb58.fbt
│   ├── 0018_20101017_014708.all.mb58.fbt
│   ├── 0019_20101017_021707.all.mb58.fbt
│   ├── 0020_20101017_024502.all.mb58.fbt
│   ├── 0021_20101017_024605.all.mb58.fbt
│   ├── 0022_20101017_031606.all.mb58.fbt
│   ├── 0027_20101017_042227.all.mb58.fbt
│   ├── 0029_20101017_043936.all.mb58.fbt
│   ├── 0030_20101017_045640.all.mb58.fbt
│   ├── 0031_20101017_050058.all.mb58.fbt
│   ├── 0032_20101017_050408.all.mb58.fbt
│   ├── 0033_20101017_050509.all.mb58.fbt
│   ├── 0034_20101017_053508.all.mb58.fbt
│   ├── 0035_20101017_060508.all.mb58.fbt
│   ├── 0036_20101017_063509.all.mb58.fbt
│   ├── 0037_20101017_070509.all.mb58.fbt
│   ├── 0038_20101017_073509.all.mb58.fbt
│   ├── 0039_20101017_080508.all.mb58.fbt
│   ├── 0040_20101017_083508.all.mb58.fbt
│   ├── 0041_20101017_090508.all.mb58.fbt
│   ├── 0042_20101017_093509.all.mb58.fbt
│   ├── 0043_20101017_100509.all.mb58.fbt
│   ├── 0044_20101017_103511.all.mb58.fbt
│   ├── 0045_20101017_110513.all.mb58.fbt
│   ├── 0046_20101017_113508.all.mb58.fbt
│   └── 0047_20101017_120509.all.mb58.fbt
├── TN264
│   ├── 0013_20110522_080050_ShipName.all.mb58.fbt
│   ├── 0014_20110522_083051_ShipName.all.mb58.fbt
│   ├── 0015_20110522_090050_ShipName.all.mb58.fbt
│   ├── 0016_20110522_093049_ShipName.all.mb58.fbt
│   ├── 0017_20110522_100050_ShipName.all.mb58.fbt
│   ├── 0019_20110522_110050_ShipName.all.mb58.fbt
│   ├── 0020_20110522_113050_ShipName.all.mb58.fbt
│   ├── 0021_20110522_120050_ShipName.all.mb58.fbt
│   ├── 0022_20110522_123050_ShipName.all.mb58.fbt
│   ├── 0023_20110522_130050_ShipName.all.mb58.fbt
│   ├── 0024_20110522_133050_ShipName.all.mb58.fbt
│   ├── 0025_20110522_140050_ShipName.all.mb58.fbt
│   ├── 0026_20110522_143050_ShipName.all.mb58.fbt
│   ├── 0027_20110522_150049_ShipName.all.mb58.fbt
│   ├── 0028_20110522_151819_ShipName.all.mb58.fbt
│   ├── 0029_20110522_154819_ShipName.all.mb58.fbt
│   ├── 0030_20110522_161819_ShipName.all.mb58.fbt
│   ├── 0031_20110522_164819_ShipName.all.mb58.fbt
│   ├── 0032_20110522_171819_ShipName.all.mb58.fbt
│   ├── 0033_20110522_174819_ShipName.all.mb58.fbt
│   ├── 0034_20110522_181707_ShipName.all.mb58.fbt
│   ├── 0035_20110522_184706_ShipName.all.mb58.fbt
│   ├── 0036_20110522_191706_ShipName.all.mb58.fbt
│   ├── 0037_20110522_194706_ShipName.all.mb58.fbt
│   ├── 0038_20110522_201706_ShipName.all.mb58.fbt
│   ├── 0039_20110522_202921_ShipName.all.mb58.fbt
│   ├── 0040_20110522_205920_ShipName.all.mb58.fbt
│   ├── 0041_20110522_212918_ShipName.all.mb58.fbt
│   ├── 0087_20110523_175750_ShipName.all.mb58.fbt
│   ├── 0088_20110523_182747_ShipName.all.mb58.fbt
│   ├── 0089_20110523_185747_ShipName.all.mb58.fbt
│   ├── 0090_20110523_192747_ShipName.all.mb58.fbt
│   ├── 0091_20110523_214701_ShipName.all.mb58.fbt
│   └── 0092_20110523_222057_ShipName.all.mb58.fbt
├── TN265
│   ├── 0009_20110529_120005_TGT.all.mb58.fbt
│   ├── 0300_20110607_164726_TGT.all.mb58.fbt
│   ├── 0301_20110607_171726_TGT.all.mb58.fbt
│   ├── 0302_20110607_174727_TGT.all.mb58.fbt
│   ├── 0303_20110607_181728_TGT.all.mb58.fbt
│   ├── 0304_20110607_184726_TGT.all.mb58.fbt
│   ├── 0305_20110607_203505_TGT.all.mb58.fbt
│   ├── 0306_20110607_210459_TGT.all.mb58.fbt
│   ├── 0307_20110607_213500_TGT.all.mb58.fbt
│   ├── 0308_20110607_220500_TGT.all.mb58.fbt
│   ├── 0309_20110607_223459_TGT.all.mb58.fbt
│   ├── 0434_20110609_222319_TGT.all.mb58.fbt
│   ├── 0435_20110609_225319_TGT.all.mb58.fbt
│   ├── 0436_20110609_232320_TGT.all.mb58.fbt
│   ├── 0437_20110609_234901_TGT.all.mb58.fbt
│   ├── 0438_20110610_001901_TGT.all.mb58.fbt
│   ├── 0439_20110610_004901_TGT.all.mb58.fbt
│   ├── 0523_20110611_160252_TGT.all.mb58.fbt
│   ├── 0524_20110611_163252_TGT.all.mb58.fbt
│   ├── 0537_20110611_221302_TGT.all.mb58.fbt
│   ├── 0538_20110611_224302_TGT.all.mb58.fbt
│   ├── 1135_20110623_232523_TGT.all.mb58.fbt
│   ├── 1136_20110623_235523_TGT.all.mb58.fbt
│   ├── 1137_20110624_002524_TGT.all.mb58.fbt
│   ├── 1138_20110624_005522_TGT.all.mb58.fbt
│   ├── 1139_20110624_012523_TGT.all.mb58.fbt
│   └── 1140_20110624_015522_TGT.all.mb58.fbt
├── TN267
│   ├── 0014_20110730_021802_TGT.all.mb58.fbt
│   ├── 0015_20110730_024802_TGT.all.mb58.fbt
│   ├── 0016_20110730_031802_TGT.all.mb58.fbt
│   ├── 0017_20110730_034802_TGT.all.mb58.fbt
│   ├── 0279_20110804_021211_TGT.all.mb58.fbt
│   ├── 0280_20110804_024211_TGT.all.mb58.fbt
│   ├── 0281_20110804_031210_TGT.all.mb58.fbt
│   ├── 0282_20110804_034211_TGT.all.mb58.fbt
│   ├── 0283_20110804_041211_TGT.all.mb58.fbt
│   ├── 0284_20110804_044211_TGT.all.mb58.fbt
│   ├── 0285_20110804_051210_TGT.all.mb58.fbt
│   ├── 0303_20110805_025656_TGT.all.mb58.fbt
│   ├── 0309_20110808_213120_TGT.all.mb58.fbt
│   ├── 0310_20110808_214059_TGT.all.mb58.fbt
│   ├── 0311_20110808_214321_TGT.all.mb58.fbt
│   ├── 0312_20110808_214345_TGT.all.mb58.fbt
│   ├── 0313_20110808_214931_TGT.all.mb58.fbt
│   └── 0314_20110808_221931_TGT.all.mb58.fbt
├── TN268
│   ├── 0008_20110813_162338_TGT.all.mb58.fbt
│   ├── 0009_20110813_165338_TGT.all.mb58.fbt
│   ├── 0010_20110813_172337_TGT.all.mb58.fbt
│   ├── 0011_20110813_175337_TGT.all.mb58.fbt
│   └── 0012_20110813_182337_TGT.all.mb58.fbt
├── TN269
│   ├── 0019_20110910_234153_TGT.all.mb58.fbt
│   ├── 0020_20110911_001153_TGT.all.mb58.fbt
│   ├── 0021_20110911_004153_TGT.all.mb58.fbt
│   ├── 0022_20110911_011153_TGT.all.mb58.fbt
│   ├── 0415_20110930_034450_TGT.all.mb58.fbt
│   ├── 0416_20110930_041450_TGT.all.mb58.fbt
│   ├── 0417_20110930_044450_TGT.all.mb58.fbt
│   ├── 0418_20110930_051450_TGT.all.mb58.fbt
│   └── 0419_20110930_054450_TGT.all.mb58.fbt
├── TN270
│   ├── 0000_20111013_231528_TGT.all.mb58.fbt
│   ├── 0001_20111013_234518_TGT.all.mb58.fbt
│   ├── 0002_20111014_025530_TGT.all.mb58.fbt
│   ├── 0003_20111014_032530_TGT.all.mb58.fbt
│   ├── 0004_20111014_035530_TGT.all.mb58.fbt
│   ├── 0005_20111014_042530_TGT.all.mb58.fbt
│   ├── 0006_20111014_045530_TGT.all.mb58.fbt
│   ├── 0007_20111014_052530_TGT.all.mb58.fbt
│   ├── 0008_20111014_055530_TGT.all.mb58.fbt
│   ├── 0009_20111014_062530_TGT.all.mb58.fbt
│   ├── 0010_20111014_065530_TGT.all.mb58.fbt
│   ├── 0011_20111014_065726_TGT.all.mb58.fbt
│   ├── 0012_20111014_072715_TGT.all.mb58.fbt
│   ├── 0013_20111014_075715_TGT.all.mb58.fbt
│   ├── 0014_20111014_082715_TGT.all.mb58.fbt
│   ├── 0015_20111014_085715_TGT.all.mb58.fbt
│   ├── 0070_20111013_172034_TGT.all.mb58.fbt
│   ├── 0071_20111013_173841_TGT.all.mb58.fbt
│   ├── 0072_20111013_174309_TGT.all.mb58.fbt
│   ├── 0073_20111013_181309_TGT.all.mb58.fbt
│   ├── 0074_20111013_184309_TGT.all.mb58.fbt
│   ├── 0075_20111013_191309_TGT.all.mb58.fbt
│   ├── 0076_20111013_194309_TGT.all.mb58.fbt
│   ├── 0077_20111013_201246_TGT.all.mb58.fbt
│   ├── 0078_20111013_204246_TGT.all.mb58.fbt
│   ├── 0079_20111013_211246_TGT.all.mb58.fbt
│   ├── 0080_20111013_214246_TGT.all.mb58.fbt
│   ├── 0081_20111013_221246_TGT.all.mb58.fbt
│   └── 0082_20111013_224246_TGT.all.mb58.fbt
├── TN279
│   ├── 0208_20120429_051026_TGT.all.mb58.fbt
│   ├── 0209_20120429_054027_TGT.all.mb58.fbt
│   ├── 0210_20120429_061026_TGT.all.mb58.fbt
│   ├── 0211_20120429_064027_TGT.all.mb58.fbt
│   ├── 0212_20120429_071026_TGT.all.mb58.fbt
│   └── 0213_20120429_074027_TGT.all.mb58.fbt
├── TN280
│   ├── 0025_20120517_002340_TGT.all.mb58.fbt
│   ├── 0026_20120517_005341_TGT.all.mb58.fbt
│   ├── 0027_20120517_012341_TGT.all.mb58.fbt
│   ├── 0028_20120517_015340_TGT.all.mb58.fbt
│   ├── 0029_20120517_022340_TGT.all.mb58.fbt
│   ├── 0030_20120517_025340_TGT.all.mb58.fbt
│   ├── 0031_20120517_032340_TGT.all.mb58.fbt
│   ├── 0032_20120517_035341_TGT.all.mb58.fbt
│   ├── 0033_20120517_042341_TGT.all.mb58.fbt
│   ├── 0034_20120517_045340_TGT.all.mb58.fbt
│   ├── 0035_20120517_052340_TGT.all.mb58.fbt
│   ├── 0036_20120517_055340_TGT.all.mb58.fbt
│   ├── 0148_20120522_000850_TGT.all.mb58.fbt
│   ├── 0149_20120522_003849_TGT.all.mb58.fbt
│   ├── 0150_20120522_010850_TGT.all.mb58.fbt
│   ├── 0151_20120522_013850_TGT.all.mb58.fbt
│   ├── 0152_20120522_020850_TGT.all.mb58.fbt
│   └── 0153_20120522_023849_TGT.all.mb58.fbt
├── TN281
│   ├── 0016_20120525_033836_TGT.all.mb58.fbt
│   ├── 0017_20120525_040836_TGT.all.mb58.fbt
│   ├── 0018_20120525_043836_TGT.all.mb58.fbt
│   ├── 0019_20120525_050836_TGT.all.mb58.fbt
│   ├── 0020_20120525_053836_TGT.all.mb58.fbt
│   ├── 0023_20120525_070837_TGT.all.mb58.fbt
│   ├── 0024_20120525_073836_TGT.all.mb58.fbt
│   ├── 0025_20120517_002340_TGT.all.mb58.fbt
│   ├── 0025_20120525_080836_TGT.all.mb58.fbt
│   ├── 0026_20120517_005341_TGT.all.mb58.fbt
│   ├── 0026_20120525_083836_TGT.all.mb58.fbt
│   ├── 0027_20120517_012341_TGT.all.mb58.fbt
│   ├── 0027_20120525_090836_TGT.all.mb58.fbt
│   ├── 0028_20120517_015340_TGT.all.mb58.fbt
│   ├── 0028_20120525_093836_TGT.all.mb58.fbt
│   ├── 0029_20120517_022340_TGT.all.mb58.fbt
│   ├── 0029_20120525_100836_TGT.all.mb58.fbt
│   ├── 0030_20120517_025340_TGT.all.mb58.fbt
│   ├── 0030_20120525_103836_TGT.all.mb58.fbt
│   ├── 0031_20120517_032340_TGT.all.mb58.fbt
│   ├── 0031_20120525_110836_TGT.all.mb58.fbt
│   ├── 0032_20120517_035341_TGT.all.mb58.fbt
│   ├── 0032_20120525_113836_TGT.all.mb58.fbt
│   ├── 0033_20120517_042341_TGT.all.mb58.fbt
│   ├── 0033_20120525_120837_TGT.all.mb58.fbt
│   ├── 0034_20120517_045340_TGT.all.mb58.fbt
│   ├── 0034_20120525_123836_TGT.all.mb58.fbt
│   ├── 0035_20120517_052340_TGT.all.mb58.fbt
│   ├── 0035_20120525_130836_TGT.all.mb58.fbt
│   ├── 0036_20120517_055340_TGT.all.mb58.fbt
│   ├── 0036_20120525_133836_TGT.all.mb58.fbt
│   ├── 0037_20120525_140836_TGT.all.mb58.fbt
│   ├── 0038_20120525_143836_TGT.all.mb58.fbt
│   ├── 0039_20120525_150836_TGT.all.mb58.fbt
│   ├── 0040_20120525_153836_TGT.all.mb58.fbt
│   ├── 0041_20120525_160837_TGT.all.mb58.fbt
│   ├── 0042_20120525_163836_TGT.all.mb58.fbt
│   ├── 0043_20120525_170836_TGT.all.mb58.fbt
│   ├── 0044_20120525_173837_TGT.all.mb58.fbt
│   ├── 0045_20120525_180837_TGT.all.mb58.fbt
│   ├── 0046_20120525_183836_TGT.all.mb58.fbt
│   ├── 0056_20120525_233837_TGT.all.mb58.fbt
│   ├── 0057_20120526_000837_TGT.all.mb58.fbt
│   ├── 0058_20120526_003837_TGT.all.mb58.fbt
│   ├── 0059_20120526_010837_TGT.all.mb58.fbt
│   ├── 0060_20120526_013837_TGT.all.mb58.fbt
│   ├── 0061_20120526_020837_TGT.all.mb58.fbt
│   ├── 0062_20120526_023837_TGT.all.mb58.fbt
│   ├── 0063_20120526_030838_TGT.all.mb58.fbt
│   ├── 0064_20120526_033838_TGT.all.mb58.fbt
│   ├── 0065_20120526_040838_TGT.all.mb58.fbt
│   ├── 0066_20120526_043838_TGT.all.mb58.fbt
│   ├── 0067_20120526_050837_TGT.all.mb58.fbt
│   ├── 0068_20120526_053838_TGT.all.mb58.fbt
│   ├── 0148_20120522_000850_TGT.all.mb58.fbt
│   ├── 0149_20120522_003849_TGT.all.mb58.fbt
│   ├── 0150_20120522_010850_TGT.all.mb58.fbt
│   ├── 0151_20120522_013850_TGT.all.mb58.fbt
│   ├── 0152_20120522_020850_TGT.all.mb58.fbt
│   ├── 0153_20120522_023849_TGT.all.mb58.fbt
│   ├── 0208_20120429_051026_TGT.all.mb58.fbt
│   ├── 0209_20120429_054027_TGT.all.mb58.fbt
│   ├── 0210_20120429_061026_TGT.all.mb58.fbt
│   ├── 0211_20120429_064027_TGT.all.mb58.fbt
│   ├── 0212_20120429_071026_TGT.all.mb58.fbt
│   └── 0213_20120429_074027_TGT.all.mb58.fbt
├── TN282
│   ├── 0010_20120528_063724_TGT.all.mb58.fbt
│   ├── 0011_20120528_070723_TGT.all.mb58.fbt
│   ├── 0012_20120528_073724_TGT.all.mb58.fbt
│   ├── 0013_20120528_080723_TGT.all.mb58.fbt
│   ├── 0322_20120626_040420_TGT.all.mb58.fbt
│   ├── 0323_20120626_043420_TGT.all.mb58.fbt
│   ├── 0324_20120626_050420_TGT.all.mb58.fbt
│   ├── 0325_20120626_053420_TGT.all.mb58.fbt
│   ├── 0326_20120626_060421_TGT.all.mb58.fbt
│   ├── 0327_20120626_063420_TGT.all.mb58.fbt
│   └── 0328_20120626_070421_TGT.all.mb58.fbt
├── TN283
│   ├── 0016_20120711_020129_TGT.all.mb58.fbt
│   ├── 0017_20120711_023130_TGT.all.mb58.fbt
│   ├── 0018_20120711_030130_TGT.all.mb58.fbt
│   ├── 0019_20120711_033129_TGT.all.mb58.fbt
│   ├── 0020_20120711_040129_TGT.all.mb58.fbt
│   ├── 0411_20120722_235926_TGT.all.mb58.fbt
│   ├── 0412_20120723_002926_TGT.all.mb58.fbt
│   ├── 0413_20120723_005926_TGT.all.mb58.fbt
│   ├── 0414_20120723_012926_TGT.all.mb58.fbt
│   ├── 0415_20120723_015926_TGT.all.mb58.fbt
│   ├── 0416_20120723_022926_TGT.all.mb58.fbt
│   ├── 0417_20120723_025926_TGT.all.mb58.fbt
│   ├── 0418_20120723_032926_TGT.all.mb58.fbt
│   ├── 0419_20120723_035926_TGT.all.mb58.fbt
│   ├── 0420_20120723_042926_TGT.all.mb58.fbt
│   ├── 0420_20120723_045139_TGT.all.mb58.fbt
│   └── 0421_20120723_052140_TGT.all.mb58.fbt
├── TN290B
│   ├── 0012_20130117_044306_TGT.all.mb58.fbt
│   ├── 0013_20130117_051307_TGT.all.mb58.fbt
│   ├── 0014_20130117_054307_TGT.all.mb58.fbt
│   ├── 0015_20130117_061306_TGT.all.mb58.fbt
│   ├── 0016_20130117_064307_TGT.all.mb58.fbt
│   ├── 0028_20130117_132630_TGT.all.mb58.fbt
│   ├── 0029_20130117_135630_TGT.all.mb58.fbt
│   ├── 0030_20130117_142630_TGT.all.mb58.fbt
│   ├── 0031_20130117_145630_TGT.all.mb58.fbt
│   ├── 0032_20130117_152630_TGT.all.mb58.fbt
│   ├── 0033_20130117_202630_TGT.all.mb58.fbt
│   ├── 0034_20130117_205630_TGT.all.mb58.fbt
│   ├── 0075_20130118_174054_TGT.all.mb58.fbt
│   ├── 0076_20130118_181053_TGT.all.mb58.fbt
│   ├── 0077_20130118_184053_TGT.all.mb58.fbt
│   ├── 0078_20130118_191053_TGT.all.mb58.fbt
│   ├── 0079_20130118_194054_TGT.all.mb58.fbt
│   ├── 0080_20130118_201054_TGT.all.mb58.fbt
│   ├── 0081_20130118_204054_TGT.all.mb58.fbt
│   ├── 0082_20130118_211054_TGT.all.mb58.fbt
│   ├── 0083_20130118_214054_TGT.all.mb58.fbt
│   ├── 0084_20130118_221054_TGT.all.mb58.fbt
│   ├── 0085_20130118_224059_TGT.all.mb58.fbt
│   ├── 0086_20130119_021748_TGT.all.mb58.fbt
│   ├── 0100_20130119_094424_TGT.all.mb58.fbt
│   ├── 0101_20130119_101424_TGT.all.mb58.fbt
│   ├── 0102_20130119_104424_TGT.all.mb58.fbt
│   ├── 0103_20130119_111425_TGT.all.mb58.fbt
│   ├── 0104_20130119_114338_TGT.all.mb58.fbt
│   ├── 0105_20130119_121337_TGT.all.mb58.fbt
│   └── 0106_20130119_124338_TGT.all.mb58.fbt
├── TN291
│   ├── 0255_20130202_052026_TGT.all.mb58.fbt
│   ├── 0256_20130202_055026_TGT.all.mb58.fbt
│   ├── 0257_20130202_062026_TGT.all.mb58.fbt
│   ├── 0258_20130202_065027_TGT.all.mb58.fbt
│   ├── 0259_20130202_072026_TGT.all.mb58.fbt
│   └── 0260_20130202_075026_TGT.all.mb58.fbt
├── TN296
│   ├── 0193_20130426_024156_TGT.all.mb58.fbt
│   ├── 0194_20130426_031156_TGT.all.mb58.fbt
│   ├── 0195_20130426_034156_TGT.all.mb58.fbt
│   ├── 0196_20130426_041157_TGT.all.mb58.fbt
│   └── 0197_20130426_044157_TGT.all.mb58.fbt
├── TN312
│   ├── 0090_20140627_051908_TGT.all.mb58.fbt
│   ├── 0091_20140627_054908_TGT.all.mb58.fbt
│   ├── 0137_20140628_114355_TGT.all.mb58.fbt
│   ├── 0138_20140628_121355_TGT.all.mb58.fbt
│   ├── 0139_20140628_124354_TGT.all.mb58.fbt
│   ├── 0140_20140628_141703_TGT.all.mb58.fbt
│   ├── 0141_20140628_144701_TGT.all.mb58.fbt
│   ├── 0142_20140628_151701_TGT.all.mb58.fbt
│   ├── 0143_20140628_154701_TGT.all.mb58.fbt
│   ├── 0144_20140628_171506_TGT.all.mb58.fbt
│   ├── 0326_20140705_101050_TGT.all.mb58.fbt
│   ├── 0327_20140705_104050_TGT.all.mb58.fbt
│   ├── 0328_20140705_111050_TGT.all.mb58.fbt
│   ├── 0329_20140705_114051_TGT.all.mb58.fbt
│   ├── 0330_20140705_121051_TGT.all.mb58.fbt
│   ├── 0331_20140705_124051_TGT.all.mb58.fbt
│   ├── 0332_20140705_131050_TGT.all.mb58.fbt
│   ├── 0333_20140705_134050_TGT.all.mb58.fbt
│   └── 0334_20140705_141050_TGT.all.mb58.fbt
├── TN313
│   ├── 0112_20141005_073528_TGT.all.mb58.fbt
│   ├── 0113_20141005_083528_TGT.all.mb58.fbt
│   └── 0114_20141005_093528_TGT.all.mb58.fbt
├── TN314
│   ├── 0002_20141011_085437_TGT.all.mb58.fbt
│   └── 0003_20141011_095437_TGT.all.mb58.fbt
├── TN323
│   ├── 0009_20150530_024807_TGT.all.mb58.fbt
│   ├── 0010_20150530_031808_TGT.all.mb58.fbt
│   ├── 0011_20150530_034807_TGT.all.mb58.fbt
│   ├── 0012_20150530_041808_TGT.all.mb58.fbt
│   └── 0013_20150530_044807_TGT.all.mb58.fbt
├── TN327
│   ├── 0086_20150815_004048_TGT.all.mb58.fbt
│   ├── 0087_20150815_011048_TGT.all.mb58.fbt
│   ├── 0088_20150815_014048_TGT.all.mb58.fbt
│   ├── 0089_20150815_021048_TGT.all.mb58.fbt
│   ├── 0090_20150815_024048_TGT.all.mb58.fbt
│   └── 0091_20150815_031048_TGT.all.mb58.fbt
├── TN330
│   ├── 0057_20150923_135041_TGT.all.mb58.fbt
│   ├── 0125_20150925_103553_TGT.all.mb58.fbt
│   ├── 0144_20150926_041451_TGT.all.mb58.fbt
│   ├── 0145_20150926_044451_TGT.all.mb58.fbt
│   ├── 0146_20150926_051452_TGT.all.mb58.fbt
│   ├── 0165_20150926_185254_TGT.all.mb58.fbt
│   ├── 0166_20150926_192255_TGT.all.mb58.fbt
│   ├── 0167_20150926_195255_TGT.all.mb58.fbt
│   ├── 0168_20150926_202255_TGT.all.mb58.fbt
│   └── 0169_20150926_205255_TGT.all.mb58.fbt
├── TN332
│   ├── 0073_20151024_044037_TGT.all.mb58.fbt
│   ├── 0074_20151024_051037_TGT.all.mb58.fbt
│   ├── 0075_20151024_054038_TGT.all.mb58.fbt
│   ├── 0076_20151024_061037_TGT.all.mb58.fbt
│   └── 0077_20151024_064038_TGT.all.mb58.fbt
├── TN333
│   ├── 0075_20151119_033830_TGT.all.mb58.fbt
│   ├── 0076_20151119_040830_TGT.all.mb58.fbt
│   ├── 0077_20151119_043830_TGT.all.mb58.fbt
│   ├── 0078_20151119_050830_TGT.all.mb58.fbt
│   ├── 0079_20151119_053830_TGT.all.mb58.fbt
│   ├── 0080_20151119_060829_TGT.all.mb58.fbt
│   ├── 0081_20151119_063830_TGT.all.mb58.fbt
│   ├── 0082_20151119_070829_TGT.all.mb58.fbt
│   ├── 0083_20151119_073830_TGT.all.mb58.fbt
│   ├── 0084_20151119_080830_TGT.all.mb58.fbt
│   ├── 0085_20151119_083830_TGT.all.mb58.fbt
│   ├── 0086_20151119_090830_TGT.all.mb58.fbt
│   ├── 0087_20151119_093831_TGT.all.mb58.fbt
│   ├── 0088_20151119_100830_TGT.all.mb58.fbt
│   ├── 0089_20151119_103831_TGT.all.mb58.fbt
│   ├── 0090_20151119_110831_TGT.all.mb58.fbt
│   ├── 0091_20151119_113831_TGT.all.mb58.fbt
│   ├── 0092_20151119_120830_TGT.all.mb58.fbt
│   ├── 0093_20151119_123831_TGT.all.mb58.fbt
│   ├── 0094_20151119_130830_TGT.all.mb58.fbt
│   ├── 0095_20151119_133830_TGT.all.mb58.fbt
│   ├── 0096_20151119_140830_TGT.all.mb58.fbt
│   ├── 0097_20151119_143830_TGT.all.mb58.fbt
│   ├── 0098_20151119_150830_TGT.all.mb58.fbt
│   ├── 0099_20151119_153830_TGT.all.mb58.fbt
│   ├── 0100_20151119_160830_TGT.all.mb58.fbt
│   ├── 0101_20151119_163830_TGT.all.mb58.fbt
│   ├── 0102_20151119_170830_TGT.all.mb58.fbt
│   ├── 0103_20151119_173830_TGT.all.mb58.fbt
│   ├── 0104_20151119_180830_TGT.all.mb58.fbt
│   ├── 0105_20151120_001110_TGT.all.mb58.fbt
│   ├── 0106_20151120_004110_TGT.all.mb58.fbt
│   ├── 0107_20151120_011110_TGT.all.mb58.fbt
│   ├── 0108_20151120_014110_TGT.all.mb58.fbt
│   ├── 0109_20151120_021111_TGT.all.mb58.fbt
│   ├── 0110_20151120_024111_TGT.all.mb58.fbt
│   ├── 0111_20151120_031111_TGT.all.mb58.fbt
│   ├── 0112_20151120_034111_TGT.all.mb58.fbt
│   ├── 0113_20151120_041110_TGT.all.mb58.fbt
│   ├── 0114_20151120_044110_TGT.all.mb58.fbt
│   ├── 0120_20151120_074111_TGT.all.mb58.fbt
│   ├── 0121_20151120_081111_TGT.all.mb58.fbt
│   ├── 0122_20151120_084111_TGT.all.mb58.fbt
│   ├── 0123_20151120_091111_TGT.all.mb58.fbt
│   ├── 0124_20151120_094112_TGT.all.mb58.fbt
│   ├── 0125_20151120_101112_TGT.all.mb58.fbt
│   └── 0126_20151120_104112_TGT.all.mb58.fbt
├── TN342
│   ├── 0024_20160505_142527_TGT.all.mb58.fbt
│   ├── 0025_20160505_145527_TGT.all.mb58.fbt
│   ├── 0026_20160505_152527_TGT.all.mb58.fbt
│   ├── 0031_20160505_222948_TGT.all.mb58.fbt
│   ├── 0032_20160505_225948_TGT.all.mb58.fbt
│   ├── 0033_20160505_232948_TGT.all.mb58.fbt
│   ├── 0034_20160505_235948_TGT.all.mb58.fbt
│   ├── 0035_20160506_192224_TGT.all.mb58.fbt
│   ├── 0036_20160506_195224_TGT.all.mb58.fbt
│   ├── 0037_20160506_225840_TGT.all.mb58.fbt
│   ├── 0038_20160506_232841_TGT.all.mb58.fbt
│   ├── 0144_20160513_080219_TGT.all.mb58.fbt
│   ├── 0145_20160513_083219_TGT.all.mb58.fbt
│   ├── 0146_20160513_090219_TGT.all.mb58.fbt
│   ├── 0147_20160513_093219_TGT.all.mb58.fbt
│   ├── 0148_20160513_100219_TGT.all.mb58.fbt
│   ├── 0149_20160513_103219_TGT.all.mb58.fbt
│   ├── 0150_20160513_110219_TGT.all.mb58.fbt
│   ├── 0151_20160513_113219_TGT.all.mb58.fbt
│   ├── 0152_20160513_120219_TGT.all.mb58.fbt
│   ├── 0153_20160513_124226_TGT.all.mb58.fbt
│   ├── 0154_20160513_131226_TGT.all.mb58.fbt
│   ├── 0155_20160513_134226_TGT.all.mb58.fbt
│   ├── 0156_20160513_141226_TGT.all.mb58.fbt
│   ├── 0157_20160513_144226_TGT.all.mb58.fbt
│   ├── 0158_20160513_153314_TGT.all.mb58.fbt
│   ├── 0159_20160513_211003_TGT.all.mb58.fbt
│   ├── 0160_20160513_214004_TGT.all.mb58.fbt
│   └── 0161_20160513_221004_TGT.all.mb58.fbt
├── TN347
│   ├── 0031_20180110_060149_TGT.all.mb58.fbt
│   ├── 0032_20180110_070139_TGT.all.mb58.fbt
│   ├── 0033_20180110_080139_TGT.all.mb58.fbt
│   ├── 0093_20180112_022428_TGT.all.mb58.fbt
│   ├── 0094_20180112_032428_TGT.all.mb58.fbt
│   └── 0095_20180112_042428_TGT.all.mb58.fbt
├── TN380
│   ├── 0023_20200704_140030_TGT.all.mb58.fbt
│   ├── 0024_20200704_143030_TGT.all.mb58.fbt
│   ├── 0025_20200704_150030_TGT.all.mb58.fbt
│   ├── 0026_20200704_153030_TGT.all.mb58.fbt
│   ├── 0027_20200704_160030_TGT.all.mb58.fbt
│   ├── 0028_20200704_163030_TGT.all.mb58.fbt
│   ├── 0029_20200704_170030_TGT.all.mb58.fbt
│   ├── 0030_20200704_195652_TGT.all.mb58.fbt
│   ├── 0031_20200704_202651_TGT.all.mb58.fbt
│   ├── 0032_20200704_205651_TGT.all.mb58.fbt
│   ├── 0033_20200704_223736_TGT.all.mb58.fbt
│   ├── 0034_20200704_230718_TGT.all.mb58.fbt
│   ├── 0035_20200704_233718_TGT.all.mb58.fbt
│   ├── 0036_20200705_000718_TGT.all.mb58.fbt
│   ├── 0037_20200705_003719_TGT.all.mb58.fbt
│   ├── 0038_20200705_010719_TGT.all.mb58.fbt
│   ├── 0126_20200707_080438_TGT.all.mb58.fbt
│   ├── 0127_20200707_083438_TGT.all.mb58.fbt
│   ├── 0128_20200707_090439_TGT.all.mb58.fbt
│   ├── 0129_20200707_093439_TGT.all.mb58.fbt
│   ├── 0130_20200707_100438_TGT.all.mb58.fbt
│   ├── 0131_20200707_103439_TGT.all.mb58.fbt
│   ├── 0132_20200707_110439_TGT.all.mb58.fbt
│   ├── 0133_20200707_113439_TGT.all.mb58.fbt
│   ├── 0134_20200707_120439_TGT.all.mb58.fbt
│   ├── 0139_20200707_143439_TGT.all.mb58.fbt
│   ├── 0140_20200707_150439_TGT.all.mb58.fbt
│   ├── 0141_20200707_153439_TGT.all.mb58.fbt
│   ├── 0142_20200707_160439_TGT.all.mb58.fbt
│   ├── 0152_20200707_210439_TGT.all.mb58.fbt
│   ├── 0153_20200707_213439_TGT.all.mb58.fbt
│   ├── 0154_20200707_220439_TGT.all.mb58.fbt
│   ├── 0155_20200707_223439_TGT.all.mb58.fbt
│   ├── 0161_20200708_013445_TGT.all.mb58.fbt
│   ├── 0162_20200708_020440_TGT.all.mb58.fbt
│   ├── 0399_20200714_153210_TGT.all.mb58.fbt
│   ├── 0400_20200714_160210_TGT.all.mb58.fbt
│   ├── 0401_20200714_163210_TGT.all.mb58.fbt
│   ├── 0402_20200714_170210_TGT.all.mb58.fbt
│   ├── 0403_20200714_173211_TGT.all.mb58.fbt
│   ├── 0404_20200714_180211_TGT.all.mb58.fbt
│   ├── 0405_20200714_183211_TGT.all.mb58.fbt
│   ├── 0406_20200714_190211_TGT.all.mb58.fbt
│   ├── 0407_20200714_193211_TGT.all.mb58.fbt
│   ├── 0408_20200714_200211_TGT.all.mb58.fbt
│   ├── 0409_20200714_203211_TGT.all.mb58.fbt
│   ├── 0410_20200714_210211_TGT.all.mb58.fbt
│   ├── 0411_20200714_213211_TGT.all.mb58.fbt
│   ├── 0412_20200714_220211_TGT.all.mb58.fbt
│   ├── 0413_20200714_223211_TGT.all.mb58.fbt
│   ├── 0414_20200714_231017_TGT.all.mb58.fbt
│   ├── 0415_20200714_234017_TGT.all.mb58.fbt
│   ├── 0416_20200715_001017_TGT.all.mb58.fbt
│   ├── 0417_20200715_004017_TGT.all.mb58.fbt
│   ├── 0434_20200715_091018_TGT.all.mb58.fbt
│   ├── 0435_20200715_094026_TGT.all.mb58.fbt
│   ├── 0438_20200715_111026_TGT.all.mb58.fbt
│   ├── 0439_20200715_114028_TGT.all.mb58.fbt
│   ├── 0440_20200715_121021_TGT.all.mb58.fbt
│   ├── 0452_20200715_181018_TGT.all.mb58.fbt
│   ├── 0453_20200715_184018_TGT.all.mb58.fbt
│   ├── 0454_20200715_191019_TGT.all.mb58.fbt
│   ├── 0455_20200715_194019_TGT.all.mb58.fbt
│   ├── 0456_20200715_201019_TGT.all.mb58.fbt
│   ├── 0457_20200715_204019_TGT.all.mb58.fbt
│   ├── 0458_20200715_211019_TGT.all.mb58.fbt
│   ├── 0459_20200715_214019_TGT.all.mb58.fbt
│   ├── 0460_20200715_221019_TGT.all.mb58.fbt
│   ├── 0461_20200715_224019_TGT.all.mb58.fbt
│   ├── 0462_20200715_231019_TGT.all.mb58.fbt
│   ├── 0463_20200715_234019_TGT.all.mb58.fbt
│   ├── 0464_20200716_001019_TGT.all.mb58.fbt
│   ├── 0465_20200716_004019_TGT.all.mb58.fbt
│   ├── 0466_20200716_011019_TGT.all.mb58.fbt
│   └── 0467_20200716_014019_TGT.all.mb58.fbt
├── TN384
│   ├── 0034_20200926_150059_TGT.all.mb58.fbt
│   ├── 0035_20200926_153059_TGT.all.mb58.fbt
│   ├── 0036_20200926_160059_TGT.all.mb58.fbt
│   ├── 0037_20200926_163100_TGT.all.mb58.fbt
│   ├── 0038_20200926_170059_TGT.all.mb58.fbt
│   ├── 0039_20200926_173059_TGT.all.mb58.fbt
│   ├── 0040_20200927_015007_TGT.all.mb58.fbt
│   ├── 0058_20200928_024311_TGT.all.mb58.fbt
│   ├── 0059_20200928_031311_TGT.all.mb58.fbt
│   ├── 0060_20200928_033159_TGT.all.mb58.fbt
│   ├── 0061_20200928_040155_TGT.all.mb58.fbt
│   ├── 0080_20200928_133155_TGT.all.mb58.fbt
│   ├── 0081_20200928_140156_TGT.all.mb58.fbt
│   ├── 0082_20200928_143156_TGT.all.mb58.fbt
│   ├── 0083_20200928_150155_TGT.all.mb58.fbt
│   ├── 0084_20200928_153156_TGT.all.mb58.fbt
│   ├── 0085_20200928_160156_TGT.all.mb58.fbt
│   ├── 0086_20200928_163156_TGT.all.mb58.fbt
│   ├── 0087_20200928_170156_TGT.all.mb58.fbt
│   ├── 0088_20200928_173156_TGT.all.mb58.fbt
│   ├── 0089_20200928_180156_TGT.all.mb58.fbt
│   ├── 0090_20200928_183156_TGT.all.mb58.fbt
│   ├── 0091_20200928_190156_TGT.all.mb58.fbt
│   ├── 0092_20200928_193156_TGT.all.mb58.fbt
│   ├── 0093_20200928_200156_TGT.all.mb58.fbt
│   ├── 0094_20200928_203156_TGT.all.mb58.fbt
│   ├── 0095_20200928_210156_TGT.all.mb58.fbt
│   ├── 0096_20200928_213157_TGT.all.mb58.fbt
│   ├── 0097_20200928_220156_TGT.all.mb58.fbt
│   ├── 0098_20200928_223157_TGT.all.mb58.fbt
│   ├── 0099_20200928_230156_TGT.all.mb58.fbt
│   ├── 0100_20200928_233157_TGT.all.mb58.fbt
│   ├── 0101_20200929_000156_TGT.all.mb58.fbt
│   ├── 0102_20200929_003157_TGT.all.mb58.fbt
│   ├── 0103_20200929_010158_TGT.all.mb58.fbt
│   ├── 0104_20200929_013158_TGT.all.mb58.fbt
│   ├── 0105_20200929_020158_TGT.all.mb58.fbt
│   ├── 0106_20200929_023157_TGT.all.mb58.fbt
│   ├── 0128_20200929_124737_TGT.all.mb58.fbt
│   ├── 0129_20200929_131737_TGT.all.mb58.fbt
│   ├── 0130_20200929_134737_TGT.all.mb58.fbt
│   ├── 0131_20200929_141737_TGT.all.mb58.fbt
│   ├── 0132_20200929_144737_TGT.all.mb58.fbt
│   ├── 0133_20200929_151737_TGT.all.mb58.fbt
│   ├── 0134_20200929_154737_TGT.all.mb58.fbt
│   ├── 0135_20200929_161737_TGT.all.mb58.fbt
│   ├── 0136_20200929_164737_TGT.all.mb58.fbt
│   ├── 0137_20200929_171738_TGT.all.mb58.fbt
│   └── 0138_20200929_174738_TGT.all.mb58.fbt
├── TN393
│   ├── 0015_20210802_023713_TGT.all.mb58.fbt
│   ├── 0016_20210802_030714_TGT.all.mb58.fbt
│   ├── 0017_20210802_033714_TGT.all.mb58.fbt
│   ├── 0018_20210802_040714_TGT.all.mb58.fbt
│   └── 0019_20210802_043714_TGT.all.mb58.fbt
├── TN394
│   ├── 0072_20210911_104326_TGT.all.mb58.fbt
│   ├── 0073_20210911_111326_TGT.all.mb58.fbt
│   ├── 0074_20210911_114327_TGT.all.mb58.fbt
│   ├── 0075_20210911_121327_TGT.all.mb58.fbt
│   ├── 0076_20210911_124327_TGT.all.mb58.fbt
│   ├── 0077_20210911_131326_TGT.all.mb58.fbt
│   ├── 0078_20210911_134326_TGT.all.mb58.fbt
│   ├── 0079_20210911_141327_TGT.all.mb58.fbt
│   ├── 0080_20210911_144327_TGT.all.mb58.fbt
│   ├── 0081_20210911_151327_TGT.all.mb58.fbt
│   ├── 0082_20210911_154327_TGT.all.mb58.fbt
│   ├── 0083_20210911_161327_TGT.all.mb58.fbt
│   ├── 0086_20210911_174327_TGT.all.mb58.fbt
│   ├── 0087_20210911_181327_TGT.all.mb58.fbt
│   ├── 0088_20210911_184327_TGT.all.mb58.fbt
│   ├── 0089_20210911_191327_TGT.all.mb58.fbt
│   ├── 0093_20210911_211327_TGT.all.mb58.fbt
│   ├── 0094_20210911_214327_TGT.all.mb58.fbt
│   ├── 0095_20210911_221327_TGT.all.mb58.fbt
│   ├── 0096_20210911_224327_TGT.all.mb58.fbt
│   ├── 0097_20210911_231327_TGT.all.mb58.fbt
│   ├── 0098_20210911_234327_TGT.all.mb58.fbt
│   ├── 0099_20210912_001327_TGT.all.mb58.fbt
│   ├── 0100_20210912_004327_TGT.all.mb58.fbt
│   ├── 0101_20210912_011319_TGT.all.mb58.fbt
│   ├── 0102_20210912_014319_TGT.all.mb58.fbt
│   ├── 0238_20210914_224229_TGT.all.mb58.fbt
│   ├── 0239_20210914_231230_TGT.all.mb58.fbt
│   ├── 0241_20210915_001230_TGT.all.mb58.fbt
│   ├── 0242_20210915_004230_TGT.all.mb58.fbt
│   ├── 0243_20210915_011230_TGT.all.mb58.fbt
│   ├── 0244_20210915_014230_TGT.all.mb58.fbt
│   ├── 0303_20210920_105050_TGT.all.mb58.fbt
│   ├── 0304_20210920_112050_TGT.all.mb58.fbt
│   ├── 0305_20210920_115050_TGT.all.mb58.fbt
│   ├── 0306_20210920_122050_TGT.all.mb58.fbt
│   ├── 0307_20210920_125050_TGT.all.mb58.fbt
│   ├── 0308_20210920_132050_TGT.all.mb58.fbt
│   ├── 0309_20210920_135050_TGT.all.mb58.fbt
│   ├── 0310_20210920_142050_TGT.all.mb58.fbt
│   └── 0311_20210920_145050_TGT.all.mb58.fbt
├── TN396
│   ├── 0019_20211112_040633_TGT.all.mb58.fbt
│   ├── 0020_20211112_043632_TGT.all.mb58.fbt
│   ├── 0021_20211112_050633_TGT.all.mb58.fbt
│   ├── 0022_20211112_053633_TGT.all.mb58.fbt
│   ├── 0023_20211112_060632_TGT.all.mb58.fbt
│   └── 0024_20211112_063632_TGT.all.mb58.fbt
└── TN409
    ├── 0034_20221010_134412_TGT.all.mb58.fbt
    ├── 0035_20221010_141412_TGT.all.mb58.fbt
    ├── 0036_20221010_144412_TGT.all.mb58.fbt
    ├── 0037_20221010_151412_TGT.all.mb58.fbt
    ├── 0038_20221010_154412_TGT.all.mb58.fbt
    ├── 0039_20221010_161412_TGT.all.mb58.fbt
    ├── 0040_20221010_164412_TGT.all.mb58.fbt
    ├── 0041_20221010_171412_TGT.all.mb58.fbt
    ├── 0042_20221010_174108_TGT.all.mb58.fbt
    ├── 0043_20221010_181101_TGT.all.mb58.fbt
    ├── 0044_20221010_184101_TGT.all.mb58.fbt
    ├── 0045_20221010_191101_TGT.all.mb58.fbt
    ├── 0046_20221010_194101_TGT.all.mb58.fbt
    ├── 0047_20221010_201101_TGT.all.mb58.fbt
    ├── 0048_20221010_204101_TGT.all.mb58.fbt
    ├── 0049_20221010_211101_TGT.all.mb58.fbt
    ├── 0050_20221010_214101_TGT.all.mb58.fbt
    ├── 0051_20221010_221101_TGT.all.mb58.fbt
    ├── 0052_20221010_224101_TGT.all.mb58.fbt
    ├── 0053_20221010_231101_TGT.all.mb58.fbt
    ├── 0054_20221010_234101_TGT.all.mb58.fbt
    ├── 0055_20221011_001101_TGT.all.mb58.fbt
    ├── 0056_20221011_004101_TGT.all.mb58.fbt
    ├── 0059_20221011_021101_TGT.all.mb58.fbt
    ├── 0060_20221011_024101_TGT.all.mb58.fbt
    ├── 0061_20221011_031056_TGT.all.mb58.fbt
    ├── 0062_20221011_034056_TGT.all.mb58.fbt
    ├── 0063_20221011_041056_TGT.all.mb58.fbt
    ├── 0064_20221011_044056_TGT.all.mb58.fbt
    ├── 0065_20221011_051056_TGT.all.mb58.fbt
    ├── 0066_20221011_054056_TGT.all.mb58.fbt
    ├── 0067_20221011_061056_TGT.all.mb58.fbt
    ├── 0068_20221011_064056_TGT.all.mb58.fbt
    ├── 0069_20221011_071056_TGT.all.mb58.fbt
    ├── 0070_20221011_074057_TGT.all.mb58.fbt
    ├── 0071_20221011_081056_TGT.all.mb58.fbt
    ├── 0072_20221011_084057_TGT.all.mb58.fbt
    ├── 0073_20221011_091057_TGT.all.mb58.fbt
    ├── 0074_20221011_094057_TGT.all.mb58.fbt
    ├── 0075_20221011_101057_TGT.all.mb58.fbt
    ├── 0076_20221011_104057_TGT.all.mb58.fbt
    ├── 0077_20221011_111057_TGT.all.mb58.fbt
    ├── 0078_20221011_113511_TGT.all.mb58.fbt
    ├── 0079_20221011_144512_TGT.all.mb58.fbt
    ├── 0080_20221011_151511_TGT.all.mb58.fbt
    ├── 0081_20221011_155216_TGT.all.mb58.fbt
    ├── 0082_20221011_162211_TGT.all.mb58.fbt
    ├── 0083_20221011_165211_TGT.all.mb58.fbt
    ├── 0084_20221011_172211_TGT.all.mb58.fbt
    ├── 0085_20221011_175211_TGT.all.mb58.fbt
    ├── 0086_20221011_182211_TGT.all.mb58.fbt
    ├── 0087_20221011_185211_TGT.all.mb58.fbt
    ├── 0088_20221011_192211_TGT.all.mb58.fbt
    ├── 0089_20221011_195211_TGT.all.mb58.fbt
    ├── 0090_20221011_202211_TGT.all.mb58.fbt
    ├── 0091_20221011_205211_TGT.all.mb58.fbt
    ├── 0092_20221011_212211_TGT.all.mb58.fbt
    ├── 0093_20221011_215211_TGT.all.mb58.fbt
    ├── 0094_20221011_222211_TGT.all.mb58.fbt
    ├── 0095_20221011_225211_TGT.all.mb58.fbt
    ├── 0096_20221011_234455_TGT.all.mb58.fbt
    ├── 0097_20221012_063907_TGT.all.mb58.fbt
    ├── 0098_20221012_070907_TGT.all.mb58.fbt
    ├── 0099_20221012_073906_TGT.all.mb58.fbt
    ├── 0100_20221012_080907_TGT.all.mb58.fbt
    ├── 0101_20221012_083907_TGT.all.mb58.fbt
    ├── 0102_20221012_090907_TGT.all.mb58.fbt
    ├── 0103_20221012_093907_TGT.all.mb58.fbt
    ├── 0104_20221012_100907_TGT.all.mb58.fbt
    ├── 0105_20221012_103907_TGT.all.mb58.fbt
    ├── 0106_20221012_110907_TGT.all.mb58.fbt
    ├── 0107_20221012_113907_TGT.all.mb58.fbt
    ├── 0108_20221012_120907_TGT.all.mb58.fbt
    ├── 0109_20221012_123907_TGT.all.mb58.fbt
    ├── 0110_20221012_130907_TGT.all.mb58.fbt
    ├── 0111_20221012_133907_TGT.all.mb58.fbt
    ├── 0112_20221012_140907_TGT.all.mb58.fbt
    ├── 0113_20221012_143907_TGT.all.mb58.fbt
    ├── 0114_20221012_150907_TGT.all.mb58.fbt
    ├── 0115_20221012_153907_TGT.all.mb58.fbt
    ├── 0116_20221012_160907_TGT.all.mb58.fbt
    ├── 0117_20221012_163907_TGT.all.mb58.fbt
    ├── 0118_20221012_170907_TGT.all.mb58.fbt
    ├── 0119_20221012_173907_TGT.all.mb58.fbt
    ├── 0120_20221012_180907_TGT.all.mb58.fbt
    ├── 0121_20221012_183908_TGT.all.mb58.fbt
    ├── 0122_20221012_190907_TGT.all.mb58.fbt
    ├── 0123_20221012_193907_TGT.all.mb58.fbt
    ├── 0124_20221012_200907_TGT.all.mb58.fbt
    ├── 0126_20221012_210908_TGT.all.mb58.fbt
    ├── 0127_20221012_213908_TGT.all.mb58.fbt
    ├── 0129_20221012_223446_TGT.all.mb58.fbt
    ├── 0130_20221012_230446_TGT.all.mb58.fbt
    ├── 0131_20221012_233446_TGT.all.mb58.fbt
    ├── 0132_20221013_000446_TGT.all.mb58.fbt
    ├── 0133_20221013_003446_TGT.all.mb58.fbt
    ├── 0134_20221013_010446_TGT.all.mb58.fbt
    ├── 0135_20221013_013446_TGT.all.mb58.fbt
    ├── 0136_20221013_020446_TGT.all.mb58.fbt
    ├── 0137_20221013_023446_TGT.all.mb58.fbt
    ├── 0138_20221013_030446_TGT.all.mb58.fbt
    ├── 0139_20221013_033446_TGT.all.mb58.fbt
    ├── 0140_20221013_040446_TGT.all.mb58.fbt
    ├── 0141_20221013_043446_TGT.all.mb58.fbt
    ├── 0142_20221013_050446_TGT.all.mb58.fbt
    ├── 0143_20221013_053446_TGT.all.mb58.fbt
    ├── 0144_20221013_060446_TGT.all.mb58.fbt
    ├── 0145_20221013_063447_TGT.all.mb58.fbt
    ├── 0146_20221013_070447_TGT.all.mb58.fbt
    ├── 0147_20221013_073446_TGT.all.mb58.fbt
    ├── 0148_20221013_080446_TGT.all.mb58.fbt
    ├── 0149_20221013_083446_TGT.all.mb58.fbt
    ├── 0150_20221013_090446_TGT.all.mb58.fbt
    ├── 0151_20221013_093447_TGT.all.mb58.fbt
    ├── 0152_20221013_100447_TGT.all.mb58.fbt
    ├── 0153_20221013_103447_TGT.all.mb58.fbt
    ├── 0154_20221013_110447_TGT.all.mb58.fbt
    ├── 0154_20221013_123134_TGT.all.mb58.fbt
    ├── 0155_20221013_130131_TGT.all.mb58.fbt
    ├── 0156_20221013_133131_TGT.all.mb58.fbt
    ├── 0157_20221013_140131_TGT.all.mb58.fbt
    ├── 0158_20221013_143131_TGT.all.mb58.fbt
    ├── 0159_20221013_150131_TGT.all.mb58.fbt
    ├── 0160_20221013_153131_TGT.all.mb58.fbt
    ├── 0161_20221013_160131_TGT.all.mb58.fbt
    ├── 0162_20221013_163131_TGT.all.mb58.fbt
    ├── 0163_20221013_170131_TGT.all.mb58.fbt
    ├── 0164_20221013_173131_TGT.all.mb58.fbt
    ├── 0165_20221013_180131_TGT.all.mb58.fbt
    ├── 0166_20221013_183131_TGT.all.mb58.fbt
    ├── 0167_20221013_190131_TGT.all.mb58.fbt
    ├── 0168_20221013_193131_TGT.all.mb58.fbt
    ├── 0169_20221013_200131_TGT.all.mb58.fbt
    ├── 0170_20221013_203131_TGT.all.mb58.fbt
    ├── 0171_20221013_210131_TGT.all.mb58.fbt
    ├── 0172_20221013_213131_TGT.all.mb58.fbt
    ├── 0173_20221013_220131_TGT.all.mb58.fbt
    ├── 0174_20221013_223131_TGT.all.mb58.fbt
    ├── 0175_20221013_230131_TGT.all.mb58.fbt
    ├── 0176_20221013_233131_TGT.all.mb58.fbt
    ├── 0177_20221014_000131_TGT.all.mb58.fbt
    ├── 0178_20221014_003131_TGT.all.mb58.fbt
    ├── 0179_20221014_010131_TGT.all.mb58.fbt
    ├── 0180_20221014_013131_TGT.all.mb58.fbt
    ├── 0181_20221014_020132_TGT.all.mb58.fbt
    ├── 0182_20221014_023131_TGT.all.mb58.fbt
    ├── 0183_20221014_030132_TGT.all.mb58.fbt
    ├── 0184_20221014_033131_TGT.all.mb58.fbt
    ├── 0185_20221014_040132_TGT.all.mb58.fbt
    ├── 0186_20221014_043132_TGT.all.mb58.fbt
    ├── 0187_20221014_050132_TGT.all.mb58.fbt
    ├── 0188_20221014_053132_TGT.all.mb58.fbt
    ├── 0189_20221014_060132_TGT.all.mb58.fbt
    ├── 0190_20221014_063132_TGT.all.mb58.fbt
    ├── 0191_20221014_070132_TGT.all.mb58.fbt
    ├── 0192_20221014_073132_TGT.all.mb58.fbt
    ├── 0193_20221014_080132_TGT.all.mb58.fbt
    ├── 0194_20221014_083132_TGT.all.mb58.fbt
    ├── 0195_20221014_090132_TGT.all.mb58.fbt
    ├── 0196_20221014_093132_TGT.all.mb58.fbt
    └── 0197_20221014_100132_TGT.all.mb58.fbt
```

#### Crowd-Sourced Bathymetry
```bash
fetches -R tiles_1_9.shp csb
```

### Topography / Near-shore Bathymetry

#### Digital Coast Lidar
```bash
fetches -R tiles_1_9.shp digital_coast:datatype=lidar
```

#### USGS Lidar

#### USGS DEMs

#### CUDEMs
```bash
fetches -R tiles_1_9.shp CUDEM
```

```
digital_coast/
├── 8483
│   ├── ncei19_n47x00_w0124x25_2018v1.tif
│   ├── ncei19_n47x25_w0124x25_2018v1.tif
│   ├── ncei19_n47x50_w0124x25_2018v1.tif
│   ├── ncei19_n47x50_w0124x50_2018v1.tif
│   ├── ncei19_n47x75_w0124x50_2018v1.tif
│   ├── ncei19_n48x00_w0124x50_2018v1.tif
│   ├── ncei19_n48x00_w0124x75_2018v1.tif
│   ├── ncei19_n48x25_w0124x75_2018v1.tif
│   ├── ncei19_n48x25_w124x25_2021v1.tif
│   ├── ncei19_n48x25_w124x50_2021v1.tif
│   ├── ncei19_n48x50_w0124x75_2018v1.tif
│   └── ncei19_n48x50_w124x50_2021v1.tif
└── 8580
    ├── ncei13_n47x25_w0124x50_2018v1.tif
    ├── ncei13_n47x75_w0124x75_2018v1.tif
    ├── ncei13_n48x25_w0125x00_2018v1.tif
    └── ncei13_n48x50_w0125x00_2018v1.tif
```


## Make a datalist

- see fetches --help for available datasets, or gather and use your own data.

## Generate a test tile

pick a tile and generate an on-the-fly DEM to see what it looks like

either, use the region dimensions of the desired tile or select and export the tile to a new vector using a GIS.

