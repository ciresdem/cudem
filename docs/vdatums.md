# vdatums

Transform data between vertical datums

## Synopsis

```
vdatums [OPTIONS] input_grid output_grid
```

## Description

Transform rasters between vertical datums.

## Options

`-i, --vdatum_in`
the input vertical datum as EPSG code

`-o, --vdatum_out`
the output vertical datum as EPSG code

`-D, --cache-dir`
CACHE Directory for storing temp data.
Default Cache Directory is ~/.cudem_cache; cache will be cleared after a vdatums session
to retain the data, use the --keep-cache flag

`-k, --keep-cache`
KEEP the cache data intact after run

`-l, --list-epsg`
List the supported EPSG codes and their names

`-q, --quiet`
Lower verbosity to a quiet.

`--help`
Print the usage text

`--version`
Print the version information

## Python API

## Examples