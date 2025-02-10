
```bash
for i in ATL24*.h5; do fetches icesat2:filename_filter=$(echo $i | awk -F_ '{print $2}'); done
for i in ATL03*.h5; do icesat2dump.py $i > $(basename $i .h5)_bathy.xyz; done
wget https://www.ngdc.noaa.gov/thredds/fileServer/crm/cudem/crm_vol9_2023.nc
wget https://www.ngdc.noaa.gov/thredds/fileServer/crm/cudem/crm_vol3_2023.nc
for i in ATL03*_bathy.xyz; do gdal_query.py crm_vol9_2023.nc $i -d_format xyzgdsp | grep -v "\-99999.0000" > ATL24_vs_CRMv9.xyzgd; done
for i in ATL03*_bathy.xyz; do gdal_query.py crm_vol3_2023.nc $i -d_format xyzgdsp | grep -v "\-99999.0000" > ATL24_vs_CRMv3.xyzgd; done
```

```python
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

diff_file = sys.argv[1]
if os.path.exists(diff_file):
    data = np.loadtxt(diff_file)
    plt.hist(data, bins=30)
    plt.xlabel('diff')
    plt.ylabel('freq')
    plt.title('{}'.format(' '.join(diff_file.split('.')[0].split('_'))))
    plt.savefig('{}.png'.format(diff_file.split('.')[0]))
    data = None
```

```python
import sys
import pygmt
from cudem import utils
from cudem import regions
from cudem import dlim
from cudem import perspecto
import numpy as np

xyd_file = sys.argv[1]
grid_file = sys.argv[2]

etopo = perspecto.generate_etopo_cpt(-5000, 50)
data = dlim.DatasetFactory(mod=xyd_file, data_format='168:xpos=0:ypos=1:zpos=2')._acquire_module().initialize()    
this_inf = data.inf()
this_region = regions.Region().from_list(this_inf.minmax)
dp = data.export_data_as_pandas()
fig = pygmt.Figure()
fig.grdimage(grid=grid_file, cmap=etopo, projection="M15c", frame=True, region=this_region.format('str'))
pygmt.makecpt(cmap="viridis", series=[abs(dp.z).min(), abs(dp.z).max()])
fig.plot(
    x=dp.x,
    y=dp.y,
    size=0.08 * np.sqrt(abs(dp.z)),
    style="cc",
    fill=abs(dp.z),
    cmap=True,
    pen="black"
)
fig.colorbar(frame="xaf+lDifference (m)")
fig.savefig('{}_fig.png'.format(utils.fn_basename2(xyd_file)))
```