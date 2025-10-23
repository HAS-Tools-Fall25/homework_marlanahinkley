import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from dask.diagnostics import ProgressBar

gpm_file = 'https://ncsa.osn.xsede.org/Pangeo/pangeo-forge/gpcp-feedstock/gpcp.zarr'
gpm = xr.open_dataset(gpm_file, engine='zarr', chunks={})
gpm = gpm['precip']

sst_file = 'https://ncsa.osn.xsede.org/Pangeo/pangeo-forge/HadISST-feedstock/hadisst.zarr'
sst = xr.open_dataset(sst_file, engine='zarr', chunks={})
sst = sst['sst']
gpm

print(sst)