import os
import sys
sys.path.append('../src')
import glob

import xarray as xr
import xvec
import numpy as np
import pandas as pd
import geopandas as gpd

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors

# import custom_basemaps as cbm

import param_nwm3
import misc

nc_gridded_data = xr.open_dataset(os.path.join(param_nwm3.out_dir,'03_postproc_stats/nc','NWM_1981_2020_summary_AZ_HUC8.nc'))
nc_huc8_data = xr.open_dataset(os.path.join(param_nwm3.out_dir,'03_postproc_huc8_stats/nc','NWM_1981_2020_summary_VEC_HUC8.nc')).xvec.decode_cf()

ds_P = nc_huc8_data['P']
ds_ET = nc_huc8_data['ET']
ds_qSfcLatRunoff = nc_huc8_data['qSfcLatRunoff']
ds_gw_inflow = nc_huc8_data['gw_inflow']
ds_qBucket = nc_huc8_data['qBucket']
ds_qlateral = nc_huc8_data['q_lateral']

ds_WB = ds_P - ds_ET - ds_qlateral
ds_WB = ds_WB.rename('WB')


fig,ax = plt.subplots()
ds_qSfcLatRunoff.sel(season='WY').mean('geometry').plot(ax=ax,label='qSfcLatRunoff')
ds_gw_inflow.sel(season='WY').mean('geometry').plot(ax=ax,label='gw_inflow')
ds_qBucket.sel(season='WY').mean('geometry').plot(ax=ax,label='qBucket')
ds_qlateral.sel(season='WY').mean('geometry').plot(ax=ax,label='qlateral')
(ds_qBucket+ds_qSfcLatRunoff).sel(season='WY').mean('geometry').plot(ax=ax,label='ds_qBucket+ds_qSfcLatRunoff')
ax.legend()

# ds_WB = ds_WB.sel(season='WY').xvec.to_geodataframe(geometry='geometry')


