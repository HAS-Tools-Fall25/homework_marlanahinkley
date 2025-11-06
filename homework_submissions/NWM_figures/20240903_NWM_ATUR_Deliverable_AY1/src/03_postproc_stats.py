import os
import glob

import xarray as xr
import xvec
import rioxarray
import numpy as np
import pandas as pd
import geopandas as gpd
from rasterio.enums import Resampling
import dask

import misc
import param_nwm3

start_wy = 1981
end_wy = 2020

savedir = os.path.join('../out',os.path.basename(__file__).split('.')[0])
savedir = misc.makedir(os.path.join(savedir,f'NWM_{start_wy}_{end_wy}'))
savedir_nc = os.path.join(savedir,'nc')
savedir_nc = misc.makedir(savedir_nc)
savedir_tif = os.path.join(savedir,'tif')
savedir_tif = misc.makedir(savedir_tif)

gdf_wbdhu8_boundary_az_lcc = gpd.read_file(param_nwm3.gp_wbdhu8_boundary_az_lcc)
gdf_wbdhu8_basins_az_lcc = gpd.read_file(param_nwm3.gp_wbdhu8_basins_az_lcc)


ds_crs = param_nwm3.crs_nwm_proj4_lcc
target_crs = param_nwm3.crs_utm12n_nad83
dir_stats = os.path.join(f'../out/02_calc_stats/NWM_{start_wy}_{end_wy}')

with dask.config.set(**{'array.slicing.split_large_chunks': True}):

    ds_stats = xr.open_mfdataset(sorted(glob.glob(os.path.join(dir_stats,'*.nc'))),decode_timedelta=False)

    # Clip and Reproject 1 km NWM data to HUC8 boundary
    ds_stats = ds_stats.rio.set_spatial_dims(x_dim='x', y_dim='y')
    ds_stats = ds_stats.rio.write_crs(ds_crs)
    ds_stats_clipped = ds_stats.rio.clip(gdf_wbdhu8_boundary_az_lcc.geometry)
    dss_stats_reproj = []
    for season in ds_stats_clipped.season:
        ds_stats_reproj = ds_stats_clipped.sel(season=season).rio.reproject(target_crs,nodata=np.nan)
        dss_stats_reproj.append(ds_stats_reproj)
    ds_stats_reproj = xr.concat(dss_stats_reproj,dim='season')
    ds_stats_reproj = ds_stats_reproj.dropna(dim='x',how='all')
    ds_stats_reproj = ds_stats_reproj.dropna(dim='y',how='all')
    ds_stats_reproj.to_netcdf(os.path.join(savedir_nc,f'NWM_{start_wy}_{end_wy}_summary_AZ_HUC8.nc'))

# Create Deliverable Tiff files
# deliv_vars = {
#               'P_LTM': 'P_LTM',
#               'P_LTSD': 'P_LTSD',
#               'P_LTCV': 'P_LTCV',
#               'P_maxD_LTM': 'P_maxD_LTM',
#               'P_maxD_LTSD': 'P_maxD_LTSD',
#               'P_maxD_LTCV': 'P_maxD_LTCV',
#               'P_maxH_LTM': 'P_maxH_LTM',
#               'P_maxH_LTSD': 'P_maxH_LTSD',
#               'P_maxH_LTCV': 'P_maxH_LTCV',
#               'P_WET_DAYS_LTM': 'P_WET_DAYS_LTM',
#               'P_WET_DAYS_LTSD': 'P_WET_DAYS_LTSD',
#               'P_WET_DAYS_LTCV': 'P_WET_DAYS_LTCV',
#               'P_DRY_DAYS_LTM': 'P_DRY_DAYS_LTM',
#               'P_DRY_DAYS_LTSD': 'P_DRY_DAYS_LTSD',
#               'P_DRY_DAYS_LTCV': 'P_DRY_DAYS_LTCV',
#               'ET_LTM': 'ET_LTM',
#               'ET_LTSD': 'ET_LTSD',
#               'ET_LTCV': 'ET_LTCV',
#               'qSfcLatRunoff_LTM': 'SR_LTM',
#               'qSfcLatRunoff_LTSD': 'SR_LTSD',
#               'qSfcLatRunoff_LTCV': 'SR_LTCV',
#               'gw_inflow_LTM': 'Re_LTM',
#               'gw_inflow_LTSD': 'Re_LTSD',
#               'gw_inflow_LTCV': 'Re_LTCV',
#               }

deliv_vars = {
              'P_LTM': 'P_M',
              'P_LTSD': 'P_SD',
              'P_LTCV': 'P_CV',
              'P_maxD_LTM': 'P_mD_M',
              'P_maxD_LTSD': 'P_mD_SD',
              'P_maxD_LTCV': 'P_mD_CV',
              'P_maxH_LTM': 'P_mH_M',
              'P_maxH_LTSD': 'P_mH_SD',
              'P_maxH_LTCV': 'P_mH_CV',
              'P_WET_DAYS_LTM': 'P_WD_M',
              'P_WET_DAYS_LTSD': 'P_WD_SD',
              'P_WET_DAYS_LTCV': 'P_WD_CV',
              'P_DRY_DAYS_LTM': 'P_DD_M',
              'P_DRY_DAYS_LTSD': 'P_DD_SD',
              'P_DRY_DAYS_LTCV': 'P_DD_CV',
            #   'T_LTM': 'T_M',
            #   'T_LTSD': 'T_SD',
            #   'T_LTCV': 'T_CV',
              'ET_LTM': 'ET_M',
              'ET_LTSD': 'ET_SD',
              'ET_LTCV': 'ET_CV',
              'qSfcLatRunoff_LTM': 'SR_M',
              'qSfcLatRunoff_LTSD': 'SR_SD',
              'qSfcLatRunoff_LTCV': 'SR_CV',
              'q_lateral_LTM': 'TR_M',
              'q_lateral_LTSD': 'TR_SD',
              'q_lateral_LTCV': 'TR_CV',
              'gw_inflow_LTM': 'Re_M',
              'gw_inflow_LTSD': 'Re_SD',
              'gw_inflow_LTCV': 'Re_CV',
              }

deliv_seasons = {'WY':'WY',
                 'JFM':'Q1',
                 'AMJ':'Q2',
                 'JAS':'Q3',
                 'OND':'Q4'}

for deliv_season in deliv_seasons:
    for deliv_var in deliv_vars:
        var = deliv_vars[deliv_var]
        season = deliv_seasons[deliv_season]
        ds_stats_reproj[deliv_var].sel(season=deliv_season).rio.to_raster(os.path.join(savedir_tif,f'{var}_{season}.tif'))


# Create HUC8 level statistics
gdf_wbdhu8_basins_az_lcc.index = gdf_wbdhu8_basins_az_lcc['huc8']
# with dask.config.set(**{'array.slicing.split_large_chunks': True}):
ds_stats_huc8 = ds_stats.xvec.zonal_stats(gdf_wbdhu8_basins_az_lcc.geometry,x_coords='x',y_coords='y',stats='mean')




