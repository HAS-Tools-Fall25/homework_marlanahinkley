import os
import glob

import xarray as xr
import xvec #need to have this installed in our environment
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
misc.makedir(savedir)
savedir_shp = misc.makedir(os.path.join(savedir,'shp'))
misc.makedir(savedir_shp)
savedir_nc = misc.makedir(os.path.join(savedir,'nc'))
misc.makedir(savedir_nc)
savedir_csv = misc.makedir(os.path.join(savedir,'csv'))
misc.makedir(savedir_csv)
savedir_parquet = misc.makedir(os.path.join(savedir,'parquet'))
misc.makedir(savedir_parquet)

ds_crs = param_nwm3.crs_nwm_proj4_lcc
target_crs = param_nwm3.crs_utm12n_nad83

gdf_wbdhu8_boundary_az_lcc = gpd.read_file(param_nwm3.gp_wbdhu8_boundary_az_lcc)
gdf_wbdhu8_basins_az_lcc = gpd.read_file(param_nwm3.gp_wbdhu8_basins_az_lcc)
gdf_wbdhu8_basins_az_lcc.index = gdf_wbdhu8_basins_az_lcc['huc8']
gdf_wbdhu8_basins_az_utm12n_nad83 = gdf_wbdhu8_basins_az_lcc.to_crs(target_crs)

gdf_az_adwr_gwbasins_utm12n_nad83 = gpd.read_file('/workspaces/homework_marlanahinkley/homework_submissions/NWM_figures/Groundwater_Basin/Groundwater_Basin.shp') #need to change to full GW basin shapefile
gdf_az_adwr_gwbasins_utm12n_nad83 = gdf_az_adwr_gwbasins_utm12n_nad83[gdf_az_adwr_gwbasins_utm12n_nad83.geometry!=None]
gdf_az_adwr_gwbasins_lcc = gdf_az_adwr_gwbasins_utm12n_nad83.to_crs(ds_crs)

dir_stats = os.path.join(f'../out/02_calc_stats/NWM_{start_wy}_{end_wy}')


# with dask.config.set(**{'array.slicing.split_large_chunks': True}):

ds_stats = xr.open_mfdataset(sorted(glob.glob(os.path.join(dir_stats,'*.nc'))),decode_timedelta=False)

# Create ADWR GWbasin level statistics
ds_stats = ds_stats.rio.write_crs(ds_crs)
ds_stats_adwr_gwbasin = ds_stats.xvec.zonal_stats(gdf_az_adwr_gwbasins_lcc.geometry,x_coords='x',y_coords='y',stats='mean')
ds_stats_adwr_gwbasin_reproj = ds_stats_adwr_gwbasin.xvec.to_crs(geometry=target_crs) # reproject to UTM12N
ds_stats_adwr_gwbasin_encoded = ds_stats_adwr_gwbasin_reproj.xvec.encode_cf()
ds_stats_adwr_gwbasin_encoded.to_netcdf(os.path.join(savedir_nc,f'NWM_{start_wy}_{end_wy}_summary_VEC_GWBasins.nc'),mode='w')

# Create HUC8 level statistics Shapefile with deliverables
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
            'ET_LTM': 'ET_M',
            'ET_LTSD': 'ET_SD',
            'ET_LTCV': 'ET_CV',
            # 'T_LTM': 'T_M',
            # 'T_LTSD': 'T_SD',
            # 'T_LTCV': 'T_CV',
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


# Create  ADWR GWBasins statistics Shapefile with deliverables
gdf_stats_adwr_gwbasins_deliv = gdf_az_adwr_gwbasins_utm12n_nad83.copy()
for var in deliv_vars.keys():
    for season in deliv_seasons.keys():
        var_season = deliv_vars[var]+'_'+deliv_seasons[season]
        gdf_stats_adwr_gwbasins_deliv[var_season] = ds_stats_adwr_gwbasin_reproj[var].sel(season=season).values
gdf_stats_adwr_gwbasins_deliv.to_file(os.path.join(savedir_shp,f'NWM_{start_wy}_{end_wy}_summary_VEC_GWBasins.shp'))
# gdf_stats_adwr_gwbasins_deliv = gdf_stats_adwr_gwbasins_deliv.drop(columns=['basin_name'])

# # Create HUC8 level statistics Shapefile with deliverables
# gdf_stats_huc8_deliv = gdf_wbdhu8_basins_az_utm12n_nad83.copy()
# gdf_stats_huc8_deliv = gdf_stats_huc8_deliv.drop(columns=['huc8'])
# for var in deliv_vars.keys():
#     for season in deliv_seasons.keys():
#         var_season = deliv_vars[var]+'_'+deliv_seasons[season]
#         gdf_stats_huc8_deliv[var_season] = ds_stats_huc8_reproj[var].sel(season=season).values
# gdf_stats_huc8_deliv.to_file(os.path.join(savedir_shp,f'NWM_{start_wy}_{end_wy}_summary_VEC_HUC8.shp'))



