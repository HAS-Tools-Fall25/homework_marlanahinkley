import os
import glob
import time
import pickle
import multiprocessing as mp

import xarray as xr
import numpy as np
import pandas as pd
import geopandas as gpd

import misc
import param_nwm3

import matplotlib.colors as colors

def read_spatialweights(comids,dir_spatialweights):
    print('Reading sptialweights...')
    dict_spatialweights = {}
    i_comid = 0
    for comid in comids:
        print(i_comid,comid)
        comid = int(comid)
        ds_spatialweights = xr.open_dataset(os.path.join(dir_spatialweights,'{}_spatialweights.nc'.format(comid)))
        df_spatialweights = ds_spatialweights.to_dataframe().reset_index().drop(columns=['feature_id'])
        df_spatialweights = df_spatialweights[df_spatialweights['weight'].notna()]
        df_spatialweights = df_spatialweights.reset_index(drop=True)
        dict_spatialweights[comid] = df_spatialweights
        ds_spatialweights.close()
        i_comid += 1
    return dict_spatialweights

def nwm_catchment_to_grid(ds,ds_grid_template,dict_spatialweights,target_resolution):
# Prepare Template Grid
    varname = list(ds_grid_template.data_vars)[0]
    ds_grid_template = ds_grid_template.rename({varname:'VAR'})
    ds_grid_template = xr.zeros_like(ds_grid_template['VAR']).isel(time=0)
    ds_grid_template.attrs = {}
    ds_grid_template = ds_grid_template.drop_vars('time')

    # NWM to Grid
    varname = list(ds.data_vars)[0]
    i_comid = 0
    for comid in ds['feature_id'].values:
        comid = int(comid)
        var_comid_value = ds[varname].sel(feature_id=comid).values[0]
        if var_comid_value != 0:
            # print(i_comid,comid,var_comid_value)
            df_spatialweights = dict_spatialweights[comid]
            for i in df_spatialweights.index:
                ds_grid_template.loc[dict(x=df_spatialweights['x'][i],y=df_spatialweights['y'][i])] += df_spatialweights['weight'][i]*var_comid_value
        i_comid += 1

    ds = (ds_grid_template/(target_resolution*target_resolution))*1000 # m3/M --> m/M --> mm/M 
    return ds

def proc_nwm_catchment_to_grid_mp(nc_file,ds_grid_template,dict_spatialweights,target_resolution,savedir):
    start_ET = time.time()
    # print(nc_file.split('/')[-1].split('.')[0])
    ds = xr.open_dataset(nc_file)
    varname = list(ds.data_vars)[0]
    ds_1km = nwm_catchment_to_grid(ds,ds_grid_template,dict_spatialweights,target_resolution)
    ds_1km = ds_1km.expand_dims('time')
    ds_1km = ds_1km.assign_coords(time=pd.to_datetime(ds['time'].values))
    ds_1km = ds_1km.to_dataset()
    ds_1km = ds_1km.rename({'VAR':varname})
    ds_1km[varname].attrs['units'] = 'mm M-1'
    ds_1km[varname].attrs['long_name'] = 'Monthly Lateral Surface Runoff'
    encoding = {varname: {'zlib': True, 'complevel': 5}}
    savedir_var = os.path.join(savedir,varname)
    misc.makedir(savedir_var)
    ds_1km.to_netcdf(os.path.join(savedir_var,nc_file.split('/')[-1]),encoding=encoding)
    end_ET = time.time()
    print(varname,nc_file.split('/')[-1].split('.')[0],'Elapsed Time: ',end_ET-start_ET)


def reindex_spatialweight(ds_spatialweight,ds_grid_template):
    ds_spatialweight_reindex = ds_spatialweight.reindex({'x': ds_grid_template['x'], 'y': ds_grid_template['y']}, method=None)
    return ds_spatialweight_reindex

if __name__ == '__main__':
    dir_nwm3_retro = '../inp/nwm3/retrospective/'
    dir_spatialweights_1km = '../inp/nwm3/spatial_weights_1km'
    target_resolution = 1000 # m
    gdf_nwm_catchments = gpd.read_parquet(param_nwm3.gp_nwm_catchments_az_lcc)
    ds_grid_template_1km = xr.open_dataset(os.path.join(dir_nwm3_retro,'forcing/monthly/precip/198010.nc')) # Any file with the same domain
    savedir = os.path.join('../out',os.path.basename(__file__).split('.')[0])
    savedir = misc.makedir(savedir)

    # Read Spatial Weights
    # print('Reading spatial weights...')
    pkl_file = os.path.join(savedir,'dict_spatialweights.pkl')

    # Check if the pickle file exists
    if os.path.exists(pkl_file):
        # If the pickle file exists, load it
        with open(pkl_file, 'rb') as f:
            dict_spatialweights = pickle.load(f)
    else:
        # If the pickle file does not exist, generate dict_spatialweights and save it to a pickle file
        dict_spatialweights = read_spatialweights(gdf_nwm_catchments['ID'].values,dir_spatialweights_1km)
        with open(pkl_file, 'wb') as f:
            pickle.dump(dict_spatialweights, f)

    # NWM Catchment to Grid


    # # qSfcLatRunoff
    # print('NWM Catchment to Grid...')
    # dir_nc = os.path.join(dir_nwm3_retro,'chrtout/monthly/qSfcLatRunoff')
    # nc_files = sorted(glob.glob(os.path.join(dir_nc,'*.nc')))

    # arg_list = []
    # for nc_file in nc_files:
    #     arg = (nc_file,ds_grid_template_1km,dict_spatialweights,target_resolution,savedir)
    #     arg_list.append(arg)
    # with mp.Pool(mp.cpu_count()) as pool:
    #     pool.starmap(proc_nwm_catchment_to_grid_mp,arg_list)

    # # qBucket
    # print('NWM Catchment to Grid...')
    # dir_nc = os.path.join(dir_nwm3_retro,'chrtout/monthly/qBucket')
    # nc_files = sorted(glob.glob(os.path.join(dir_nc,'*.nc')))

    # arg_list = []
    # for nc_file in nc_files:
    #     arg = (nc_file,ds_grid_template_1km,dict_spatialweights,target_resolution,savedir)
    #     arg_list.append(arg)
    # with mp.Pool(mp.cpu_count()) as pool:
    #     pool.starmap(proc_nwm_catchment_to_grid_mp,arg_list)

    # # q_lateral
    # print('NWM Catchment to Grid...')
    # dir_nc = os.path.join(dir_nwm3_retro,'chrtout/monthly/q_lateral')
    # nc_files = sorted(glob.glob(os.path.join(dir_nc,'*.nc')))

    # arg_list = []
    # for nc_file in nc_files:
    #     arg = (nc_file,ds_grid_template_1km,dict_spatialweights,target_resolution,savedir)
    #     arg_list.append(arg)
    # with mp.Pool(mp.cpu_count()) as pool:
    #     pool.starmap(proc_nwm_catchment_to_grid_mp,arg_list)


    # gwout - inflow
    print('NWM Catchment to Grid...')
    dir_nc = os.path.join(dir_nwm3_retro,'gwout/monthly/inflow')
    nc_files = sorted(glob.glob(os.path.join(dir_nc,'*.nc')))

    # arg_list = []
    # for nc_file in nc_files:
    #     arg = (nc_file,ds_grid_template_1km,dict_spatialweights,target_resolution,savedir)
    #     arg_list.append(arg)
    # with mp.Pool(mp.cpu_count()) as pool:
    #     pool.starmap(proc_nwm_catchment_to_grid_mp,arg_list)

    for nc_file in nc_files[-40:]:
        proc_nwm_catchment_to_grid_mp(nc_file,ds_grid_template_1km,dict_spatialweights,target_resolution,savedir)


