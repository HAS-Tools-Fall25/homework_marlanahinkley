import os
import glob

import xarray as xr
import numpy as np
import pandas as pd
import geopandas as gpd

import misc
import param_nwm3


def calc_summary(ds,start_wy,end_wy,seasons_dict,P_threshold=1):
    ds = ds.sel(time=slice(str(start_wy-1)+'-10-01',str(end_wy)+'-09-30'))

    # Resample by Water Year
    ds_WY = ds.resample(time='YE-SEP').sum('time')
    ds_WY = ds_WY.expand_dims('season')
    ds_WY = ds_WY.assign_coords(season=['WY'])
    ds_WY = ds_WY.sel(time=slice(str(start_wy),str(end_wy))) 

    # Resample by Season
    ds_seasons = []
    for season in seasons_dict.keys():
        ds_season = ds.where(ds['time.month'].isin(seasons_dict[season]), drop=True)
        ds_season = ds_season.resample(time='YE-SEP').sum('time')
        ds_season = ds_season.expand_dims('season')
        ds_season = ds_season.assign_coords(season=[season])
        ds_seasons.append(ds_season)
    ds_seasons = xr.merge(ds_seasons)

    ds_summary = xr.merge([ds_WY,ds_seasons])
    varname = list(ds_summary.data_vars)[0]
    ds_summary['{}_LTM'.format(varname)] = ds_summary[varname].mean('time')
    ds_summary['{}_LTSD'.format(varname)] = ds_summary[varname].std('time')
    ds_summary['{}_LTCV'.format(varname)] = ds_summary['{}_LTSD'.format(varname)]/ds_summary['{}_LTM'.format(varname)]

    return ds_summary

def calc_summary_mean(ds,start_wy,end_wy,seasons_dict,P_threshold=1):
    ds = ds.sel(time=slice(str(start_wy-1)+'-10-01',str(end_wy)+'-09-30'))

    # Resample by Water Year
    ds_WY = ds.resample(time='YE-SEP').mean('time')
    ds_WY = ds_WY.expand_dims('season')
    ds_WY = ds_WY.assign_coords(season=['WY'])
    ds_WY = ds_WY.sel(time=slice(str(start_wy),str(end_wy))) 

    # Resample by Season
    ds_seasons = []
    for season in seasons_dict.keys():
        ds_season = ds.where(ds['time.month'].isin(seasons_dict[season]), drop=True)
        ds_season = ds_season.resample(time='YE-SEP').mean('time')
        ds_season = ds_season.expand_dims('season')
        ds_season = ds_season.assign_coords(season=[season])
        ds_seasons.append(ds_season)
    ds_seasons = xr.merge(ds_seasons)

    ds_summary = xr.merge([ds_WY,ds_seasons])
    varname = list(ds_summary.data_vars)[0]
    ds_summary['{}_LTM'.format(varname)] = ds_summary[varname].mean('time')
    ds_summary['{}_LTSD'.format(varname)] = ds_summary[varname].std('time')
    ds_summary['{}_LTCV'.format(varname)] = ds_summary['{}_LTSD'.format(varname)]/ds_summary['{}_LTM'.format(varname)]

    return ds_summary

def calc_summary_max(ds,start_wy,end_wy,seasons_dict):
    ds = ds.sel(time=slice(str(start_wy-1)+'-10-01',str(end_wy)+'-09-30'))

    # Resample by Water Year
    ds_WY = ds.resample(time='YE-SEP').max('time')
    ds_WY = ds_WY.expand_dims('season')
    ds_WY = ds_WY.assign_coords(season=['WY'])
    ds_WY = ds_WY.sel(time=slice(str(start_wy),str(end_wy))) 

    # Resample by Season
    ds_seasons = []
    for season in seasons_dict.keys():
        ds_season = ds.where(ds['time.month'].isin(seasons_dict[season]), drop=True)
        ds_season = ds_season.resample(time='YE-SEP').max('time')
        ds_season = ds_season.expand_dims('season')
        ds_season = ds_season.assign_coords(season=[season])
        ds_seasons.append(ds_season)
    ds_seasons = xr.merge(ds_seasons)

    ds_summary = xr.merge([ds_WY,ds_seasons])
    varname = list(ds_summary.data_vars)[0]
    ds_summary['{}_LTM'.format(varname)] = ds_summary[varname].mean('time')
    ds_summary['{}_LTSD'.format(varname)] = ds_summary[varname].std('time')
    ds_summary['{}_LTCV'.format(varname)] = ds_summary['{}_LTSD'.format(varname)]/ds_summary['{}_LTM'.format(varname)]

    return ds_summary


def calc_summary_max_P_WET_DRY_days(ds,start_wy,end_wy,seasons_dict,P_type,P_threshold):
    varname = list(ds.data_vars)[0]
    
    if P_type == 'DRY':
        ds_P_type_days = ds.where(ds['P']<P_threshold)
    elif P_type == 'WET':
        ds_P_type_days = ds.where(ds['P']>=P_threshold)

    # Count by Water Year
    ds_P_type_days_WY = ds_P_type_days.resample(time='YE-SEP').count('time')
    ds_P_type_days_WY = ds_P_type_days_WY.expand_dims('season')
    ds_P_type_days_WY = ds_P_type_days_WY.assign_coords(season=['WY'])
    ds_P_type_days_WY = ds_P_type_days_WY.sel(time=slice(str(start_wy),str(end_wy)))
    ds_P_type_days_WY = ds_P_type_days_WY.astype('timedelta64[D]').astype('float32')

    # Count by Season
    ds_P_type_days_seasons = []
    for season in seasons_dict.keys():
        ds_P_type_days_season = ds_P_type_days.where(ds_P_type_days['time.month'].isin(seasons_dict[season]), drop=True)
        ds_P_type_days_season = ds_P_type_days_season.resample(time='YE-SEP').count('time')
        ds_P_type_days_season = ds_P_type_days_season.expand_dims('season')
        ds_P_type_days_season = ds_P_type_days_season.assign_coords(season=[season])
        ds_P_type_days_season = ds_P_type_days_season.astype('timedelta64[D]').astype('float32')
        ds_P_type_days_seasons.append(ds_P_type_days_season)
    ds_P_type_days_seasons = xr.merge(ds_P_type_days_seasons)

    ds_summary = xr.merge([ds_P_type_days_WY,ds_P_type_days_seasons])
    ds_summary = ds_summary.rename({varname:'{}_{}_DAYS'.format(varname,P_type)})
    varname = list(ds_summary.data_vars)[0]
    ds_summary['{}_LTM'.format(varname)] = ds_summary[varname].mean('time')
    ds_summary['{}_LTSD'.format(varname)] = ds_summary[varname].std('time')
    ds_summary['{}_LTCV'.format(varname)] = ds_summary['{}_LTSD'.format(varname)]/ds_summary['{}_LTM'.format(varname)]
    return ds_summary



start_wy = 1991#1981
end_wy = 2020

seasons_dict = {'JFM':[1,2,3],'AMJ':[4,5,6],'JAS':[7,8,9],'OND':[10,11,12]}


nwm3_retro_dir = '../inp/nwm3/retrospective/'
savedir = os.path.join('../out',os.path.basename(__file__).split('.')[0])
savedir = misc.makedir(os.path.join(savedir,'NWM_{}_{}'.format(start_wy,end_wy)))

# NetCDF Attribute

model_name = 'NWMv3.0 Retrospective'
domain_name = 'AZ-HUC8'
seasons_def = 'WY: Oct 1 - Sep 30\n'+\
               'JFM: Jan 1 - Mar 31\n'+\
                'AMJ: Apr 1 - Jun 30\n'+\
                'JAS: Jul 1 - Sep 30\n'+\
                'OND: Oct 1 - Dec 31\n'
resolution = '1 km'
ds_attrs = {'P':{'vars_attrs':{'P':{'long_name':'Precipitation',
                                   'units':'mm'},
                              'P_LTM':{'long_name':'Long-term mean precipitation',
                                       'units':'mm'},
                              'P_LTSD':{'long_name':'Long-term standard deviation precipitation',
                                         'units':'mm'},
                              'P_LTCV':{'long_name':'Long-term coefficient of variation precipitation',
                                         'units':'-'},
                              'P_maxD':{'long_name':'Maximum daily precipitation',
                                        'units':'mm'},
                              'P_maxD_LTM':{'long_name':'Long-term mean maximum daily precipitation',
                                            'units':'mm'},
                              'P_maxD_LTSD':{'long_name':'Long-term standard deviation maximum daily precipitation',
                                              'units':'mm'},
                              'P_maxD_LTCV':{'long_name':'Long-term coefficient of variation maximum daily precipitation',
                                                'units':'-'},   
                                'P_maxH':{'long_name':'Maximum hourly precipitation',
                                            'units':'mm'},
                                'P_maxH_LTM':{'long_name':'Long-term mean maximum hourly precipitation',
                                              'units':'mm'},
                                'P_maxH_LTSD':{'long_name':'Long-term standard deviation maximum hourly precipitation',
                                                'units':'mm'},
                                'P_maxH_LTCV':{'long_name':'Long-term coefficient of variation maximum hourly precipitation',
                                                  'units':'-'},
                                'P_WET_DAYS':{'long_name':'Number of wet days',
                                                'units':'days'},
                                'P_WET_DAYS_LTM':{'long_name':'Long-term mean number of wet days',
                                                  'units':'days'},
                                'P_WET_DAYS_LTSD':{'long_name':'Long-term standard deviation number of wet days',
                                                    'units':'days'},
                                'P_WET_DAYS_LTCV':{'long_name':'Long-term coefficient of variation number of wet days',
                                                      'units':'-'},
                                'P_DRY_DAYS':{'long_name':'Number of dry days',
                                                'units':'days'},
                                'P_DRY_DAYS_LTM':{'long_name':'Long-term mean number of dry days',
                                                  'units':'days'},
                                'P_DRY_DAYS_LTSD':{'long_name':'Long-term standard deviation number of dry days',
                                                    'units':'days'},
                                'P_DRY_DAYS_LTCV':{'long_name':'Long-term coefficient of variation number of dry days',
                                                      'units':'-'}}},
            'ET':{'vars_attrs':{'ET':{'long_name':'Evapotranspiration',
                                      'units':'mm'},
                                'ET_LTM':{'long_name':'Long-term mean evapotranspiration',
                                          'units':'mm'},
                                'ET_LTSD':{'long_name':'Long-term standard deviation evapotranspiration',
                                            'units':'mm'},
                                'ET_LTCV':{'long_name':'Long-term coefficient of variation evapotranspiration',
                                            'units':'-'}}},
            'T':{'vars_attrs':{'T':{'long_name':'Air temperature',
                                    'units':'C'},
                               'T_LTM':{'long_name':'Long-term mean air temperature',
                                    'units':'C'},
                               'T_LTSD':{'long_name':'Long-term standard deviation air temperature',
                                    'units':'C'},
                                'T_LTCV':{'long_name':'Long-term coefficient of variation air temperature',
                                    'units':'-'}}}, 
            'qSfcLatRunoff':{'vars_attrs':{'qSfcLatRunoff':{'long_name':'Surface lateral runoff',
                                                        'units':'mm'},
                                        'qSfcLatRunoff_LTM':{'long_name':'Long-term mean surface lateral runoff',
                                                            'units':'mm'},
                                        'qSfcLatRunoff_LTSD':{'long_name':'Long-term standard deviation surface lateral runoff',
                                                              'units':'mm'},
                                        'qSfcLatRunoff_LTCV':{'long_name':'Long-term coefficient of variation surface lateral runoff',
                                                              'units':'-'}}},
            'qBucket':{'vars_attrs':{'qBucket':{'long_name':'Bucket runoff',
                                                'units':'mm'},
                                    'qBucket_LTM':{'long_name':'Long-term mean bucket runoff',
                                                    'units':'mm'},
                                    'qBucket_LTSD':{'long_name':'Long-term standard deviation bucket runoff',
                                                      'units':'mm'},
                                    'qBucket_LTCV':{'long_name':'Long-term coefficient of variation bucket runoff',
                                                      'units':'-'}}},
            'q_lateral':{'vars_attrs':{'q_lateral':{'long_name':'Lateral runoff',
                                                'units':'mm'},
                                    'q_lateral_LTM':{'long_name':'Long-term mean lateral runoff',
                                                    'units':'mm'},
                                    'q_lateral_LTSD':{'long_name':'Long-term standard deviation lateral runoff',
                                                      'units':'mm'},
                                    'q_lateral_LTCV':{'long_name':'Long-term coefficient of variation lateral runoff',
                                                      'units':'-'}}},
            'gw_inflow':{'vars_attrs':{'gw_inflow':{'long_name':'Groundwater inflow',
                                                        'units':'mm'},
                                        'gw_inflow_LTM':{'long_name':'Long-term mean groundwater inflow',
                                                            'units':'mm'},
                                        'gw_inflow_LTSD':{'long_name':'Long-term standard deviation groundwater inflow',
                                                              'units':'mm'},
                                        'gw_inflow_LTCV':{'long_name':'Long-term coefficient of variation groundwater inflow',
                                                              'units':'-'}
                                   }},
   
            'model':model_name,
            'domain':domain_name,
            'resolution':resolution,
            'start_water_year':start_wy,
            'end_water_year':end_wy,
            'seasons':seasons_def}


# # P Summary
# P_threshold = 1
# print('Calculating P Summary')
# dir_P_daily_nc = os.path.join(nwm3_retro_dir,'forcing/daily/precip/*.nc')
# dir_P_monthly_nc = os.path.join(nwm3_retro_dir,'forcing/monthly/precip/*.nc')
# dir_P_maxD_nc = os.path.join(nwm3_retro_dir,'forcing/maxD/precip/*.nc')
# dir_P_maxH_nc = os.path.join(nwm3_retro_dir,'forcing/maxH/precip/*.nc')

# ds_P_daily = xr.open_mfdataset(sorted(glob.glob(dir_P_daily_nc))).rename({'RAINRATE':'P'})
# ds_P_monthly = xr.open_mfdataset(sorted(glob.glob(dir_P_monthly_nc))).rename({'RAINRATE':'P'})
# ds_P_maxD = xr.open_mfdataset(sorted(glob.glob(dir_P_maxD_nc))).rename({'RAINRATE':'P_maxD'})
# ds_P_maxH = xr.open_mfdataset(sorted(glob.glob(dir_P_maxH_nc))).rename({'RAINRATE':'P_maxH'})

# P_summary = calc_summary(ds_P_monthly,start_wy,end_wy,seasons_dict)
# P_summary_maxD = calc_summary_max(ds_P_maxD,start_wy,end_wy,seasons_dict)
# P_summary_maxH = calc_summary_max(ds_P_maxH,start_wy,end_wy,seasons_dict)
# P_summary_WET_days = calc_summary_max_P_WET_DRY_days(ds_P_daily,start_wy,end_wy,seasons_dict,'WET',P_threshold)
# P_summary_DRY_days = calc_summary_max_P_WET_DRY_days(ds_P_daily,start_wy,end_wy,seasons_dict,'DRY',P_threshold)
# P_summary = xr.merge([P_summary,P_summary_maxD,P_summary_maxH,P_summary_WET_days,P_summary_DRY_days])
# P_summary.attrs = {}
# P_summary.attrs.update(dict((k, ds_attrs[k]) for k in ('model','domain','resolution','start_water_year','end_water_year','seasons')))
# for var in P_summary.data_vars:
#     P_summary[var].attrs = {}
#     P_summary[var].attrs.update(ds_attrs['P']['vars_attrs'][var])
# P_summary.to_netcdf(os.path.join(savedir,'NWM_{}_{}_P_summary.nc'.format(start_wy,end_wy)))

# # ET Summary
# print('Calculating ET Summary')
# dir_ET_monthly_nc = os.path.join(nwm3_retro_dir,'ldasout/monthly/ET/*.nc')
# ds_ET_monthly = xr.open_mfdataset(dir_ET_monthly_nc)
# ET_summary = calc_summary(ds_ET_monthly,start_wy,end_wy,seasons_dict)
# ET_summary.attrs = {}
# ET_summary.attrs.update(dict((k, ds_attrs[k]) for k in ('model','domain','resolution','start_water_year','end_water_year','seasons')))
# for var in ET_summary.data_vars:
#     ET_summary[var].attrs = {}
#     ET_summary[var].attrs.update(ds_attrs['ET']['vars_attrs'][var])
# ET_summary.to_netcdf(os.path.join(savedir,'NWM_{}_{}_ET_summary.nc'.format(start_wy,end_wy)))

# # T Summary
# print('Calculating T Summary')
# dir_T_monthly_nc = os.path.join(nwm3_retro_dir,'forcing/monthly/t2d/*.nc')
# ds_T_monthly = xr.open_mfdataset(dir_T_monthly_nc)
# ds_T_monthly = ds_T_monthly.rename({'T2D':'T'})

# T_summary = calc_summary_mean(ds_T_monthly,start_wy,end_wy,seasons_dict)
# T_summary.attrs = {}
# T_summary.attrs.update(dict((k, ds_attrs[k]) for k in ('model','domain','resolution','start_water_year','end_water_year','seasons')))
# for var in T_summary.data_vars:
#     T_summary[var].attrs = {}
#     T_summary[var].attrs.update(ds_attrs['T']['vars_attrs'][var])
# T_summary.to_netcdf(os.path.join(savedir,'NWM_{}_{}_T_summary.nc'.format(start_wy,end_wy)))

# # qSfcLatRunoff Summary
# print('Calculating qSfcLatRunoff Summary')
# dir_qSfcLatRunoff_monthly_nc = os.path.join(nwm3_retro_dir,'chrtout_1km/monthly/qSfcLatRunoff/*.nc')
# ds_qSfcLatRunoff_monthly = xr.open_mfdataset(dir_qSfcLatRunoff_monthly_nc)
# qSfcLatRunoff_summary = calc_summary(ds_qSfcLatRunoff_monthly,start_wy,end_wy,seasons_dict)
# qSfcLatRunoff_summary.attrs = {}
# qSfcLatRunoff_summary.attrs.update(dict((k, ds_attrs[k]) for k in ('model','domain','resolution','start_water_year','end_water_year','seasons')))
# for var in qSfcLatRunoff_summary.data_vars:
#     qSfcLatRunoff_summary[var].attrs = {}
#     qSfcLatRunoff_summary[var].attrs.update(ds_attrs['qSfcLatRunoff']['vars_attrs'][var])
# qSfcLatRunoff_summary.to_netcdf(os.path.join(savedir,'NWM_{}_{}_qSfcLatRunoff_summary.nc'.format(start_wy,end_wy)))

# # qBucket Summary
# print('Calculating qBucket Summary')
# dir_qBucket_monthly_nc = os.path.join(nwm3_retro_dir,'chrtout_1km/monthly/qBucket/*.nc')
# ds_qBucket_monthly = xr.open_mfdataset(dir_qBucket_monthly_nc)
# qBucket_summary = calc_summary(ds_qBucket_monthly,start_wy,end_wy,seasons_dict)
# qBucket_summary.attrs = {}
# qBucket_summary.attrs.update(dict((k, ds_attrs[k]) for k in ('model','domain','resolution','start_water_year','end_water_year','seasons')))
# for var in qBucket_summary.data_vars:
#     qBucket_summary[var].attrs = {}
#     qBucket_summary[var].attrs.update(ds_attrs['qBucket']['vars_attrs'][var])
# qBucket_summary.to_netcdf(os.path.join(savedir,'NWM_{}_{}_qBucket_summary.nc'.format(start_wy,end_wy)))

# # q_lateral Summary
# print('Calculating q_lateral Summary')
# dir_q_lateral_monthly_nc = os.path.join(nwm3_retro_dir,'chrtout_1km/monthly/q_lateral/*.nc')
# ds_q_lateral_monthly = xr.open_mfdataset(dir_q_lateral_monthly_nc)
# q_lateral_summary = calc_summary(ds_q_lateral_monthly,start_wy,end_wy,seasons_dict)
# q_lateral_summary.attrs = {}
# q_lateral_summary.attrs.update(dict((k, ds_attrs[k]) for k in ('model','domain','resolution','start_water_year','end_water_year','seasons')))
# for var in q_lateral_summary.data_vars:
#     q_lateral_summary[var].attrs = {}
#     q_lateral_summary[var].attrs.update(ds_attrs['q_lateral']['vars_attrs'][var])
# q_lateral_summary.to_netcdf(os.path.join(savedir,'NWM_{}_{}_q_lateral_summary.nc'.format(start_wy,end_wy)))

# # gwout - inflow summary
# print('Calculating gwout inflow Summary')
# dir_gwout_inflow_monthly_nc = os.path.join(nwm3_retro_dir,'gwout_1km/monthly/inflow/*.nc')
# ds_gwout_inflow_monthly = xr.open_mfdataset(dir_gwout_inflow_monthly_nc)
# gwout_inflow_summary = calc_summary(ds_gwout_inflow_monthly,start_wy,end_wy,seasons_dict)
# for varname in gwout_inflow_summary.data_vars:
#     gwout_inflow_summary = gwout_inflow_summary.rename({varname:'gw_{}'.format(varname)})
# gwout_inflow_summary.attrs = {}
# gwout_inflow_summary.attrs.update(dict((k, ds_attrs[k]) for k in ('model','domain','resolution','start_water_year','end_water_year','seasons')))
# for var in gwout_inflow_summary.data_vars:
#     gwout_inflow_summary[var].attrs = {}
#     gwout_inflow_summary[var].attrs.update(ds_attrs['gw_inflow']['vars_attrs'][var])
# gwout_inflow_summary.to_netcdf(os.path.join(savedir,'NWM_{}_{}_gwout_inflow_summary.nc'.format(start_wy,end_wy)))


