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
from matplotlib.cm import ScalarMappable

import custom_basemaps as cbm

import param_nwm3
import misc

ds_nwm_1km_grid = xr.open_dataset(os.path.join(param_nwm3.out_dir,'03_postproc_stats/NWM_1991_2020/nc','NWM_1991_2020_summary_AZ_HUC8.nc'))
ds_nwm_huc8 = xr.open_dataset(os.path.join(param_nwm3.out_dir,'03_postproc_huc8_stats/NWM_1991_2020/nc','NWM_1991_2020_summary_VEC_HUC8.nc')).xvec.decode_cf()

savedir = os.path.join(param_nwm3.out_dir,'figures',os.path.basename(__file__).split('.')[0])
misc.makedir(savedir)

# P_LTM_WY
fig,ax = plt.subplots(1,1,subplot_kw={'projection':param_nwm3.cartopy_crs_atur_utm12n},figsize=(4,4))
ax=cbm.az_huc8_boundaries(ax=ax,gridlines_flag=False)
bounds = np.array([40,50,60,70,80,90,100,200,300,400,500,600,700,800,900,1000,1300])
norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)
cmap = 'turbo_r'
img = ds_nwm_1km_grid['P_LTM'].sel(season='WY').plot(ax=ax,norm=norm,cmap=cmap,add_colorbar=False)
cbar = plt.colorbar(img, ax=ax,shrink=0.8)
cbar.set_ticks(bounds)
cbar.set_label('Long-term mean\nprecipitation (mm)',fontsize=12)
ax.axis('off')
ax.set_title('P_M_WY',fontsize=14,weight='bold')
plt.tight_layout()
fig.savefig(os.path.join(savedir,'P_M_WY.png'),dpi=300,bbox_inches='tight')

# P_LTSD_WY
fig,ax = plt.subplots(1,1,subplot_kw={'projection':param_nwm3.cartopy_crs_atur_utm12n},figsize=(4,4))
ax=cbm.az_huc8_boundaries(ax=ax,gridlines_flag=False)
bounds = np.array([20,30,40,50,60,70,80,90,100,200,300,400])
norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)
cmap = 'turbo_r'
img = ds_nwm_1km_grid['P_LTSD'].sel(season='WY').plot(ax=ax,norm=norm,cmap=cmap,add_colorbar=False)
cbar = plt.colorbar(img, ax=ax,shrink=0.8)
cbar.set_ticks(bounds)
cbar.set_label('Long-term\nstandard deviation\nprecipitation (mm)',fontsize=12)
ax.axis('off')
ax.set_title('P_SD_WY',fontsize=14,weight='bold')
plt.tight_layout()
fig.savefig(os.path.join(savedir,'P_SD_WY.png'),dpi=300,bbox_inches='tight')

# P_LTCV_WY
fig,ax = plt.subplots(1,1,subplot_kw={'projection':param_nwm3.cartopy_crs_atur_utm12n},figsize=(4,4))
ax=cbm.az_huc8_boundaries(ax=ax,gridlines_flag=False)
bounds = np.array([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1])
norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)
cmap = 'turbo_r'
img = ds_nwm_1km_grid['P_LTCV'].sel(season='WY').plot(ax=ax,norm=norm,cmap=cmap,add_colorbar=False)
cbar = plt.colorbar(img, ax=ax,shrink=0.8)
cbar.set_ticks(bounds)
cbar.set_label('Long-term\ncoefficient of variation\n precipitation (-)',fontsize=12)
ax.axis('off')
ax.set_title('P_CV_WY',fontsize=14,weight='bold')
plt.tight_layout()
fig.savefig(os.path.join(savedir,'P_CV_WY.png'),dpi=300,bbox_inches='tight')

# P_maxD_LTM_WY
fig,ax = plt.subplots(1,1,subplot_kw={'projection':param_nwm3.cartopy_crs_atur_utm12n},figsize=(4,4))
ax=cbm.az_huc8_boundaries(ax=ax,gridlines_flag=False)
bounds = np.array([0,10,20,30,40,50,60,70,80,90,100])
norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)
cmap = 'turbo_r'
img = ds_nwm_1km_grid['P_maxD_LTM'].sel(season='WY').plot(ax=ax,norm=norm,cmap=cmap,add_colorbar=False)
cbar = plt.colorbar(img, ax=ax,shrink=0.8)
cbar.set_ticks(bounds)
cbar.set_label('Long-term mean\nmaximum daily\n precipitation (mm)',fontsize=12)
ax.axis('off')
ax.set_title('P_mD_M_WY',fontsize=14,weight='bold')
plt.tight_layout()
fig.savefig(os.path.join(savedir,'P_mD_M_WY.png'),dpi=300,bbox_inches='tight')

# P_maxH_LTM_WY
fig,ax = plt.subplots(1,1,subplot_kw={'projection':param_nwm3.cartopy_crs_atur_utm12n},figsize=(4,4))
ax=cbm.az_huc8_boundaries(ax=ax,gridlines_flag=False)
bounds = np.array([3,4,5,6,7,8,9,10,11,12,13,14,15,20,25,30])
norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)
cmap = 'turbo_r'
img = ds_nwm_1km_grid['P_maxH_LTM'].sel(season='WY').plot(ax=ax,norm=norm,cmap=cmap,add_colorbar=False)
cbar = plt.colorbar(img, ax=ax,shrink=0.8)
cbar.set_ticks(bounds)
cbar.set_label('Long-term mean\nmaximum hourly\n precipitation (mm)',fontsize=12)
ax.axis('off')
ax.set_title('P_mH_M_WY',fontsize=14,weight='bold')
plt.tight_layout()
fig.savefig(os.path.join(savedir,'P_mH_M_WY.png'),dpi=300,bbox_inches='tight')

# P_WET_DAYS_LTM_WY
fig,ax = plt.subplots(1,1,subplot_kw={'projection':param_nwm3.cartopy_crs_atur_utm12n},figsize=(4,4))
ax=cbm.az_huc8_boundaries(ax=ax,gridlines_flag=False)
bounds = np.array([10,20,30,40,50,60,70,80,90,100,110,120,130])
norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)
cmap = 'turbo_r'
img = ds_nwm_1km_grid['P_WET_DAYS_LTM'].sel(season='WY').dt.days.plot(ax=ax,norm=norm,cmap=cmap,add_colorbar=False)
cbar = plt.colorbar(img, ax=ax,shrink=0.8)
cbar.set_ticks(bounds)
cbar.set_label('Long-term mean\nnumber of wet days (days)',fontsize=12)
ax.axis('off')
ax.set_title('P_WD_M_WY',fontsize=14,weight='bold')
plt.tight_layout()
fig.savefig(os.path.join(savedir,'P_WD_M_WY.png'),dpi=300,bbox_inches='tight')

# P_DRY_DAYS_LTM_WY
fig,ax = plt.subplots(1,1,subplot_kw={'projection':param_nwm3.cartopy_crs_atur_utm12n},figsize=(4,4))
ax=cbm.az_huc8_boundaries(ax=ax,gridlines_flag=False)
bounds = np.array([230,240,250,260,270,280,290,300,310,320,330,340,350,360])
norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)
cmap = 'turbo'
img = ds_nwm_1km_grid['P_DRY_DAYS_LTM'].sel(season='WY').dt.days.plot(ax=ax,norm=norm,cmap=cmap,add_colorbar=False)
cbar = plt.colorbar(img, ax=ax,shrink=0.8)
cbar.set_ticks(bounds)
cbar.set_label('Long-term mean\nnumber of dry days (days)',fontsize=12)
ax.axis('off')
ax.set_title('P_DD_M_WY',fontsize=14,weight='bold')
plt.tight_layout()
fig.savefig(os.path.join(savedir,'P_DD_M_WY.png'),dpi=300,bbox_inches='tight')

# P_LTM_JFM
fig,ax = plt.subplots(1,1,subplot_kw={'projection':param_nwm3.cartopy_crs_atur_utm12n},figsize=(4,4))
ax=cbm.az_huc8_boundaries(ax=ax,gridlines_flag=False)
bounds = np.array([0,10,20,30,40,50,60,70,80,90,100,150,200,300,400,500,600])
norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)
cmap = 'turbo_r'
img = ds_nwm_1km_grid['P_LTM'].sel(season='JFM').plot(ax=ax,norm=norm,cmap=cmap,add_colorbar=False)
cbar = plt.colorbar(img, ax=ax,shrink=0.8)
cbar.set_ticks(bounds)
cbar.set_label('Long-term mean\nprecipitation (mm)',fontsize=12)
ax.axis('off')
ax.set_title('P_M_Q1',fontsize=14,weight='bold')
plt.tight_layout()
fig.savefig(os.path.join(savedir,'P_M_Q1.png'),dpi=300,bbox_inches='tight')

# P_LTM_AMJ
fig,ax = plt.subplots(1,1,subplot_kw={'projection':param_nwm3.cartopy_crs_atur_utm12n},figsize=(4,4))
ax=cbm.az_huc8_boundaries(ax=ax,gridlines_flag=False)
bounds = np.array([0,10,20,30,40,50,60,70,80,90,100,150,200,300,400,500,600])
norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)
cmap = 'turbo_r'
img = ds_nwm_1km_grid['P_LTM'].sel(season='AMJ').plot(ax=ax,norm=norm,cmap=cmap,add_colorbar=False)
cbar = plt.colorbar(img, ax=ax,shrink=0.8)
cbar.set_ticks(bounds)
cbar.set_label('Long-term mean\nprecipitation (mm)',fontsize=12)
ax.axis('off')
ax.set_title('P_M_Q2',fontsize=14,weight='bold')
plt.tight_layout()
fig.savefig(os.path.join(savedir,'P_M_Q2.png'),dpi=300,bbox_inches='tight')

# P_LTM_JAS
fig,ax = plt.subplots(1,1,subplot_kw={'projection':param_nwm3.cartopy_crs_atur_utm12n},figsize=(4,4))
ax=cbm.az_huc8_boundaries(ax=ax,gridlines_flag=False)
bounds = np.array([0,10,20,30,40,50,60,70,80,90,100,150,200,300,400,500,600])
norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)
cmap = 'turbo_r'
img = ds_nwm_1km_grid['P_LTM'].sel(season='JAS').plot(ax=ax,norm=norm,cmap=cmap,add_colorbar=False)
cbar = plt.colorbar(img, ax=ax,shrink=0.8)
cbar.set_ticks(bounds)
cbar.set_label('Long-term mean\nprecipitation (mm)',fontsize=12)
ax.axis('off')
ax.set_title('P_M_Q3',fontsize=14,weight='bold')
plt.tight_layout()
fig.savefig(os.path.join(savedir,'P_M_Q3.png'),dpi=300,bbox_inches='tight')

# P_LTM_OND
fig,ax = plt.subplots(1,1,subplot_kw={'projection':param_nwm3.cartopy_crs_atur_utm12n},figsize=(4,4))
ax=cbm.az_huc8_boundaries(ax=ax,gridlines_flag=False)
bounds = np.array([0,10,20,30,40,50,60,70,80,90,100,150,200,300,400,500,600])
norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)
cmap = 'turbo_r'
img = ds_nwm_1km_grid['P_LTM'].sel(season='OND').plot(ax=ax,norm=norm,cmap=cmap,add_colorbar=False)
cbar = plt.colorbar(img, ax=ax,shrink=0.8)
cbar.set_ticks(bounds)
cbar.set_label('Long-term mean\nprecipitation (mm)',fontsize=12)
ax.axis('off')
ax.set_title('P_M_Q4',fontsize=14,weight='bold')
plt.tight_layout()
fig.savefig(os.path.join(savedir,'P_M_Q4.png'),dpi=300,bbox_inches='tight')


# ET_LTM_WY
fig,ax = plt.subplots(1,1,subplot_kw={'projection':param_nwm3.cartopy_crs_atur_utm12n},figsize=(4,4))
ax=cbm.az_huc8_boundaries(ax=ax,gridlines_flag=False)
bounds = np.array([0,10,20,30,40,50,60,70,80,90,100,200,300,400,500,600,700,800,900,1000,1200])
norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)
cmap = 'turbo_r'
img = ds_nwm_1km_grid['ET_LTM'].sel(season='WY').plot(ax=ax,norm=norm,cmap=cmap,add_colorbar=False)
cbar = plt.colorbar(img, ax=ax,shrink=0.8)
cbar.set_ticks(bounds)
cbar.set_label('Long-term mean\nevapotranspiration (mm)',fontsize=12)
ax.axis('off')
ax.set_title('ET_M_WY',fontsize=14,weight='bold')
plt.tight_layout()
fig.savefig(os.path.join(savedir,'ET_M_WY.png'),dpi=300,bbox_inches='tight')

# SR_LTM_WY
fig,ax = plt.subplots(1,1,subplot_kw={'projection':param_nwm3.cartopy_crs_atur_utm12n},figsize=(4,4))
ax=cbm.az_huc8_boundaries(ax=ax,gridlines_flag=False)
bounds = np.array([0,1,2,3,4,5,10,20,30,40,50,60,70,80,90,100,200,300,400,500])
norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)
cmap = 'turbo_r'
img = ds_nwm_1km_grid['qSfcLatRunoff_LTM'].sel(season='WY').plot(ax=ax,norm=norm,cmap=cmap,add_colorbar=False)
cbar = plt.colorbar(img, ax=ax,shrink=0.8,extend='max')
cbar.set_ticks(bounds)
cbar.set_label('Long-term mean\nsurface runoff (mm)',fontsize=12)
ax.axis('off')
ax.set_title('SR_M_WY',fontsize=14,weight='bold')
plt.tight_layout()
fig.savefig(os.path.join(savedir,'SR_M_WY.png'),dpi=300,bbox_inches='tight')

# Re_LTM_WY
fig,ax = plt.subplots(1,1,subplot_kw={'projection':param_nwm3.cartopy_crs_atur_utm12n},figsize=(4,4))
ax=cbm.az_huc8_boundaries(ax=ax,gridlines_flag=False)
bounds = np.array([0,1,2,3,4,5,10,20,30,40,50,60,70,80,90,100,200,300,400,500,600,700,800])
norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)
cmap = 'turbo_r'
img = ds_nwm_1km_grid['gw_inflow_LTM'].sel(season='WY').plot(ax=ax,norm=norm,cmap=cmap,add_colorbar=False)
cbar = plt.colorbar(img, ax=ax,shrink=0.8)
cbar.set_ticks(bounds)
cbar.set_label('Long-term mean\ngroundwater recharge (mm)',fontsize=12)
ax.axis('off')
ax.set_title('Re_M_WY',fontsize=14,weight='bold')
plt.tight_layout()
fig.savefig(os.path.join(savedir,'Re_M_WY.png'),dpi=300,bbox_inches='tight')

# P_LTM_WY (HUC8)
fig,ax = plt.subplots(1,1,subplot_kw={'projection':param_nwm3.cartopy_crs_atur_utm12n},figsize=(4,4))
ax=cbm.az_huc8_boundaries(ax=ax,gridlines_flag=False)
bounds = np.array([60,70,80,90,100,150,200,250,300,350,400,450,500,600,700])
norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)
cmap = 'turbo_r'
img = ds_nwm_huc8['P_LTM'].sel(season='WY').xvec.to_geodataframe(geometry='geometry').plot(ax=ax,column='P_LTM',norm=norm,cmap=cmap,legend=False)

# Create a ScalarMappable with the norm and cmap
sm = ScalarMappable(norm=norm, cmap=cmap)
# The ScalarMappable does not have any data associated with it, so we give it some
sm.set_array([])
cbar = plt.colorbar(sm, ax=ax, shrink=0.8, boundaries=bounds)
cbar.set_ticks(bounds)
cbar.set_label('Long-term mean\nprecipitation (mm)',fontsize=12)
ax.axis('off')
ax.set_title('P_M_WY [HUC8]',fontsize=14,weight='bold')
plt.tight_layout()
fig.savefig(os.path.join(savedir,'P_M_WY_HUC8.png'),dpi=300,bbox_inches='tight')
# cbar = plt.colorbar(img, ax=ax,shrink=0.8)
# cbar.set_ticks(bounds)
# cbar.set_label('Long-term mean\nprecipitation (mm)',fontsize=12)
# ax.axis('off')
# ax.set_title('P_LTM_WY [HUC8]',fontsize=14,weight='bold')

# P_Re_WY (HUC8)
fig,ax = plt.subplots(1,1,subplot_kw={'projection':param_nwm3.cartopy_crs_atur_utm12n},figsize=(4,4))
ax=cbm.az_huc8_boundaries(ax=ax,gridlines_flag=False)
bounds = np.array([0,1,2,3,4,5,10,25,50,100,150,200])
norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)
cmap = 'turbo_r'
img = ds_nwm_huc8['gw_inflow_LTM'].sel(season='WY').xvec.to_geodataframe(geometry='geometry').plot(ax=ax,column='gw_inflow_LTM',norm=norm,cmap=cmap,legend=False)

# Create a ScalarMappable with the norm and cmap
sm = ScalarMappable(norm=norm, cmap=cmap)
# The ScalarMappable does not have any data associated with it, so we give it some
sm.set_array([])
cbar = plt.colorbar(sm, ax=ax, shrink=0.8, boundaries=bounds)
cbar.set_ticks(bounds)
cbar.set_label('Long-term mean\nrecharge (mm)',fontsize=12)
ax.axis('off')
ax.set_title('Re_M_WY [HUC8]',fontsize=14,weight='bold')
plt.tight_layout()
fig.savefig(os.path.join(savedir,'Re_M_WY_HUC8.png'),dpi=300,bbox_inches='tight')
# cbar = plt.colorbar(img, ax=ax,shrink=0.8)
# cbar.set_ticks(bounds)
# cbar.set_label('Long-term mean\nprecipitation (mm)',fontsize=12)
# ax.axis('off')
# ax.set_title('P_LTM_WY [HUC8]',fontsize=14,weight='bold')



