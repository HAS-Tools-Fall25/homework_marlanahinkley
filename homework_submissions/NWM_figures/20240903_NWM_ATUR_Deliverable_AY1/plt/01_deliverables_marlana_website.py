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

ds_nwm_1km_grid = xr.open_dataset(os.path.join(param_nwm3.out_dir,'03_postproc_stats/NWM_1981_2020/nc','NWM_1981_2020_summary_AZ_HUC8.nc'))
ds_nwm_huc8 = xr.open_dataset(os.path.join(param_nwm3.out_dir,'03_postproc_huc8_stats/NWM_1981_2020/nc','NWM_1981_2020_summary_VEC_HUC8.nc')).xvec.decode_cf()

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
# cbar.set_label('Long-term mean\nprecipitation (mm)',fontsize=12)
cbar.set_label('Precipitation (mm)',fontsize=12)
ax.axis('off')
# ax.set_title('P_M_WY',fontsize=14,weight='bold')
ax.set_title('Long-term mean (1981-2020) \nannual precipitation',fontsize=14)
plt.tight_layout()
fig.savefig(os.path.join(savedir,'P_M_WY.png'),dpi=600,bbox_inches='tight')
fig.savefig(os.path.join(savedir,'P_M_WY.pdf'),dpi=600,bbox_inches='tight')

# ET_LTM_WY
fig,ax = plt.subplots(1,1,subplot_kw={'projection':param_nwm3.cartopy_crs_atur_utm12n},figsize=(4,4))
ax=cbm.az_huc8_boundaries(ax=ax,gridlines_flag=False)
bounds = np.array([0,10,20,30,40,50,60,70,80,90,100,200,300,400,500,600,700,800,900,1000,1200])
norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)
cmap = 'turbo'
img = ds_nwm_1km_grid['ET_LTM'].sel(season='WY').plot(ax=ax,norm=norm,cmap=cmap,add_colorbar=False)
cbar = plt.colorbar(img, ax=ax,shrink=0.8)
cbar.set_ticks(bounds)
# cbar.set_label('Long-term mean\nevapotranspiration (mm)',fontsize=12)
cbar.set_label('Evapotranspiration (mm)',fontsize=12)
ax.axis('off')
ax.set_title('ET_M_WY',fontsize=14,weight='bold')
ax.set_title('Long-term mean (1981-2020) \nannual evapotranspiration',fontsize=14)
plt.tight_layout()
fig.savefig(os.path.join(savedir,'ET_M_WY.png'),dpi=600,bbox_inches='tight')
fig.savefig(os.path.join(savedir,'ET_M_WY.pdf'),dpi=600,bbox_inches='tight')

# SR_LTM_WY
fig,ax = plt.subplots(1,1,subplot_kw={'projection':param_nwm3.cartopy_crs_atur_utm12n},figsize=(4,4))
ax=cbm.az_huc8_boundaries(ax=ax,gridlines_flag=False)
bounds = np.array([0,1,2,3,4,5,10,20,30,40,50,60,70,80,90,100,200,300,400,500])
norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)
cmap = 'turbo_r'
img = ds_nwm_1km_grid['qSfcLatRunoff_LTM'].sel(season='WY').plot(ax=ax,norm=norm,cmap=cmap,add_colorbar=False)
cbar = plt.colorbar(img, ax=ax,shrink=0.8,extend='max')
cbar.set_ticks(bounds)
cbar.set_label('Surface runoff (mm)',fontsize=12)
ax.axis('off')
ax.set_title('Long-term mean (1981-2020) \nannual surface runoff',fontsize=14)
plt.tight_layout()
fig.savefig(os.path.join(savedir,'SR_M_WY.png'),dpi=300,bbox_inches='tight')
fig.savefig(os.path.join(savedir,'SR_M_WY.pdf'),dpi=600,bbox_inches='tight')

# Re_LTM_WY
fig,ax = plt.subplots(1,1,subplot_kw={'projection':param_nwm3.cartopy_crs_atur_utm12n},figsize=(4,4))
ax=cbm.az_huc8_boundaries(ax=ax,gridlines_flag=False)
bounds = np.array([0,1,2,3,4,5,10,20,30,40,50,60,70,80,90,100,200,300,400,500,600,700,800])
norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)
cmap = 'turbo_r'
img = ds_nwm_1km_grid['gw_inflow_LTM'].sel(season='WY').plot(ax=ax,norm=norm,cmap=cmap,add_colorbar=False)
cbar = plt.colorbar(img, ax=ax,shrink=0.8)
cbar.set_ticks(bounds)
cbar.set_label('Groundwater recharge (mm)',fontsize=12)
ax.axis('off')
ax.set_title('Long-term mean (1981-2020) \nannual groundwater recharge',fontsize=14)
plt.tight_layout()
fig.savefig(os.path.join(savedir,'Re_M_WY.png'),dpi=300,bbox_inches='tight')
fig.savefig(os.path.join(savedir,'Re_M_WY.pdf'),dpi=600,bbox_inches='tight')

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
cbar.set_label('Precipitation (mm)',fontsize=12)
ax.axis('off')
ax.set_title('Long-term mean (1981-2020) \nannual precipitation (HUC8)',fontsize=14)
plt.tight_layout()
fig.savefig(os.path.join(savedir,'P_M_WY_HUC8.png'),dpi=300,bbox_inches='tight')
fig.savefig(os.path.join(savedir,'P_M_WY_HUC8.pdf'),dpi=300,bbox_inches='tight')

# cbar = plt.colorbar(img, ax=ax,shrink=0.8)
# cbar.set_ticks(bounds)
# cbar.set_label('Long-term mean\nprecipitation (mm)',fontsize=12)
# ax.axis('off')
# ax.set_title('P_LTM_WY [HUC8]',fontsize=14,weight='bold')

# ET_LTM_WY (HUC8)
fig,ax = plt.subplots(1,1,subplot_kw={'projection':param_nwm3.cartopy_crs_atur_utm12n},figsize=(4,4))
ax=cbm.az_huc8_boundaries(ax=ax,gridlines_flag=False)
bounds = np.array([60,70,80,90,100,150,200,250,300,350,400,450,500,600])
norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)
cmap = 'turbo_r'
img = ds_nwm_huc8['ET_LTM'].sel(season='WY').xvec.to_geodataframe(geometry='geometry').plot(ax=ax,column='ET_LTM',norm=norm,cmap=cmap,legend=False)

# Create a ScalarMappable with the norm and cmap
sm = ScalarMappable(norm=norm, cmap=cmap)
# The ScalarMappable does not have any data associated with it, so we give it some
sm.set_array([])
cbar = plt.colorbar(sm, ax=ax, shrink=0.8, boundaries=bounds)
cbar.set_ticks(bounds)
cbar.set_label('Evapotranspiration (mm)',fontsize=12)
ax.axis('off')
ax.set_title('Long-term mean (1981-2020) \nannual evapotranspiration (HUC8)',fontsize=14)
plt.tight_layout()
fig.savefig(os.path.join(savedir,'ET_M_WY_HUC8.png'),dpi=300,bbox_inches='tight')
fig.savefig(os.path.join(savedir,'ET_M_WY_HUC8.pdf'),dpi=300,bbox_inches='tight')

# cbar = plt.colorbar(img, ax=ax,shrink=0.8)
# cbar.set_ticks(bounds)
# cbar.set_label('Long-term mean\nprecipitation (mm)',fontsize=12)
# ax.axis('off')
# ax.set_title('P_LTM_WY [HUC8]',fontsize=14,weight='bold')

# SR_LTM_WY (HUC8)
fig,ax = plt.subplots(1,1,subplot_kw={'projection':param_nwm3.cartopy_crs_atur_utm12n},figsize=(4,4))
ax=cbm.az_huc8_boundaries(ax=ax,gridlines_flag=False)
bounds = np.array([0,1,2,3,4,5,10,25,50,60,70,80,90])
norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)
cmap = 'turbo_r'
img = ds_nwm_huc8['qSfcLatRunoff_LTM'].sel(season='WY').xvec.to_geodataframe(geometry='geometry').plot(ax=ax,column='qSfcLatRunoff_LTM',norm=norm,cmap=cmap,legend=False)

# Create a ScalarMappable with the norm and cmap
sm = ScalarMappable(norm=norm, cmap=cmap)
# The ScalarMappable does not have any data associated with it, so we give it some
sm.set_array([])
cbar = plt.colorbar(sm, ax=ax, shrink=0.8, boundaries=bounds)
cbar.set_ticks(bounds)
cbar.set_label('Surface runoff (mm)',fontsize=12)
ax.axis('off')
ax.set_title('Long-term mean (1981-2020) \nsurface runoff (HUC8)',fontsize=14)
plt.tight_layout()
fig.savefig(os.path.join(savedir,'SR_M_WY_HUC8.png'),dpi=300,bbox_inches='tight')
# cbar = plt.colorbar(img, ax=ax,shrink=0.8)
# cbar.set_ticks(bounds)
# cbar.set_label('Long-term mean\nprecipitation (mm)',fontsize=12)
# ax.axis('off')
# ax.set_title('P_LTM_WY [HUC8]',fontsize=14,weight='bold')


# Re_LTM_WY (HUC8)
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
cbar.set_label('Groundwater recharge (mm)',fontsize=12)
ax.axis('off')
ax.set_title('Long-term mean (1981-2020) \ngroundwater recharge (HUC8)',fontsize=14)
plt.tight_layout()
fig.savefig(os.path.join(savedir,'Re_M_WY_HUC8.png'),dpi=300,bbox_inches='tight')
# cbar = plt.colorbar(img, ax=ax,shrink=0.8)
# cbar.set_ticks(bounds)
# cbar.set_label('Long-term mean\nprecipitation (mm)',fontsize=12)
# ax.axis('off')
# ax.set_title('P_LTM_WY [HUC8]',fontsize=14,weight='bold')



