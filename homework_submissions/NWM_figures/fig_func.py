import rioxarray as rxr
import pandas as pd
import geopandas as gpd
import numpy as np  
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.patches import Rectangle
import re



#Loading in the NWM mean annual recharge raster (Re_WY) and the ADWR groundwater basin shapfile (GW_basin)
Re_WY = rxr.open_rasterio("/workspaces/homework_marlanahinkley/homework_submissions/NWM_figures/Re/Re_M_WY.tif")
GW_basin = gpd.read_file("/workspaces/homework_marlanahinkley/homework_submissions/NWM_figures/Groundwater_Basin/Groundwater_Basin.shp")


#Function to create basin recharge map
def create_basin_Re_map(basin): #input the geodataframe of the basin to plot
    #Creating a figure
    fig, axes = plt.subplots(figsize=(10,10))

    #Clipping the Re_WY raster to the GW basin shapefile (AZ Boundary)
    Re_WY_clipped = Re_WY.rio.clip(GW_basin.geometry, GW_basin.crs, drop=False)

    #Plotting the clipped Re_WY raster and the GW basin polygons
    Re_WY_clipped.plot(ax=axes, cmap="Blues", vmin=0, vmax=20, cbar_kwargs={"label": "Groundwater recharge (mm)"})
    GW_basin.plot(ax=axes, facecolor="none", edgecolor="black" , linewidth=1)
    axes.axis('off')
    plt.title("Long-term Mean (1981-2020) Annual Groundwater Recharge in Arizona")
        
    # create inset axes and plot zoomed raster + feature
    axins = inset_axes(axes, width="35%", height="35%", loc="upper right", borderpad=2)

    #Highlighting the specified groundwater basin
    feature = basin
    #if feature.empty:
        # fallback: try alternative column names or choose first feature
        #feature = GW_basin.iloc[[0]]

    #Clipping the raster to the specified groundwater basin
    Re_WY_basin_clip = Re_WY.rio.clip(feature.geometry, GW_basin.crs, drop=False)
    Re_WY_basin_clip.plot(ax=axins, cmap="Blues", vmin=0, vmax=20, add_colorbar=False)
    feature.plot(ax=axins, facecolor="none", edgecolor="black", linewidth=1.5)

    # draw a red rectangle on the main map showing inset area
    minx, miny, maxx, maxy = feature.total_bounds
    pad = 0.05 * max(maxx - minx, maxy - miny)
    rect = Rectangle((minx - pad, miny - pad),
                     (maxx - minx) + 2*pad,
                     (maxy - miny) + 2*pad,
                     edgecolor="red", facecolor="none", linewidth=3)
    axes.add_patch(rect)
    
    axins.set_xlim(minx - pad, maxx + pad)
    axins.set_ylim(miny - pad, maxy + pad)
    axins.set_xticks([])
    axins.set_yticks([])
    axins.set_title("")
    axins.tick_params(labelleft=False, labelright=False, labeltop=False, labelbottom=False)
    basin_name = feature.loc[0, 'BASIN_NAME']
    plt.savefig(f"{basin_name}.png")
    


#Function to create basin recharge pie chart
def plot_divP_pie(row, label_col='SUBBASIN_N', outfn=None, figsize=(4,4)):
    # require columns
    required = ['Re_divP','ET_divP','SR_divP']
    missing = [c for c in required if c not in row.index]
    if missing:
        print("missing columns:", missing); return

    vals = [float(row.get(c, 0) or 0) for c in required]
    total = sum(vals)
    if total <= 0:
        print("no positive values to plot for", row.get(label_col, 'unknown')); return

    sizes = [v/total for v in vals]  # normalize to sum==1
    labels = ['Re (recharge)','ET (evapotranspiration)','SR (surface runoff)']
    colors = ['#4c72b0','#dd8452','#55a868']

    fig, ax = plt.subplots(figsize=figsize)
    wedges, texts, autotexts = ax.pie(
        sizes,
        labels=labels,
        colors=colors,
        autopct=lambda pct: f"{pct:.1f}%",
        startangle=90,
        counterclock=False
    )
    ax.set_aspect('equal')
    title = str(row.get(label_col, '')).strip()
    ax.set_title(title if title else "Basin")
    if outfn is None:
        safe = re.sub(r'[^A-Za-z0-9_-]+', '_', title).strip('_')[:80] or 'basin'
        outfn = f"pie_{safe}.png"
    plt.savefig(outfn, bbox_inches='tight', dpi=150)
    plt.close(fig)
    print("saved:", outfn)



#Function to create basin recharge map V2 
#Switching the figure so that the inset is the full state of Arizona and the main plot is the specific basin. Note, AI made this. 

def create_basin_Re_map_V2(basin):  # input: one-row GeoDataFrame for the basin
    fig, ax = plt.subplots(figsize=(10,10))

    feature = basin.copy()
    if feature.empty:
        print("empty feature"); return

    # MAIN: zoomed basin raster + basin polygon
    try:
        basin_clip = Re_WY.rio.clip(feature.geometry, GW_basin.crs, drop=False)
    except Exception as e:
        print("basin clip failed:", e); return

    basin_clip.plot(ax=ax, cmap="Blues", vmin=0, vmax=20, cbar_kwargs={"label": "Groundwater recharge (mm)"})
    feature.plot(ax=ax, facecolor="none", edgecolor="black", linewidth=1.5)

    # zoom main axes to basin bounds
    minx, miny, maxx, maxy = feature.total_bounds
    pad = 0.05 * max(maxx - minx, maxy - miny)
    ax.set_xlim(minx - pad, maxx + pad)
    ax.set_ylim(miny - pad, maxy + pad)
    ax.axis("off")
    ax.set_title(str(feature.iloc[0].get('BASIN_NAME', 'Basin')))

    # INSET: full state context
    axins = inset_axes(ax, width="35%", height="35%", loc="upper right", borderpad=1)
    try:
        state_clip = Re_WY.rio.clip(GW_basin.geometry, GW_basin.crs, drop=False)
    except Exception as e:
        print("state clip failed:", e); state_clip = Re_WY  # fallback

    # plot state raster & basin outlines in inset
    state_clip.plot(ax=axins, cmap="Blues", vmin=0, vmax=20, add_colorbar=False)
    GW_basin.plot(ax=axins, facecolor="none", edgecolor="black", linewidth=0.6)

    # draw red rectangle on inset showing basin location
    rect = Rectangle((minx - pad, miny - pad),
                     (maxx - minx) + 2 * pad,
                     (maxy - miny) + 2 * pad,
                     edgecolor="red", facecolor="none", linewidth=2)
    axins.add_patch(rect)
    axins.set_xticks([]); axins.set_yticks([])
    axins.set_title("Arizona (context)")

    # save and close
    basin_name = str(feature.iloc[0].get('BASIN_NAME', 'basin')).strip()
    safe = re.sub(r'[^A-Za-z0-9_-]+', '_', basin_name).strip('_')[:80] or "basin"
    plt.savefig(f"{safe}.png", bbox_inches="tight", dpi=200)
    plt.close(fig)