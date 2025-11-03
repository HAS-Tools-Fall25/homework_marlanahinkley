import xarray as xr
import matplotlib.pyplot as plt
from dask.diagnostics import ProgressBar

def calculate_climatology(ds, latitude, longitude):
    location_ds = ds.sel(latitude=latitude, longitude=longitude, method='nearest')
    clim_ds = location_ds.groupby(location_ds['time'].dt.dayofyear).mean(skipna=True)
    return clim_ds

cities_to_plot = {
    'Tucson': (32.22, -110.97),
    'Chicago': (41.8832, -87),
    'New York': (40.7128, -74.0060)
    }


def plot_climatology(ds):
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), dpi=150)
    ds['temperature_2m'].plot(ax=axes[0,0], color='salmon')
    ds['precipitation_surface'].plot(ax=axes[0,1], color='skyblue')
    ds['wind_u_10m'].plot(ax=axes[1,0], color='green')
    ds['wind_v_10m'].plot(ax=axes[1,1], color='limegreen')
    plt.tight_layout()
    return fig, axes

def run_pipeline(city, ds, lat, long):
    print(f'Running climatology for {city}')
    ds = calculate_climatology(ds, lat, long)
    fig, axes = plot_climatology(ds)
    plt.savefig(f'{city}_climatology.png')

if __name__ == "__main__":
    ds = xr.open_zarr("https://data.dynamical.org/noaa/gfs/analysis-hourly/latest.zarr?email=optional@email.com")
    sub_ds = ds.isel(time=slice(-8760, None))  # Last year of hourly data
    for city, (lat,lon) in cities_to_plot.items():
        run_pipeline(city, sub_ds, lat, lon)