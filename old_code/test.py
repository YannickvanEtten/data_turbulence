import xarray as xr

ds = xr.open_dataset(
    r"C:\Users\yanni\OneDrive\Documenten\Universiteit\Turbulence project\Data\climate_data_02_07_2016.grib",
    engine="cfgrib"
)
# Keep 2 timesteps, a 15°×15° window — still real, ~2 MB as NetCDF
small = ds.isel(time=slice(0, 2)).sel(latitude=slice(55, 40), longitude=slice(-60, -45))
small.to_netcdf("era5_validation_subset.nc")