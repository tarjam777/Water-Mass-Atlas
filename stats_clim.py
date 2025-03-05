####description...

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from scipy.stats import binned_statistic_2d
import cmocean
from matplotlib.cm import get_cmap 

# Load the Argo dataset
name = 'SASTMW_26.6'
path_argo = '/Users/tarjaluna/UHM_Ocean_BGC_Group Dropbox/Datasets/Homemade_Products/Argo_glodap_per_water_mass/Argo_watermasses/'
file_argo = 'Argo_glodap_final_snap022024_%s.nc' % name
glo = xr.open_dataset(path_argo + file_argo)
layer = float(name[-4:])

# Load climatology dataset
path = '/Users/tarjaluna/UHM_Ocean_BGC_Group Dropbox/Datasets/Homemade_Products/Argo_glodap_per_water_mass/Clim_per_water_mass/'
dsclim = xr.open_dataset(path + 'clim_neutraldens_glodap.nc')

# Extract Date_layer, latitudes, and longitudes
date_layer = glo['Date_layer'].values
lats = glo['Lats_layer'].values
lons = glo['Lons_layer'].values

# Adjust longitudes to -180 to 180 range
lons_adjusted = np.where(lons > 180, lons - 360, lons)
dsclim['lon'] = np.where(dsclim.lon > 180, dsclim.lon - 360, dsclim.lon)
dsclim = dsclim.sortby('lon')

variables = ['Nitrate_layer', 'Temperature_layer', 'Doxy_layer', 'Salinity_layer']
#'DIC_Layer', 'Phosphate_layer'

names = {'Doxy_layer': 'Oxygen',
          'Nitrate_layer': 'Nitrate',
          'Salinity_layer': 'Sal',
          'Temperature_layer': 'Temp',
          #'DIC_Layer': 'DIC',
          'Phosphate_layer': 'Phosphate'}

# Define colormaps and limits
color_dict = {
    'Nitrate_layer': get_cmap('cubehelix'),
    'Temperature_layer': get_cmap('RdPu'),
    'Doxy_layer': get_cmap('turbo'),
    'Salinity_layer': cmocean.cm.haline,
    #'DIC_layer': cmocean.cm.deep,
    'Phosphate_layer': cmocean.cm.matter
}

unit_dict = {
    "Doxy_layer": "µmol/kg",
    "Nitrate_layer": "µmol/kg",
    "Temperature_layer": "°C",
    "Salinity_layer": "PSU",
    "pH_layer": "pH units",
    #"DIC_layer": "µmol/kg",
    "Alkalinity_layer": "µmol/kg"
}

# Define min and max values for 'pure' data
vmin_dict = {var: np.nanmin(glo[var].values) for var in variables if var in glo}
vmax_dict = {var: np.nanmax(glo[var].values) for var in variables if var in glo}

vlim_pure_dict = {
    'Nitrate_layer': 3,
    'Temperature_layer': 1,  
    'Doxy_layer': 10,
    'Salinity_layer': 0.1,     
    #'DIC_Layer': 2.5,
    'Phosphate_layer': 2.5
}

# Use vlim_pure_dict directly without scale factor
adjusted_vmin_dict = {
    var: vmin - vlim_pure_dict.get(var, 0)  # Default to 0 if var not in vlim_pure_dict
    for var, vmin in vmin_dict.items()
}
adjusted_vmax_dict = {
    var: vmax + vlim_pure_dict.get(var, 0)  
    for var, vmax in vmax_dict.items()
}

# Reverse the anomaly colormap
reversed_anomaly_cmap = plt.colormaps['RdBu_r']  # '_r' reverses the colormap

# Include watermass name in the plot title
watermass_title = name

# Define limits for anomaly data - avoid saturated grids
vlim_dict = {
    'Nitrate_layer': 2,
    'Temperature_layer': 1,  
    'Doxy_layer': 15,
    'Salinity_layer': 0.1,     
    #'DIC_Layer': 2.5,
    'Phosphate_layer': 2.5
}

# Set opacity for anomaly plots
anomaly_opacity = 0.8

# Time periods for plotting (10 years except for last)
# Time periods for plotting (10 years except for last)
time_periods = [
    (2012, 2022),
    (2001, 2011),
    (1992, 2002),
    (1983, 1993),
    (1974, 1984),
    #(1970, 1975),
    #(1988,1992),
]

lon_margin = 5
lat_margin = 5
lon_bins = np.arange(np.nanmax([-180, np.floor(lons_adjusted.min()) - lon_margin]), np.nanmin([180, np.ceil(lons_adjusted.max()) + lon_margin + 1]), 1)
lat_bins = np.arange(np.floor(lats.min()) - lat_margin, np.ceil(lats.max()) + lat_margin + 1, 1)

# Assume glo, date_layer, lons_adjusted, lats, dsclim, and layer are defined datasets/arrays
for var in variables:  # Process all variables in the list
    clean_name = names.get(var, var)  # Get the cleaned name from the dictionary
    unit = unit_dict.get(var, "")  # Get the unit from the dictionary

    # Step 1: Determine valid years
    valid_years = []
    for start_year, end_year in time_periods:
        time_filtered_indices = (date_layer >= start_year) & (date_layer <= end_year)
        if time_filtered_indices.sum() != 0:  # Check if there is data for this period
            valid_years.append((start_year, end_year))

    # Step 2: Adjust subplot grid dynamically
    nrows = (len(valid_years) // 2) + 1  # Adjust number of rows for subplots (2 columns for example)

    #if lon_bins[0] < -170 and lon_bins[-1] > 170:
        #central_longitude = 180
    #else:
        #central_longitude = 0

    if ('Pac' in name) | ('NPCMW' in name):
        central_longitude = 220
    else:
        central_longitude = 0

    fig, axs = plt.subplots(len(valid_years), 2, figsize=(10, 12), subplot_kw={'projection': ccrs.PlateCarree(central_longitude=central_longitude)}, sharex=True, sharey=True)  # PlateCarree projection

    fig.suptitle(f'Water Mass: {watermass_title}', fontsize=10)

    for i, (start_year, end_year) in enumerate(valid_years):
        ax_pure = axs[i, 0]  # For pure data
        ax_anom = axs[i, 1]  # For anomalies

        # Filter the dataset by time period
        time_filtered_indices = (date_layer >= start_year) & (date_layer <= end_year)
        lons_filtered = lons_adjusted[time_filtered_indices]
        lats_filtered = lats[time_filtered_indices]
        var_values = glo[var].values[time_filtered_indices]

        if time_filtered_indices.sum() != 0:
            # Compute the binned data
            statistic, x_edges, y_edges, binnumber = binned_statistic_2d(
                lons_filtered, lats_filtered, var_values,
                statistic=np.nanmean, bins=[lon_bins, lat_bins]
            )

            # Resample or interpolate clim data to match statistic dimensions
            clim_resampled = dsclim[names[var]].sel(gamma=layer).interp(
                lon=(x_edges[:-1] + x_edges[1:]) / 2,
                lat=(y_edges[:-1] + y_edges[1:]) / 2
            )
            clim = dsclim[names[var]].sel(gamma=layer).sel(lon=slice(lon_bins[0], lon_bins[-1]), lat=slice(lat_bins[0], lat_bins[-1]))

            # Compute anomaly
            anom = statistic.T - clim

            # Set extent and plot coastlines for pure data
            ax_pure.coastlines(resolution='50m')
            ax_pure.add_feature(cfeature.BORDERS, linestyle=':')
            ax_pure.add_feature(cfeature.LAND, facecolor='lightgray')
            ax_pure.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)

            # Plot the pure data with specified colormap and limits
            mesh_pure = ax_pure.pcolormesh(x_edges, y_edges, statistic.T, cmap=color_dict[var],
                                        shading='auto', vmin=adjusted_vmin_dict[var], vmax=adjusted_vmax_dict[var],
                                        transform=ccrs.PlateCarree())
            ax_pure.set_title(f'{clean_name} Pure ({start_year} to {end_year}) {unit}', fontsize=10)

            # Set extent for lons and lats
            #ax_pure.set_extent([lon_bins[0] - 2, lon_bins[-1] + 2, lat_bins[0] - 2, lat_bins[-1] + 2])
            ax_pure.set_extent([lon_bins[0] - 2, lon_bins[-1] + 2, lat_bins[0] - 2, lat_bins[-1] + 2], crs=ccrs.PlateCarree())
            ax_anom.set_extent([lon_bins[0] - 2, lon_bins[-1] + 2, lat_bins[0] - 2, lat_bins[-1] + 2], crs=ccrs.PlateCarree())

            # Set extent and plot coastlines for anomaly
            ax_anom.coastlines(resolution='50m')
            ax_anom.add_feature(cfeature.BORDERS, linestyle=':')
            ax_anom.add_feature(cfeature.LAND, facecolor='lightgray')
            ax_anom.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)

            # Plot the anomaly with specified opacity and RdBu colormap
            mesh_anom = ax_anom.pcolormesh(x_edges, y_edges, anom, cmap=reversed_anomaly_cmap, shading='auto',
                                        vmin=-vlim_dict[var], vmax=vlim_dict[var],
                                        alpha=anomaly_opacity, transform=ccrs.PlateCarree())
            ax_anom.set_title(f'{clean_name} Anomaly ({start_year} to {end_year}) {unit}', fontsize=10)

    # Add more space between each subplot
    plt.subplots_adjust(hspace=0.8)
    plt.subplots_adjust(wspace=0.00001)

    # Pure data colorbar on the left
    if 'mesh_pure' in locals():  # Check if mesh_pure exists
        cbar_pure = fig.colorbar(mesh_pure, ax=axs[:, 0], orientation='vertical', location='left', label=f'{clean_name} {unit} (Pure Data)', pad=0.1)
        cbar_pure.ax.tick_params(labelsize= 8)
        cbar_pure.set_label(f'{clean_name} Pure {unit}', fontsize=10)

        # Move the pure data colorbar further to the left
        cbar_pure.ax.set_position([cbar_pure.ax.get_position().x0 - 0.05,
                                cbar_pure.ax.get_position().y0,
                                cbar_pure.ax.get_position().width,
                                cbar_pure.ax.get_position().height])

    # Anomaly colorbar on the right
    if 'mesh_anom' in locals():  # Check if mesh_anom exists
        cbar_anom = fig.colorbar(mesh_anom, ax=axs[:, 1], orientation='vertical', location='right', label=f'{clean_name} {unit} Anomaly', pad=0.1)
        cbar_anom.ax.tick_params(labelsize=8)
        cbar_anom.set_label(f'{clean_name} Anomaly {unit}', fontsize=10)

        # Move the anomaly colorbar further to the right
        cbar_anom.ax.set_position([cbar_anom.ax.get_position().x0 + 0.04,
                                cbar_anom.ax.get_position().y0,
                                cbar_anom.ax.get_position().width,
                                cbar_anom.ax.get_position().height])

    # Save and display the plot
    plt.savefig(f'/Users/tarjaluna/UHM_Ocean_BGC_Group Dropbox/Collaborative_projects/griddedmaps/figures/stats_clim/{name}/10yr/{var}_pure_anomaly_gridded_.png')
