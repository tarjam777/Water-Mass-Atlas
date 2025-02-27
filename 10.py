import xarray as xr
import numpy as np
import plotly.graph_objects as go
import os

# Load bathymetry data
path = '/Users/tarjaluna/UHM_Ocean_BGC_Group Dropbox/Datasets/Data_Products/Bathymetry/'
file = 'coarsened_bathy.nc'
bath = xr.open_dataset(path + file, engine='netcdf4')

# Load Argo data for density layers
name = 'Ind'
path_argo = '/Users/tarjaluna/UHM_Ocean_BGC_Group Dropbox/Datasets/Homemade_Products/Argo_glodap_per_water_mass/'
list_of_files = os.listdir(path_argo)
list_of_files = [path_argo + each for each in list_of_files if name in each]

# Extract densities from file names (assuming they follow a specific pattern)
densities = [file[-7:-2] for file in list_of_files]

# Filter the list to include only files that match the extracted densities
dens_layers = [file for file in list_of_files if any(density in file for density in densities)]

print(f"Name: {name}")
print(f"Densities: {densities}")
print(f"Density Layer Files: {dens_layers}")

# Variables to plot
variables = ['Temperature_layer', 
'Salinity_layer', 
'Nitrate_layer', 
'Doxy_layer', 
'Phosphate_layer', 
'DIC_layer', 
'pH_layer', 
'alk_layer']

alldata = {var: [] for var in variables}
lons, lats, pressure = [], [], []

# Loop through density layer files and extract relevant data
for file in dens_layers:
    ds = xr.open_dataset(file)
    for var in variables:
        alldata[var].extend(ds[var].values.flatten())  # Flatten to ensure 1D
    lons.extend(ds.Lons_layer.values.flatten())
    lats.extend(ds.Lats_layer.values.flatten())
    pressure.extend(ds.Pressure_layer.values.flatten())

# Convert to numpy arrays
lons = np.array(lons)
lats = np.array(lats)
pressure = np.array(pressure)

# Correct longitudes
lons_cor = np.where(lons > 180, lons - 360, lons)

# Handle longitude wraparound in bathymetry data
if np.min(lons_cor) < -170 and np.max(lons_cor) > 170:
    bath['lon'] = np.where(bath.lon < 0, bath.lon + 360, bath.lon)
    bath = bath.sortby(bath.lon)
    lons_cor = lons  # No correction needed

# Find data range
lonmin, lonmax = np.min(lons_cor), np.max(lons_cor)
latmin, latmax = np.min(lats), np.max(lats)

# Subset bathymetry data
bathy = bath.sel(lat=slice(latmin, latmax), lon=slice(lonmin, lonmax))

# Create meshgrids
lon_mesh, lat_mesh = np.meshgrid(bathy.lon, bathy.lat)
elevation_flipped = -bathy['elevation']

# Directory to save screenshots
save_dir = f'/Users/tarjaluna/UHM_Ocean_BGC_Group Dropbox/UHMBGC Lab/Luna/Watermasses/{name}/frames/'

# Color dictionary
color_dict = {
    'Temperature_layer': 'magma',
    'Salinity_layer': 'viridis',
    'Nitrate_layer': 'cividis',
    'Doxy_layer': 'plasma',
    'Phosphate_layer': 'inferno',
    'DIC_layer': 'thermal',
    'pH_layer': 'tealrose',
    'alk_layer': 'icefire'
}

# Unit dictionary 
varunits = {
    'Temperature_layer': '°C',
    'Salinity_layer': 'psu',
    'Nitrate_layer': 'µmol/kg',
    'Doxy_layer': 'µmol/kg',
    'Phosphate_layer': 'µmol/kg',
    'DIC_layer': 'µmol/kg',
    'pH_layer': '',
    'alk_layer': 'µmol/kg'
}

varname = {
    'Temperature_layer': 'Temperature',
    'Salinity_layer': 'Salinity',
    'Nitrate_layer': 'Nitrate',
    'Doxy_layer': 'Oxygen',
    'Phosphate_layer': 'Phosphate',
    'DIC_layer': 'DIC',
    'pH_layer': 'pH',
    'alk_layer': 'Alkalinity'
}

# Calculate the axis ranges based on the bathymetry data and Argo data
x_range = [lonmin, lonmax]
y_range = [latmin, latmax]
z_range = [np.min(elevation_flipped.values), np.max(elevation_flipped.values)]
print('Calculating axis ranges...')

# Define specific vertical extents for initial and zoomed-in views
initial_z_range = [np.max(elevation_flipped.values), 0]  # viewing bottom of bathy
zoomed_in_z_range = [np.max(pressure) + 200, 0]  # zooming into Argo data 

# Function to create content of plot and save screenshot
def plot_and_save(variable, angle, elevation, zoom, frame_number, z_axis_range):
    fig = go.Figure()

    # Plot bathymetry surface
    fig.add_trace(go.Surface(z=elevation_flipped, x=lon_mesh, y=lat_mesh, colorscale='Greys', opacity=0.7, showscale=False))

    # Extract variable data
    variable_data = np.array(alldata[variable])

    # Filter out invalid data points
    valid_indices = ~np.isnan(variable_data)
    variable_data = variable_data[valid_indices]
    valid_lons = lons_cor[valid_indices]
    valid_lats = lats[valid_indices]
    valid_pressure = pressure[valid_indices]

    # Subsample data (every 10th point)
    print(f"Subsampling data for {variable}...")
    variable_data = variable_data[::10]
    valid_lons = valid_lons[::10]
    valid_lats = valid_lats[::10]
    valid_pressure = valid_pressure[::10]

    # Determine if the axes should be hidden based on the angle and elevation
    hide_axes = (elevation > np.pi/6) and (angle > np.pi/2) and (angle < 3*np.pi/2)
    
    # Hide axes if the condition is met
    if hide_axes:
        fig.update_scenes(xaxis_visible=False, yaxis_visible=False)

    # Plot Argo scatter for merged density layer
    fig.add_trace(go.Scatter3d(
        x=valid_lons, y=valid_lats, z=valid_pressure,
        mode='markers',
        marker=dict(
            size=5, 
            color=variable_data, 
            colorscale=color_dict[variable], 
            colorbar=dict(
                title=dict(
                    text=f'{varname[variable]} ({varunits[variable]})',  # Add units here
                    side='right'
                )
            )
        ),
        name=f'{name}-{densities} kg/m³'
    ))

    # Control layout and camera changes (fine-tuning of above)
    fig.update_layout(
        template = 'none', # to try to get rid of flipping axis
        margin=dict(l=100, b=100, r=100, t=100),
        scene=dict(
            xaxis=dict(title='Longitude', range=x_range),
            yaxis=dict(title='Latitude', range=y_range),
            zaxis=dict(title='Pressure (db)', range=z_axis_range),
            camera=dict(
                eye=dict(
                    x=zoom * np.cos(angle),
                    y=zoom * np.sin(angle),
                    z=zoom * np.sin(elevation)
                ),
                up=dict(x=0, y=0, z=1)
            )
        ),
        title=f'{name} {"-".join(densities)} kg/m³'
        #title=f'{name}-{densities} kg/m³'
    )

    # Save screenshot
    print(f"Saving screenshot for {variable}, frame {frame_number}...")
    fig.write_image(f"{save_dir}{variable}_frame_{frame_number:04d}.png")

# Parameters for the animation sequence
zoom = 2.0
print('Starting zoom...')
elevations = np.linspace(0, np.pi/4, num=25)
angles = np.linspace(0, 2 * np.pi, num=50)

# Outer loop to iterate through variables
for variable in variables:
    frame_number = 0

    # Full rotation at elevation 0
    for angle in angles:
        print(f'Rotating {variable} at angle {angle}...')
        plot_and_save(variable, angle, 0, zoom, frame_number, initial_z_range)
        frame_number += 1

    # Gradually tilt the camera to 45 degrees
    for elevation in elevations:
        print(f'Tilting {variable} to elevation {elevation}...')
        plot_and_save(variable, 0, elevation, zoom, frame_number, initial_z_range)
        frame_number += 1

    # Full rotation at elevation 45 degrees
    for angle in angles:
        print(f'Rotating {variable} at elevation 45 degrees and angle {angle}...')
        plot_and_save(variable, angle, np.pi/4, zoom, frame_number, initial_z_range)
        frame_number += 1

    # Gradually tilt the camera back to 0 degrees
    for elevation in reversed(elevations):
        print(f'Returning {variable} to origin from elevation {elevation}...')
        plot_and_save(variable, 0, elevation, zoom, frame_number, initial_z_range)
        frame_number += 1

    # Zoom in and full rotation at elevation 0
    for angle in angles:
        print(f'Zooming and rotating {variable} at angle {angle}...')
        plot_and_save(variable, angle, 0, zoom, frame_number, zoomed_in_z_range)
        frame_number += 1

    # Gradually tilt the camera to 45 degrees and zoom in further
    for elevation in elevations:
        print(f'Tilting up and zooming in {variable} to elevation {elevation}...')
        plot_and_save(variable, 0, elevation, zoom, frame_number, zoomed_in_z_range)
        frame_number += 1

    # Full rotation at elevation 45 degrees with further zoom
    for angle in angles:
        print(f'Rotating at {variable} at elevation 45 degrees and angle {angle}...')
        plot_and_save(variable, angle, np.pi/4, zoom, frame_number, zoomed_in_z_range)
        frame_number += 1

    # Gradually tilt the camera back to 0 degrees
    for elevation in reversed(elevations):
        print(f'Returning {variable} to origin from elevation {elevation}...')
        plot_and_save(variable, 0, elevation, zoom, frame_number, zoomed_in_z_range)
        frame_number += 1

print("All screenshots saved successfully.")
