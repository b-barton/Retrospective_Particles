import numpy as np
from netCDF4 import Dataset
import cartopy.crs as ccrs
import os
import matplotlib.pyplot as plt
import haversine as ca
import shapely
from shapely.ops import unary_union

from pylag.processing.plot import FVCOMPlotter
from pylag.processing.plot import create_figure, colourmap
from pylag.processing.coordinate import utm_from_lonlat, lonlat_from_utm
from pylag.processing.release_zone import create_release_zone, ReleaseZone
from pylag.processing.input import create_initial_positions_file_single_group, create_initial_positions_file_multi_group



# Root directory for PyLag example input files
data_dir = '/dssgfs01/scratch/benbar/SSW_RS/SSW_RS_v1.2_2016_12_26/Turb/'
amm15_dir = '/dssgfs01/scratch/benbar/Processed_Data/'

simulation_dir = '../Simulation'
# Create input sub-directory
input_dir = '{}/input'.format(simulation_dir)

# Grid metrics file
grid_metrics_file_name = '{}/grid_metrics_cart.nc'.format(input_dir)
part_init = '{}/initial_pos_cart_coast_full_484096.dat'.format(input_dir)

# Read in the bathymetry
ds = Dataset(grid_metrics_file_name, 'r')
bathy = -ds.variables['h'][:]
x_m = ds.variables['x'][:]
y_m = ds.variables['y'][:]
nv = ds.variables['nv'][:]
ds.close()
del(ds)

triangles = nv.transpose()


# initial positions for the grid
if 1:
  with open(part_init) as f:
    lines = f.readlines()
  n_seed_particles = int(float(lines[0]))

  seed_x = np.zeros((n_seed_particles))
  seed_y = np.zeros((n_seed_particles))
  seed_z = np.zeros((n_seed_particles))
  for i in range(n_seed_particles):
    parts = lines[i + 1].split()
    seed_x[i] = float(parts[1])
    seed_y[i] = float(parts[2])
    seed_z[i] = float(parts[3])


print(len(seed_x))

lons = []
lats = []
lons1, lats1 = lonlat_from_utm(seed_x,
                             seed_y,
                             epsg_code='32630')
lons.extend(lons1)
lats.extend(lats1)


# Create figure
font_size = 15
cmap = colourmap('h_r')
#fig, ax = create_figure(figure_size=(14, 10), axis_position=[0.15, 0.15, 0.65, 0.75], projection=ccrs.PlateCarree(), font_size=font_size, bg_color='gray')
mrc = ccrs.PlateCarree()
fig = plt.figure(figsize=(10, 8))
ax = fig.add_axes([0.1, 0.1, 0.75, 0.8], projection=mrc)

ax.set_extent([lons[0], lons[1], lats[0], lats[1]], crs=mrc)

# Configure plotter
plotter = FVCOMPlotter(grid_metrics_file_name,
                       geographic_coords=True,
                       font_size=font_size)
#print(np.min(plotter.x), np.max(plotter.x))
plotter.x[plotter.x > 180] = plotter.x[plotter.x > 180] - 360
plotter.xc[plotter.xc > 180] = plotter.xc[plotter.xc > 180] - 360

# Plot bathymetry
extents = np.array([-8, 6.0, 54, 62.0], dtype=float)
ax, plot = plotter.plot_field(ax, bathy, extents=extents, add_colour_bar=True, cb_label='Depth (m)',
                              vmin=-200., vmax=0., cmap=cmap)
#ax.pcolormesh(nc_lon, nc_lat, nc_bathy, vmin=-200, vmax=0, cmap=cmap)
#ax.contour(nc_lon, nc_lat, nc_bathy, [-200])

# Overlay grid
#plotter.draw_grid(ax, linewidth=1.0)

# Plot particle initial positions
plotter.scatter(ax, lons, lats, s=8, color='#e50000', edgecolors='none')

if 0: # turn on for contours
  for i in range(len(poly1.geoms)):
    x_p, y_p = poly1.geoms[i].exterior.xy
    lon_p, lat_p = lonlat_from_utm(x_p.tolist(), y_p.tolist(), epsg_code='32630')
  #  ax.plot(lon_p, lat_p)
  x_p, y_p = poly2.exterior.xy
  lon_p, lat_p = lonlat_from_utm(x_p.tolist(), y_p.tolist(), epsg_code='32630')
  ax.plot(lon_p, lat_p)
  for i in range(len(poly3.geoms)):
    x_p, y_p = poly3.geoms[i].exterior.xy
    lon_p, lat_p = lonlat_from_utm(x_p.tolist(), y_p.tolist(), epsg_code='32630')
    ax.plot(lon_p, lat_p)


plt.savefig('./Figures/particle_start_coast.png', dpi=300)
