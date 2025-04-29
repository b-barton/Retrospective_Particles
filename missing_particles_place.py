import numpy as np
from netCDF4 import Dataset
import cartopy.crs as ccrs
import os
import matplotlib.pyplot as plt
import scipy.interpolate as interp
import haversine as ca
import shapely
from shapely.ops import unary_union
import glob

from pylag.processing.plot import FVCOMPlotter
from pylag.processing.plot import create_figure, colourmap
from pylag.processing.coordinate import utm_from_lonlat, lonlat_from_utm
from pylag.processing.release_zone import create_release_zone, ReleaseZone
from pylag.processing.input import create_initial_positions_file_single_group, create_initial_positions_file_multi_group


# NOTE problem is with remove_particles() function in other script
# it performs the division check on locations in each layer not
# total particle number with groups

# Root directory for PyLag example input files
data_dir = '/dssgfs01/scratch/benbar/SSW_RS/SSW_RS_v1.2_2016_12_26/Turb/'
amm15_dir = '/dssgfs01/scratch/benbar/Processed_Data/'


# Keep a copy of the cwd
cwd = os.getcwd()

# Create run directory
simulation_dir = '{}/Simulation'.format(cwd)
try:
    os.makedirs(simulation_dir)
except FileExistsError:
    pass

# Create input sub-directory
input_dir = '{}/input'.format(simulation_dir)
try:
    os.makedirs(input_dir)
except FileExistsError:
    pass


# Grid metrics file
grid_metrics_file_name = '{}/grid_metrics_cart.nc'.format(input_dir)

part_init1 = '{}/initial_pos_cart_coast_full_430080.dat'.format(input_dir)

bathy_mod = amm15_dir + 'bathy_meter_ORIGINAL_AMM15.nc'

# number of processors
nproc = 128


# initial positions for the grid with missing particles
if 1:
  with open(part_init1) as f:
    lines = f.readlines()
  n_seed_particles = int(float(lines[0]))
  half_rg = np.zeros((n_seed_particles))
  half_x = np.zeros((n_seed_particles))
  half_y = np.zeros((n_seed_particles))
  half_z = np.zeros((n_seed_particles))
  for i in range(n_seed_particles):
    parts = lines[i + 1].split()
    half_rg[i] = int(parts[0])
    half_x[i] = float(parts[1])
    half_y[i] = float(parts[2])
    half_z[i] = float(parts[3])

# make grid like the origonal but populate all locations

array_gxy = np.array([half_rg, half_x, half_y])
array_uni = np.unique(array_gxy, axis=1)

new_z = (np.arange(0, -200, -10)[:, np.newaxis])
new_zg = np.tile(new_z, array_uni.shape[1]).flatten()
array_gxyz = np.tile(array_uni, len(new_z))
array_gxyz = np.append(array_gxyz, new_zg[np.newaxis, :], axis=0)
print(array_gxyz.shape)

# Convert utm coords to degrees
lons, lats = lonlat_from_utm(array_gxyz[1, :], array_gxyz[2, :],
                             epsg_code='32630')

# Remove particles below sea floor

with Dataset(bathy_mod, 'r') as nc_fid:
  nc_lon = nc_fid.variables['lon'][:]
  nc_lat = nc_fid.variables['lat'][:]
  nc_bathy = -nc_fid.variables['Bathymetry'][:]

# Mask points on land including in bathymetry
# interpolate bathymetry to the particle grid and mask particles deeper than bathy

interp_func = interp.LinearNDInterpolator(list(zip(nc_lon.flatten(), nc_lat.flatten())), nc_bathy.flatten())
particle_bathy = interp_func(lons, lats)
off_bathy = array_gxyz[3, :] >= (particle_bathy + 5) # true if outside bathy
sub_gxyz = array_gxyz[:, off_bathy]
print(sub_gxyz.shape)


# Make comparison between particle sets

half_gxyz = np.array([half_rg, half_x, half_y, half_z])
print(half_gxyz.shape)
lons, lats = lonlat_from_utm(half_gxyz[1, :], half_gxyz[2, :],
                             epsg_code='32630')
particle_bathy = interp_func(lons, lats)
off_bathy_h = half_gxyz[3, :] >= (particle_bathy + 5) # true if outside bathy
half_sub_gxyz = half_gxyz[:, off_bathy_h]
print(half_sub_gxyz.shape)

both_gxyz = np.append(sub_gxyz, half_sub_gxyz, axis=1)
uni_xyz = np.unique(both_gxyz[1:, :], axis=1)

keep_bool = np.zeros(both_gxyz.shape[1], dtype=bool)
c_miss = 0
for i in range(uni_xyz.shape[1]):
  bool_loc = ((both_gxyz[1, :] == uni_xyz[0, i])
                & (both_gxyz[2, :] == uni_xyz[1, i])
                & (both_gxyz[3, :] == uni_xyz[2, i]))
  n_part_in_location = np.sum(bool_loc)
  if n_part_in_location == 24: 
    # 24 if only in one set or 48 if in both
    keep_bool = keep_bool | bool_loc
    c_miss = c_miss + 1
print(c_miss, sum(keep_bool))

keep_gxyz = both_gxyz[:, keep_bool]
#plt.scatter(keep_gxyz[1, :], keep_gxyz[2, :], c=keep_gxyz[3, :])
#plt.show()

# Remove particles so they are divisable by processors

if keep_gxyz.shape[1] % nproc != 0:
  rem = keep_gxyz.shape[1] % nproc
  keep_gxyz = keep_gxyz[:, :-rem]
print(rem)
print(keep_gxyz.shape)

# Write particles to file

file_name = '{}/initial_pos_cart_coast_missing.dat'.format(input_dir)
with open(file_name, 'w') as f_out:
  f_out.write('{}\n'.format(int(keep_gxyz.shape[1])))
  for i in range(keep_gxyz.shape[1]):
    f_out.write('{} '.format(int(keep_gxyz[0, i])) 
        + '{} '.format(keep_gxyz[1, i]) 
        + '{} '.format(keep_gxyz[2, i])
        + '{}\n'.format(keep_gxyz[3, i]))

# Figure

lons, lats = lonlat_from_utm(keep_gxyz[1, :], keep_gxyz[2, :],
                             epsg_code='32630')

# Read in the bathymetry
ds = Dataset(grid_metrics_file_name, 'r')
bathy = -ds.variables['h'][:]
x_m = ds.variables['x'][:]
y_m = ds.variables['y'][:]
nv = ds.variables['nv'][:]
ds.close()
del(ds)

out_dir = '/dssgfs01/scratch/benbar/Particles/Main_Run/Store_Output/'
fnames = sorted(glob.glob(out_dir + 'connect_2d_2014*.npz'))
data = np.load(fnames[0], allow_pickle=True)
x = data['x']
y = data['y']
dx = data['dx']
dy = data['dy']
x_bound = data['x_bound']
y_bound = data['y_bound']
mask = data['mask']
data.close()


# calculate particle indices
s = 2
x_ind = ((keep_gxyz[1, :] - x_bound[0, 0]) / (dx * 2 * s)).astype(int)
y_ind = ((keep_gxyz[2, :] - y_bound[0, 0]) / (dy * 2 * s)).astype(int)
z_ind = ((keep_gxyz[3, :] - 5) / -10).astype(int)

lon_c, lat_c = lonlat_from_utm(x_bound[::s, ::s], y_bound[::s, ::s], epsg_code='32630')

num_init = np.ma.zeros((lon_c.shape[0], lon_c.shape[1], 20))
for i in range(len(x_ind)):
  num_init[y_ind[i], x_ind[i], z_ind[i]] = num_init[y_ind[i], x_ind[i], z_ind[i]] + 1
sum_init = np.ma.sum(num_init, axis=2)

# Create figure
font_size = 15
cmap = colourmap('h_r')
mrc = ccrs.PlateCarree()
fig = plt.figure(figsize=(10, 8))
ax = fig.add_axes([0.1, 0.55, 0.75, 0.4], projection=mrc)
ax2 = fig.add_axes([0.1, 0.05, 0.75, 0.4])
ax2c = fig.add_axes([0.9, 0.05, 0.02, 0.4])

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

# Plot particle initial positions
plotter.scatter(ax, lons, lats, s=8, color='#e50000', edgecolors='none')

c2 = ax2.pcolormesh(lon_c, lat_c, sum_init)
plt.colorbar(c2, cax=ax2c)

plt.savefig('./Figures/particle_start_missing.png', dpi=150)
