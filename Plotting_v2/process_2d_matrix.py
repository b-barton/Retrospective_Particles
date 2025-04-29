import datetime as dt
import numpy as np
from netCDF4 import Dataset
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import glob
import sys
import copy
from pylag.processing.ncview import Viewer
from pylag.processing.plot import create_figure, colourmap
from pylag.processing.plot import FVCOMPlotter
from pylag.processing.coordinate import utm_from_lonlat, lonlat_from_utm
import shapely.geometry as geom

cwd = '../' 

simulation_dir = '{}/Simulation'.format(cwd)
input_dir = '{}/input'.format(simulation_dir)

out_dir = '/dssgfs01/scratch/benbar/Particles/Main_Run/Merged_Output'

file_names = sorted(glob.glob(out_dir + '/retrospective_particle_tracks_2017-*.nc'))

grid_metrics_file_name = '{}/grid_metrics_cart.nc'.format(input_dir)
part_init = '{}/initial_pos_cart_coast_full_484096.dat'.format(input_dir)
in_dir = '/dssgfs01/scratch/benbar/SSW_RS/SSW_RS_v1.2_2016_12_26/Turb'

# initial positions for the grid
if 0:
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

# make a grid

lat = 58.7
lon = -2.0

easting, northing, _ = utm_from_lonlat([lon], [lat], epsg_code='32630')
radius = 350000.0
z_vals = np.arange(0, -200, -10)
n_particles = 3000

area = np.pi * radius * radius
delta = np.sqrt(area / float(n_particles))
n_coords_xy = int(2.0 * radius / delta)

# Form a regular square grid of particle positions centered on centre
x_vals = np.linspace(easting - radius, easting + radius, n_coords_xy)
y_vals = np.linspace(northing - radius, northing + radius, n_coords_xy)
x_coords, y_coords = np.meshgrid(x_vals, y_vals)


# define the bounding boxes of the grid
print(x_coords.shape, x_vals.shape, y_vals.shape, z_vals.shape)
x_delta = np.mean((x_coords[:, 1:] - x_coords[:, :-1]) / 2)
x_bound = x_coords - x_delta
x_bound = np.append(x_bound, x_coords[:, -2:] + (x_delta * 2), axis=1)
x_bound = np.append(x_bound, x_bound[-2:, :], axis=0)
y_delta = np.mean((y_coords[1:, :] - y_coords[:-1, :]) / 2)
y_bound = y_coords - y_delta
y_bound = np.append(y_bound, y_coords[-2:, :] + (y_delta * 2), axis=0)
y_bound = np.append(y_bound, y_bound[:, -2:], axis=1)
print(x_delta, y_delta)
print(x_bound.shape, x_coords.shape)

# start small, scale up later
s = 2
x_sub_bound = x_bound[::s, ::s]
y_sub_bound = y_bound[::s, ::s]
x_csub = x_coords[::s, ::s]
y_csub = y_coords[::s, ::s]

# Mask for initial grid of compressed file

with Dataset(file_names[0], 'r') as ds:
  x_part = ds.variables['x'][0, :] # time, particle
  y_part = ds.variables['y'][0, :]

x_ind = ((x_part - x_bound[0, 0]) / (x_delta * 2 * s)).astype(int)
y_ind = ((y_part - y_bound[0, 0]) / (y_delta * 2 * s)).astype(int)
mask_seed = np.zeros(x_csub.shape, dtype=bool)
mask_seed[y_ind, x_ind] = 1 # true where particles start

flat_mask = (mask_seed.flatten()) # true for particle
flat_size = np.sum(flat_mask)
undo_mask = np.nonzero(flat_mask)[0]

x_save = x_csub.flatten()[flat_mask]
y_save = y_csub.flatten()[flat_mask]

print(x_sub_bound.shape, mask_seed.shape)

# Grid the particle locations

min_age = np.arange(80 * 24, 0 * 24, -10 * 24).astype(int) # hours
min_tmp = 0 # C
count_map = np.ma.zeros(mask_seed.shape)
date_list = np.zeros((len(file_names)), dtype=object)
date_list[:] = dt.datetime(1960, 1, 1)
connect_full = np.zeros((len(file_names), flat_size, flat_size)) 
settle_full = np.zeros((len(file_names), flat_size, flat_size)) 
print(connect_full.shape)

mn_index = 11 # change for differnt month start
ref = dt.datetime(1960, 1, 1)
for i in range(304, len(file_names)):
  now = dt.datetime.now()
  with Dataset(file_names[i], 'r') as ds:
    x_part = ds.variables['x'][:] # time, particle
    y_part = ds.variables['y'][:]
    z_part = ds.variables['z'][:]
    time = ds.variables['time'][:]

  #ind = z_part[0, :] > -10
  #plt.scatter(x_part[-1, ind], y_part[-1, ind])
  #plt.show()

  date_file = np.zeros((len(time)), dtype=object)
  for t in range(len(time)):
    date_file[t] = ref + dt.timedelta(seconds=int(time[t]))   

  if date_file[0].month != mn_index:
#  if i == 1:
    mn_index = date_file[0].month
    mn = np.array([d.month for d in date_list])
    yr = np.array([d.year for d in date_list])
    this_mn = (mn == date_list[i -1].month) & (yr == date_list[i -1].year)

    date_save = date_list[this_mn]
    connect_save = connect_full[this_mn, :, :]
    settle_save = settle_full[this_mn, :, :]
    print(date_save) 

    np.savez(out_dir + '/connect_2d_' + 
        date_save[-1].strftime('%Y-%m') + '.npz', 
        connect=connect_save, settle=settle_save, 
        x=x_save, y=y_save, x_bound=x_bound, y_bound=y_bound, 
        dx=x_delta, dy=y_delta, x_coords=x_csub, y_coords=y_csub, 
        mask=np.invert(mask_seed), undo_mask=undo_mask, date_list=date_save)

  date_list[i] = copy.deepcopy(date_file[0])
  print(date_file[0])

  connect_part = np.ma.zeros((mask_seed.flatten().shape[0], mask_seed.flatten().shape[0])) # time, source, sink
  settle_part = np.ma.zeros((mask_seed.flatten().shape[0], mask_seed.flatten().shape[0])) # time, source, sink


  # calculate particle indices
  x_ind = ((x_part - x_bound[0, 0]) / (x_delta * 2 * s)).astype(int)
  y_ind = ((y_part - y_bound[0, 0]) / (y_delta * 2 * s)).astype(int)

  condition = ((x_ind < 0) | (x_ind >= mask_seed.shape[1]) 
            | (y_ind < 0) | (y_ind >= mask_seed.shape[0]))
  x_ind = np.ma.masked_where(condition, x_ind)
  y_ind = np.ma.masked_where(condition, y_ind)
  x_ind = x_ind.filled(0)
  y_ind = y_ind.filled(0)

  # calculate particle source
  
  ind_arr = np.ma.array([y_ind[0, :], x_ind[0, :]]) # first timestep (source)
  source_ind = np.ravel_multi_index(ind_arr, mask_seed.shape, mode='clip') # flatten the index to the connectivity array

  # each in-domain cell will be referenced 25 times
  # source_ind will have length 430080

  for j in range(0, x_ind.shape[1], 1): # particles loop
    # calculate particle sink indices
    # could use particle age and/or temperature threshold before sink is logged
    if i == 0:
      count_map[y_ind[0, j], x_ind[0, j]] = count_map[
          y_ind[0, j], x_ind[0, j]] + 1

    for k in range(len(min_age)):
      ind_arr = np.ma.array([y_ind[min_age[k]:min_age[k] + 10 * 24, j], 
                              x_ind[min_age[k]:min_age[k] + 10 * 24, j]]) 
      sink_ind = np.ravel_multi_index(ind_arr, mask_seed.shape, mode='clip')
      sink_ind = np.unique(sink_ind) 
      
      # Each particle can only reference a cell once within the 10 day window
      # Sinks have potential for more or less particles than a 
      # source (which is fixed)

      for l in range(len(sink_ind)):
        if (sink_ind[l] == 0 | source_ind[l] == 0):
          # Skip out of domain
          # Howevere if all particles from a source go out of domain,
          # a source has no sink
          continue
#        if min_age[k] == 10 * 24: # only use 10 day for connection count
        connect_part[source_ind[j], sink_ind[l]] = (
              connect_part[source_ind[j], sink_ind[l]] + 1)
        settle_part[source_ind[j], sink_ind[l]] = min_age[k]

  # mask and remove empty grid cells

  connect_mask = np.invert(np.outer(flat_mask, flat_mask)) # true for no particle
  connect_part_m = np.ma.masked_where(connect_mask, connect_part)
  settle_part_m = np.ma.masked_where(connect_mask, settle_part)

  connect_compress = []
  settle_compress = []
  for t1 in range(connect_part_m.shape[0]):
    connect_compress.append(connect_part_m[t1, :].compressed())
    settle_compress.append(settle_part_m[t1, :].compressed())

  c = -1
  for t1 in range(0, len(connect_compress)):
    if len(connect_compress[t1]) != 0:
      c = c + 1
      connect_full[i, c, :] = connect_compress[t1][:]
      settle_full[i, c, :] = settle_compress[t1][:]


  print('Day:', i, 'Time:', dt.datetime.now() - now)

mn = np.array([d.month for d in date_list])
yr = np.array([d.year for d in date_list])
this_mn = (mn == date_list[i -1].month) & (yr == date_list[i -1].year)

date_save = date_list[this_mn]
connect_save = connect_full[this_mn, :, :]
settle_save = settle_full[this_mn, :, :]
 
np.savez(out_dir + '/connect_2d_' + 
    date_save[-1].strftime('%Y-%m') + '.npz', 
    connect=connect_save, settle=settle_save, 
    x=x_save, y=y_save, x_bound=x_bound, y_bound=y_bound, 
    dx=x_delta, dy=y_delta, x_coords=x_csub, y_coords=y_csub, 
    mask=np.invert(mask_seed), undo_mask=undo_mask, date_list=date_save)



# map of particle connections

count_map = np.ma.masked_where(np.invert(mask_seed), count_map)
count_map1 = count_map * 1
count_map[count_map == 0] = -200
count_map[count_map > 0] = 0
count_map[count_map == -200] = 1

print(connect_full.shape)
print(np.ma.min(connect_full), np.ma.max(connect_full), np.ma.count_masked(connect_full))
print(np.sum((connect_full > 1).astype(int)), np.sum((connect_full == 0).astype(int)))


np.savez(out_dir + '/connect_2d.npz', connect=connect_full, settle=settle_full, x=x_save, y=y_save, x_bound=x_bound, y_bound=y_bound, dx=x_delta, dy=y_delta, x_coords=x_csub, y_coords=y_csub, count_map=count_map1, mask=np.invert(mask_seed), undo_mask=undo_mask, date_list=date_list)

# Simple plots

fig1 = plt.figure(figsize=(8, 8))
ax1 = fig1.add_axes([0.1, 0.1, 0.8, 0.8])
ax1.pcolormesh((connect_full[0, :, :] > 1).astype(int), vmin=0, vmax=1)
fig1.savefig('./Figures/basic_connect_2d.png', dpi=250)

fig2 = plt.figure(figsize=(8, 8))
ax1 = fig2.add_axes([0.1, 0.1, 0.8, 0.8])
ax1.pcolormesh(count_map)
fig2.savefig('./Figures/count_map_2d.png')

#plt.plot(np.ma.sum(connect_part, axis=1))
#plt.savefig('./Figures/particles_in_source.png')



