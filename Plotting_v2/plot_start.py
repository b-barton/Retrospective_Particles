import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from netCDF4 import Dataset
import networkx as nx
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from pylag.processing.coordinate import lonlat_from_utm, utm_from_lonlat
import glob
import seaborn as sns
import datetime as dt
import shapely

out_dir = '/dssgfs01/scratch/benbar/Particles/Main_Run/Store_Output/'

fnames = sorted(glob.glob(out_dir + 'connect_2d_2014*.npz'))
fnames.extend(sorted(glob.glob(out_dir + 'connect_2d_2017*.npz')))
print(fnames)

input_dir = '../Simulation/input'
part_init = '{}/initial_pos_cart_coast_full_484096.dat'.format(input_dir)
 

file_names = sorted(glob.glob(out_dir + '/pylag_2014-*.nc'))

for i in range(len(fnames)):

  data = np.load(fnames[i], allow_pickle=True)
  connect = data['connect'] # source, sink
  settle = data['settle'] / 24 # days
  x = data['x']
  y = data['y']
  x_bound = data['x_bound']
  y_bound = data['y_bound']
  x_coord = data['x_coords']
  y_coord = data['y_coords']
  dx = data['dx']
  dy = data['dy']
  mask = data['mask'] # map, true for no particle
  u_mask = data['undo_mask'] # 1d, index for no particle
  date_list = data['date_list']
  data.close()
  print(connect.shape, x.shape)

  if i == 0:
    connect_full = connect * 1
    settle_full = settle * 1
    date_full = date_list[:]
  else:
    connect_full = np.append(connect_full, connect, axis=0)
    settle_full = np.append(settle_full, settle, axis=0)
    date_full = np.append(date_full, date_list, axis=0)

# Load the initial particle location and count how many per grid location

#with Dataset(file_names[0], 'r') as ds:
#  x_part = ds.variables['x'][0] # time, particle
#  y_part = ds.variables['y'][0]
#  z_part = ds.variables['z'][0]
#  time = ds.variables['time'][:]

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

# calculate particle indices
s = 2
x_ind = ((seed_x - x_bound[0, 0]) / (dx * 2 * s)).astype(int)
y_ind = ((seed_y - y_bound[0, 0]) / (dy * 2 * s)).astype(int)
z_ind = ((seed_z - 5) / -10).astype(int)

condition = ((x_ind < 0) | (x_ind >= mask.shape[1])
    | (y_ind < 0) | (y_ind >= mask.shape[0]))
#x_ind = np.ma.masked_where(condition, x_ind)
#y_ind = np.ma.masked_where(condition, y_ind)
#x_ind = x_ind.filled(0)
#y_ind = y_ind.filled(0)

num_init = np.ma.zeros((mask.shape[0], mask.shape[1], 20))
#num_init = np.ma.masked_where(mask, num_init)

for i in range(len(x_ind)):
  num_init[y_ind[i], x_ind[i], z_ind[i]] = num_init[y_ind[i], x_ind[i], z_ind[i]] + 1
print(np.ma.max(num_init))
sum_init = np.ma.sum(num_init, axis=2)
print(np.ma.sum(num_init))

x, y, _ = utm_from_lonlat([-1], [61], epsg_code='32630')
print(np.sum((seed_y > y).astype(int)))

plt.pcolormesh(sum_init[:, :])
plt.plot([18, 18], [0, 30], '-r')
plt.colorbar()
plt.show()

plt.pcolormesh(num_init[:, 18, :].T)
plt.gca().invert_yaxis()
plt.colorbar()
plt.show()

