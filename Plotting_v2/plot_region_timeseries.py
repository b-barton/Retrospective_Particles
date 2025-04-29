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

out_dir = '/dssgfs01/scratch/benbar/Particles/Main_Run/Merged_Output/'

fnames = sorted(glob.glob(out_dir + 'connect_2d_2014*.npz'))
fnames.extend(sorted(glob.glob(out_dir + 'connect_2d_2017*.npz')))
print(fnames)

input_dir = '../Simulation/input'
part_init = '{}/initial_pos_cart_coast_full_484096.dat'.format(input_dir)
#part_init = '{}/initial_pos_cart_coast_full_fix.dat'.format(input_dir)
 

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

data = np.load(out_dir + 'installation_mask_martins.npz')
#lon = data['lon']
#lat = data['lat']
i_mask = data['i_mask'] # true where installation
data.close()

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

# mask and remove blank cells

i_mask = np.ma.masked_where(mask, i_mask)
i_flat = i_mask.flatten().compressed()


lon, lat = lonlat_from_utm(x, y, epsg_code='32630')
lon = np.array(lon)
lat = np.array(lat)
print(np.min(lon), np.max(lon))
print(lon.shape)

s = 2
lon_a, lat_a = lonlat_from_utm(x_coord, y_coord, epsg_code='32630')
print(lon_a.shape)

names = ['Edge', 'Coast', 'Central', 'North', 'South']
#lonc = [np.array([-4, -0.7]), 
#        np.array([-5, -1]),
#        np.array([-0.7, 4]),
#        np.array([-0.7, 4])]
#latc = [np.array([55, 59.5]), 
#        np.array([59.5, 62]), 
#        np.array([59, 62]), 
#        np.array([55, 59])]

bx = []
by = []
bx.append([-0.7, -1.5, -2, -3.8, -5.4, -1.5, 0.5, 0.5, -0.7]) # 1
by.append([59.6, 58.8, 59.5, 59.5, 60.1, 61.7, 61.7, 61.2, 59.6])
bx.append([-4.2, -3.8, -2, -1.5, -0.4, -1, -4.2]) # 2
by.append([55.5, 59.5, 59.5, 58.8, 57.8, 55.5, 55.5])
bx.append([-0.4, -1.5, -0.7, 1.2, 2.5, 0.5, -0.4]) # 3
by.append([57.8, 58.8, 59.6, 59.1, 58.5, 57.5, 57.8])
bx.append([0.5, 0.5, 3.2, 3.9, 3.9, 2.5, 1.2, -0.7, 0.5]) # 4
by.append([61.2, 61.7, 61.1, 59.5, 58, 58.5, 59.1, 59.6, 61.2])
bx.append([-0.4, 0.5, 2.5, 3.9, 3.3, 2, 0, -1, -0.4]) # 5
by.append([57.8, 57.5, 58.5, 58, 57, 56.1, 55.5, 55.5, 57.8])


poly = []
mask_region = []
for i in range(len(names)):
  point_pair = list(zip(bx[i], by[i]))
  poly.append(shapely.geometry.Polygon(point_pair))
  mask_region.append(np.ma.zeros(i_mask.shape, dtype=bool))

flat_region = []
for i in range(lon_a.shape[0]):
  for j in range(lon_a.shape[1]):
    ll = np.array([[lon_a[i, j], lat_a[i, j]]])
    pp = shapely.geometry.Point(ll[0].tolist())
    for n in range(len(names)):
      mask_region[n][i, j] = poly[n].contains(pp)

for n in range(len(names)):
  mask_region[n] = np.ma.masked_where(i_mask.mask, mask_region[n])
  mask_region[n] = i_mask & mask_region[n]
  flat_region.append(mask_region[n].flatten().compressed())

#  tmp_mask = ((lon_a >= lonc[i][0]) & (lon_a < lonc[i][1]) 
#        & (lat_a >= latc[i][0]) & (lat_a < latc[i][1]))
#  mask_region.append(tmp_mask & i_mask)
#  flat_region.append(mask_region[i].flatten().compressed())


def undo_mask(u_mask, mask, var):
  var_m = np.ma.zeros(mask.shape)
  c = -1
  for i in range(var.shape[0]):
    for j in range(var.shape[1]):
      c = c + 1
      var_m[u_mask[0][c], u_mask[1][c]] = var[i, j]
  var_a = np.ma.masked_where(np.invert(mask), var_m)
  return var_a

def compress_mask(connect_part_m, flat_size):
  connect_compress = []
  for i in range(connect_part_m.shape[0]):
    connect_compress.append(connect_part_m[i, :].compressed())

  connect_full = np.zeros((flat_size, flat_size))
  c = -1
  for i in range(0, len(connect_compress)):
    if len(connect_compress[i]) != 0:
      c = c + 1
      connect_full[c, :] = connect_compress[i][:]
  return connect_full

def re_map(u_mask, mask, var):
  #var_m = np.ma.zeros(mask.flatten().shape)
  #var_m[u_mask] = var
  var_m = var
  var_a = np.reshape(var_m, mask.shape)
  var_a = np.ma.masked_where(mask, var_a)
  return var_a

# redo the masking

flat_mask = np.invert(mask).flatten() # true for particle
connect_mask = np.outer(flat_mask, flat_mask)
undo_con_m = np.nonzero(connect_mask)

print(len(date_full))
connect_m = np.ma.zeros((len(date_full), connect_mask.shape[0], connect_mask.shape[1]))
settle_m = np.ma.zeros((len(date_full), connect_mask.shape[0], connect_mask.shape[1]))

for i in range(len(date_full)):
  connect_m[i, :, :] = undo_mask(undo_con_m, connect_mask, connect_full[i, :, :])
  settle_m[i, :, :] = undo_mask(undo_con_m, connect_mask, settle_full[i, :, :])

flat_size = np.sum(flat_mask)
print(connect_m.shape, connect_full.shape, flat_region[0].shape)

# Plot

def set_map(ax1):
  ax1.set_extent(extents, crs=mrc)
  ax1.add_feature(land_10m, edgecolor='None', facecolor='0.5', zorder=100)
#  ax1.set_facecolor('0.5')
  gl = ax1.gridlines(draw_labels=True)
  gl.top_labels = False
  gl.right_labels = False
  gl.xlines = False
  gl.ylines = False

land_10m = cfeature.NaturalEarthFeature('physical', 'land', '10m',
                                        edgecolor='face',
                                        facecolor=cfeature.COLORS['land'])

my_cm = plt.cm.plasma
mrc = ccrs.PlateCarree()
fig1 = plt.figure(figsize=(8, 6))
ax1 = fig1.add_axes([0.08, 0.08, 0.8, 0.84], projection=mrc)
cax1 = fig1.add_axes([0.9, 0.3, 0.02, 0.4])
extents = np.array([-8.5, 4.0, 55, 62])
set_map(ax1)

map_region = mask_region[0] * 1
map_region[mask_region[1]] = 2
map_region[mask_region[2]] = 3
map_region[mask_region[3]] = 4
map_region[mask_region[4]] = 5
#map_region[mask_region[5]] = 6

np.savez(out_dir + 'map_region.npz', map_region=map_region.filled(-999), lat_a=lat_a, lon_a=lon_a)

cs1 = ax1.pcolormesh(lon_a, lat_a, map_region, cmap=my_cm, vmin=0, vmax=6)

cbar = plt.colorbar(cs1, cax=cax1)
cax1.set_xlabel('Installation Region')

fig1.savefig('./Figures/region_map.png', dpi=300)



fig2 = plt.figure(figsize=(12, 8))
ax1 = [None] * 5
ax1[0] = fig2.add_axes([0.08, 0.72, 0.19, 0.24])
ax1[1] = fig2.add_axes([0.56, 0.72, 0.19, 0.24])
ax1[2] = fig2.add_axes([0.08, 0.40, 0.19, 0.24])
ax1[3] = fig2.add_axes([0.56, 0.40, 0.19, 0.24])
ax1[4] = fig2.add_axes([0.08, 0.08, 0.19, 0.24])
#ax1[5] = fig2.add_axes([0.56, 0.08, 0.19, 0.24])
ax2 = [None] * 5
ax2[0] = fig2.add_axes([0.27, 0.72, 0.19, 0.24])
ax2[1] = fig2.add_axes([0.75, 0.72, 0.19, 0.24])
ax2[2] = fig2.add_axes([0.27, 0.40, 0.19, 0.24])
ax2[3] = fig2.add_axes([0.75, 0.40, 0.19, 0.24])
ax2[4] = fig2.add_axes([0.27, 0.08, 0.19, 0.24])
#ax2[5] = fig2.add_axes([0.75, 0.08, 0.19, 0.24])



c_list = sns.color_palette('tab10', n_colors=6)
pl_name = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)']

for i in range(len(flat_region)): # sink
#  site_source = connect_m[source_ind[i], :] # from one source
  for j in range(len(flat_region)): # source
    if i == j:
      if i == 3:
        ax2[i].plot([date_full[0], date_full[1]], [0, 0], '-', lw=2, color=c_list[i], label=names[i])
      continue
    print(connect_full[:, flat_region[j], :][:, :, flat_region[i]].shape, sum(flat_region[i]))
    
    source_init = np.ma.sum(sum_init[map_region == j])
    site_sink = np.ma.sum(connect_full[:, flat_region[j], :][:, :, flat_region[i]], axis=(1, 2)) / (source_init * sum(flat_region[i]))
    #site_sink = np.ma.sum(connect_full[:, flat_region[j], :][:, :, flat_region[i]], axis=(1, 2)) / sum(flat_region[i])
    # divide by the number of grid points in the sink

    ax1[i].plot(date_full, site_sink*100, '-', lw=2, color=c_list[j], label=names[j])
    ax2[i].plot(date_full, site_sink*100, '-', lw=2, color=c_list[j], label=names[j])

  if i == 3:
    ax2[i].legend(bbox_to_anchor=(1.05, 1.05), loc='upper right')

#  ax1[i].set_yscale('log')
#  ax2[i].set_yscale('log')
  ax1[i].set_ylabel('Particles\nin Region (%)')

  ax1[i].annotate(pl_name[i] + ' Destination ' + names[i], (0.05, 0.9), xycoords='axes fraction', fontsize=12, bbox=dict(boxstyle="round", fc="w"), zorder=105)

  ax1[i].xaxis.set_major_formatter(mdates.DateFormatter('%b-%d'))
  ax2[i].xaxis.set_major_formatter(mdates.DateFormatter('%b-%d'))
  # Rotates and right-aligns the x labels so they don't crowd each other.
  for label in ax1[i].get_xticklabels(which='major'):
    label.set(rotation=20, horizontalalignment='right')

  for label in ax2[i].get_xticklabels(which='major'):
    label.set(rotation=20, horizontalalignment='right')

  ax1[i].set_xlim([dt.datetime(2014, 1, 1), dt.datetime(2015, 1, 1)])
  ax2[i].set_xlim([dt.datetime(2017, 1, 1), dt.datetime(2018, 1, 1)])
  ax2[i].set_ylim(ax1[i].get_ylim())
  ax2[i].set_yticks([])
  ax1[i].set_xlabel('2014')
  ax2[i].set_xlabel('2017')

fig2.savefig('./Figures/timeseries_regions.png', dpi=300)




