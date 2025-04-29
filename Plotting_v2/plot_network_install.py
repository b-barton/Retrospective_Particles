import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from pylag.processing.coordinate import lonlat_from_utm, utm_from_lonlat
import glob

out_dir = '/dssgfs01/scratch/benbar/Particles/Main_Run/Merged_Output/'

yr = '2014'
fnames = sorted(glob.glob(out_dir + 'connect_2d_' + yr + '*.npz'))
yr = '2017'
fnames.extend(sorted(glob.glob(out_dir + 'connect_2d_' + yr + '*.npz')))
print(fnames)

for i in range(len(fnames)):

  data = np.load(fnames[i], allow_pickle=True)
  connect = data['connect'][:, :, :] # time, source, sink
  #count_map = data['count_map']
  x = data['x']
  y = data['y']
  x_bound = data['x_bound']
  y_bound = data['y_bound']
  x_coord = data['x_coords']
  y_coord = data['y_coords']
  dx = data['dx']
  dy = data['dy']
  mask = data['mask']
  u_mask = data['undo_mask']
  date_list = data['date_list']
  print(connect.shape, x.shape)

  if i == 0:
    connect_full = connect * 1
    date_full = date_list[:]
  else:
    connect_full = np.append(connect_full, connect, axis=0)
    date_full = np.append(date_full, date_list, axis=0)

data = np.load(out_dir + 'installation_mask_martins.npz')
lon = data['lon']
lat = data['lat']
i_mask = data['i_mask'] # true where installation
data.close()

# mask and remove blank cells

i_mask = np.ma.masked_where(mask, i_mask)
i_flat = i_mask.flatten().compressed()

print(connect_full.shape)
print(i_flat.shape)

connect_source = connect_full * 1
connect_sink = connect_full * 1
connect_source[:, np.invert(i_flat), :] = 0 # installations are only source
connect_sink[:, :, np.invert(i_flat)] = 0 # installations are only sink
connect_inst = connect_source * 1
connect_inst[:, :, np.invert(i_flat)] = 0 # installations only sink and source

mn = np.array([d.month for d in date_full])
connect_sink = np.ma.mean(connect_sink, axis=0) # full time mean
connect_source = np.ma.mean(connect_source, axis=0) # full time mean
connect_inst = np.ma.mean(connect_inst, axis=0) # full time mean
connect_all = np.ma.mean(connect_full, axis=0) # full time mean

# Make it so installations are the only source and sink
connect_source = connect_inst * 1
connect_sink = connect_inst * 1

#count_map = np.ma.masked_where(mask, count_map)

lon, lat = lonlat_from_utm(x, y, epsg_code='32630')
lon = np.array(lon)
lat = np.array(lat)
print(np.min(lon), np.max(lon))

s = 2
lon_a, lat_a = lonlat_from_utm(x_coord, y_coord, epsg_code='32630')

lon_bound, lat_bound = lonlat_from_utm(x_bound, y_bound, epsg_code='32630')
lon_bound = lon_bound[::s, ::s]
lat_bound = lat_bound[::s, ::s]

n_point = connect_sink.shape[0]

dg1 = nx.DiGraph()
dg2 = nx.DiGraph()
dg3 = nx.DiGraph()
pos = {}
for i in range(n_point):
  dg1.add_node(i, x=x[i], y=y[i]) 
  dg2.add_node(i, x=x[i], y=y[i]) 
  dg3.add_node(i, x=x[i], y=y[i]) 
  pos[i] = (lon[i], lat[i])

for i in range(n_point):
  for j in range(n_point):
    if (i == j) :
      continue
    if connect_sink[i, j] > 100:
      dg1.add_edge(i, j)
    if connect_source[i, j] > 100:
      dg2.add_edge(i, j)
    if connect_all[i, j] > 100:
      dg3.add_edge(i, j)




def plot_box(ax, box_list_x, box_list_y):
  for i in range(len(box_list_x)):
    ax.plot(box_list_x[i], box_list_y[i], '-k', lw=1, zorder=105)


out_degree = np.array(list(zip(*dg1.out_degree()))[1])
in_degree = np.array(list(zip(*dg2.in_degree()))[1])

#bet = np.array(list(nx.bridges(dg1.to_undirected())))
bet = np.array(list(nx.betweenness_centrality(dg1).values()))
close_in = np.array(list(nx.closeness_centrality(dg2).values())) # incoming
close_out = np.array(list(nx.closeness_centrality(dg1.reverse()).values())) # outgoing
#eig = nx.eigenvector_centrality(dg)
cluster = np.array(list(nx.clustering(dg2).values()))

def re_map(u_mask, mask, var):
  var_m = np.ma.zeros(mask.flatten().shape)
  var_m[u_mask] = var
  var_a = np.reshape(var_m, mask.shape)
  var_a = np.ma.masked_where(mask, var_a)
  return var_a

in_d = (re_map(u_mask, mask, in_degree) / (n_point - 1)) * 100
out_d = (re_map(u_mask, mask, out_degree) / (n_point - 1)) * 100
bet_d = re_map(u_mask, mask, bet) / ((n_point - 1) * (n_point - 2)) # rescale
bet_d = ((bet_d - np.min(bet_d)) / (np.max(bet_d) - np.min(bet_d))) * 100
clo_i = re_map(u_mask, mask, close_in) * 100 # %
clo_o = re_map(u_mask, mask, close_out) * 100 # %
clust = re_map(u_mask, mask, cluster) * 100 # %

def set_map(ax1):
  ax1.set_extent(extents, crs=mrc)
  ax1.pcolormesh(lon_bound, lat_bound, np.ma.masked_where(i_mask, i_mask), cmap=plt.cm.Greys, vmin=-1, vmax=2, zorder=99)
  ax1.add_feature(cfeature.LAND, edgecolor='None', facecolor='w', zorder=100)
  ax1.set_facecolor('0.5')
  gl = ax1.gridlines(draw_labels=True)
  gl.top_labels = False
  gl.right_labels = False
  gl.xlines = False
  gl.ylines = False


mrc = ccrs.PlateCarree()
fig1 = plt.figure(figsize=(8, 6))

ax1 = fig1.add_axes([0.08, 0.69, 0.36, 0.27], projection=mrc)
ax2 = fig1.add_axes([0.53, 0.69, 0.36, 0.27], projection=mrc)
ax3 = fig1.add_axes([0.08, 0.37, 0.36, 0.27], projection=mrc)
ax4 = fig1.add_axes([0.53, 0.37, 0.36, 0.27], projection=mrc)
ax5 = fig1.add_axes([0.08, 0.05, 0.33, 0.27], projection=mrc)
ax6 = fig1.add_axes([0.56, 0.05, 0.33, 0.27], projection=mrc)

cax1 = fig1.add_axes([0.91, 0.73, 0.01, 0.19])
cax2 = fig1.add_axes([0.91, 0.41, 0.01, 0.19])
cax3 = fig1.add_axes([0.42, 0.09, 0.01, 0.19])
cax4 = fig1.add_axes([0.91, 0.09, 0.01, 0.19])
extents = np.array([-8.5, 4.0, 55, 62])
my_cm = plt.cm.plasma

set_map(ax1)
set_map(ax2)
set_map(ax3)
set_map(ax4)
set_map(ax5)
set_map(ax6)

cs1 = ax1.pcolormesh(lon_bound, lat_bound, out_d, cmap=my_cm, vmin=0, vmax=8)
#cbar = plt.colorbar(cs1, cax=cax1)
#cax1.set_xlabel('Out (Source)')

cs2 = ax2.pcolormesh(lon_bound, lat_bound, in_d, cmap=my_cm, vmin=0, vmax=8)
cbar = plt.colorbar(cs2, cax=cax1)
cax1.set_ylabel('In/Out Degree (%)')

cs3 = ax3.pcolormesh(lon_bound, lat_bound, clo_o, cmap=my_cm, vmin=0, vmax=8)
#cbar = plt.colorbar(cs3, cax=cax3, orientation='horizontal')
#cax3.set_xlabel('Out-Closeness')

cs4 = ax4.pcolormesh(lon_bound, lat_bound, clo_i, cmap=my_cm, vmin=0, vmax=8)
cbar = plt.colorbar(cs4, cax=cax2)
cax2.set_ylabel('In/Out Closeness (%)')

cs5 = ax5.pcolormesh(lon_bound, lat_bound, bet_d, cmap=my_cm, vmin=0, vmax=50)
cbar = plt.colorbar(cs5, cax=cax3, extend='max')
cax3.set_ylabel('Betweenness (Norm. %)')

#cs6 = ax6.contourf(lon_a, lat_a, clust, cmap=my_cm, vmin=0, vmax=100)
cs6 = ax6.pcolormesh(lon_bound, lat_bound, clust, cmap=my_cm, vmin=0, vmax=100)
cbar = plt.colorbar(cs6, cax=cax4, extend='max')
cax4.set_ylabel('Clustering (%)')

#def add_regions(ax, ann=0):
#  b1_x = [-0.7, -1.5, -2, -3.8, -5, -1.5, 0.5, 0.5, -0.7]
#  b1_y = [59.6, 58.8, 59.5, 59.5, 60, 61.7, 61.7, 61.2, 59.6]
#  b2_x = [-4.2, -3.8, -2, -1.5, -0.4, -1, -4.2]
#  b2_y = [57.5, 59.5, 59.5, 58.8, 57.8, 57.5, 57.5]
#  b3_x = [-0.4, -1.5, -0.7, 1.2, 1.5, 0.5, -0.4]
#  b3_y = [57.8, 58.8, 59.6, 59.4, 58.5, 57.5, 57.8]
#  b4_x = [0.5, 0.5, 3.2, 3.9, 3.9, 1.5, 1.2, -0.7, 0.5]
#  b4_y = [61.2, 61.7, 61.1, 59.5, 58, 58.5, 59.4, 59.6, 61.2]
#  b5_x = [-1, -0.4, 0.5, 1.5, 3.9, 3.3, 2, 0, -3.5, -3.5, -1]
#  b5_y = [57.5, 57.8, 57.5, 58.5, 58, 57, 56.1, 55.5, 55.5, 57.5, 57.5]

def add_regions(ax, ann=0):
  b1_x = [-0.7, -1.5, -2, -3.8, -5.4, -1.5, 0.5, 0.5, -0.7]
  b1_y = [59.6, 58.8, 59.5, 59.5, 60.1, 61.7, 61.7, 61.2, 59.6]
  b2_x = [-4.2, -3.8, -2, -1.5, -0.4, -1, -4.2]
  b2_y = [55.5, 59.5, 59.5, 58.8, 57.8, 55.5, 55.5]
  b3_x = [-0.4, -1.5, -0.7, 1.2, 2.5, 0.5, -0.4]
  b3_y = [57.8, 58.8, 59.6, 59.1, 58.5, 57.5, 57.8]
  b4_x = [0.5, 0.5, 3.2, 3.9, 3.9, 2.5, 1.2, -0.7, 0.5]
  b4_y = [61.2, 61.7, 61.1, 59.5, 58, 58.5, 59.1, 59.6, 61.2]
  b5_x = [-0.4, 0.5, 2.5, 3.9, 3.3, 2, 0, -1, -0.4]
  b5_y = [57.8, 57.5, 58.5, 58, 57, 56.1, 55.5, 55.5, 57.8]


  ax.plot(b1_x, b1_y, 'k-', zorder=105)
  ax.plot(b2_x, b2_y, 'k-', zorder=105)
  ax.plot(b3_x, b3_y, 'k-', zorder=105)
  ax.plot(b4_x, b4_y, 'k-', zorder=105)
  ax.plot(b5_x, b5_y, 'k-', zorder=105)

  if ann:
    ax.annotate('1', (np.mean(b1_x), np.mean(b1_y)+0.4), xycoords='data', zorder=105)
    ax.annotate('2', (np.mean(b2_x), np.mean(b2_y)), xycoords='data', zorder=105)
    ax.annotate('3', (np.mean(b3_x), np.mean(b3_y)-0.3), xycoords='data', zorder=105)
    ax.annotate('4', (np.mean(b4_x), np.mean(b4_y)), xycoords='data', zorder=105)
    ax.annotate('5', (np.mean(b5_x), np.mean(b5_y)-1), xycoords='data', zorder=105)


#add_regions(ax5)
add_regions(ax6, 1)


ax1.annotate('(a) Out (Source)', (0.05, 0.95), xycoords='axes fraction', bbox=dict(boxstyle="round", fc="w"), zorder=105)
ax2.annotate('(b) In (Sink)', (0.05, 0.95), xycoords='axes fraction', bbox=dict(boxstyle="round", fc="w"), zorder=105)
ax3.annotate('(c) Out (Source)', (0.05, 0.95), xycoords='axes fraction', bbox=dict(boxstyle="round", fc="w"), zorder=105)
ax4.annotate('(d) In (Sink)', (0.05, 0.95), xycoords='axes fraction', bbox=dict(boxstyle="round", fc="w"), zorder=105)
ax5.annotate('(e)', (0.05, 0.95), xycoords='axes fraction', bbox=dict(boxstyle="round", fc="w"), zorder=105)
ax6.annotate('(f) (Sink)', (0.05, 0.95), xycoords='axes fraction', bbox=dict(boxstyle="round", fc="w"), zorder=105)


fig1.savefig('./Figures/between_install.png', dpi=300)



