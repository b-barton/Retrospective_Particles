import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from pylag.processing.coordinate import lonlat_from_utm, utm_from_lonlat
import glob

out_dir = '/dssgfs01/scratch/benbar/Particles/Main_Run/Merged_Output/'

yr = '2017'

thresh = 100

summer = 1

fnames = sorted(glob.glob(out_dir + 'connect_2d_' + yr + '*.npz'))
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

print(x_bound.shape)
mn = np.array([d.month for d in date_full])
if summer:
  connect_m = np.ma.mean(connect_full[(mn > 6) & (mn <= 9), :, :], axis=0)
else:
  connect_m = np.ma.mean(connect_full[mn <= 3, :, :], axis=0)

#count_map = np.ma.masked_where(mask, count_map)

lon, lat = lonlat_from_utm(x, y, epsg_code='32630')
lon = np.array(lon)
lat = np.array(lat)
print(np.min(lon), np.max(lon))

s = 2
lon_a, lat_a = lonlat_from_utm(x_coord, y_coord, epsg_code='32630')

n_point = connect_m.shape[0]

dg = nx.DiGraph()
pos = {}
for i in range(n_point):
  dg.add_node(i, x=x[i], y=y[i]) 
  pos[i] = (lon[i], lat[i])

for i in range(n_point):
  for j in range(n_point):
    if (i == j) :
      continue
    if connect_m[i, j] > thresh:
      dg.add_edge(i, j)


names = ['1', '2', '3', '4', '5', '6']
lonc = np.array([-3.111, -2.395, -1.175, 1.642, 1.403, 1.265])
latc = np.array([58.025, 58.137, 58.180, 57.599, 58.047, 58.478])
xc, yc, _ = utm_from_lonlat(lonc, latc, epsg_code='32630')

x_ind = ((xc - x_bound[0, 0]) / (dx * 2 * s)).astype(int)
y_ind = ((yc - y_bound[0, 0]) / (dy * 2 * s)).astype(int)

lon_bound, lat_bound = lonlat_from_utm(x_bound, y_bound, epsg_code='32630')
lon_bound = lon_bound[::s, ::s]
lat_bound = lat_bound[::s, ::s]
box_list_x = []
box_list_y = []
for i in range(len(x_ind)):
  box_list_x.append(np.array([lon_bound[y_ind[i], x_ind[i]], 
                              lon_bound[y_ind[i] + 1, x_ind[i]], 
                              lon_bound[y_ind[i] + 1, x_ind[i] + 1], 
                              lon_bound[y_ind[i], x_ind[i] + 1], 
                              lon_bound[y_ind[i], x_ind[i]]]))
  box_list_y.append(np.array([lat_bound[y_ind[i], x_ind[i]], 
                              lat_bound[y_ind[i] + 1, x_ind[i]], 
                              lat_bound[y_ind[i] + 1, x_ind[i] + 1], 
                              lat_bound[y_ind[i], x_ind[i] + 1], 
                              lat_bound[y_ind[i], x_ind[i]]]))

def plot_box(ax, box_list_x, box_list_y):
  for i in range(len(box_list_x)):
    ax.plot(box_list_x[i], box_list_y[i], '-k', lw=1, zorder=105)


out_degree = np.array(list(zip(*dg.out_degree()))[1])
in_degree = np.array(list(zip(*dg.in_degree()))[1])

bet = np.array(list(nx.betweenness_centrality(dg, normalized=True).values()))
close_in = np.array(list(nx.closeness_centrality(dg).values())) # incoming
close_out = np.array(list(nx.closeness_centrality(dg.reverse()).values())) # outgoing
#eig = nx.eigenvector_centrality(dg)

def re_map(u_mask, mask, var):
  var_m = np.ma.zeros(mask.flatten().shape)
  var_m[u_mask] = var
  var_a = np.reshape(var_m, mask.shape)
  var_a = np.ma.masked_where(mask, var_a)
  return var_a

in_d = (re_map(u_mask, mask, in_degree) / (n_point - 1)) * 100
out_d = (re_map(u_mask, mask, out_degree) / (n_point - 1)) * 100
bet_d = re_map(u_mask, mask, bet) *100 # Already normalised #/ ((n_point - 1) * (n_point - 2)) # rescale
#bet_d = ((bet_d - np.min(bet_d)) / (np.max(bet_d) - np.min(bet_d))) * 100
clo_i = re_map(u_mask, mask, close_in)
clo_o = re_map(u_mask, mask, close_out)

def set_map(ax1):
  ax1.set_extent(extents, crs=mrc)
  ax1.add_feature(cfeature.LAND, edgecolor='None', facecolor='w', zorder=100)
  ax1.set_facecolor('0.5')
  gl = ax1.gridlines(draw_labels=True)
  gl.top_labels = False
  gl.right_labels = False
  gl.xlines = False
  gl.ylines = False

mrc = ccrs.PlateCarree()
fig1 = plt.figure(figsize=(8, 8))
ax1 = fig1.add_axes([0.1, 0.1, 0.75, 0.8], projection=mrc)
cax1 = fig1.add_axes([0.9, 0.3, 0.01, 0.4])
extents = np.array([-8.5, 4.0, 55, 62])
set_map(ax1)

my_cm = plt.cm.plasma
nx.draw(dg.to_undirected(), pos, ax=ax1, node_size=10, width=0.1, alpha=0.5)
#nx.draw_networkx_nodes(dg)
#nx.draw_networkx_edges(dg)
nx.draw_networkx_nodes(dg, pos, ax=ax1, node_color=close_in, node_size=30, vmin=0, vmax=0.3, cmap=my_cm)
sm = plt.cm.ScalarMappable(cmap=my_cm, norm=plt.Normalize(vmin=0, vmax=0.3))
sm._A = []
cbar = plt.colorbar(sm, cax=cax1)
cax1.set_ylabel('Closeness')

plot_box(ax1, box_list_x, box_list_y)


mrc = ccrs.PlateCarree()
fig2 = plt.figure(figsize=(8, 6))

ax1 = fig2.add_axes([0.08, 0.69, 0.36, 0.27], projection=mrc)
ax2 = fig2.add_axes([0.53, 0.69, 0.36, 0.27], projection=mrc)
ax3 = fig2.add_axes([0.08, 0.37, 0.36, 0.27], projection=mrc)
ax4 = fig2.add_axes([0.53, 0.37, 0.36, 0.27], projection=mrc)
ax5 = fig2.add_axes([0.08, 0.05, 0.36, 0.27], projection=mrc)
#ax6 = fig2.add_axes([0.53, 0.05, 0.36, 0.27], projection=mrc)

cax1 = fig2.add_axes([0.91, 0.73, 0.01, 0.19])
cax2 = fig2.add_axes([0.91, 0.41, 0.01, 0.19])
cax3 = fig2.add_axes([0.46, 0.09, 0.01, 0.19])

set_map(ax1)
set_map(ax2)
set_map(ax3)
set_map(ax4)
set_map(ax5)

#plot_box(ax1, box_list_x, box_list_y)
#plot_box(ax2, box_list_x, box_list_y)
#plot_box(ax3, box_list_x, box_list_y)
#plot_box(ax4, box_list_x, box_list_y)
#plot_box(ax5, box_list_x, box_list_y)

cs1 = ax1.pcolormesh(lon_bound, lat_bound, out_d, cmap=my_cm, vmin=0, vmax=20)
#cbar = plt.colorbar(cs1, cax=cax1)
#cax1.set_xlabel('Out (Source)')

cs2 = ax2.pcolormesh(lon_bound, lat_bound, in_d, cmap=my_cm, vmin=0, vmax=20)
cbar = plt.colorbar(cs2, cax=cax1)
cax1.set_ylabel('In/Out Degree (%)')

cs3 = ax3.pcolormesh(lon_bound, lat_bound, clo_o, cmap=my_cm, vmin=0, vmax=0.4)
#cbar = plt.colorbar(cs3, cax=cax3, orientation='horizontal')
#cax3.set_xlabel('Out-Closeness')

cs4 = ax4.pcolormesh(lon_bound, lat_bound, clo_i, cmap=my_cm, vmin=0, vmax=0.4)
cbar = plt.colorbar(cs4, cax=cax2)
cax2.set_ylabel('In/Out Closeness (degree)')

cs5 = ax5.pcolormesh(lon_bound, lat_bound, bet_d, cmap=my_cm, vmin=0, vmax=3)
cbar = plt.colorbar(cs5, cax=cax3, extend='max')
cax3.set_ylabel('Betweenness Norm. (%)')

ax1.annotate('(a) Out-Degree', (0.05, 0.95), xycoords='axes fraction', bbox=dict(boxstyle="round", fc="w"), zorder=105)
ax2.annotate('(b) In-Degree', (0.05, 0.95), xycoords='axes fraction', bbox=dict(boxstyle="round", fc="w"), zorder=105)
ax3.annotate('(c) Out-Closeness', (0.05, 0.95), xycoords='axes fraction', bbox=dict(boxstyle="round", fc="w"), zorder=105)
ax4.annotate('(d) In-Closeness', (0.05, 0.95), xycoords='axes fraction', bbox=dict(boxstyle="round", fc="w"), zorder=105)
ax5.annotate('(e) Betweenness', (0.05, 0.95), xycoords='axes fraction', bbox=dict(boxstyle="round", fc="w"), zorder=105)


fig3 = plt.figure(figsize=(8, 8))
ax1 = fig3.add_axes([0.1, 0.1, 0.75, 0.8], projection=mrc)
cax1 = fig3.add_axes([0.86, 0.3, 0.01, 0.4])
set_map(ax1)

#cs1 = ax1.pcolormesh(lon_bound, lat_bound, count_map, cmap=my_cm)
#cbar = plt.colorbar(cs1, cax=cax1)
cax1.set_ylabel('Number of Particles')

plot_box(ax1, box_list_x, box_list_y)


fig1.savefig('./Figures/network_graph_' + yr + '.png', dpi=300)
if summer:
  fig2.savefig('./Figures/between_' + yr + '_JAS' + '_' + str(thresh) + '.png', dpi=300)
else:
  fig2.savefig('./Figures/between_' + yr + '_' + str(thresh) + '.png', dpi=300)
#fig3.savefig('./Figures/start_particles.png')




