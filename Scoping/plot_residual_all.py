#!/usr/bin/env python3

""" Plot a surface from an FVCOM model output.

"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits.basemap import Basemap
import datetime as dt
import netCDF4 as nc
import glob
import matplotlib.tri as tri
from PyFVCOM.read import ncread as readFVCOM
from PyFVCOM.read import FileReader
from PyFVCOM.plot import Plotter
from PyFVCOM.grid import elems2nodes
from PyFVCOM.preproc import Model

in_dir = '/scratch/benbar/Processed_Data_V3.02/'
mjas = '/scratch/benbar/JASMIN/Model_Output_V3.02/Hourly/'
fn = sorted(glob.glob(mjas + '*/SSWRS*V3.02*hr*' + 'RE.nc'))

fin = '/projectsa/SSW_RS/FVCOM/input/SSW_Reanalysis/'
fgrd = fin + 'SSW_Hindcast_riv_xe3_grd.dat'
in_nao = '/scratch/benbar/NAO/nao_data.csv'

data = np.load(in_dir + 'res_vel_all.npz', allow_pickle=True)

mag_r = data['mag_d']# year, season, node
dir_r = data['drc_d']
lat = data['latc']
lon = data['lonc']
data.close()

year = np.arange(1993, 2020)

st_date = dt.datetime(1992, 12, 30, 0, 0, 0)
en_date = dt.datetime(1993, 1, 31, 23, 0, 0)

fvg = Model(st_date, en_date, grid=fgrd, 
                      native_coordinates='spherical', zone='30N')
x = fvg.grid.x
y = fvg.grid.y
xc = fvg.grid.xc
yc = fvg.grid.yc

with open(in_nao) as f:
  lines = f.readlines()

nao_index = np.zeros((len(lines) - 3))
nao_date = np.zeros((len(lines) - 3), dtype=object)
for i in range(3, len(lines)):
  part = lines[i].split(',')
  nao_date[i - 3] = dt.datetime.strptime(part[0], '%Y%m')
  nao_index[i - 3] = float(part[1])

nao_mn = np.array([d.month for d in nao_date])
nao_yr = np.array([d.year for d in nao_date])

nao_mn_index = np.zeros((len(np.unique(nao_yr))))
nao_mn_yr = np.zeros((len(np.unique(nao_yr))))
for yr in range(len(nao_mn_yr)):
  ind = ((nao_mn <= 3) & (nao_yr == yr + nao_yr[0])) 
  nao_mn_index[yr] = np.mean(nao_index[ind])
  nao_mn_yr[yr] = yr + nao_yr[0]

# Extract only the first 24 time steps.
dims = {'time':[0, -1]}
#dims = {'time':range(9)}
# List of the variables to extract.
vars = ['lonc', 'latc', 'nv', 'zeta', 
    'ua', 'va', 'Times']
fvfile = FileReader(fn[0], vars)#, dims=dims)
ind = np.where(fvfile.grid.lon > 180) # for Scottish shelf domain
fvfile.grid.lon[ind]=fvfile.grid.lon[ind]-360 # for Scottish shelf domain
ind = np.where(fvfile.grid.lonc > 180) # for Scottish shelf domain
fvfile.grid.lonc[ind]=fvfile.grid.lonc[ind]-360 # for Scottish shelf domain


mag_r = np.ma.masked_where(mag_r == -999, mag_r)
dir_r = np.ma.masked_where(dir_r == -999, dir_r)

mean_mag = np.ma.mean(mag_r, axis=0)


dims = {'time':':24'}
# List of the variables to extract.
vars = ('nv')
FVCOM = readFVCOM(fn[-1], vars, dims=dims)
nv = FVCOM['nv']
triangles = nv.transpose() -1

u_f = elems2nodes(mag_r * np.ma.cos(dir_r), triangles)
v_f = elems2nodes(mag_r * np.ma.sin(dir_r), triangles)
u_f = np.ma.mean(u_f, axis=0)
v_f = np.ma.mean(v_f, axis=0)

def plot_streamlines(ax, x, y, lon, lat, m, triangles, u, v, c='w', lw=0.5):
    lon_g, lat_g, xg, yg = m.makegrid(100, 100, returnxy=True)
    print(np.min(lon_g), np.max(lon_g), np.min(lat_g), np.max(lat_g))
    trio = tri.Triangulation(lon, lat, triangles=np.asarray(triangles))
    interpolator_u = tri.LinearTriInterpolator(trio, u)
    interpolator_v = tri.LinearTriInterpolator(trio, v)

    grid_u = interpolator_u(lon_g, lat_g)
    grid_v = interpolator_v(lon_g, lat_g)

    ax.streamplot(xg, yg, grid_u, grid_v, density=(1, 2), color=c, linewidth=lw, zorder=102)


def add_map(ax):

  extents = np.array((-7, 3, 56, 61)) # lon1, lon2, lat1, lat2

  m1 = Basemap(llcrnrlon=extents[:2].min(),
          llcrnrlat=extents[-2:].min(),
          urcrnrlon=extents[:2].max(),
          urcrnrlat=extents[-2:].max(),
          rsphere=(6378137.00, 6356752.3142),
          resolution='h',
          projection='merc',
          lat_0=extents[-2:].mean(),
          lon_0=extents[:2].mean(),
          lat_ts=extents[-2:].mean(),
          ax=ax)  

  parallels = np.arange(np.floor(extents[2]), np.ceil(extents[3]), 5)  
  meridians = np.arange(np.floor(extents[0]), np.ceil(extents[1]), 5) 
#  m1.drawmapboundary()

#  m1.drawparallels(parallels, labels=[1, 0, 0, 0],
#                  fontsize=10, linewidth=0)
#  m1.drawmeridians(meridians, labels=[0, 0, 0, 1],
#                  fontsize=10, linewidth=0)

  return m1



fig1, axs1 = plt.subplots(5, 6, figsize=(12, 9))  # size in inches
#axs1 = [None] * 30
#axs1[0] = fig1.add_axes([0.08, 0.69, 0.36, 0.27])
cax1 = fig1.add_axes([0.68, 0.17, 0.29, 0.02])
cax2 = fig1.add_axes([0.68, 0.08, 0.29, 0.02])

print(len(axs1))

season = 0 # 0-4 : JFM, AMJ, JAS, OND
my_cm = plt.cm.plasma

c = -1
for i in range(axs1.shape[0]):
  for j in range(axs1.shape[1]):
    c = c + 1
    if c >= 28:
      axs1[i, j].set_visible(False)
      continue
    m1 = add_map(axs1[i, j])
    fvplt = Plotter(fvfile, figure=fig1, axes=axs1[i, j], m=m1, extend='max')

    mx, my = fvplt.mxc, fvplt.myc
    x1, y1 = m1(fvfile.grid.lon, fvfile.grid.lat)

    if c == 0: # First figure is mean
      cs1 = axs1[i, j].tripcolor(x1, y1, triangles, mean_mag[season, :], vmin=0, vmax=0.1, zorder=100)
      plot_streamlines(axs1[i, j], x, y, fvfile.grid.lon, fvfile.grid.lat, m1, triangles, u_f[season, :], v_f[season, :], c='k', lw=1.6)
      plot_streamlines(axs1[i, j], x, y, fvfile.grid.lon, fvfile.grid.lat, m1, triangles, u_f[season, :], v_f[season, :], lw=1)

      plt.colorbar(cs1, cax=cax1, orientation='horizontal', extend='max')
      cax1.set_xlabel('Mean Current (m/s)')

      axs1[i, j].annotate('Mean Current' , (0.05, 0.95), xycoords='axes fraction', bbox=dict(boxstyle="round", fc="w"), zorder=105)

    else:
      print(year[c-1])
      cs1 = axs1[i, j].tripcolor(x1, y1, triangles, mag_r[c-1, season, :] - mean_mag[season, :], cmap=my_cm, vmin=-0.04, vmax=0.04, zorder=100)
    
      if c == 1:
        plt.colorbar(cs1, cax=cax2, orientation='horizontal', extend='both')
        cax2.set_xlabel('Residual Current Anom. (m/s)')
    #ua_unit = 1 * np.ma.cos(dir_r)
    #va_unit = 1 * np.ma.sin(dir_r)
      n = 100

    #u_f = elems2nodes(mag_r * np.ma.cos(dir_r), triangles)
    #v_f = elems2nodes(mag_r * np.ma.sin(dir_r), triangles)

    #plot_streamlines(axs1[i, j], x, y, fvfile.grid.lon, fvfile.grid.lat, 
    #    m1, triangles,
    #    u_f[0, 0, :], v_f[0, 0, :])

      cmap = matplotlib.cm.get_cmap('bwr')
      norm = matplotlib.colors.Normalize(
          vmin=np.min(nao_mn_index), vmax=np.max(nao_mn_index))
      m = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)

      nao_str = '{:.2f}'.format(nao_mn_index[nao_mn_yr == year[c-1]][0])
      axs1[i, j].annotate(str(year[c-1]) , (0.05, 0.95), xycoords='axes fraction', bbox=dict(boxstyle="round", fc="w"), zorder=105)
      axs1[i, j].annotate(nao_str , (0.5, 0.95), xycoords='axes fraction', bbox=dict(boxstyle="round", fc=m.to_rgba(nao_mn_index[nao_mn_yr == year[c-1]][0])), zorder=105)

    axs1[i, j].axes.get_xaxis().set_ticks([])
    axs1[i, j].axes.get_yaxis().set_ticks([])


plt.draw()
plt.tight_layout(pad=0.5, w_pad=0.1, h_pad=0.1)

if season == 0:
  s_string = 'winter'
elif season == 2:
  s_string = 'summer'
fig1.savefig('./Figures/residual_vel_V3.02_all_' + s_string + '.png', dpi=300)

