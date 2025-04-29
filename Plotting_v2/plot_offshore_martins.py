import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import glob
import fiona
import datetime as dt
from shapely.geometry import shape
from pylag.processing.coordinate import lonlat_from_utm, utm_from_lonlat
#import geopandas as gpd
#import shapefile as shp

fdata = '/dssgfs01/scratch/benbar/Offshore_Structures/'

f_fix = fdata + 'Martins_et_al_2023_Dataset/Data_283897180/Fixed_Platforms_230321_NorthSeaonly.shp'
f_float = fdata + 'Martins_et_al_2023_Dataset/Data_283897180/Floating_Platforms_230321_NorthSeaonly.shp'
f_renew = fdata + 'Martins_et_al_2023_Dataset/Data_283897180/Windturbines_230321_NorthSeaonly.shp'

out_dir = '/dssgfs01/scratch/benbar/Particles/Main_Run/Merged_Output/'


pres_abs_fix = []
geom_fix = []

with fiona.open(f_fix) as shapefile:
  # Iterate over the records
  print(shapefile[0])
  print(shapefile[0]['properties']['PRESENT_AB'])
  for record in shapefile:
    # Get the geometry from the record
    pres_abs_fix.append(shapefile[0]['properties']['PRESENT_AB'] == 'PRESENT')
    geom_fix.append(shape(record['geometry']))

pres_abs_float = []
geom_float = []

with fiona.open(f_float) as shapefile:
  # Iterate over the records
  print(shapefile[0])
  print(shapefile[0]['properties']['PRESENT_AB'])
  for record in shapefile:
    # Get the geometry from the record
    pres_abs_float.append(shapefile[0]['properties']['PRESENT_AB'] == 'PRESENT')
    geom_float.append(shape(record['geometry']))

pres_abs_renew = []
geom_renew = []

yr_now = dt.datetime.now().year

with fiona.open(f_renew) as shapefile:
  # Iterate over the records
  print(shapefile[0])
  print(shapefile[0]['properties']['YEAR_CONST'])
  for record in shapefile:
    # Get the geometry from the record
    pres_abs_renew.append(shapefile[0]['properties']['YEAR_CONST'] <= yr_now)
    geom_renew.append(shape(record['geometry']))


yr = '2017'

fnames = sorted(glob.glob(out_dir + 'connect_2d_' + yr + '*.npz'))
print(fnames)

data = np.load(fnames[0], allow_pickle=True)
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
data.close()


data = np.load(out_dir + 'map_region.npz')
map_region = data['map_region']
lon_a = data['lon_a']
lat_a = data['lat_a']
data.close()

map_region = np.ma.masked_where(map_region == -999, map_region)

def set_map(ax1):
  ax1.set_extent(extents, crs=mrc)
  ax1.add_feature(land_10m, edgecolor='0.5', facecolor=cfeature.COLORS['land'], zorder=100)
#  ax1.set_facecolor('0.5')
  gl = ax1.gridlines(draw_labels=True)
  gl.top_labels = False
  gl.right_labels = False
  gl.xlines = False
  gl.ylines = False

land_10m = cfeature.NaturalEarthFeature('physical', 'land', '10m',
                                        edgecolor='face',
                                        facecolor=cfeature.COLORS['land'])

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




mrc = ccrs.PlateCarree()
fig1 = plt.figure(figsize=(8, 4))
ax1 = fig1.add_axes([0.1, 0.1, 0.75, 0.8], projection=mrc)
#cax1 = fig1.add_axes([0.9, 0.3, 0.02, 0.4])
extents = np.array([-8.5, 6.0, 55, 62])
set_map(ax1)

lon_inst = []
lat_inst = []

for i in range(len(geom_fix)):
  if pres_abs_fix[i]:
    x, y = geom_fix[i].xy
    lon_inst.extend(x)
    lat_inst.extend(y)
    ax1.plot(x, y, 'ok', ms=3)

for i in range(len(geom_float)):
  if pres_abs_float[i]:
    x, y = geom_float[i].xy
    lon_inst.extend(x)
    lat_inst.extend(y)
    ax1.plot(x, y, 'ok', ms=3)

for i in range(len(geom_renew)):
  if pres_abs_renew[i]:
    x, y = geom_renew[i].xy
    lon_inst.extend(x)
    lat_inst.extend(y)
    ax1.plot(x, y, 'ok', ms=3)


#add_regions(ax1)

lon_inst = np.array(lon_inst)
lat_inst = np.array(lat_inst)
print(lon_inst.shape)

keep = np.invert((lon_inst > -7) & (lon_inst < -5) 
    & (lat_inst > 55) & (lat_inst < 57))
lon_inst = lon_inst[keep]
lat_inst = lat_inst[keep]


s = 2
xc, yc, _ = utm_from_lonlat(lon_inst, lat_inst, epsg_code='32630')

x_ind = ((xc - x_bound[0, 0]) / (dx * 2 * s)).astype(int)
y_ind = ((yc - y_bound[0, 0]) / (dy * 2 * s)).astype(int)

print(mask.shape)
condition = (x_ind < 0) | (x_ind >= mask.shape[0]) | (y_ind < 0) | (y_ind >= mask.shape[0])
x_ind[condition] = -1
y_ind[condition] = -1

install = np.zeros_like(mask)
for i in range(len(x_ind)):
  if x_ind[i] == -1:
    continue
  install[y_ind[i], x_ind[i]] = 1

install = np.ma.masked_where(mask, install)
lon_bound, lat_bound = lonlat_from_utm(x_bound, y_bound, epsg_code='32630')
lon_bound = lon_bound[::s, ::s]
lat_bound = lat_bound[::s, ::s]
#ax1.pcolormesh(lon_bound, lat_bound, install)

my_cm = plt.cm.plasma
cmaplist = [my_cm(i) for i in range(my_cm.N)]
# force the first color entry to be grey
cmaplist[0] = (.75, .75, .75, 1.0)
# create the new map
my_cm = colors.LinearSegmentedColormap.from_list(
    'Custom cmap', cmaplist, my_cm.N)

map_region1 = map_region*1
map_region1[map_region >= 1] = 4
norm = colors.BoundaryNorm(np.arange(-0.5,6), my_cm.N) 
cs1 = ax1.pcolormesh(lon_bound, lat_bound, map_region1, cmap=my_cm, norm=norm)
#cs1 = ax1.pcolormesh(lon_bound, lat_bound, map_region, cmap=my_cm, norm=norm)

#cbar = plt.colorbar(cs1, cax=cax1, ticks=np.arange(0, 6))
#cax1.set_ylabel('Installation Module')



np.savez(out_dir + 'installation_mask_martins.npz', lon=lon_bound, lat=lat_bound, i_mask=install)

fig1.savefig('./Figures/installations_martins.png', dpi=300)



