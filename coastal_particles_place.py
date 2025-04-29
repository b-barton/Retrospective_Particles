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


# Read in the bathymetry
ds = Dataset(grid_metrics_file_name, 'r')
bathy = -ds.variables['h'][:]
x_m = ds.variables['x'][:]
y_m = ds.variables['y'][:]
nv = ds.variables['nv'][:]
ds.close()
del(ds)

triangles = nv.transpose()

bathy_mod = amm15_dir + 'bathy_meter_ORIGINAL_AMM15.nc'

with Dataset(bathy_mod, 'r') as nc_fid:
  nc_lon = nc_fid.variables['lon'][:] 
  nc_lat = nc_fid.variables['lat'][:]
  nc_bathy = -nc_fid.variables['Bathymetry'][:] 

# Coastal distance file
in_data = '/dssgfs01/scratch/benbar/Processed_Data/'
ds = np.load(in_data + 'coast_distance.npz')
dist = ds['dist_c']
lat_grid = ds['lat']
lon_grid = ds['lon']
ds.close()

coast = (dist == 0)
c_lat = lat_grid[coast]
c_lon = lon_grid[coast]

# The group ID of this particle set
#group_id = 1
group_id = np.arange(1, 50)

# Lat and lon coordiantes for the centre of the release zone
lat = 58.7
lon = -2.0

# Convert to UTM coordinates
easting, northing, _ = utm_from_lonlat([lon], [lat], epsg_code='32630')
nc_x, nc_y, _ = utm_from_lonlat(nc_lon.flatten(), nc_lat.flatten(), 
                epsg_code='32630')
nc_x = nc_x.reshape(nc_lon.shape)
nc_y = nc_y.reshape(nc_lat.shape)

# Release zone radius (m)
radius = 350000.0

# Target number of particles to be released. Only a target,
# since we are evenly distributing particles in the release
# zone, which has no unique solution.
n_particles_target = 1800 # 476672 
#n_particles_target = 2000 # 514304
#n_particles_target = 3000 # 878080

# number of processors
n_proc = 256

# Release depths
depth_below_surface = np.arange(0, -200, -10)
#depth_below_surface = np.arange(0, -200, -20)

release_zones = []

for i in range(len(depth_below_surface)):
  for j in range(len(group_id)):
    # Create the release zone
    release_zone1 = create_release_zone(group_id = group_id[j],
                                           radius = radius,
                                           centre = [easting, northing],
                                           n_particles = n_particles_target,
                                           depth = depth_below_surface[i],
                                           random = False)

    release_zones.append(release_zone1)

# Convert utm coords to degrees
lons = []
lats = []
for rz in range(len(release_zones)):
  lons1, lats1 = lonlat_from_utm(release_zones[rz].get_eastings(),
                             release_zones[rz].get_northings(),
                             epsg_code='32630')
  lons.extend(lons1)
  lats.extend(lats1)


# Get the actual number of particles
n_particles = np.sum([r.get_number_of_particles() for r in release_zones])
print(n_particles)
print('North61', np.sum((np.array(lats) > 61).astype(int)))

# Mask points on land including in bathymetry

# loop over 10 m depth intervals in bathy and making polygons

step = 5 #depth_below_surface[0] - depth_below_surface[1]
max_d = 200

off_bathy = []

for i in range(int(max_d / step)):
  # Get contours from bathy
  #cs = plt.tricontour(x_m, y_m, triangles, alt_bathy, [i * -step])
  cs = plt.contour(nc_x, nc_y, nc_bathy + 5, [i * -step]) 
  poly_list = []
  for p in cs.collections[0].get_paths():
    v = p.vertices
    #c_x = v[:,0]
    #c_y = v[:,1]
    poly_list.append(shapely.geometry.Polygon(v).buffer(0))
  poly1 = unary_union(poly_list)
  #for geom in poly1.geoms:
  #  plt.plot(*geom.exterior.xy)
  #plt.show()
  #continue
  for rz in range(len(release_zones)):
    pt_layer = ((np.asarray(release_zones[rz].get_depths()) > ((i+1) * -step))
                  & (np.asarray(release_zones[rz].get_depths()) <= (i * -step)))
    print(i * -step, np.sum(pt_layer), len(poly1.geoms))

    part_points = shapely.geometry.MultiPoint(np.vstack((
                      np.asarray(release_zones[rz].get_eastings())[pt_layer],
                      np.asarray(release_zones[rz].get_northings())[pt_layer])).T)

    zone_bool = np.zeros((release_zones[rz].get_number_of_particles()), dtype=bool)
    off_bathy.append(zone_bool)
    off_tmp = np.zeros((len(part_points.geoms)), dtype=bool)
    for j in range(np.sum(pt_layer)):
      for p in range(len(poly1.geoms)):
        p_in = poly1.geoms[p].contains(part_points.geoms[j])
        if p_in:
          off_tmp[j] = off_tmp[j] | p_in
          #print(p_in, j)
          break
    
    off_bathy[rz][pt_layer] = off_tmp[:]

print(np.sum(np.array(off_bathy)))


# Mask points that are outside the model grid

points = shapely.geometry.MultiPoint(np.asarray((x_m, y_m)).T)
poly2 = points.convex_hull

off_grid = []
for rz in range(len(release_zones)):
  n_particles = release_zones[rz].get_number_of_particles()
  zone_bool = np.zeros((n_particles), dtype=bool)
  off_grid.append(zone_bool)
 
  part_points = shapely.geometry.MultiPoint(np.asarray((
                            release_zones[rz].get_eastings(),
                            release_zones[rz].get_northings())).T)

  for i in range(n_particles):
    off_grid[rz][i] = poly2.contains(part_points.geoms[i])

# Mask points not on the shelf

land_x, land_y, _ = utm_from_lonlat([17, 17, 17], [62, 55, 49], 
                epsg_code='32630')
fland = np.array([land_x, land_y]).T[::-1, :]
print(fland)

cs = plt.contour(nc_x, nc_y, nc_bathy, [-200]) 
poly_list = []

n_cont = []
for p in cs.collections[0].get_paths():
  n_cont.append(len(p.vertices))
max_cont = np.max(np.array(n_cont))

for p in cs.collections[0].get_paths():
  v = p.vertices
  if len(v) == max_cont:
    v = np.vstack((v, fland))
    print(v[:2, :])
    print(v[-5:, :])
  poly_list.append(shapely.geometry.Polygon(v).buffer(0))

poly3 = unary_union(poly_list)

off_shelf = []
del_mask = []
for rz in range(len(release_zones)):
  n_particles = release_zones[rz].get_number_of_particles()
  zone_bool = np.zeros((n_particles), dtype=bool)
  off_shelf.append(zone_bool)
  del_mask.append(zone_bool)

  part_points = shapely.geometry.MultiPoint(np.asarray((
                            release_zones[rz].get_eastings(),
                            release_zones[rz].get_northings())).T)

  for j in range(n_particles):
    for p in range(len(poly3.geoms)):
      p_in = poly3.geoms[p].contains(part_points.geoms[j])
      if p_in:
        off_shelf[rz][j] = off_shelf[rz][j] | p_in
        break
  rise_x, rise_y, _ = utm_from_lonlat([-5], [60],
                  epsg_code='32630')
  off_shelf[rz] = (off_shelf[rz] 
            & (((np.asarray(release_zones[rz].get_eastings()) < rise_x)
            & (np.asarray(release_zones[rz].get_northings()) > rise_y)) == 0))
 
  del_mask[rz] = ((off_bathy[rz] == 1) | (off_grid[rz] == 0) 
                  | (off_shelf[rz] == 0))

def remove_particles(release_zones, del_mask, nproc):
  release_zone_sub = []
  for rz in range(len(release_zones)):
    parts = release_zones[rz].get_particle_set()
    p_id = release_zones[rz].get_group_id()
    p_radius = release_zones[rz].get_radius()
    p_centre = release_zones[rz].get_centre()
    parts_array = np.array(parts)
    sub_parts = parts_array[del_mask[rz] == 0].tolist()
    if len(sub_parts) % nproc != 0:
      rem = len(sub_parts) % nproc
      sub_parts = sub_parts[:-rem]
    release_zone_sub.append(ReleaseZone(group_id=p_id, 
                                  radius=p_radius, 
                                  centre=p_centre))
    for i in range(len(sub_parts)):
      release_zone_sub[rz].add_particle(
          sub_parts[i][0], sub_parts[i][1], sub_parts[i][2])
  return release_zone_sub

test1 = remove_particles(release_zones, off_bathy, n_proc)
lons = []
lats = []
for rz in range(len(test1)):
  lons1, lats1 = lonlat_from_utm(test1[rz].get_eastings(),
                             test1[rz].get_northings(),
                             epsg_code='32630')
  lons.extend(lons1)
  lats.extend(lats1)
print('North61', np.sum((np.array(lats) > 61).astype(int)))
test2 = remove_particles(release_zones, off_grid, n_proc)
for rz in range(len(test2)):
  lons1, lats1 = lonlat_from_utm(test2[rz].get_eastings(),
                             test2[rz].get_northings(),
                             epsg_code='32630')
  lons.extend(lons1)
  lats.extend(lats1)
print('North61', np.sum((np.array(lats) > 61).astype(int)))
test3 = remove_particles(release_zones, off_shelf, n_proc)
for rz in range(len(test3)):
  lons1, lats1 = lonlat_from_utm(test3[rz].get_eastings(),
                             test3[rz].get_northings(),
                             epsg_code='32630')
  lons.extend(lons1)
  lats.extend(lats1)
print('North61', np.sum((np.array(lats) > 61).astype(int)))

release_zones = remove_particles(release_zones, del_mask, n_proc)
n_particles = np.sum([r.get_number_of_particles() for r in release_zones])

print(n_particles)

lons = []
lats = []
for rz in range(len(release_zones)):
  lons1, lats1 = lonlat_from_utm(release_zones[rz].get_eastings(),
                             release_zones[rz].get_northings(),
                             epsg_code='32630')
  lons.extend(lons1)
  lats.extend(lats1)

# Output filename
file_name = '{}/initial_pos_cart_coast_full.dat'.format(input_dir)

# Write data to file
#create_initial_positions_file_single_group(file_name,
#                                           n_particles,
#                                           group_id,
#                                           lons,
#                                           lats,
#                                           surface_release_zone.get_depths())
create_initial_positions_file_multi_group(file_name, release_zones)


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
