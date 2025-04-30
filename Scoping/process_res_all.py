#!/usr/bin/env python3

""" Plot a surface from an FVCOM model output.

"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from PyFVCOM.read import ncread as readFVCOM
import datetime as dt
import netCDF4 as nc
import glob
import matplotlib.tri as tri
import PyFVCOM as fvcom

if __name__ == '__main__':

  mjas = '/scratch/benbar/JASMIN/Model_Output_V3.02/Daily/'
  fn = []
  for yr in range(1993, 2020):
    fn.append(sorted(glob.glob(mjas + '*/SSWRS*V3.02*dy*' + str(yr) + '*RE.nc')))
  print(len(fn))
  out_dir = '/scratch/benbar/Processed_Data_V3.02/'


  fin = '/projectsa/SSW_RS/FVCOM/input/SSW_Reanalysis/'
  fgrd = fin + 'SSW_Hindcast_riv_xe3_grd.dat'

  # Extract only the first 24 time steps.
  dims = {'time':':10'}
  # List of the variables to extract.
  vars = ('lon', 'lat', 'latc', 'lonc', 'nv', 'zeta', 'temp', 'salinity', 'ua', 'va', 'Times')
  FVCOM = readFVCOM(fn[-1][0], vars, dims=dims)

  # Create the triangulation table array (with Python indexing
  # [zero-based])
  triangles = FVCOM['nv'].transpose() - 1
  # Find the domain extents.
  I=np.where(FVCOM['lon'] > 180) # MICDOM: for Scottish shelf domain
  FVCOM['lon'][I]=FVCOM['lon'][I]-360 # MICDOM: for Scottish shelf domain
  I=np.where(FVCOM['lonc'] > 180) # MICDOM: for Scottish shelf domain
  FVCOM['lonc'][I]=FVCOM['lonc'][I]-360 # MICDOM: for Scottish shelf domain
  extents = np.array((FVCOM['lon'].min(),
                    FVCOM['lon'].max(),
                     FVCOM['lat'].min(),
                     FVCOM['lat'].max()))

  lon = FVCOM['lon']
  lat = FVCOM['lat']
  lonc = FVCOM['lonc']
  latc = FVCOM['latc']

  st_date = dt.datetime(1992, 12, 30, 0, 0, 0)
  en_date = dt.datetime(1993, 1, 31, 23, 0, 0)

  fvg = fvcom.preproc.Model(st_date, en_date, grid=fgrd, 
                      native_coordinates='spherical', zone='30N')
  x = fvg.grid.x
  y = fvg.grid.y
  xc = fvg.grid.xc
  yc = fvg.grid.yc

  #trio = tri.Triangulation(x, y, triangles=np.asarray(triangles))

  vars = ('ua', 'va', 'Times', 'Itime')

  u_sum = np.ma.zeros((4, FVCOM['ua'].shape[1]))
  v_sum = np.ma.zeros((4, FVCOM['ua'].shape[1]))
  u_mean = np.ma.zeros((len(fn), 4, FVCOM['ua'].shape[1]))
  v_mean = np.ma.zeros((len(fn), 4, FVCOM['ua'].shape[1]))

  ntime1 = 0
  ntime2 = 0
  ntime3 = 0
  ntime4 = 0
  ref = dt.datetime(1858, 11, 17)

  for j in range(len(fn)):
    for i in range(len(fn[j])):
      print(j / len(fn) *100, '%')
      yr = 1993 + j
      FVCOM = readFVCOM(fn[j][i], vars, dims=dims)

      date = np.zeros((len(FVCOM['Itime'][:])), dtype=object)
      mn = np.zeros((len(FVCOM['Itime'][:])))
      for d in range(len(date)):
        date[d] = ref + dt.timedelta(days=int(FVCOM['Itime'][d]))
        mn[d] = date[d].month
      #print(FVCOM['ua'][:].shape)

      #u_app = FVCOM['u'][:, :, :] # time, depth, elem
      #v_app = FVCOM['v'][:, :, :]
      dind1 = (mn >= 1) & (mn <= 3)
      dind2 = (mn >= 4) & (mn <= 6)
      dind3 = (mn >= 7) & (mn <= 9)
      dind4 = (mn >= 10) & (mn <= 12)
      u_sum[0, :] = u_sum[0, :] + np.ma.sum(FVCOM['ua'][dind1], axis = 0)
      v_sum[0, :] = v_sum[0, :] + np.ma.sum(FVCOM['va'][dind1], axis = 0)
      u_sum[1, :] = u_sum[1, :] + np.ma.sum(FVCOM['ua'][dind2], axis = 0)
      v_sum[1, :] = v_sum[1, :] + np.ma.sum(FVCOM['va'][dind2], axis = 0)
      u_sum[2, :] = u_sum[2, :] + np.ma.sum(FVCOM['ua'][dind3], axis = 0)
      v_sum[2, :] = v_sum[2, :] + np.ma.sum(FVCOM['va'][dind3], axis = 0)
      u_sum[3, :] = u_sum[3, :] + np.ma.sum(FVCOM['ua'][dind4], axis = 0)
      v_sum[3, :] = v_sum[3, :] + np.ma.sum(FVCOM['va'][dind4], axis = 0)
      ntime1 = ntime1 + np.sum(dind1)
      ntime2 = ntime2 + np.sum(dind2)
      ntime3 = ntime3 + np.sum(dind3)
      ntime4 = ntime4 + np.sum(dind4)

    u_mean[j, 0, :] = u_sum[0, :] / ntime1
    v_mean[j, 0, :] = v_sum[0, :] / ntime1
    u_mean[j, 1, :] = u_sum[1, :] / ntime2
    v_mean[j, 1, :] = v_sum[1, :] / ntime2
    u_mean[j, 2, :] = u_sum[2, :] / ntime3
    v_mean[j, 2, :] = v_sum[2, :] / ntime3
    u_mean[j, 3, :] = u_sum[3, :] / ntime4
    v_mean[j, 3, :] = v_sum[3, :] / ntime4
    u_sum = np.ma.zeros((4, FVCOM['ua'].shape[1]))
    v_sum = np.ma.zeros((4, FVCOM['ua'].shape[1]))
    ntime1 = 0
    ntime2 = 0
    ntime3 = 0
    ntime4 = 0


  mag_d = (u_mean ** 2 + v_mean ** 2) ** 0.5
  drc_d = np.ma.arctan2(v_mean, u_mean) # y, x

  np.savez(out_dir + 'res_vel_all.npz', lonc=lonc, latc=latc, mag_d=mag_d, drc_d=drc_d)


