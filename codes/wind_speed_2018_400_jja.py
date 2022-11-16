from netCDF4 import Dataset
import numpy as np
import os
os.environ['PROJ_LIB'] = '/usr/local/python/2.7/share/proj'
from mpl_toolkits.basemap import Basemap, shiftgrid
import matplotlib.pyplot as plt
from sys import exit
import matplotlib as mpl
import scipy.ndimage as ndimage
import matplotlib.colors
import metpy.calc as mpcalc

lon = Dataset('../../../data/ERA5/2018/ERA5.100_500hPa.201812.00.nc').variables['g0_lon_3'][:]
lat = Dataset('../../../data/ERA5/2018/ERA5.100_500hPa.201812.00.nc').variables['g0_lat_2'][:]
lv = Dataset('../../../data/ERA5/2018/ERA5.100_500hPa.201812.00.nc').variables['lv_ISBL1'][:]

u_jja = np.zeros((1,721,1440))
v_jja = np.zeros((1,721,1440))

jja = ['06', '07', '08']
time = ['00', '06', '12', '18']

for i in range(3):
	print i
	for j in range(4):
		u_jja = np.append(u_jja, Dataset('../../../data/ERA5/2018/ERA5.100_500hPa.2018'+jja[i]+'.'+time[j]+'.nc').variables['U_GDS0_ISBL'][:,9,:,:], axis = 0)
		v_jja = np.append(v_jja, Dataset('../../../data/ERA5/2018/ERA5.100_500hPa.2018'+jja[i]+'.'+time[j]+'.nc').variables['V_GDS0_ISBL'][:,9,:,:], axis = 0)
		print j 

print 'end reading variables'

wind_jja = ( (u_jja**(2.)) + (v_jja**(2.)) )**(1./2.)


dtr = np.pi/180.
theta = lat*dtr
dthe = 0.25
Re = 6371000.
dx = Re*np.cos(theta)*dthe*dtr
dy = Re*dthe*dtr


wind_jja = np.sum(wind_jja, axis = 0 )/(92. * 4.)

print 'end calculating meean wind speed'

titles= ['(c) 2018 JJA']
fig = plt.figure(figsize=(15,10))

m = Basemap(lon_0=270, boundinglat=0, projection='npstere',round=True)#,resolution='c')
plt.title(titles[0], fontsize=12, loc='left')
m.drawcoastlines()
m.drawcountries()
m.drawparallels(np.arange(0,81,20))
m.drawmeridians(np.arange(0,360,60))
m.fillcontinents(alpha=0)
m.drawmapboundary(fill_color='w')
lon_2d1, lat_2d1 = np.meshgrid(lon, lat)
x, y = m(lon_2d1, lat_2d1)
cmap = mpl.cm.Blues
clevs2 = np.array([20,25,30,40])
cs1 = m.contour(x,y,wind_jja,clevs2, colors='black', linewidths = 0.7, linestyles='solid')
plt.clabel(cs1,fmt='%3.0f',colors='black', fontsize = 7)

fig = plt.savefig('./figures/wind_speed__2018_400_jja')


