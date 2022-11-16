from netCDF4 import Dataset
import numpy as np
import os
os.environ['PROJ_LIB'] = '/usr/local/python/2.7/share/proj'
from mpl_toolkits.basemap import Basemap, shiftgrid
import matplotlib.pyplot as plt
from sys import exit
import matplotlib as mpl
import metpy.calc as mpcalc
from metpy.units import units
import scipy.ndimage as ndimage


lon = Dataset('data/2018-10-30_06:00:00.nc').variables['g0_lon_2'][:]
lat = Dataset('data/2018-10-30_06:00:00.nc').variables['g0_lat_1'][:]
gph500 = np.zeros((11,721,1440))
u = np.zeros((11,721,1440))
v = np.zeros((11,721,1440))

k = 0
dd = ['28','29','30']
tt = ['00','06','12','18']

for i  in range(3):
        for j in range(4):
                        if dd[i] == '30' and tt[j] == '18':
                                break
                        else :
                                gph500[k,:,:] = Dataset('data/2018-10-'+dd[i]+'_'+tt[j]+':00:00.nc').variables['Z_GDS0_ISBL'][16,:,:]
				u[k,:,:] = Dataset('data/2018-10-'+dd[i]+'_'+tt[j]+':00:00.nc').variables['U_GDS0_ISBL'][16,:,:]
				v[k,:,:] = Dataset('data/2018-10-'+dd[i]+'_'+tt[j]+':00:00.nc').variables['V_GDS0_ISBL'][16,:,:]
                                k = k+1


gph500 = gph500/9.806

wind = ((u)**2 + (v)**2 )**(1./2.)

dx = np.zeros((721,1440))
dy = np.zeros((721,1440))

for i in range(len(lat)):
        for j in range(1,len(lon)-1):
                dx[i,j] = 2*np.pi*6371000*np.cos(lat[i]*2*np.pi/360.)*(lon[j+1]-lon[j-1])/360.

for i in range(len(lat)):
        dx[i,-1] = dx[i,-2]
        dx[i,0] = dx[i,1]

for i in range(1,len(lat)-1):
        for j in range(len(lon)):
                dy[i,j] = 2*np.pi*6371000*(lat[i+1]-lat[i-1])/360.


for j in range(len(lon)):
        dy[-1,j] = 2*np.pi*6371000*(0.5)/360.
        dy[0,j] = 2*np.pi*6371000*(0.5)/360.

dudx = np.zeros((11,721,1440))
dvdx = np.zeros((11,721,1440))
dudy = np.zeros((11,721,1440))
dvdy = np.zeros((11,721,1440))

for i in range(11):
	for j in range(len(lat)):
        	for k in range(1,len(lon)-1):
               		dudx[i,j,k] = (u[i,j,k+1]-u[i,j,k-1])/(dx[j,k])
                	dvdx[i,j,k] = (v[i,j,k+1]-v[i,j,k-1])/(dx[j,k])

for i in range(11):
	for j in range(len(lat)):
        	dudx[i,j,-1] = (u[i,j,0]-u[i,j,-2])/(dx[j,-1])
        	dvdx[i,j,-1] = (v[i,j,0]-v[i,j,-2])/(dx[j,-1])
        	dudx[i,j,0] = (u[i,j,-1]-u[i,j,1])/(dx[j,0])
        	dvdx[i,j,0] = (v[i,j,-1]-v[i,j,1])/(dx[j,0])

for i in range(11):
	for j in range(1,len(lat)-1):
        	for k in range(len(lon)):
                	dudy[i,j,k] = (u[i,j+1,k]-u[i,j-1,k])/(dy[j,k])
                	dvdy[i,j,k] = (v[i,j+1,k]-v[i,j-1,k])/(dy[j,k])

for i in range(11):
	for j in range(len(lon)):
        	dudy[i,-1,j] = dudy[i,-2,j]
        	dvdy[i,-1,j] = dvdy[i,-2,j]
        	dudy[i,0,j] = dudy[i,1,j]
        	dvdy[i,0,j] = dvdy[i,1,j]

vo = np.zeros((11,721,1440))

for i in range(11):
	for j in range(len(lat)):
        	for k in range(len(lon)):
			vo[i,j,k] = (dvdx[i,j,k]-dudy[i,j,k])+((2.)*(2.*np.pi/86400.)*(np.sin(lat[j]*np.pi/180.)))

for i in range(11):
        for j in range(1,len(lat)-1):
                for k in range(1,len(lon)-1):
                        vo[i,j,k] = (((vo[i,j-1,k] + ((2.)*vo[i,j,k]) + vo[i,j+1,k])/(4.))+((vo[i,j,k-1] + ((2.)*vo[i,j,k]) + vo[i,j,k+1])/(4.)))/(2.)

for i in range(11):
        for j in range(len(lon)):
                vo[i,0,j] = vo[i,1,j]
                vo[i,-1,j] = vo[i,-2,j]

for i in range(11):
        for j in range(len(lat)):
                vo[i,j,0] = vo[i,j,1]
                vo[i,j,-1] = vo[i,j,-2]


vo = (100000)*(vo)
gph500 = np.around(gph500)
gph500 = gph500.astype(int)

vo = ndimage.gaussian_filter(vo, sigma=3, order=0)


m = Basemap(lon_0=270, boundinglat=0, projection='npstere',round=True)#,resolution='c')
m.drawcoastlines()
m.drawcountries()
m.drawparallels(np.arange(0,81,20))
m.drawmeridians(np.arange(0,360,60))
m.fillcontinents(alpha=0)
m.drawmapboundary(fill_color='w')
lon_2d1, lat_2d1 = np.meshgrid(lon, lat)
x, y = m(lon_2d1, lat_2d1)
cmap = mpl.cm.Blues
clevs2 = np.arange(40, 120, 20)
cs=m.contourf(x,y,vo[0,:,:], cmap=cmap)
cs1 = m.contour(x,y,wind[0,:,:],clevs2, colors='red', linewidths = 1, linestyles='solid')
cbar=m.colorbar(cs)
plt.clabel(cs1,fmt='%3.0f',colors='black', fontsize = 7)
cbar.set_label('PVU > 1.5 frequency', fontsize = 7)

plt.show()
