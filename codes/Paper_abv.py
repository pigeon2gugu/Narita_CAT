from netCDF4 import Dataset
import numpy as np
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


n = 2
titles= ['(a) 20181028 12UTC', '(b) 20181029 12UTC', '(c) 20181030 12UTC']

fig = plt.figure(figsize=(15,10))
for i in range(3):
        plt.subplot(3, 4, i+1)
        m = Basemap(llcrnrlat=20,urcrnrlat=60,\
        llcrnrlon=100,urcrnrlon=160, resolution='i')#,resolution='c')
        plt.title(titles[i], fontsize=12, loc='left')
        m.drawcoastlines()
        m.drawcountries()
        m.drawparallels(np.arange(20,70,10),labels=[1,0,0,0], fontsize = 7)
        m.drawmeridians(np.arange(100,170,10),labels=[0,0,0,1], fontsize = 7)
        m.fillcontinents(alpha=0)
        m.drawmapboundary(fill_color='w')
        lon_2d1, lat_2d1 = np.meshgrid(lon, lat)
        x, y = m(lon_2d1, lat_2d1)
        x1, y1 = m(lon[562], lat[215]) #CAT incounter
	clevs1= np.arange(int(4720),int(5960+60),60)
	gph500[n,:,:][ gph500[n,:,:] > max(clevs1) ] = max(clevs1)
	gph500[n,:,:][ gph500[n,:,:] < min(clevs1) ] = min(clevs1)
	m.plot(x1, y1, color='r', marker='*')
	cmap = mpl.cm.Blues
	cs=m.contour(x,y,gph500[n,:,:],clevs1,colors='red', linewidths = 0.3, linestyles='solid')
        colormesh = m.pcolormesh(x, y, vo[i,:,:], vmin = 6, vmax = 30, cmap=cmap)
	cbar=m.colorbar(colormesh)
	plt.clabel(cs, [int(5380), int(5860)], fmt='%3.0f', colors='red')
        cbar.set_label('Absolute Vorticity $(10^{-5}/s)$', fontsize = 7)
	n = n+4

plt.subplots_adjust(wspace = 0.42, hspace = 0.30)
plt.show()
#fig = plt.savefig('./figures/Paper_abv')
