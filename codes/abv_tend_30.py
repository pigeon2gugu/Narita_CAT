from netCDF4 import Dataset
import os
import numpy as np
os.environ['PROJ_LIB'] = '/usr/local/python/2.7/share/proj'
from mpl_toolkits.basemap import Basemap, shiftgrid
import matplotlib.pyplot as plt
from sys import exit
import matplotlib as mpl
import metpy.calc as mpcalc
from metpy.units import units
import scipy.ndimage as ndimage

lon = Dataset('data/2018_10_30_12:00:00.nc').variables['g0_lon_2'][:]
lat = Dataset('data/2018_10_30_12:00:00.nc').variables['g0_lat_1'][:]
gph500 = np.zeros((3,721,1440))
o500 = np.zeros((3,721,1440))
u450 = np.zeros((3,721,1440))
u500 = np.zeros((3,721,1440))
u550 = np.zeros((3,721,1440))
v450 = np.zeros((3,721,1440))
v500 = np.zeros((3,721,1440))
v550 = np.zeros((3,721,1440))
vo450 = np.zeros((3,721,1440))
vo500 = np.zeros((3,721,1440))
vo550 = np.zeros((3,721,1440))

dd = ['28','29','30']

for i  in range(3):
	gph500[i,:,:] = Dataset('data/2018_10_'+dd[i]+'_12:00:00.nc').variables['Z_GDS0_ISBL'][21,:,:]
	u450[i,:,:] = Dataset('data/2018_10_'+dd[i]+'_12:00:00.nc').variables['U_GDS0_ISBL'][20,:,:]
	u500[i,:,:] = Dataset('data/2018_10_'+dd[i]+'_12:00:00.nc').variables['U_GDS0_ISBL'][21,:,:]
	u550[i,:,:] = Dataset('data/2018_10_'+dd[i]+'_12:00:00.nc').variables['U_GDS0_ISBL'][22,:,:]
	v450[i,:,:] = Dataset('data/2018_10_'+dd[i]+'_12:00:00.nc').variables['V_GDS0_ISBL'][20,:,:]
	v500[i,:,:] = Dataset('data/2018_10_'+dd[i]+'_12:00:00.nc').variables['V_GDS0_ISBL'][21,:,:]
	v550[i,:,:] = Dataset('data/2018_10_'+dd[i]+'_12:00:00.nc').variables['V_GDS0_ISBL'][22,:,:]
	o500[i,:,:] = Dataset('data/2018_10_'+dd[i]+'_12:00:00.nc').variables['W_GDS0_ISBL'][21,:,:]
	vo450[i,:,:] = Dataset('data/2018_10_'+dd[i]+'_12:00:00.nc').variables['VO_GDS0_ISBL'][20,:,:]
	vo500[i,:,:] = Dataset('data/2018_10_'+dd[i]+'_12:00:00.nc').variables['VO_GDS0_ISBL'][21,:,:]
	vo550[i,:,:] = Dataset('data/2018_10_'+dd[i]+'_12:00:00.nc').variables['VO_GDS0_ISBL'][22,:,:]

gph500 = gph500/9.806

hadv = np.zeros((3,721,1440))
vadv = np.zeros((3,721,1440))
stre = np.zeros((3,721,1440))
tilt = np.zeros((3,721,1440))
tend = np.zeros((3,721,1440))

dtr = np.pi/180.
theta = lat*dtr
dthe = 0.25
Re = 6371000.
dx = Re*np.cos(theta)*dthe*dtr
dy = Re*dthe*dtr

for i in range(3):
	for j in range(1,len(lat)-1):
        	for k in range(1,len(lon)-1):
			vo450[i,j,k] = vo450[i,j,k] + ((2.)*(2.*np.pi/86400.)*(np.sin(lat[j]*dtr)))
                        vo500[i,j,k] = vo500[i,j,k] + ((2.)*(2.*np.pi/86400.)*(np.sin(lat[j]*dtr)))
                        vo550[i,j,k] = vo550[i,j,k] + ((2.)*(2.*np.pi/86400.)*(np.sin(lat[j]*dtr)))


for i in range(3):
        for j in range(1,len(lat)-1):
                for k in range(1,len(lon)-1):
			hadv[i,j,k] = -( ((u500[i,j,k]) * ((vo500[i,j,k+1]-vo500[i,j,k-1])/(2*dx[j]))) + ((v500[i,j,k]) * ((vo500[i,j+1,k]-vo500[i,j-1,k])/(-2*dy))) )
			vadv[i,j,k] = -((o500[i,j,k]) * ((vo450[i,j,k]-vo550[i,j,k]) / (-10000)))
			stre[i,j,k] = -((vo500[i,j,k]) * (((u500[i,j,k+1]-u500[i,j,k-1]) / (2*dx[j])) + ((v500[i,j+1,k]-v500[i,j-1,k]) / (-2*dy))))
			tilt[i,j,k] = - ((((o500[i,j,k+1]-o500[i,j,k-1]) / (2*dx[j])) * ((v450[i,j,k]-v550[i,j,k]) / (-10000))) - (((o500[i,j+1,k]-o500[i,j-1,k]) / (-2*dy)) * ((u450[i,j,k]-u550[i,j,k]) / (-10000))))

tend = hadv + vadv + stre + tilt


lat1 = lat[::-1]

su, newlon1 = shiftgrid(180., u500, lon, start=False)
sv, newlon1 = shiftgrid(180., v500, lon, start=False)
su = su[:,::-1,:]
sv = sv[:,::-1,:]

hadv = (hadv)*(1000000000.)
vadv = (vadv)*(1000000000.)
stre = (stre)*(1000000000.)
tilt = (tilt)*(1000000000.)
tend = (tend)*(1000000000.)
vo500 = (vo500)*(100000.)

vo500 = ndimage.gaussian_filter(vo500, sigma=3, order=0)
hadv = ndimage.gaussian_filter(hadv, sigma=3, order=0)
vadv = ndimage.gaussian_filter(vadv, sigma=3, order=0)
stre = ndimage.gaussian_filter(stre, sigma=3, order=0)
tilt = ndimage.gaussian_filter(tilt, sigma=3, order=0)
tend = ndimage.gaussian_filter(tend, sigma=3, order=0)

titles= ['(a) Vorticity tendency','(b) Horizontal advection','(c) Vertical advection', '(d) Vortex streatching', '(e) Vortex tilting']

fig = plt.figure(figsize=(15,10))

plt.subplot(3, 2, 1)
m = Basemap(llcrnrlat=20,urcrnrlat=60,llcrnrlon=100,urcrnrlon=160, resolution='c')#,resolution='c')
plt.title(titles[0], fontsize=12, loc='left')
m.drawcoastlines()
m.drawcountries()
m.drawparallels(np.arange(20,70,10),labels=[1,0,0,0], fontsize = 7)
m.drawmeridians(np.arange(100,170,10),labels=[0,0,0,1], fontsize = 7)
m.fillcontinents(alpha=0)
m.drawmapboundary(fill_color='w')
lon_2d1, lat_2d1 = np.meshgrid(lon, lat)
x, y = m(lon_2d1, lat_2d1)
x1, y1 = m(lon[562], lat[215]) #CAT incounter
m.plot(x1, y1, color='r', marker='*')
cmap = mpl.cm.bwr
clevs1= np.arange(20,40,4)
colormesh = m.pcolormesh(x, y, tend[2,:,:], vmin = -10, vmax = 10, cmap=cmap)
cs=m.contour(x,y,vo500[2,:,:], clevs1, colors='black', linewidths = 1.5, linestyles='--')
cbar=m.colorbar(colormesh)
uproj,vproj,xx,yy = m.transform_vector(su[2,:,:],sv[2,:,:],newlon1,lat1,10,10,returnxy=True,masked=True)
x2 , y2 = m(xx, yy)
barbs = m.barbs(x2, y2, uproj, vproj, length=6, pivot='middle', linewidth=0.5)
plt.clabel(cs, fmt='%3.0f', colors='black', fontsize = 7)
cbar.set_label('Vorticity tendency $(10^{-9}/s^{2})$', fontsize = 7)

plt.subplot(3, 2, 2)
m = Basemap(llcrnrlat=20,urcrnrlat=60,llcrnrlon=100,urcrnrlon=160, resolution='c')#,resolution='c')
plt.title(titles[1], fontsize=12, loc='left')
m.drawcoastlines()
m.drawcountries()
m.drawparallels(np.arange(20,70,10),labels=[1,0,0,0], fontsize = 7)
m.drawmeridians(np.arange(100,170,10),labels=[0,0,0,1], fontsize = 7)
m.fillcontinents(alpha=0)
m.drawmapboundary(fill_color='w')
lon_2d1, lat_2d1 = np.meshgrid(lon, lat)
x, y = m(lon_2d1, lat_2d1)
x1, y1 = m(lon[562], lat[215]) #CAT incounter
clevs1= np.arange(20,40,4)
m.plot(x1, y1, color='r', marker='*')
cmap = mpl.cm.bwr
colormesh = m.pcolormesh(x, y, hadv[2,:,:], vmin = -14, vmax = 14, cmap=cmap)
cs=m.contour(x,y,vo500[2,:,:], clevs1, colors='black', linewidths = 1.5, linestyles='--')
cbar=m.colorbar(colormesh)
plt.clabel(cs, fmt='%3.0f', colors='black', fontsize = 7)
cbar.set_label('Horizontal advection $(10^{-9}/s^{2})$', fontsize = 7)

plt.subplot(3, 2, 3)
m = Basemap(llcrnrlat=20,urcrnrlat=60,llcrnrlon=100,urcrnrlon=160, resolution='c')#,resolution='c')
plt.title(titles[2], fontsize=12, loc='left')
m.drawcoastlines()
m.drawcountries()
m.drawparallels(np.arange(20,70,10),labels=[1,0,0,0], fontsize = 7)
m.drawmeridians(np.arange(100,170,10),labels=[0,0,0,1], fontsize = 7)
m.fillcontinents(alpha=0)
m.drawmapboundary(fill_color='w')
lon_2d1, lat_2d1 = np.meshgrid(lon, lat)
x, y = m(lon_2d1, lat_2d1)
x1, y1 = m(lon[562], lat[215]) #CAT incounter
clevs1= np.arange(20,40,4)
m.plot(x1, y1, color='r', marker='*')
cmap = mpl.cm.bwr
colormesh = m.pcolormesh(x, y, vadv[2,:,:], vmin = -6, vmax = 6, cmap=cmap)
cs=m.contour(x,y,vo500[2,:,:], clevs1, colors='black', linewidths = 1.5, linestyles='--')
cbar=m.colorbar(colormesh)
plt.clabel(cs, fmt='%3.0f', colors='black', fontsize = 7)
cbar.set_label('Vertical advection $(10^{-9}/s^{2})$', fontsize = 7)

plt.subplot(3, 2, 4)
m = Basemap(llcrnrlat=20,urcrnrlat=60,llcrnrlon=100,urcrnrlon=160, resolution='c')#,resolution='c')
plt.title(titles[3], fontsize=12, loc='left')
m.drawcoastlines()
m.drawcountries()
m.drawparallels(np.arange(20,70,10),labels=[1,0,0,0], fontsize = 7)
m.drawmeridians(np.arange(100,170,10),labels=[0,0,0,1], fontsize = 7)
m.fillcontinents(alpha=0)
m.drawmapboundary(fill_color='w')
lon_2d1, lat_2d1 = np.meshgrid(lon, lat)
x, y = m(lon_2d1, lat_2d1)
x1, y1 = m(lon[562], lat[215]) #CAT incounter
clevs1= np.arange(20,40,4)
m.plot(x1, y1, color='r', marker='*')
cmap = mpl.cm.bwr
colormesh = m.pcolormesh(x, y, stre[2,:,:], vmin = -6, vmax = 6, cmap=cmap)
cs=m.contour(x,y,vo500[2,:,:], clevs1, colors='black', linewidths = 1.5, linestyles='--')
cbar=m.colorbar(colormesh)
plt.clabel(cs, fmt='%3.0f', colors='black', fontsize = 7)
cbar.set_label('Vortex stretching $(10^{-9}/s^{2})$', fontsize = 7)

plt.subplot(3, 2, 5)
m = Basemap(llcrnrlat=20,urcrnrlat=60,llcrnrlon=100,urcrnrlon=160, resolution='c')#,resolution='c')
plt.title(titles[4], fontsize=12, loc='left')
m.drawcoastlines()
m.drawcountries()
m.drawparallels(np.arange(20,70,10),labels=[1,0,0,0], fontsize = 7)
m.drawmeridians(np.arange(100,170,10),labels=[0,0,0,1], fontsize = 7)
m.fillcontinents(alpha=0)
m.drawmapboundary(fill_color='w')
lon_2d1, lat_2d1 = np.meshgrid(lon, lat)
x, y = m(lon_2d1, lat_2d1)
x1, y1 = m(lon[562], lat[215]) #CAT incounter
clevs1= np.arange(20,40,4)
m.plot(x1, y1, color='r', marker='*')
cmap = mpl.cm.bwr
colormesh = m.pcolormesh(x, y, tilt[2,:,:], vmin = -6, vmax = 6, cmap=cmap)
cs=m.contour(x,y,vo500[2,:,:], clevs1, colors='black', linewidths = 1.5, linestyles='--')
cbar=m.colorbar(colormesh)
plt.clabel(cs, fmt='%3.0f', colors='black', fontsize = 7)
cbar.set_label('Vortex tilting $(10^{-9}/s^{2})$', fontsize = 7)

plt.subplots_adjust(wspace = 0.42, hspace = 0.30)
#plt.show()
fig = plt.savefig('./figures/vorticity tendency_30')
