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
gph400 = np.zeros((3,721,1440))
o400 = np.zeros((3,721,1440))
u400 = np.zeros((3,721,1440))
v400 = np.zeros((3,721,1440))
the350 = np.zeros((3,721,1440))
the400 = np.zeros((3,721,1440))
the450 = np.zeros((3,721,1440))
hthe350 = np.zeros((3,721,1440))
hthe400 = np.zeros((3,721,1440))
hthe450 = np.zeros((3,721,1440))


dd = ['28','29','30']

for i  in range(3):
	gph400[i,:,:] = Dataset('data/2018_10_'+dd[i]+'_12:00:00.nc').variables['Z_GDS0_ISBL'][19,:,:]
	u400[i,:,:] = Dataset('data/2018_10_'+dd[i]+'_12:00:00.nc').variables['U_GDS0_ISBL'][19,:,:]
	v400[i,:,:] = Dataset('data/2018_10_'+dd[i]+'_12:00:00.nc').variables['V_GDS0_ISBL'][19,:,:]
	o400[i,:,:] = Dataset('data/2018_10_'+dd[i]+'_12:00:00.nc').variables['W_GDS0_ISBL'][19,:,:]
        the350[i,:,:] = Dataset('data/2018_10_'+dd[i]+'_12:00:00.nc').variables['T_GDS0_ISBL'][18,:,:]
	the400[i,:,:] = Dataset('data/2018_10_'+dd[i]+'_12:00:00.nc').variables['T_GDS0_ISBL'][19,:,:]
	the450[i,:,:] = Dataset('data/2018_10_'+dd[i]+'_12:00:00.nc').variables['T_GDS0_ISBL'][20,:,:]

gph400 = gph400/9.806

hadv = np.zeros((3,721,1440))
vadv = np.zeros((3,721,1440))
hd = np.zeros((3,721,1440))
tilt = np.zeros((3,721,1440))
tend = np.zeros((3,721,1440))
wind = np.zeros((3,721,1440))

dtr = np.pi/180.
theta = lat*dtr
dthe = 0.25
Re = 6371000.
dx = Re*np.cos(theta)*dthe*dtr
dy = Re*dthe*dtr

the350 = ((the350)*(((1000./350.))**0.286))
the400 = ((the400)*(((1000./400.))**0.286))
the450 = ((the450)*(((1000./450.))**0.286))

for i in range(3):
	for j in range(1,len(lat)-1):
        	for k in range(1,len(lon)-1):
                        hthe350[i,j,k] = ( (((the350[i,j,k+1]-the350[i,j,k-1])/(2*dx[j]))**(2.)) + (((the350[i,j+1,k]-the350[i,j-1,k])/(-2*dy))**(2.)) )**(1./2.)
			hthe400[i,j,k] = ( (((the400[i,j,k+1]-the400[i,j,k-1])/(2*dx[j]))**(2.)) + (((the400[i,j+1,k]-the400[i,j-1,k])/(-2*dy))**(2.)) )**(1./2.)
                        hthe450[i,j,k] = ( (((the450[i,j,k+1]-the450[i,j,k-1])/(2*dx[j]))**(2.)) + (((the450[i,j+1,k]-the450[i,j-1,k])/(-2*dy))**(2.)) )**(1./2.)



for i in range(3):
        for j in range(1,len(lat)-1):
                for k in range(1,len(lon)-1):
			hadv[i,j,k] = -( ((u400[i,j,k]) * ((hthe400[i,j,k+1]-hthe400[i,j,k-1])/(2*dx[j]))) + ((v400[i,j,k]) * ((hthe400[i,j+1,k]-hthe400[i,j-1,k])/(-2*dy))))
			vadv[i,j,k] = -((o400[i,j,k]) * ((hthe350[i,j,k]-hthe450[i,j,k]) / (-10000)))
			hd[i,j,k] = (-1./hthe400[i,j,k]) * ( ((((the400[i,j,k+1]-the400[i,j,k-1])/(2*dx[j]))**(2.)) * ((u400[i,j,k+1]-u400[i,j,k-1])/(2*dx[j]))) + ( ((the400[i,j,k+1]-the400[i,j,k-1])/(2*dx[j])) * ((the400[i,j+1,k]-the400[i,j-1,k])/(-2*dy)) * (((v400[i,j,k+1]-v400[i,j,k-1])/(2*dx[j])) + ((u400[i,j+1,k]-u400[i,j-1,k])/(-2*dy))) ) + ((((the400[i,j+1,k]-the400[i,j-1,k])/(2*dy))**(2.)) * ((v400[i,j+1,k]-v400[i,j-1,k])/(-2*dy))))
			tilt[i,j,k] = (-1./hthe400[i,j,k]) * ((((the400[i,j,k+1]-the400[i,j,k-1])/(2*dx[j])) * ((the350[i,j,k]-the450[i,j,k]) / (-10000)) * ((o400[i,j,k+1]-o400[i,j,k-1])/(2*dx[j]))) + (((the400[i,j+1,k]-the400[i,j-1,k])/(-2*dy)) * ((the350[i,j,k]-the450[i,j,k])/(-10000)) * ((o400[i,j+1,k]-o400[i,j-1,k])/(-2*dy))))
			tend[i,j,k] = hadv[i,j,k] + vadv[i,j,k] + hd[i,j,k] + tilt[i,j,k]


lat1 = lat[::-1]

su, newlon1 = shiftgrid(180., u400, lon, start=False)
sv, newlon1 = shiftgrid(180., v400, lon, start=False)
su = su[:,::-1,:]
sv = sv[:,::-1,:]


hadv = (hadv)*(1000000000.)
vadv = (vadv)*(1000000000.)
hd = (hd)*(1000000000.)
tilt = (tilt)*(1000000000.)
tend = (tend)*(1000000000.)
hthe400 = (hthe400)*(100000.)

hthe400 = ndimage.gaussian_filter(hthe400, sigma=3, order=0)
hadv = ndimage.gaussian_filter(hadv, sigma=3, order=0)
vadv = ndimage.gaussian_filter(vadv, sigma=3, order=0)
hd = ndimage.gaussian_filter(hd, sigma=3, order=0)
tilt = ndimage.gaussian_filter(tilt, sigma=3, order=0)
tend = ndimage.gaussian_filter(tend, sigma=3, order=0)

titles= ['(a) Frontogenesis tendency','(b) Horizontal advection','(c) Vertical advection', '(d) Horizontal divergence and deformation', '(e) Tilting']

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
clevs1= np.arange(2,8,1)
colormesh = m.pcolormesh(x, y, tend[2,:,:], vmin = -4, vmax = 4, cmap=cmap)
cs=m.contour(x,y,hthe400[2,:,:], clevs1, colors='black', linewidths = 1.5, linestyles='--')
cbar=m.colorbar(colormesh)
uproj,vproj,xx,yy = m.transform_vector(su[2,:,:],sv[2,:,:],newlon1,lat1,10,10,returnxy=True,masked=True)
x2 , y2 = m(xx, yy)
barbs = m.barbs(x2, y2, uproj, vproj, length=6, pivot='middle', linewidth=0.5)
plt.clabel(cs, fmt='%3.0f', colors='black', fontsize = 7)
cbar.set_label('Frontogenesis tendency $(10^{-9}K/m^{1}s^{1})$', fontsize = 7)

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
clevs1= np.arange(2,8,1)
m.plot(x1, y1, color='r', marker='*')
cmap = mpl.cm.bwr
colormesh = m.pcolormesh(x, y, hadv[2,:,:], vmin = -4, vmax = 4, cmap=cmap)
cs=m.contour(x,y,hthe400[2,:,:], clevs1, colors='black', linewidths = 1.5, linestyles='--')
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
clevs1= np.arange(2,8,1)
m.plot(x1, y1, color='r', marker='*')
cmap = mpl.cm.bwr
colormesh = m.pcolormesh(x, y, vadv[2,:,:], vmin = -4, vmax = 4, cmap=cmap)
cs=m.contour(x,y,hthe400[2,:,:], clevs1, colors='black', linewidths = 1.5, linestyles='--')
cbar=m.colorbar(colormesh)
plt.clabel(cs, fmt='%3.0f', colors='black', fontsize = 7)
cbar.set_label('Vertical advection $(10^{-9}K/m^{1}s^{1})$', fontsize = 7)

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
clevs1= np.arange(2,8,1)
m.plot(x1, y1, color='r', marker='*')
cmap = mpl.cm.bwr
colormesh = m.pcolormesh(x, y, hd[2,:,:], vmin = -4, vmax = 4, cmap=cmap)
cs=m.contour(x,y,hthe400[2,:,:], clevs1, colors='black', linewidths = 1.5, linestyles='--')
cbar=m.colorbar(colormesh)
plt.clabel(cs, fmt='%3.0f', colors='black', fontsize = 7)
cbar.set_label('Divergence and deformation $(10^{-9}K/m^{1}s^{1})$', fontsize = 7)

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
clevs1= np.arange(2,8,1)
m.plot(x1, y1, color='r', marker='*')
cmap = mpl.cm.bwr
colormesh = m.pcolormesh(x, y, tilt[2,:,:], vmin = -4, vmax = 4, cmap=cmap)
cs=m.contour(x,y,hthe400[2,:,:], clevs1, colors='black', linewidths = 1.5, linestyles='--')
cbar=m.colorbar(colormesh)
plt.clabel(cs, fmt='%3.0f', colors='black', fontsize = 7)
cbar.set_label('Tilting $(10^{-9}K/m^{1}s^{1})$', fontsize = 7)

plt.subplots_adjust(wspace = 0.42, hspace = 0.30)

#plt.show()
fig = plt.savefig('./figures/frontogenesis tendency_30')
