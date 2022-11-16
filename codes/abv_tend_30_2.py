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

olon = Dataset('data/2018_10_30_12:00:00.nc').variables['g0_lon_2'][:]
olat = Dataset('data/2018_10_30_12:00:00.nc').variables['g0_lat_1'][:]

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



no500 = np.zeros((3,41,61))
nu450 = np.zeros((3,41,61))
nu500 = np.zeros((3,41,61))
nu550 = np.zeros((3,41,61))
nv450 = np.zeros((3,41,61))
nv500 = np.zeros((3,41,61))
nv550 = np.zeros((3,41,61))
nvo450 = np.zeros((3,41,61))
nvo500 = np.zeros((3,41,61))
nvo550 = np.zeros((3,41,61))


for i in range(3):
        for j in range(41):
                for k in range(61):
                        for l in range(0,4):
                                for m in range(0,4):
                                        nu450[i,j,k] = nu450[i,j,k] + u450[i,(4*j)+120+l,(4*k)+400+m]
                                        nu500[i,j,k] = nu500[i,j,k] + u500[i,(4*j)+120+l,(4*k)+400+m]
                                        nu550[i,j,k] = nu550[i,j,k] + u550[i,(4*j)+120+l,(4*k)+400+m]
                                        nv450[i,j,k] = nv450[i,j,k] + v450[i,(4*j)+120+l,(4*k)+400+m]
                                        nv500[i,j,k] = nv500[i,j,k] + v500[i,(4*j)+120+l,(4*k)+400+m]
                                        nv550[i,j,k] = nv550[i,j,k] + v550[i,(4*j)+120+l,(4*k)+400+m]
                                        no500[i,j,k] = no500[i,j,k] + o500[i,(4*j)+120+l,(4*k)+400+m]
                                        nvo450[i,j,k] = nvo450[i,j,k] + vo450[i,(4*j)+120+l,(4*k)+400+m]
                                        nvo500[i,j,k] = nvo500[i,j,k] + vo500[i,(4*j)+120+l,(4*k)+400+m]
                                        nvo550[i,j,k] = nvo550[i,j,k] + vo550[i,(4*j)+120+l,(4*k)+400+m]


nu450 = nu450/16.
nu500 = nu500/16.
nu550 = nu550/16.
nv450 = nv450/16.
nv500 = nv500/16.
nv550 = nv550/16.
no500 = no500/16.
nvo450 = nvo450/16.
nvo500 = nvo500/16.
nvo550 = nvo550/16.

lat = np.arange(60,19,-1)
lon = np.arange(100,161,1)


hadv = np.zeros((3,41,61))
vadv = np.zeros((3,41,61))
stre = np.zeros((3,41,61))
tilt = np.zeros((3,41,61))
tend = np.zeros((3,41,61))



dtr = np.pi/180.
theta = lat*dtr
dthe = 1.0
Re = 6371000.
dx = Re*np.cos(theta)*dthe*dtr
dy = Re*dthe*dtr


f = 2.*(2.*np.pi/86400.)*(np.sin(lat*dtr))
f = np.full((61,41), f)
f= np.transpose(f, (1,0))
f = np.full((3,41,61), f)


nvo450 = nvo450 + f
nvo500 = nvo500 + f
nvo550 = nvo550 + f

for i in range(3):
        for j in range(1,len(lat)-1):
                for k in range(1,len(lon)-1):
                        hadv[i,j,k] = -( ((nu500[i,j,k]) * ((nvo500[i,j,k+1]-nvo500[i,j,k-1])/(2*dx[j]))) + ((nv500[i,j,k]) * ((nvo500[i,j+1,k]-nvo500[i,j-1,k])/(-2*dy))) )
                        vadv[i,j,k] = -((no500[i,j,k]) * ((nvo450[i,j,k]-nvo550[i,j,k]) / (-10000)))
                        stre[i,j,k] = -((nvo500[i,j,k]) * (((nu500[i,j,k+1]-nu500[i,j,k-1]) / (2*dx[j])) + ((nv500[i,j+1,k]-nv500[i,j-1,k]) / (-2*dy))))
                        tilt[i,j,k] = - ((((no500[i,j,k+1]-no500[i,j,k-1]) / (2*dx[j])) * ((nv450[i,j,k]-nv550[i,j,k]) / (-10000))) - (((no500[i,j+1,k]-no500[i,j-1,k]) / (-2*dy)) * ((nu450[i,j,k]-nu550[i,j,k]) / (-10000))))



tend = hadv + vadv + stre + tilt


lat1 = olat[::-1]

su, newlon1 = shiftgrid(180., u500, olon, start=False)
sv, newlon1 = shiftgrid(180., v500, olon, start=False)
su = su[:,::-1,:]
sv = sv[:,::-1,:]

hadv = (hadv)*(1000000000.)
vadv = (vadv)*(1000000000.)
stre = (stre)*(1000000000.)
tilt = (tilt)*(1000000000.)
tend = (tend)*(1000000000.)
nvo500 = (nvo500)*(100000.)

nvo500 = ndimage.gaussian_filter(nvo500, sigma=1, order=0)
hadv = ndimage.gaussian_filter(hadv, sigma=1, order=0)
vadv = ndimage.gaussian_filter(vadv, sigma=1, order=0)
stre = ndimage.gaussian_filter(stre, sigma=1, order=0)
tilt = ndimage.gaussian_filter(tilt, sigma=1, order=0)
tend = ndimage.gaussian_filter(tend, sigma=1, order=0)


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
x1, y1 = m(lon[41], lat[24]) #CAT incounter
m.plot(x1, y1, color='r', marker='*')
cmap = mpl.cm.bwr
clevs1= np.arange(16,40,4)
colormesh = m.pcolormesh(x, y, tend[2,:,:], vmin = -10, vmax = 10, cmap=cmap)
cs=m.contour(x,y,nvo500[2,:,:], clevs1, colors='black', linewidths = 1.5, linestyles='--')
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
x1, y1 = m(lon[41], lat[24]) #CAT incounter
clevs1= np.arange(16,40,4)
m.plot(x1, y1, color='r', marker='*')
cmap = mpl.cm.bwr
colormesh = m.pcolormesh(x, y, hadv[2,:,:], vmin = -14, vmax = 14, cmap=cmap)
cs=m.contour(x,y,nvo500[2,:,:], clevs1, colors='black', linewidths = 1.5, linestyles='--')
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
x1, y1 = m(lon[41], lat[24]) #CAT incounter
clevs1= np.arange(16,40,4)
m.plot(x1, y1, color='r', marker='*')
cmap = mpl.cm.bwr
colormesh = m.pcolormesh(x, y, vadv[2,:,:], vmin = -6, vmax = 6, cmap=cmap)
cs=m.contour(x,y,nvo500[2,:,:], clevs1, colors='black', linewidths = 1.5, linestyles='--')
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
x1, y1 = m(lon[41], lat[24]) #CAT incounter
clevs1= np.arange(16,40,4)
m.plot(x1, y1, color='r', marker='*')
cmap = mpl.cm.bwr
colormesh = m.pcolormesh(x, y, stre[2,:,:], vmin = -6, vmax = 6, cmap=cmap)
cs=m.contour(x,y,nvo500[2,:,:], clevs1, colors='black', linewidths = 1.5, linestyles='--')
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
x1, y1 = m(lon[41], lat[24]) #CAT incounter
clevs1= np.arange(16,40,4)
m.plot(x1, y1, color='r', marker='*')
cmap = mpl.cm.bwr
colormesh = m.pcolormesh(x, y, tilt[2,:,:], vmin = -6, vmax = 6, cmap=cmap)
cs=m.contour(x,y,nvo500[2,:,:], clevs1, colors='black', linewidths = 1.5, linestyles='--')
cbar=m.colorbar(colormesh)
plt.clabel(cs, fmt='%3.0f', colors='black', fontsize = 7)
cbar.set_label('Vortex tilting $(10^{-9}/s^{2})$', fontsize = 7)

plt.subplots_adjust(wspace = 0.42, hspace = 0.30)
fig = plt.savefig('./figures/vorticity tendency_30_2')
