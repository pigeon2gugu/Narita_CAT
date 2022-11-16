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

lon = Dataset('data/2018_10_30_12:00:00.nc').variables['g0_lon_2'][:]
lat = Dataset('data/2018_10_30_12:00:00.nc').variables['g0_lat_1'][:]
lv = Dataset('data/2018_10_30_12:00:00.nc').variables['lv_ISBL0'][:]
gph = np.zeros((3,37,721,1440))
u = np.zeros((3,37,721,1440))
v = np.zeros((3,37,721,1440))
T = np.zeros((3,37,721,1440))
wind = np.zeros((3,37,721,1440))
PT = np.zeros((3,37,721,1440))

k = 0
dd = ['28','29','30']

for i  in range(3):
	gph[i,:,:,:] = Dataset('data/2018_10_'+dd[i]+'_12:00:00.nc').variables['Z_GDS0_ISBL'][:,:,:]
	u[i,:,:,:] = Dataset('data/2018_10_'+dd[i]+'_12:00:00.nc').variables['U_GDS0_ISBL'][:,:,:]
	v[i,:,:,:] = Dataset('data/2018_10_'+dd[i]+'_12:00:00.nc').variables['V_GDS0_ISBL'][:,:,:]
	T[i,:,:,:] = Dataset('data/2018_10_'+dd[i]+'_12:00:00.nc').variables['T_GDS0_ISBL'][:,:,:]


gph = gph/9.806

dtr = np.pi/180.
theta = lat*dtr
dthe = 0.25
Re = 6371000.
dx = Re*np.cos(theta)*dthe*dtr
dy = Re*dthe*dtr

dvdx = np.zeros((3,37,721,1440))
dudy = np.zeros((3,37,721,1440))

for i in range(3):
        for l in range(37):
                for j in range(118,283):
                        for k in range(398,643):
                                #wind[i,l,j,k] = (((u[i,l,j,k])**(2.))+((v[i,l,j,k])**(2.)))**(1./2.)
                                PT[i,l,j,k] = ((T[i,l,j,k])*((1000./float(int(lv[l])))**0.286))
                		dvdx[i,l,j,k] = (v[i,l,j,k+1]-v[i,l,j,k-1])/(2*dx[j])
                		dudy[i,l,j,k] = (u[i,l,j+1,k]-u[i,l,j-1,k])/(-2*dy)

PVU = np.zeros((3,37,721,1440))
PVUT = np.zeros((3,721,1440))
PVUP = np.zeros((3,721,1440))
FPVU = np.zeros((3,721,1440))

for i in range(3):
	for l in range(1,36):
		for j in range(118,283):
			for k in range(398,643):
				PVU[i,l,j,k] = (10**6)*(-9.806)*((PT[i,l-1,j,k]-PT[i,l+1,j,k])/(100.*(float(lv[l-1])-float(lv[l+1]))))*((dvdx[i,l,j,k]-dudy[i,l,j,k])+((2.)*(2.*np.pi/86400.)*(np.sin(lat[j]*dtr))))
				if 1.5<=PVU[i,l,j,k]:
					PVUT[i,j,k] = PT[i,l,j,k]
					PVUP[i,j,k] = lv[l]
					FPVU[i,j,k] = PVU[i,l,j,k]

PVUP = ndimage.gaussian_filter(PVUP, sigma=3, order=0)


titles= ['(a) 20181028 12UTC','(b) 20181029 12UTC', '(c) 20181030 12UTC']
fig = plt.figure(figsize=(15,10))

for i in range(3):
        plt.subplot(1, 3, i+1)
	plt.title(titles[i], fontsize=12, loc='left')
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
	cmap = matplotlib.colors.ListedColormap(['navy', 'blue', 'royalblue', 'deepskyblue','darkgreen','lime','greenyellow','yellow','gold','orange'])
        clevs1 = np.arange(150,700,50)
	colormesh = m.pcolormesh(x, y, PVUP[i,:,:], vmin = 150, vmax = 650, cmap=cmap)
	clevs1 = np.arange(150,750,100)
	cs2=m.contour(x,y,PVUP[i,:,:], clevs1, colors='white', linewidths = 1, linestyles='solid')
	cbar=m.colorbar(colormesh)
	plt.clabel(cs2,fmt='%3.0f',colors='white', fontsize = 7)
	cbar.set_label('Pressure $(hPa)$', fontsize = 7)

plt.subplots_adjust(wspace = 0.42, hspace = 0.30)
#plt.show()
fig = plt.savefig('./figures/PVU_pressure')
