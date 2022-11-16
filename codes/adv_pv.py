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
lv = Dataset('data/2018-10-30_06:00:00.nc').variables['lv_ISBL0'][:]
gph300 = np.zeros((3,721,1440))
gph500 = np.zeros((3,721,1440))
gph600 = np.zeros((3,721,1440))
gph850 = np.zeros((3,721,1440))
T300 = np.zeros((3,721,1440))
T500 = np.zeros((3,721,1440))
T600 = np.zeros((3,721,1440))
T850 = np.zeros((3,721,1440))
u300 = np.zeros((3,721,1440))
u500 = np.zeros((3,721,1440))
u600 = np.zeros((3,721,1440))
u850 = np.zeros((3,721,1440))
v300 = np.zeros((3,721,1440))
v500 = np.zeros((3,721,1440))
v600 = np.zeros((3,721,1440))
v850 = np.zeros((3,721,1440))

gph = np.zeros((12,721,1440))
T = np.zeros((12,721,1440))
u = np.zeros((12,721,1440))
v = np.zeros((12,721,1440))
Tall = np.zeros((3,32,721,1440))
PT = np.zeros((3,32,721,1440))

dd = ['28','29','30']

for i  in range(3):
	gph300[i,:,:] = Dataset('data/2018-10-'+dd[i]+'_12:00:00.nc').variables['Z_GDS0_ISBL'][12,:,:]
	gph500[i,:,:] = Dataset('data/2018-10-'+dd[i]+'_12:00:00.nc').variables['Z_GDS0_ISBL'][16,:,:]
	gph600[i,:,:] = Dataset('data/2018-10-'+dd[i]+'_12:00:00.nc').variables['Z_GDS0_ISBL'][18,:,:]
	gph850[i,:,:] = Dataset('data/2018-10-'+dd[i]+'_12:00:00.nc').variables['Z_GDS0_ISBL'][25,:,:]
	u300[i,:,:] = Dataset('data/2018-10-'+dd[i]+'_12:00:00.nc').variables['U_GDS0_ISBL'][12,:,:]
	u500[i,:,:] = Dataset('data/2018-10-'+dd[i]+'_12:00:00.nc').variables['U_GDS0_ISBL'][16,:,:]
	u600[i,:,:] = Dataset('data/2018-10-'+dd[i]+'_12:00:00.nc').variables['U_GDS0_ISBL'][18,:,:]
	u850[i,:,:] = Dataset('data/2018-10-'+dd[i]+'_12:00:00.nc').variables['U_GDS0_ISBL'][25,:,:]
	v300[i,:,:] = Dataset('data/2018-10-'+dd[i]+'_12:00:00.nc').variables['V_GDS0_ISBL'][12,:,:]
	v500[i,:,:] = Dataset('data/2018-10-'+dd[i]+'_12:00:00.nc').variables['V_GDS0_ISBL'][16,:,:]
	v600[i,:,:] = Dataset('data/2018-10-'+dd[i]+'_12:00:00.nc').variables['V_GDS0_ISBL'][18,:,:]
	v850[i,:,:] = Dataset('data/2018-10-'+dd[i]+'_12:00:00.nc').variables['V_GDS0_ISBL'][25,:,:]
	T300[i,:,:] = Dataset('data/2018-10-'+dd[i]+'_12:00:00.nc').variables['T_GDS0_ISBL'][12,:,:]
	T500[i,:,:] = Dataset('data/2018-10-'+dd[i]+'_12:00:00.nc').variables['T_GDS0_ISBL'][16,:,:]
	T600[i,:,:] = Dataset('data/2018-10-'+dd[i]+'_12:00:00.nc').variables['T_GDS0_ISBL'][18,:,:]
	T850[i,:,:] = Dataset('data/2018-10-'+dd[i]+'_12:00:00.nc').variables['T_GDS0_ISBL'][25,:,:]
	Tall[i,:,:,:] = Dataset('data/2018-10-'+dd[i]+'_12:00:00.nc').variables['T_GDS0_ISBL'][:,:,:]

for i in range(12):
	if i in range(0,3):
		gph[i,:,:] = gph300[i,:,:]
		u[i,:,:] = u300[i,:,:]
		v[i,:,:] = v300[i,:,:]
		T[i,:,:] = T300[i,:,:]
	elif i in range(3,6):
		gph[i,:,:] = gph500[i-3,:,:]
		u[i,:,:] = u500[i-3,:,:]
		v[i,:,:] = v500[i-3,:,:]
		T[i,:,:] = T500[i-3,:,:]
	elif i in range(6,9):
		gph[i,:,:] = gph600[i-6,:,:]
		u[i,:,:] = u600[i-6,:,:]
		v[i,:,:] = v600[i-6,:,:]
		T[i,:,:] = T600[i-6,:,:]
	elif i in range(9,12):
		gph[i,:,:] = gph850[i-12,:,:]
		u[i,:,:] = u850[i-12,:,:]
		v[i,:,:] = v850[i-12,:,:]
		T[i,:,:] = T850[i-12,:,:]

gph = gph/9.806

dx = np.zeros((721,1440))
dy = np.zeros((721,1440))
dvdx = np.zeros((12,721,1440))
dudy = np.zeros((12,721,1440))

for i in range(1,len(lat)-1):
        for j in range(1,len(lon)-1):
                dx[i,j] = 2*np.pi*6371000*np.cos(lat[i]*2*np.pi/360.)*(lon[j+1]-lon[j-1])/360.
                dy[i,j] = 2*np.pi*6371000*(lat[i+1]-lat[i-1])/360.

dy *= -1

abv = np.zeros((12,721,1440))
dvdx = np.zeros((12,721,1440))
dudy = np.zeros((12,721,1440))

for i in range(3):
	for l in range(len(lv)):
        	for j in range(len(lat)):
                	for k in range(len(lon)):
				PT[i,l,j,k] = ((Tall[i,l,j,k])*((1000./float(int(lv[l])))**0.286))

PV = np.zeros((12,721,1440))
PV_adv = np.zeros((12,721,1440))

for i in range(12):
        for j in range(1,len(lat)-1):
                for k in range(1,len(lon)-1):
			if i in range(0,3):
				dvdx[i,j,k] = (v[i,j,k+1]-v[i,j,k-1])/(dx[j,k])
	                        dudy[i,j,k] = (u[i,j+1,k]-u[i,j-1,k])/(dy[j,k])
				PV[i,j,k] = (-9.806)*((PT[i,11,j,k]-PT[i,13,j,k])/(100.*(float(lv[11])-float(lv[13]))))*((dvdx[i,j,k]-dudy[i,j,k])+((2.)*(2.*np.pi/86400.)*(np.sin(lat[j]*np.pi/180.))))
			elif i in range(3,6):
				PV[i,j,k] = (-9.806)*((PT[i-3,15,j,k]-PT[i-3,17,j,k])/(100.*(float(lv[15])-float(lv[17]))))*((dvdx[i,j,k]-dudy[i,j,k])+((2.)*(2.*np.pi/86400.)*(np.sin(lat[j]*np.pi/180.))))
			elif i in range(6,9):
				PV[i,j,k] = (-9.806)*((PT[i-6,17,j,k]-PT[i-6,19,j,k])/(100.*(float(lv[17])-float(lv[19]))))*((dvdx[i,j,k]-dudy[i,j,k])+((2.)*(2.*np.pi/86400.)*(np.sin(lat[j]*np.pi/180.))))
			elif i in range(9,12):
				PV[i,j,k] = (-9.806)*((PT[i-9,24,j,k]-PT[i-9,26,j,k])/(100.*(float(lv[24])-float(lv[26]))))*((dvdx[i,j,k]-dudy[i,j,k])+((2.)*(2.*np.pi/86400.)*(np.sin(lat[j]*np.pi/180.))))


dx, dy = mpcalc.lat_lon_grid_deltas(lon, lat)
dy *= -1

for i in range(12):
        PV_adv[i,:,:] = mpcalc.advection(PV[i,:,:], [u[i,:,:], v[i,:,:]], (dx, dy), dim_order='yx')

PV_adv = (PV_adv)*(3600.)*(10000000.)
PV_adv = ndimage.gaussian_filter(PV_adv, sigma=3, order=0)


titles= ['(a) 20181028 12UTC','(b) 20181029 12UTC', '(c) 20181030 12UTC', '(d) 20181028 12UTC', '(e) 20181029 12UTC', '(f) 20181030 12UTC', '(g) 20181028 12UTC', '(h) 20181029 12UTC', '(i) 20181030 12UTC', '(j) 20181028 12UTC', '(k) 20181029 12UTC', '(l) 20181030 12UTC']

fig = plt.figure(figsize=(15,10))
for i in range(12):
        plt.subplot(4, 3, i+1)
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
	clevs1= np.arange(np.min(gph[i,120:281,400:641]),np.max(gph[i,120:281,400:641])+60,60)
	m.plot(x1, y1, color='r', marker='*')
	cmap = mpl.cm.bwr
	colormesh = m.pcolormesh(x, y, PV_adv[i,:,:], vmin = -10, vmax = 10, cmap=cmap)
	cs=m.contour(x,y,gph[i,:,:], clevs1, colors='black', linewidths = 0.3, linestyles='solid')
	cbar=m.colorbar(colormesh)
	plt.clabel(cs, fmt='%3.0f', colors='black', fontsize = 7)
        cbar.set_label('PV advection $(10^{-7}K*m^{2}/kg*s/hr)$', fontsize = 7)

plt.subplots_adjust(wspace = 0.42, hspace = 0.30)
fig = plt.savefig('./figures/advection_pv')
