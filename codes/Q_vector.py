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

for i in range(12):
	if i in range(0,3):
		gph[i,:,:] = gph300[i,:,:,]
		u[i,:,:] = u300[i,:,:,]
		v[i,:,:] = v300[i,:,:,]
		T[i,:,:] = T300[i,:,:,]
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
		gph[i,:,:] = gph850[i-9,:,:]
		u[i,:,:] = u850[i-9,:,:]
		v[i,:,:] = v850[i-9,:,:]
		T[i,:,:] = T850[i-9,:,:]

gph = gph/9.806

dx, dy = mpcalc.lat_lon_grid_deltas(lon, lat)

uqvect = np.zeros((12,721,1440))
vqvect = np.zeros((12,721,1440))
q_div = np.zeros((12,721,1440))

lv = [300, 300, 300, 500, 500, 500, 600, 600, 600, 850, 850, 850]
lv = np.array(lv)

for i in range(12):
	uqvect[i,:,:], vqvect[i,:,:] = mpcalc.q_vector(u[i,:,:], v[i,:,:], T[i,:,:], lv[i]*units.hPa, dx, dy)
	q_div[i,:,:] = mpcalc.divergence(uqvect[i,:,:], vqvect[i,:,:], dx, dy, dim_order='yx')

q_div = q_div*(10.**(18.))
q_div = ndimage.gaussian_filter(q_div, sigma=2.5, order=0)

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
	clevs_qdiv = list(range(-30, -4, 5))+list(range(5, 35, 5))
	q_div[ q_div > max(clevs_qdiv) ] = max(clevs_qdiv)
	cs1 = m.contourf(x, y, q_div[i,:,:], clevs_qdiv, cmap=cmap)
	cs=m.contour(x,y,gph[i,:,:], clevs1, colors='black', linewidths = 0.3, linestyles='solid')
	cbar=m.colorbar(cs1)
	#wind_slice = (slice(None, None, 5), slice(None, None, 5))
	#m.quiver(lon[wind_slice[0]], lat[wind_slice[1]],
         # uqvect[i][wind_slice],
          #vqvect[i][wind_slice],
          #pivot='mid', color='black', scale=10**(-5), scale_units='inches')
	plt.clabel(cs, fmt='%3.0f', colors='black', fontsize = 7)
        cbar.set_label('Q vector divergence $(10^{18} m s^{-1} kg^{-1})$', fontsize = 7)

plt.subplots_adjust(wspace = 0.42, hspace = 0.30)
fig = plt.savefig('./figures/Q_vector')
