from netCDF4 import Dataset
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from sys import exit

data1=Dataset('data/2018-10-30_06:00:00_sfc.nc')

lon = data1.variables['g0_lon_1'][:]
lat = data1.variables['g0_lat_0'][:]
SLP = np.zeros((11,721,1440))

k = 0
dd = ['28','29','30']
tt = ['00','06','12','18']

for i in range(3):
	for j in range(4):
			if dd[i] == '30' and tt[j] == '18':
				break
			else :
				SLP[k,:,:] = Dataset('data/2018-10-'+dd[i]+'_'+tt[j]+':00:00_sfc.nc').variables['MSL_GDS0_SFC'][:][:]
                		k = k+1


print 'end reading variables'

SLP = SLP/100.

for i in range(11):
	for j in range(1,len(lat)-1):
		for k in range(1,len(lon)-1):
			SLP[i,j,k] = (((SLP[i,j-1,k] + ((2.)*SLP[i,j,k]) + SLP[i,j+1,k])/(4.))+((SLP[i,j,k-1] + ((2.)*SLP[i,j,k]) + SLP[i,j,k+1])/(4.)))/(2.)

for i in range(11):
	for j in range(len(lon)):
		SLP[i,0,j] = SLP[i,1,j]
		SLP[i,-1,j] = SLP[i,-2,j]

for i in range(11):
	for j in range(len(lat)):
		SLP[i,j,0] = SLP[i,j,1]
		SLP[i,j,-1] = SLP[i,j,-2]

SLP = np.around(SLP)
SLP = SLP.astype(int)

titles= ['(a) 20181028 00UTC','(b) 20181028 06UTC', '(c) 20181028 12UTC', '(d) 20181028 18UTC', '(e) 20181029 00UTC', '(f) 20181029 06UTC', '(g) 20181029 12UTC', '(h) 20181029 18UTC', '(i) 20181030 00UTC', '(j) 20181030 06UTC', '(k) 20181030 12UTC']
fig = plt.figure(figsize=(15,10))

for i in range(11):
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
	x1, y1 = m(lon[562], lat[215]) #CAT incounter
	x, y = m(lon_2d1, lat_2d1)
	clevs= np.arange(956,1032,4)
	m.plot(x1, y1, color='r', marker='*')
	cs=m.contour(x,y,SLP[i,:,:], clevs.astype(int), colors='red', linewidths = 0.5, linestyles='solid')
	plt.clabel(cs,fmt='%3.0f',colors='black')

plt.subplots_adjust(wspace = 0.42, hspace = 0.30)

fig = plt.savefig('./figures/SLP_20181028_29_30')


