from netCDF4 import Dataset
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib.colorbar
from matplotlib.colorbar import make_axes
from sys import exit

file = './Narita_lon.txt'
f= open(file)
lines=f.readlines()

lon = []

for line in lines:
        lon.append(float(line.strip()))

f.close()

file = './Narita_lat.txt'
g = open(file)

lines=g.readlines()

lat = []

for line in lines:
        lat.append(float(line.strip()))

lon = np.array(lon)
lat = np.array(lat)


data1=Dataset('../../../data/ERA5/2018/ERA5.100_500hPa.201810.06.nc')
data2=Dataset('../../../data/ERA5/2018/ERA5.100_500hPa.201810.12.nc')

lon1 = data1.variables['g0_lon_3'][:]
lat1 = data1.variables['g0_lat_2'][:]
ggp1 = data1.variables['Z_GDS0_ISBL'][29,:,:,:]
uu1 = data1.variables['U_GDS0_ISBL'][29,:,:,:]
vv1 = data1.variables['V_GDS0_ISBL'][29,:,:,:]
ww1 = data1.variables['W_GDS0_ISBL'][29,:,:,:]

lon2 = data2.variables['g0_lon_3'][:]
lat2 = data2.variables['g0_lat_2'][:]
ggp2 = data2.variables['Z_GDS0_ISBL'][29,:,:,:]
uu2 = data2.variables['U_GDS0_ISBL'][29,:,:,:]
vv2 = data2.variables['V_GDS0_ISBL'][29,:,:,:]
ww2 = data2.variables['W_GDS0_ISBL'][29,:,:,:]

WS1 = np.zeros((721,1440))
WS2 = np.zeros((721,1440))

for i in range(len(lat)):
	for j in range(len(lon)):
		WS1[i,j] = (((uu1[10,i,j])**(2.))+(vv1[10,i,j])**(2.))**(1./2.)
		WS2[i,j] = (((uu2[10,i,j])**(2.))+(vv2[10,i,j])**(2.))**(1./2.)

exit()

titles= ['(a) 20181030_06UTC','(b) 20181030_12UTC']
fig = plt.figure(figsize=(10,10))

plt.subplot(2, 1, 1)
m = Basemap(projection='merc',llcrnrlat=33,urcrnrlat=37,\
                        llcrnrlon=138,urcrnrlon=142,resolution='c')
plt.title(titles[0], fontsize=12, loc='left')
m.drawcoastlines()
m.drawcountries()
m.drawparallels(np.arange(33,37,1),labels=[1,1,0,0])
m.drawmeridians(np.arange(138,142,1),labels=[0,0,0,1])
m.fillcontinents(alpha=0)
m.drawmapboundary(fill_color='w')
lon_2d, lat_2d = np.meshgrid(lon1, lat1)
x1, y1 = m(lon_2d, lat_2d)
x, y= m(lon, lat)
clevs1= np.arange(np.min(ggp1[10]),np.max(ggp1[10])+1,200)
clevs2= np.arange(44,64,0.1)
WS1[ WS1 < min(clevs2) ] = min(clevs2)
cs1=m.contour(x1,y1,ggp1[10],clevs1,colors='black', linestyles='dashed')
cs2=m.contourf(x1,y1,WS1,clevs2)
plt.clabel(cs1,colors='black')
cbar2 = plt.colorbar(cs2, shrink=0.9, pad = 0.05)
cbar2.set_label('m/s')
for i in range(len(lon)):
        m.plot(x[i], y[i], color='r',marker=',')

plt.subplot(2, 1, 2)
m = Basemap(projection='merc',llcrnrlat=33,urcrnrlat=37,\
                        llcrnrlon=138,urcrnrlon=142,resolution='c')
plt.title(titles[1], fontsize=12, loc='left')
m.drawcoastlines()
m.drawcountries()
m.drawparallels(np.arange(33,37,1),labels=[1,1,0,0])
m.drawmeridians(np.arange(138,142,1),labels=[0,0,0,1])
m.fillcontinents(alpha=0)
m.drawmapboundary(fill_color='w')
lon_2d, lat_2d = np.meshgrid(lon1, lat1)
x1, y1 = m(lon_2d, lat_2d)
x, y= m(lon, lat)
clevs1= np.arange(np.min(ggp2[10]),np.max(ggp2[10])+1,200)
clevs2= np.arange(44,64,0.1)
WS2[ WS2 < min(clevs2) ] = min(clevs2)
cs3=m.contour(x1,y1,ggp2[10],clevs1,colors='black', linestyles='dashed')
cs4=m.contourf(x1,y1,WS2,clevs2)
plt.clabel(cs3, colors='black')
cbar3 = plt.colorbar(cs4, shrink=0.9, pad = 0.05)
cbar3.set_label('m/s')
for i in range(len(lon)):
        m.plot(x[i], y[i], color='r',marker=',')

plt.show()
