from netCDF4 import Dataset
import numpy as np
from mpl_toolkits.basemap import Basemap, shiftgrid
import matplotlib.pyplot as plt
from sys import exit


data1=Dataset('data/2018-10-30_06:00:00.nc')
data2=Dataset('data/2018-10-30_12:00:00.nc')

lon1 = data1.variables['g0_lon_2'][:]
lat1 = data1.variables['g0_lat_1'][:]
gp500 = data1.variables['Z_GDS0_ISBL'][16,:,:]
gp1000 = data1.variables['Z_GDS0_ISBL'][31,:,:]
uu1 = data1.variables['U_GDS0_ISBL'][:,:,:]
vv1 = data1.variables['V_GDS0_ISBL'][:,:,:]

lon2 = data2.variables['g0_lon_2'][:]
lat2 = data2.variables['g0_lat_1'][:]
ggp500 = data2.variables['Z_GDS0_ISBL'][16,:,:]
ggp1000 = data2.variables['Z_GDS0_ISBL'][31,:,:]
uu2 = data2.variables['U_GDS0_ISBL'][:,:,:]
vv2 = data2.variables['V_GDS0_ISBL'][:,:,:]

gp500 = gp500/9.806
gp1000 = gp1000/9.806
ggp500 = ggp500/9.806
ggp1000 = ggp1000/9.806

gpt1 = gp500-gp1000
gpt2 = ggp500-ggp1000

wind1 = np.zeros((721,1440))
wind2 = np.zeros((721,1440))

for i in range(len(lat1)):
	for j in range(len(lon1)):
		wind1[i,j] = (((uu1[12,i,j])**(2.))+((vv1[12,i,j])**(2.)))**(1./2.)
		wind2[i,j] = (((uu2[12,i,j])**(2.))+((vv2[12,i,j])**(2.)))**(1./2.)

lat11 = lat1[::-1]
lat22 = lat2[::-1]

su1, newlon1 = shiftgrid(180., uu1, lon1, start=False)
sv1, newlon1 = shiftgrid(180., vv1, lon1, start=False)
su1 = su1[:,::-1,:]
sv1 = sv1[:,::-1,:]

su2, newlon2 = shiftgrid(180., uu2, lon2, start=False)
sv2, newlon2 = shiftgrid(180., vv2, lon2, start=False)
su2 = su2[:,::-1,:]
sv2 = sv2[:,::-1,:]


titles= ['(a) 20181030 06UTC geopotential height and thickness','(b) 20181030 12UTC geopotential height and thickness','(c) 20181030 06UTC wind speed','(d) 20181030 12UTC wind speed']
fig = plt.figure(figsize=(12,12))

plt.subplot(2, 2, 1)
m = Basemap(llcrnrlat=25,urcrnrlat=55,\
llcrnrlon=125,urcrnrlon=155, resolution='i')#,resolution='c')
plt.title(titles[0], fontsize=12, loc='left')
m.drawcoastlines()
m.drawcountries()
m.drawparallels(np.arange(25,55,10),labels=[1,0,0,0])
m.drawmeridians(np.arange(125,155,10),labels=[0,0,0,1])
m.fillcontinents(alpha=0)
m.drawmapboundary(fill_color='w')
lon_2d1, lat_2d1 = np.meshgrid(lon1, lat1)
x, y = m(lon_2d1, lat_2d1)
clevs1= np.arange(np.min(gp500[140:260,500:620]),np.max(gp500[140:260,500:620])+80,80)
clevs2= np.arange(np.min(gpt1[140:260,500:620]),np.max(gpt1[140:260,500:620])+80,80)
cs1=m.contour(x,y,gp500,clevs1,colors='red', linestyles='solid')
cs2=m.contour(x,y,gpt1,clevs2,colors='blue', linestyles='dashed')
plt.clabel(cs1, colors='red')
plt.clabel(cs2, colors='blue')

plt.subplot(2, 2, 2)
m = Basemap(llcrnrlat=25,urcrnrlat=55,\
llcrnrlon=125,urcrnrlon=155, resolution='i')#,resolution='c')
plt.title(titles[1], fontsize=12, loc='left')
m.drawcoastlines()
m.drawcountries()
m.drawparallels(np.arange(25,55,10),labels=[1,0,0,0])
m.drawmeridians(np.arange(125,155,10),labels=[0,0,0,1])
m.fillcontinents(alpha=0)
m.drawmapboundary(fill_color='w')
lon_2d1, lat_2d1 = np.meshgrid(lon1, lat1)
x, y = m(lon_2d1, lat_2d1)
clevs1= np.arange(np.min(ggp500[140:260,500:620]),np.max(ggp500[140:260,500:620])+80,80)
clevs2= np.arange(np.min(gpt2[140:260,500:620]),np.max(gpt2[140:260,500:620])+80,80)
cs1=m.contour(x,y,ggp500,clevs1,colors='red', linestyles='solid')
cs2=m.contour(x,y,gpt2,clevs2,colors='blue', linestyles='dashed')
plt.clabel(cs1, colors='red', inline=1)
plt.clabel(cs2, colors='blue', inline=1)

plt.subplot(2, 2, 3)
m = Basemap(llcrnrlat=25,urcrnrlat=55,\
llcrnrlon=125,urcrnrlon=155, resolution='i')#,resolution='c')
plt.title(titles[2], fontsize=12, loc='left')
m.drawcoastlines()
m.drawcountries()
m.drawparallels(np.arange(25,55,10),labels=[1,0,0,0])
m.drawmeridians(np.arange(125,155,10),labels=[0,0,0,1])
m.fillcontinents(alpha=0)
m.drawmapboundary(fill_color='w')
lon_2d1, lat_2d1 = np.meshgrid(lon1, lat1)
x, y = m(lon_2d1, lat_2d1)
clevs1= np.arange(np.min(wind1[140:260,500:620]),np.max(wind1[140:260,500:620])+10,10)
cs1=m.contourf(x,y,wind1,clevs1)
uproj,vproj,xx,yy = m.transform_vector(su1[12,:,:],sv1[12,:,:],newlon1,lat11,12,12,returnxy=True,masked=True)
x1 , y1 = m(xx, yy)
barbs = m.barbs(x1, y1, uproj, vproj, length=6, pivot='middle', linewidth=0.5)
cbar=m.colorbar(cs1)#, cax=cbaxes)
cbar.set_label('wind speed')


plt.subplot(2, 2, 4)
m = Basemap(llcrnrlat=25,urcrnrlat=55,\
llcrnrlon=125,urcrnrlon=155, resolution='i')#,resolution='c')
plt.title(titles[3], fontsize=12, loc='left')
m.drawcoastlines()
m.drawcountries()
m.drawparallels(np.arange(25,55,10),labels=[1,0,0,0])
m.drawmeridians(np.arange(125,155,10),labels=[0,0,0,1])
m.fillcontinents(alpha=0)
m.drawmapboundary(fill_color='w')
lon_2d1, lat_2d1 = np.meshgrid(lon1, lat1)
x, y = m(lon_2d1, lat_2d1)
clevs1= np.arange(np.min(wind2[140:260,500:620]),np.max(wind2[140:260,500:620])+10,10)
cs1=m.contourf(x,y,wind2,clevs1)
uproj,vproj,xx,yy = m.transform_vector(su2[12,:,:],sv2[12,:,:],newlon2,lat22,12,12,returnxy=True,masked=True)
x1 , y1 = m(xx, yy)
barbs = m.barbs(x1, y1, uproj, vproj, length=6, pivot='middle', linewidth=0.5)
cbar=m.colorbar(cs1)#, cax=cbaxes)
cbar.set_label('wind speed')

plt.show()


