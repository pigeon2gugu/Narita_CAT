from netCDF4 import Dataset
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from sys import exit


data1=Dataset('ERA5.100_500hPa.201810.06.nc')
data2=Dataset('ERA5.100_500hPa.201810.12.nc')

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


dx = np.zeros((721,1440))
dy = np.zeros((721,1440))

for i in range(len(lat1)):
	for j in range(1,len(lon1)-1):
		dx[i,j] = 2*np.pi*6371000*np.cos(lat1[i]*2*np.pi/360.)*(lon1[j+1]-lon1[j-1])/360.
		
for i in range(len(lat1)):
	dx[i,-1] = 2*np.pi*6371000*np.cos(lat1[i]*np.pi/360.)*(360-lon1[j-1])/360.
	dx[i,0] = 2*np.pi*6371000*np.cos(lat1[i]*np.pi/360.)*(lon1[1]-lon1[j]+360)/360.

for i in range(1,len(lat1)-1):
	for j in range(len(lon1)):
		dy[i,j] = 2*np.pi*6371000*(lat1[i+1]-lat1[i-1])/360.


for j in range(len(lon1)):
	dy[-1,j] = 2*np.pi*6371000*(0.5)/360.
	dy[0,j] = 2*np.pi*6371000*(0.5)/360.


dudx = np.zeros((721,1440))
dudx2 = np.zeros((721,1440))

dvdx = np.zeros((721,1440))
dvdx2 = np.zeros((721,1440))

dudy = np.zeros((721,1440))
dudy2 = np.zeros((721,1440))

dvdy = np.zeros((721,1440))
dvdy2 = np.zeros((721,1440))

for j in range(len(lat1)):
	for k in range(1,len(lon1)-1):
		dudx[j,k] = (uu1[10,j,k+1]-uu1[10,j,k-1])/(dx[j,k]*((6371000+((ggp1[10,j,k+1]+ggp1[10,j,k-1])/2))/6371000))
		dudx2[j,k] = (uu2[10,j,k+1]-uu2[10,j,k-1])/(dx[j,k]*((6371000+((ggp2[10,j,k+1]+ggp2[10,j,k-1])/2))/6371000))
		dvdx[j,k] = (vv1[10,j,k+1]-vv1[10,j,k-1])/(dx[j,k]*((6371000+((ggp1[10,j,k+1]+ggp1[10,j,k-1])/2))/6371000))
		dvdx2[j,k] = (vv2[10,j,k+1]-vv2[10,j,k-1])/(dx[j,k]*((6371000+((ggp2[10,j,k+1]+ggp2[10,j,k-1])/2))/6371000))

	  		
for j in range(len(lat1)):
	dudx[j,-1] = (uu1[10,j,0]-uu1[10,j,-2])/(dx[j,-1]*((6371000+((ggp1[10,j,0]+ggp1[10,j,2])/2))/6371000))
	dudx2[j,-1] = (uu2[10,j,0]-uu2[10,j,-2])/(dx[j,-1]*((6371000+((ggp2[10,j,0]+ggp2[10,j,2])/2))/6371000))
	dvdx[j,-1] = (vv1[10,j,0]-vv1[10,j,-2])/(dx[j,-1]*((6371000+((ggp1[10,j,0]+ggp1[10,j,2])/2))/6371000))
	dvdx2[j,-1] = (vv2[10,j,0]-vv2[10,j,-2])/(dx[j,-1]*((6371000+((ggp2[10,j,0]+ggp2[10,j,2])/2))/6371000))
	dudx[j,0] = (uu1[10,j,-1]-uu1[10,j,1])/(dx[j,0]*((6371000+((ggp1[10,j,0]+ggp1[10,j,2])/2))/6371000))
	dudx2[j,0] = (uu2[10,j,-1]-uu2[10,j,1])/(dx[j,0]*((6371000+((ggp2[10,j,0]+ggp2[10,j,2])/2))/6371000))
	dvdx[j,0] = (vv1[10,j,-1]-vv1[10,j,1])/(dx[j,0]*((6371000+((ggp1[10,j,0]+ggp1[10,j,2])/2))/6371000))
	dvdx2[j,0] = (vv2[10,j,-1]-vv2[10,j,1])/(dx[j,0]*((6371000+((ggp2[10,j,0]+ggp2[10,j,2])/2))/6371000)) 

for j in range(1,len(lat1)-1):
	for k in range(len(lon1)):
		dudy[j,k] = (uu1[10,j+1,k]-uu1[10,j-1,k])/(dy[j,k]*((6371000+((ggp1[10,j+1,k]+ggp1[10,j-1,k])/2))/6371000))
		dudy2[j,k] = (uu2[10,j+1,k]-uu2[10,j-1,k])/(dy[j,k]*((6371000+((ggp2[10,j+1,k]+ggp2[10,j-1,k])/2))/6371000))

		dvdy[j,k] = (vv1[10,j+1,k]-vv1[10,j-1,k])/(dy[j,k]*((6371000+((ggp1[10,j+1,k]+ggp1[10,j-1,k])/2))/6371000))
		dvdy2[j,k] = (vv2[10,j+1,k]-vv2[10,j-1,k])/(dy[j,k]*((6371000+((ggp2[10,j+1,k]+ggp2[10,j-1,k])/2))/6371000))


for j in range(len(lon1)):
	dudy[-1,j] = dudy[-2,j]
	dudy2[-1,j] = dudy2[-2,j]
	dvdy[-1,j] = dvdy[-2,j]
	dvdy2[-1,j] = dvdy2[-2,j]
	dudy[0,j] = dudy[0,j]
	dudy2[0,j] = dudy2[0,j]
	dvdy[0,j] = dvdy[0,j]
	dvdy2[0,j] = dvdy2[0,j]

VWS = np.zeros((721,1440))
VWS2 = np.zeros((721,1440))
CVG = np.zeros((721,1440))
CVG2 = np.zeros((721,1440))
DEF = np.zeros((721,1440))
DEF2 = np.zeros((721,1440))
TI1 = np.zeros((721,1440))
TI12 = np.zeros((721,1440))
TI2 = np.zeros((721,1440))
TI22 = np.zeros((721,1440))

for i in range(len(lat1)):
	for j in range(len(lon1)):
		VWS[i,j] = abs((((uu1[11,i,j]-uu1[9,i,j])**(2.)+(vv1[11,i,j]-vv1[9,i,j])**(2.))**(1./2.))/(ggp1[11,i,j]-ggp1[9,i,j]))
		VWS2[i,j] = abs((((uu2[11,i,j]-uu2[9,i,j])**(2.)+(vv2[11,i,j]-vv2[9,i,j])**(2.))**(1./2.))/(ggp2[11,i,j]-ggp2[9,i,j]))
		CVG[i,j] = -(dudx[i,j]+dvdy[i,j])
		CVG2[i,j] = -(dudx2[i,j]+dvdy2[i,j])
		DEF[i,j] = ((dudx[i,j]-dvdy[i,j])**(2.)+(dvdx[i,j]+dudy[i,j])**(2.))**(1./2.)
		DEF2[i,j] = ((dudx2[i,j]-dvdy2[i,j])**(2.)+(dvdx2[i,j]+dudy2[i,j])**(2.))**(1./2.)
		TI1[i,j] = (VWS[i,j])*(DEF[i,j])*(10**(7.))
		TI12[i,j] = (VWS2[i,j])*(DEF2[i,j])*(10**(7.))
		TI2[i,j] = (VWS[i,j])*(DEF[i,j]+CVG[i,j])*(10**(7.))
		TI22[i,j] = (VWS2[i,j])*(DEF2[i,j]+CVG2[i,j])*(10**(7.))


titles= ['(a) 20181031_06UTC TI1 Global','(b) 20181031_06UTC TI1 Narita','(c) 20181031_06UTC TI2 Global','(d) 20181031_06UTC TI2 Narita', '(e) 20181031_12UTC TI1 Global','(f) 20181031_12UTC TI1 Narita', '(g) 20181031_12UTC TI2 Global','(h) 20181031_12UTC TI2 Narita']
fig = plt.figure(figsize=(40,20))
plt.subplot(4, 2, 1)
m = Basemap(llcrnrlat=-90,urcrnrlat=90,\
llcrnrlon=0,urcrnrlon=360)#,resolution='c')
plt.title(titles[0], fontsize=12, loc='left')
m.drawcoastlines()
m.drawcountries()
m.drawparallels(np.arange(-90,90,10),labels=[1,0,0,0])
m.drawmeridians(np.arange(0,360,30),labels=[0,0,0,1])
m.fillcontinents(alpha=0)
m.drawmapboundary(fill_color='w')
lon_2d1, lat_2d1 = np.meshgrid(lon1, lat1)
x, y = m(lon_2d1, lat_2d1)
clevs= np.arange(0.5,11.5,0.5)
TI1[ TI1 > max(clevs) ] = max(clevs)
cs=m.contourf(x,y,TI1,clevs,cmap='Blues')
cbar=m.colorbar(cs)#, cax=cbaxes)
cbar.set_label('10^(-7)sec^(-2)')

plt.subplot(4, 2, 2)
m = Basemap(llcrnrlat=33,urcrnrlat=37,\
llcrnrlon=138,urcrnrlon=142)#,resolution='c')
plt.title(titles[1], fontsize=12, loc='left')
m.drawcoastlines()
m.drawcountries()
m.drawparallels(np.arange(33,37,1),labels=[1,0,0,0])
m.drawmeridians(np.arange(138,142,1),labels=[0,0,0,1])
m.fillcontinents(alpha=0)
m.drawmapboundary(fill_color='w')
lon_2d1, lat_2d1 = np.meshgrid(lon1, lat1)
x, y = m(lon_2d1, lat_2d1)
clevs= np.arange(0.5,11.5,0.5)
TI1[ TI1 > max(clevs) ] = max(clevs)
cs=m.contourf(x,y,TI1,clevs, cmap='Blues')
cbar=m.colorbar(cs)#, cax=cbaxes)
cbar.set_label('10^(-7)sec^(-2)')

plt.subplot(4, 2, 3)
m = Basemap(llcrnrlat=-90,urcrnrlat=90,\
llcrnrlon=0,urcrnrlon=360)#,resolution='c')
plt.title(titles[2], fontsize=12, loc='left')
m.drawcoastlines()
m.drawcountries()
m.drawparallels(np.arange(-90,90,10),labels=[1,0,0,0])
m.drawmeridians(np.arange(0,360,30),labels=[0,0,0,1])
m.drawmapboundary(fill_color='w')
lon_2d1, lat_2d1 = np.meshgrid(lon1, lat1)
x, y = m(lon_2d1, lat_2d1)
clevs= np.arange(0.5,11.5,0.5)
TI2[ TI2 > max(clevs) ] = max(clevs)
cs=m.contourf(x,y,TI2,clevs, cmap='Blues')
cbar=m.colorbar(cs)#, cax=cbaxes)
cbar.set_label('10^(-7)sec^(-2)')

plt.subplot(4, 2, 4)
m = Basemap(llcrnrlat=33,urcrnrlat=37,\
llcrnrlon=138,urcrnrlon=142)#,resolution='c')
plt.title(titles[3], fontsize=12, loc='left')
m.drawcoastlines()
m.drawcountries()
m.drawparallels(np.arange(33,37,1),labels=[1,0,0,0])
m.drawmeridians(np.arange(138,142,1),labels=[0,0,0,1])
m.fillcontinents(alpha=0)
m.drawmapboundary(fill_color='w')
lon_2d1, lat_2d1 = np.meshgrid(lon1, lat1)
x, y = m(lon_2d1, lat_2d1)
clevs= np.arange(0.5,11.5,0.5)
TI2[ TI2 > max(clevs) ] = max(clevs)
cs=m.contourf(x,y,TI2,clevs, cmap='Blues')
cbar=m.colorbar(cs)#, cax=cbaxes)
cbar.set_label('10^(-7)sec^(-2)')

plt.subplot(4, 2, 5)
m = Basemap(llcrnrlat=-90,urcrnrlat=90,\
llcrnrlon=0,urcrnrlon=360)#,resolution='c')
plt.title(titles[4], fontsize=12, loc='left')
m.drawcoastlines()
m.drawcountries()
m.drawparallels(np.arange(-90,90,10),labels=[1,0,0,0])
m.drawmeridians(np.arange(0,360,30),labels=[0,0,0,1])
m.drawmapboundary(fill_color='w')
lon_2d1, lat_2d1 = np.meshgrid(lon1, lat1)
x, y = m(lon_2d1, lat_2d1)
clevs= np.arange(0.5,11.5,0.5)
TI12[ TI12 > max(clevs) ] = max(clevs)
cs=m.contourf(x,y,TI12,clevs, cmap='Blues')
cbar=m.colorbar(cs)#, cax=cbaxes)
cbar.set_label('10^(-7)sec^(-2)')

plt.subplot(4, 2, 6)
m = Basemap(llcrnrlat=33,urcrnrlat=37,\
llcrnrlon=138,urcrnrlon=142)#,resolution='c')
plt.title(titles[5], fontsize=12, loc='left')
m.drawcoastlines()
m.drawcountries()
m.drawparallels(np.arange(33,37,1),labels=[1,0,0,0])
m.drawmeridians(np.arange(138,142,1),labels=[0,0,0,1])
m.fillcontinents(alpha=0)
m.drawmapboundary(fill_color='w')
lon_2d1, lat_2d1 = np.meshgrid(lon1, lat1)
x, y = m(lon_2d1, lat_2d1)
clevs= np.arange(0.5,11.5,0.5)
TI12[ TI12 > max(clevs) ] = max(clevs)
cs=m.contourf(x,y,TI12,clevs,cmap='Blues')
cbar=m.colorbar(cs)#, cax=cbaxes)
cbar.set_label('10^(-7)sec^(-2)')

plt.subplot(4, 2, 7)
m = Basemap(llcrnrlat=-90,urcrnrlat=90,\
llcrnrlon=0,urcrnrlon=360)#,resolution='c')
plt.title(titles[6], fontsize=12, loc='left')
m.drawcoastlines()
m.drawcountries()
m.drawparallels(np.arange(-90,90,10),labels=[1,0,0,0])
m.drawmeridians(np.arange(0,360,30),labels=[0,0,0,1])
m.drawmapboundary(fill_color='w')
lon_2d1, lat_2d1 = np.meshgrid(lon1, lat1)
x, y = m(lon_2d1, lat_2d1)
clevs= np.arange(0.5,11.5,0.5)
TI22[ TI22 > max(clevs) ] = max(clevs)
cs=m.contourf(x,y,TI22,clevs, cmap='Blues')
cbar=m.colorbar(cs)#, cax=cbaxes)
cbar.set_label('10^(-7)sec^(-2)')

plt.subplot(4, 2, 8)
m = Basemap(llcrnrlat=33,urcrnrlat=37,\
llcrnrlon=138,urcrnrlon=142)#,resolution='c')
plt.title(titles[7], fontsize=12, loc='left')
m.drawcoastlines()
m.drawcountries()
m.drawparallels(np.arange(33,37,1),labels=[1,0,0,0])
m.drawmeridians(np.arange(138,142,1),labels=[0,0,0,1])
m.fillcontinents(alpha=0)
m.drawmapboundary(fill_color='w')
lon_2d1, lat_2d1 = np.meshgrid(lon1, lat1)
x, y = m(lon_2d1, lat_2d1)
clevs= np.arange(0.5,11.5,0.5)
TI22[ TI22 > max(clevs) ] = max(clevs)
cs=m.contourf(x,y,TI22,clevs, cmap='Blues')
cbar=m.colorbar(cs)#, cax=cbaxes)
cbar.set_label('10^(-7)sec^(-2)')

fig = plt.savefig('./figures/T1')


		
