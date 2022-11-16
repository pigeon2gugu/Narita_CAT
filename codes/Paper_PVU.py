from netCDF4 import Dataset
import numpy as np
import os
os.environ['PROJ_LIB'] = '/usr/local/python/2.7/share/proj'
from mpl_toolkits.basemap import Basemap, shiftgrid
import matplotlib.pyplot as plt
from sys import exit
import matplotlib as mpl

lon = Dataset('data/2018-10-30_06:00:00.nc').variables['g0_lon_2'][:]
lat = Dataset('data/2018-10-30_06:00:00.nc').variables['g0_lat_1'][:]
lv = Dataset('data/2018-10-30_06:00:00.nc').variables['lv_ISBL0'][:]
gph = np.zeros((6,32,721,3))
u = np.zeros((6,32,721,3))
v = np.zeros((6,32,721,3))
T = np.zeros((6,32,721,3))
wind = np.zeros((6,32,721,3))
PT = np.zeros((6,32,721,3))

k = 0
dd = ['25','26','27','28','29','30']

for i  in range(6):
	gph[i,:,:,:] = Dataset('data/2018-10-'+dd[i]+'_12:00:00.nc').variables['Z_GDS0_ISBL'][:,:,561:564]
	u[i,:,:,:] = Dataset('data/2018-10-'+dd[i]+'_12:00:00.nc').variables['U_GDS0_ISBL'][:,:,561:564]
	v[i,:,:,:] = Dataset('data/2018-10-'+dd[i]+'_12:00:00.nc').variables['V_GDS0_ISBL'][:,:,561:564]
	T[i,:,:,:] = Dataset('data/2018-10-'+dd[i]+'_12:00:00.nc').variables['T_GDS0_ISBL'][:,:,561:564]


gph = gph/9.806

for i in range(6):
	for l in range(32):
        	for j in range(len(lat)):
                        wind[i,l,j,1] = (((u[i,l,j,1])**(2.))+((v[i,l,j,1])**(2.)))**(1./2.)
			PT[i,l,j,1] = ((T[i,l,j,1])*((1000./float(int(lv[l])))**0.286))


dtr = np.pi/180.
theta = lat*dtr
dthe = 0.25
Re = 6371000.
dx = Re*np.cos(theta)*dthe*dtr
dy = Re*dthe*dtr


dvdx = np.zeros((6,32,721,3))
dudy = np.zeros((6,32,721,3))

for i in range(6):
	for l in range(32):
		for j in range(len(lat)):
                	dvdx[i,l,j,1] = (v[i,l,j,2]-v[i,l,j,0])/(2*dx[j])

for i in range(6):
	for l in range(32):
		for j in range(1,len(lat)-1):
                	dudy[i,l,j,1] = (u[i,l,j+1,1]-u[i,l,j-1,1])/(-2*dy)

for i in range(6):
	for l in range(32):
        	dudy[i,l,-1,1] = dudy[i,l,-2,1]
        	dudy[i,l,0,1] = dudy[i,l,1,1]

PV = np.zeros((6,32,721,3))
abv = np.zeros((6,32,721,3))

for i in range(6):
	for l in range(1,31):
		for j in range(len(lat)):
			PV[i,l,j,1] = (-9.806)*((PT[i,l-1,j,1]-PT[i,l+1,j,1])/(100.*(float(lv[l-1])-float(lv[l+1]))))*((dvdx[i,l,j,1]-dudy[i,l,j,1])+((2.)*(2.*np.pi/86400.)*(np.sin(lat[j]*np.pi/180.))))

'''
for i in range(11):
        for l in range(1,31):
                for j in range(1,len(lat)-1):
                        PV[i,l,j,1] = (((PV[i,l,j-1,1] + ((2.)*PV[i,l,j,1]) + PV[i,l,j+1,1])/(4.))+((PV[i,l-1,j,1] + ((2.)*PV[i,l,j,1]) + PV[i,l+1,j,2])/(4.)))/(2.)

for i in range(11):
        for l in range(1,31):
                PV[i,l,0,1] = PV[i,l,1,1]
                PV[i,l,-1,1] = PV[i,l,-2,1]
'''

PV = np.transpose(PV, (0,2,1,3))
wind = np.transpose(wind, (0,2,1,3))
PT = np.transpose(PT, (0,2,1,3))
u = np.transpose(u, (0,2,1,3))

lv = lv[::-1]
lat = lat[::-1]
PV = PV[:,::-1,::-1,:]
wind = wind[:,::-1,::-1,:]
PT = PT[:,::-1,::-1,:]
u = u[:,::-1,::-1,:]
PVU = PV*(1000000)
v = v.astype(int)

titles= ['(a) 20181025 12UTC','(b) 20181026 12UTC', '(c) 20181027 12UTC', '(d) 20181028 12UTC', '(e) 20181029 12UTC', '(f) 20181030 12UTC']
fig = plt.figure(figsize=(15,10))

for i in range(6):
        plt.subplot(2, 3, i+1)
	plt.title(titles[i], fontsize=12, loc='left')
	cmap = mpl.cm.coolwarm
	x, y = np.meshgrid(lv[1:31], lat[480:541])
	clevs1= np.arange(1,16,1)
	PVU[ PVU > max(clevs1) ] = max(clevs1)
	clevs2 = np.arange(150,400,5)
	clevs3 = np.arange(-50,0,10)
	clevs4 = np.arange(0,90,10)
	clevs5 = 1.5
	plt.plot(lat[505],lv[13], color='r', marker='*')
	plt.contour(y,x, PVU[i,480:541,1:31,1], clevs5, alpha = 1 , colors = 'blue' , linewidths = 1.5)
	cs = plt.contour(y, x, PT[i,480:541,1:31,1], clevs2, alpha = 0.7 , colors='black', linewidths=0.7)
	cs1 = plt.contourf(y, x, PVU[i,480:541,1:31,1], clevs1, cmap = cmap)
	cs2 = plt.contour(y, x, u[i,480:541,1:31,1], clevs3, alpha = 0.7, colors = 'red', linestyles= 'dashed',linewidths=0.7)
	cs3 = plt.contour(y, x, u[i,480:541,1:31,1], clevs4, alpha = 0.7, colors = 'red', linestyles = 'solid',linewidths=0.7)
	plt.ylabel('Pressure $(hPa)$',size=10)
	plt.xlabel('Latitude', size=10)
	cbar=plt.colorbar(cs1)
	cbar.ax.set_yticklabels(['1','3','5','7','9','11','13','>15'])
	plt.clabel(cs, fmt='%3.0f', colors = 'black')
	plt.clabel(cs2, fmt='%3.0f', colors = 'red')
	plt.clabel(cs3, fmt='%3.0f', colors = 'red')
	cbar.set_label('Potential Vorticity Unit $(10^{-6}K*m^{2}/kg*s)$', fontsize = 7.5)
	plt.gca().invert_yaxis()

plt.subplots_adjust(wspace = 0.42, hspace = 0.30)
fig = plt.savefig('./figures/Paper_PVU')
