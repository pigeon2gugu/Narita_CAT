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
gph = np.zeros((3,37,721,3))
u = np.zeros((3,37,721,3))
v = np.zeros((3,37,721,3))
T = np.zeros((3,37,721,3))
wind = np.zeros((3,37,721,3))
o = np.zeros((3,37,721,3))
PT = np.zeros((3,37,721,3))

k = 0
dd = ['28','29','30']

for i  in range(3):
	gph[i,:,:,:] = Dataset('data/2018_10_'+dd[i]+'_12:00:00.nc').variables['Z_GDS0_ISBL'][:,:,561:564]
	u[i,:,:,:] = Dataset('data/2018_10_'+dd[i]+'_12:00:00.nc').variables['U_GDS0_ISBL'][:,:,561:564]
        o[i,:,:,:] = Dataset('data/2018_10_'+dd[i]+'_12:00:00.nc').variables['W_GDS0_ISBL'][:,:,561:564]
	v[i,:,:,:] = Dataset('data/2018_10_'+dd[i]+'_12:00:00.nc').variables['V_GDS0_ISBL'][:,:,561:564]
	T[i,:,:,:] = Dataset('data/2018_10_'+dd[i]+'_12:00:00.nc').variables['T_GDS0_ISBL'][:,:,561:564]


gph = gph/9.806

dtr = np.pi/180.
theta = lat*dtr
dthe = 0.25
Re = 6371000.
dx = Re*np.cos(theta)*dthe*dtr
dy = Re*dthe*dtr

dvdx = np.zeros((3,37,721,3))
dudy = np.zeros((3,37,721,3))
ua = np.zeros((3,37,721,3))
va = np.zeros((3,37,721,3))


for i in range(3):
        for l in range(37):
                for j in range(110,290):
			PT[i,l,j,1] = ((T[i,l,j,1])*((1000./float(int(lv[l])))**0.286))
        		dvdx[i,l,j,1] = (v[i,l,j,2]-v[i,l,j,0])/(2*dx[j])
        		dudy[i,l,j,1] = (u[i,l,j+1,1]-u[i,l,j-1,1])/(-2*dy)

PVU = np.zeros((3,37,721,3))

for i in range(3):
	for l in range(1,36):
		for j in range(110,290):
			PVU[i,l,j,1] = (10**6)*(-9.806)*((PT[i,l-1,j,1]-PT[i,l+1,j,1])/(100.*(float(lv[l-1])-float(lv[l+1]))))*((dvdx[i,l,j,1]-dudy[i,l,j,1])+((2.)*(2.*np.pi/86400.)*(np.sin(lat[j]*dtr))))
			ua[i,l,j,1] = u[i,l,j,1] - ((-1./(2*(2.*np.pi/86400.)*(np.sin(lat[j]*dtr))))*((gph[i,l,j+1,1]-gph[i,l,j-1,1])/(-2*dy)))
			va[i,l,j,1] = v[i,l,j,1] - ((1./(2*(2.*np.pi/86400.)*(np.sin(lat[j]*dtr))))*((gph[i,l,j,2]-gph[i,l,j,0])/(2*dx[j])))

'''
PVU = np.transpose(PVU, (0,2,1,3))
PT = np.transpose(PT, (0,2,1,3))
va = np.transpose(va, (0,2,1,3))
o = np.transpose(o, (0,2,1,3))
lv = lv[::-1]
lat = lat[::-1]
PVU = PVU[:,::-1,::-1,:]
PT = PT[:,::-1,::-1,:]
va = va[:,::-1,::-1,:]
o = o[:,::-1,::-1,:]
'''

titles= ['(a) 20181028 12UTC','(b) 20181029 12UTC', '(c) 20181030 12UTC']
fig = plt.figure(figsize=(8,8))
#for i in range(3):
	#plt.subplot(1, 3, i+1)
plt.title(titles[0], fontsize=12, loc='left')
x, y = np.meshgrid(lat, lv)
plt.plot(lat[215],lv[23], color='r', marker='*')
cmap = matplotlib.colors.ListedColormap(['gold','orange'])
clevs1 = np.arange(1.5,3.0,0.5)
clevs2 = np.arange(150,450,10)
PVU[ PVU > max(clevs1) ] = max(clevs1)
cs = plt.contour(x, y, PT[0,:,:,1], clevs2, alpha = 0.7 , colors='black', linewidths=0.7)
cs1 = plt.contourf(x, y, PVU[0,:,:,1], clevs1, cmap = cmap)
plt.quiver(x[:,::8], y[:,::8],va[0,:,::8,1],15.*o[0,:,::8,1],pivot='middle', scale = 1, units = 'xy',  headwidth=2, headlength=2, linewidth=0.02)
plt.ylabel('Pressure $(hPa)$',size=10)
plt.xlabel('Latitude', size=10)
axes = plt.gca()
axes.set_xlim([20,60])
axes.set_ylim([10,1000])
plt.clabel(cs, fmt='%3.0f', colors = 'black')
plt.gca().invert_yaxis()

plt.subplots_adjust(wspace = 0.42, hspace = 0.30)
plt.show()

#fig = plt.savefig('./figures/ageostrophic_PVU')
