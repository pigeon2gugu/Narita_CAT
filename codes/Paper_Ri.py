from netCDF4 import Dataset
import numpy as np
from mpl_toolkits.basemap import Basemap, shiftgrid
import matplotlib.pyplot as plt
from sys import exit
import matplotlib as mpl

lon = Dataset('data/2018-10-30_06:00:00.nc').variables['g0_lon_2'][:]
lat = Dataset('data/2018-10-30_06:00:00.nc').variables['g0_lat_1'][:]
lv = Dataset('data/2018-10-30_06:00:00.nc').variables['lv_ISBL0'][:]
gph = np.zeros((11,32,721,3))
u = np.zeros((11,32,721,3))
v = np.zeros((11,32,721,3))
T = np.zeros((11,32,721,3))
wind = np.zeros((11,32,721,3))
PT = np.zeros((11,32,721,3))

k = 0
dd = ['28','29','30']
tt = ['00','06','12','18']

for i  in range(3):
        for j in range(4):
                        if dd[i] == '30' and tt[j] == '18':
                                break
                        else :
                                gph[k,:,:,:] = Dataset('data/2018-10-'+dd[i]+'_'+tt[j]+':00:00.nc').variables['Z_GDS0_ISBL'][:,:,561:564]
				u[k,:,:,:] = Dataset('data/2018-10-'+dd[i]+'_'+tt[j]+':00:00.nc').variables['U_GDS0_ISBL'][:,:,561:564]
				v[k,:,:,:] = Dataset('data/2018-10-'+dd[i]+'_'+tt[j]+':00:00.nc').variables['V_GDS0_ISBL'][:,:,561:564]
				T[k,:,:,:] = Dataset('data/2018-10-'+dd[i]+'_'+tt[j]+':00:00.nc').variables['T_GDS0_ISBL'][:,:,561:564]
                                k = k+1


gph = gph/9.806

for i in range(11):
	for l in range(32):
        	for j in range(len(lat)):
                        wind[i,l,j,1] = (((u[i,l,j,1])**(2.))+((v[i,l,j,1])**(2.)))**(1./2.)
			PT[i,l,j,1] = ((T[i,l,j,1])*((1000./float(int(lv[l])))**0.286))

dx = np.zeros((721,1440))
dy = np.zeros((721,1440))


for i in range(len(lat)):
	dx[i,562] = 2*np.pi*6371000*np.cos(lat[i]*2*np.pi/360.)*(lon[563]-lon[561])/360.

for i in range(len(lat)):
        dx[i,-1] = dx[i,-2]
        dx[i,0] = dx[i,1]

for i in range(1,len(lat)-1):
                dy[i,562] = 2*np.pi*6371000*(lat[i+1]-lat[i-1])/360.


for j in range(len(lon)):
        dy[-1,j] = 2*np.pi*6371000*(0.5)/360.
        dy[0,j] = 2*np.pi*6371000*(0.5)/360.


dvdx = np.zeros((11,32,721,3))
dudy = np.zeros((11,32,721,3))

for i in range(11):
	for l in range(32):
		for j in range(len(lat)):
                	dvdx[i,l,j,1] = (v[i,l,j,2]-v[i,l,j,0])/(dx[j,562])

for i in range(11):
	for l in range(32):
		for j in range(1,len(lat)-1):
                	dudy[i,l,j,1] = (u[i,l,j+1,1]-u[i,l,j-1,1])/(dy[j,562])

for i in range(11):
	for l in range(32):
        	dudy[i,l,-1,1] = dudy[i,l,-2,1]
        	dudy[i,l,0,1] = dudy[i,l,1,1]

Ns = np.zeros((11,32,721,3))
VWS = np.zeros((11,32,721,3))
Ri = np.zeros((11,32,721,3))

for i in range(11):
	for l in range(1,31):
		for j in range(len(lat)):
			Ns[i,l,j,1] = ((9.806*(PT[i,l-1,j,1]-PT[i,l+1,j,1]))/(PT[i,l,j,1]*(gph[i,l-1,j,1]-gph[i,l+1,j,1])))
			VWS[i,l,j,1] = abs((((u[i,l-1,j,1]-u[i,l+1,j,1])**(2.)+(v[i,l-1,j,1]-v[i,l+1,j,1])**(2.))**(1./2.))/(gph[i,l-1,j,1]-gph[i,l+1,j,1]))
			Ri[i,l,j,1] = (Ns[i,l,j,1])/((VWS[i,l,j,1])**(2.))		

'''
for i in range(11):
        for l in range(1,31):
                for j in range(1,len(lat)-1):
                        Ri[i,l,j,1] = (((Ri[i,l,j-1,1] + ((2.)*Ri[i,l,j,1]) + Ri[i,l,j+1,1])/(4.))+((Ri[i,l-1,j,1] + ((2.)*Ri[i,l,j,1]) + Ri[i,l+1,j,2])/(4.)))/(2.)

for i in range(11):
        for l in range(1,31):
                Ri[i,l,0,1] = Ri[i,l,1,1]
                Ri[i,l,-1,1] = Ri[i,l,-2,1]
'''
Ri = np.transpose(Ri, (0,2,1,3))
wind = np.transpose(wind, (0,2,1,3))
PT = np.transpose(PT, (0,2,1,3))
u = np.transpose(u, (0,2,1,3))

lv = lv[::-1]
lat = lat[::-1]
Ri = Ri[:,::-1,::-1,:]
wind = wind[:,::-1,::-1,:]
PT = PT[:,::-1,::-1,:]
u = u[:,::-1,::-1,:]
v = v.astype(int)

n = 2
titles= ['(a) 20181028 12UTC', '(b) 20181029 12UTC', '(c) 20181030 12UTC']
fig = plt.figure(figsize=(15,10))

for i in range(3):
        plt.subplot(3, 4, i+1)
	plt.title(titles[i], fontsize=12, loc='left')
	cmap = mpl.cm.coolwarm
	x, y = np.meshgrid(lv[1:31], lat[480:541])
	clevs1= np.arange(0,1.25,0.25)
	clevs2 = np.arange(150,450,10)
	clevs3 = np.arange(-50,0,10)
	clevs4 = np.arange(0,90,10)
	plt.plot(lat[505],lv[13], color='r', marker='*')
	cs = plt.contour(y, x, PT[n,480:541,1:31,1], clevs2, alpha = 0.7 , colors='black', linewidths=0.7)
	cs1 = plt.contourf(y, x, Ri[n,480:541,1:31,1], clevs1,  cmap = cmap)
	cs2 = plt.contour(y, x, u[n,480:541,1:31,1], clevs3, alpha = 0.7, colors = 'red', linestyles= 'dashed',linewidths=0.7)
	cs3 = plt.contour(y, x, u[n,480:541,1:31,1], clevs4, alpha = 0.7, colors = 'red', linestyles = 'solid',linewidths=0.7)
	plt.ylabel('Pressure $(hPa)$',size=10)
	plt.xlabel('Latitude', size=10)
	cbar=plt.colorbar(cs1)
	#cbar.ax.set_yticklabels(['-2','-1.5','-1','-0.5','0','0.5','1','1.5','2'])
	plt.clabel(cs, fmt='%3.0f', colors = 'black')
	plt.clabel(cs2, fmt='%3.0f', colors = 'red')
	plt.clabel(cs3, fmt='%3.0f', colors = 'red')
	cbar.set_label('Richardson Number', fontsize = 7)
	plt.gca().invert_yaxis()
	n = n+4

plt.subplots_adjust(wspace = 0.42, hspace = 0.30)
fig = plt.savefig('./figures/Paper_Ri')
