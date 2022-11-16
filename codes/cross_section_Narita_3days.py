from netCDF4 import Dataset
import numpy as np
from mpl_toolkits.basemap import Basemap, shiftgrid
import matplotlib.pyplot as plt
from sys import exit
import matplotlib as mpl
import matplotlib as mpl

lon = Dataset('data/2018-10-30_06:00:00.nc').variables['g0_lon_2'][:]
lat = Dataset('data/2018-10-30_06:00:00.nc').variables['g0_lat_1'][:]
lv = Dataset('data/2018-10-30_06:00:00.nc').variables['lv_ISBL0'][:]
gph500 = np.zeros((11,32,721,1440))
gph550 = np.zeros((11,32,721,1440))
gph450 = np.zeros((11,32,721,1440))
u = np.zeros((11,32,721,1440))
v = np.zeros((11,32,721,1440))
u550 = np.zeros((11,32,721,1440))
u450 = np.zeros((11,32,721,1440))
v550 = np.zeros((11,32,721,1440))
v450 = np.zeros((11,32,721,1440))

k = 0
dd = ['28','29','30']
tt = ['00','06','12','18']

for i  in range(3):
        for j in range(4):
                        if dd[i] == '30' and tt[j] == '18':
                                break
                        else :
                                gph[k,:,:,:] = Dataset('data/2018-10-'+dd[i]+'_'+tt[j]+':00:00.nc').variables['Z_GDS0_ISBL'][:,:,:]
				u[k,:,:,:] = Dataset('data/2018-10-'+dd[i]+'_'+tt[j]+':00:00.nc').variables['U_GDS0_ISBL'][:,:,:]
				v[k,:,:,:] = Dataset('data/2018-10-'+dd[i]+'_'+tt[j]+':00:00.nc').variables['V_GDS0_ISBL'][:,:,:]
                                k = k+1


gph = gph/9.806

dx = np.zeros((721,1440))
dy = np.zeros((721,1440))

for i in range(len(lat)):
        for j in range(1,len(lon)-1):
                dx[i,j] = 2*np.pi*6371000*np.cos(lat[i]*2*np.pi/360.)*(lon[j+1]-lon[j-1])/360.

for i in range(len(lat)):
        dx[i,-1] = dx[i,-2]
        dx[i,0] = dx[i,1]

for i in range(1,len(lat)-1):
        for j in range(len(lon)):
                dy[i,j] = 2*np.pi*6371000*(lat[i+1]-lat[i-1])/360.


for j in range(len(lon)):
        dy[-1,j] = 2*np.pi*6371000*(0.5)/360.
        dy[0,j] = 2*np.pi*6371000*(0.5)/360.



dudx = np.zeros((11,32,721,1440))
dvdx = np.zeros((11,32,721,1440))
dudy = np.zeros((11,32,721,1440))
dvdy = np.zeros((11,32,721,1440))

for i in range(11):
	for l in range(32):
		for j in range(len(lat)):
        		for k in range(1,len(lon)-1):
               			dudx[i,l,j,k] = (u[i,l,j,k+1]-u[i,l,j,k-1])/(dx[j,k])
                		dvdx[i,l,j,k] = (v[i,l,j,k+1]-v[i,l,j,k-1])/(dx[j,k])

for i in range(11):
	for l in range(32):
		for j in range(len(lat)):
        		dudx[i,l,j,-1] = (u[i,l,j,0]-u[i,l,j,-2])/(dx[j,-1])
        		dvdx[i,l,j,-1] = (v[i,l,j,0]-v[i,l,j,-2])/(dx[j,-1])
        		dudx[i,l,j,0] = (u[i,l,j,-1]-u[i,l,j,1])/(dx[j,0])
        		dvdx[i,l,j,0] = (v[i,l,j,-1]-v[i,l,j,1])/(dx[j,0])

for i in range(11):
	for l in range(32):
		for j in range(1,len(lat)-1):
        		for k in range(len(lon)):
                		dudy[i,l,j,k] = (u[i,l,j+1,k]-u[i,l,j-1,k])/(dy[j,k])
                		dvdy[i,l,j,k] = (v[i,l,j+1,k]-v[i,l,j-1,k])/(dy[j,k])

for i in range(11):
	for l in range(32):
		for j in range(len(lon)):
        		dudy[i,l,-1,j] = dudy[i,l,-2,j]
        		dvdy[i,l,-1,j] = dvdy[i,l,-2,j]
        		dudy[i,l,0,j] = dudy[i,l,1,j]
        		dvdy[i,l,0,j] = dvdy[i,l,1,j]

VWS = np.zeros((11,32,721,1440))
DEF = np.zeros((11,32,721,1440))
TI1 = np.zeros((11,32,721,1440))

for i in range(11):
	for l in range(1,31):
		for j in range(len(lat)):
        		for k in range(len(lon)):
				if (gph[i,l-1,j,k]-gph550[i,l+1,j,k]) == 0 :
					VWS[i,l,j,k] = np.nan
				else :
					VWS[i,l,j,k] = abs((((u[i,l-1,j,k]-u550[i,l+1,j,k])**(2.)+(v[i,l-1,j,k]-v550[i,l+1,j,k])**(2.))**(1./2.))/(gph[i,l-1,j,k]-gph550[i,l+1,j,k]))
				
for i in range(11):
	for l in range(32):
		for j in range(len(lat)):
			for k in range(len(lon)):
				DEF[i,l,j,k] = ((dudx[i,l,j,k]-dvdy[i,l,j,k])**(2.)+(dvdx[i,l,j,k]+dudy[i,l,j,k])**(2.))**(1./2.)
				TI1[i,l,j,k] = (VWS[i,l,j,k])*(DEF[i,l,,j,k])

TI1 = TI1*(1000000)
TI1 = np.around(TI1, decimals = 0)

titles= ['(a) 20181028 00UTC TI','(b) 20181028 06UTC TI', '(c) 20181028 12UTC TI', '(d) 20181028 18UTC TI', '(e) 20181029 00UTC TI', '(f) 20181029 06UTC TI', '(g) 20181029 12UTC TI', '(h) 20181029 18UTC TI', '(i) 20181030 00UTC TI', '(j) 20181030 06UTC TI', '(k) 20181030 12UTC TI']

fig = plt.figure(figsize=(15,10))
for i in range(11):
        plt.subplot(3, 4, i+1)
        m = Basemap(llcrnrlat=20,urcrnrlat=60,\
        llcrnrlon=100,urcrnrlon=160, resolution='i')#,resolution='c')
        plt.title(titles[i], fontsize=12, loc='left')
        m.drawcoastlines()
        m.drawcountries()
        m.drawparallels(np.arange(20,60,10),labels=[1,0,0,0])
        m.drawmeridians(np.arange(100,160,10),labels=[0,0,0,1])
        m.fillcontinents(alpha=0)
        m.drawmapboundary(fill_color='w')
        lon_2d1, lat_2d1 = np.meshgrid(lon, lat)
        x, y = m(lon_2d1, lat_2d1)
        x1, y1 = m(lon[562], lat[215]) #CAT incounter
	clevs1= np.arange(0,11,1)
	m.plot(x1, y1, color='r', marker='*')
	cmap = mpl.cm.Blues
	cs2=m.contourf(x,y,TI1[i,:,:],clevs1, cmap = cmap)
	cbar=m.colorbar(cs2, vmin=0, vmax=11)
        cbar.set_label('TI $(10^-6/s^-2)$', fontsize = 7)

plt.show()
plt.subplots_adjust(wspace = 0.47, hspace = 0.3)

