from netCDF4 import Dataset
import numpy as np
from mpl_toolkits.basemap import Basemap, shiftgrid
import matplotlib.pyplot as plt
from sys import exit
import matplotlib as mpl
import matplotlib as mpl

lon = Dataset('data/2018-10-30_06:00:00.nc').variables['g0_lon_2'][:]
lat = Dataset('data/2018-10-30_06:00:00.nc').variables['g0_lat_1'][:]
gph600 = np.zeros((11,721,1440))
gph650 = np.zeros((11,721,1440))
gph550 = np.zeros((11,721,1440))
u = np.zeros((11,721,1440))
v = np.zeros((11,721,1440))
u650 = np.zeros((11,721,1440))
u550 = np.zeros((11,721,1440))
v650 = np.zeros((11,721,1440))
v550 = np.zeros((11,721,1440))

k = 0
dd = ['28','29','30']
tt = ['00','06','12','18']

for i  in range(3):
        for j in range(4):
                        if dd[i] == '30' and tt[j] == '18':
                                break
                        else :
                                gph600[k,:,:] = Dataset('data/2018-10-'+dd[i]+'_'+tt[j]+':00:00.nc').variables['Z_GDS0_ISBL'][18,:,:]
				gph650[k,:,:] = Dataset('data/2018-10-'+dd[i]+'_'+tt[j]+':00:00.nc').variables['Z_GDS0_ISBL'][19,:,:]
                                gph550[k,:,:] = Dataset('data/2018-10-'+dd[i]+'_'+tt[j]+':00:00.nc').variables['Z_GDS0_ISBL'][17,:,:]

				u[k,:,:] = Dataset('data/2018-10-'+dd[i]+'_'+tt[j]+':00:00.nc').variables['U_GDS0_ISBL'][18,:,:]
				v[k,:,:] = Dataset('data/2018-10-'+dd[i]+'_'+tt[j]+':00:00.nc').variables['V_GDS0_ISBL'][18,:,:]
				u650[k,:,:] = Dataset('data/2018-10-'+dd[i]+'_'+tt[j]+':00:00.nc').variables['U_GDS0_ISBL'][19,:,:]
                                u550[k,:,:] = Dataset('data/2018-10-'+dd[i]+'_'+tt[j]+':00:00.nc').variables['U_GDS0_ISBL'][17,:,:]
                                v650[k,:,:] = Dataset('data/2018-10-'+dd[i]+'_'+tt[j]+':00:00.nc').variables['V_GDS0_ISBL'][19,:,:]
                                v550[k,:,:] = Dataset('data/2018-10-'+dd[i]+'_'+tt[j]+':00:00.nc').variables['V_GDS0_ISBL'][17,:,:]

                                k = k+1


gph600 = gph600/9.806
gph550 = gph550/9.806
gph650 = gph650/9.806

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



dudx = np.zeros((11,721,1440))
dvdx = np.zeros((11,721,1440))
dudy = np.zeros((11,721,1440))
dvdy = np.zeros((11,721,1440))

for i in range(11):
	for j in range(len(lat)):
        	for k in range(1,len(lon)-1):
               		dudx[i,j,k] = (u[i,j,k+1]-u[i,j,k-1])/(dx[j,k])
                	dvdx[i,j,k] = (v[i,j,k+1]-v[i,j,k-1])/(dx[j,k])

for i in range(11):
	for j in range(len(lat)):
        	dudx[i,j,-1] = (u[i,j,0]-u[i,j,-2])/(dx[j,-1])
        	dvdx[i,j,-1] = (v[i,j,0]-v[i,j,-2])/(dx[j,-1])
        	dudx[i,j,0] = (u[i,j,-1]-u[i,j,1])/(dx[j,0])
        	dvdx[i,j,0] = (v[i,j,-1]-v[i,j,1])/(dx[j,0])

for i in range(11):
	for j in range(1,len(lat)-1):
        	for k in range(len(lon)):
                	dudy[i,j,k] = (u[i,j+1,k]-u[i,j-1,k])/(dy[j,k])
                	dvdy[i,j,k] = (v[i,j+1,k]-v[i,j-1,k])/(dy[j,k])

for i in range(11):
	for j in range(len(lon)):
        	dudy[i,-1,j] = dudy[i,-2,j]
        	dvdy[i,-1,j] = dvdy[i,-2,j]
        	dudy[i,0,j] = dudy[i,1,j]
        	dvdy[i,0,j] = dvdy[i,1,j]

VWS = np.zeros((11,721,1440))
DEF = np.zeros((11,721,1440))
TI1 = np.zeros((11,721,1440))

for i in range(11):
	for j in range(len(lat)):
        	for k in range(len(lon)):
			if (gph550[i,j,k]-gph650[i,j,k]) == 0 :
				VWS[i,j,k] = np.nan
			else :
				VWS[i,j,k] = abs((((u550[i,j,k]-u650[i,j,k])**(2.)+(v550[i,j,k]-v650[i,j,k])**(2.))**(1./2.))/(gph550[i,j,k]-gph650[i,j,k]))
			DEF[i,j,k] = ((dudx[i,j,k]-dvdy[i,j,k])**(2.)+(dvdx[i,j,k]+dudy[i,j,k])**(2.))**(1./2.)
			TI1[i,j,k] = (VWS[i,j,k])*(DEF[i,j,k])

TI1 = TI1*(10000000)
TI1 = np.around(TI1, decimals = 0)

n = 2
titles= ['(a) 20181028 12UTC', '(b) 20181029 12UTC', '(c) 20181030 12UTC']

fig = plt.figure(figsize=(15,10))
for i in range(3):
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
        x, y = m(lon_2d1, lat_2d1)
        x1, y1 = m(lon[562], lat[215]) #CAT incounter
	clevs1= np.arange(0,13,1)
	TI1[ TI1 > max(clevs1) ] = max(clevs1)
	m.plot(x1, y1, color='r', marker='*')
	cmap = mpl.cm.Blues
	cs2=m.contourf(x,y,TI1[n,:,:],clevs1, cmap = cmap)
	cbar=m.colorbar(cs2)
        cbar.set_label('Turbulence Index 1 $(10^{-7}/s^{-2})$', fontsize = 7)
	n = n+4

plt.subplots_adjust(wspace = 0.42, hspace = 0.3)
fig = plt.savefig('./figures/Paper_TI1')
