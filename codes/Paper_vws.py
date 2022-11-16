from netCDF4 import Dataset
import numpy as np
from mpl_toolkits.basemap import Basemap, shiftgrid
import matplotlib.pyplot as plt
from sys import exit
import matplotlib as mpl


lon = Dataset('data/2018-10-30_06:00:00.nc').variables['g0_lon_2'][:]
lat = Dataset('data/2018-10-30_06:00:00.nc').variables['g0_lat_1'][:]
gph550 = np.zeros((11,721,1440))
gph650 = np.zeros((11,721,1440))
u550 = np.zeros((11,721,1440))
u650 = np.zeros((11,721,1440))
v550 = np.zeros((11,721,1440))
v650 = np.zeros((11,721,1440))

k = 0
dd = ['28','29','30']
tt = ['00','06','12','18']

for i  in range(3):
        for j in range(4):
                        if dd[i] == '30' and tt[j] == '18':
                                break
                        else :
                                gph550[k,:,:] = Dataset('data/2018-10-'+dd[i]+'_'+tt[j]+':00:00.nc').variables['Z_GDS0_ISBL'][17,:,:]
				gph650[k,:,:] = Dataset('data/2018-10-'+dd[i]+'_'+tt[j]+':00:00.nc').variables['Z_GDS0_ISBL'][19,:,:]
				u550[k,:,:] = Dataset('data/2018-10-'+dd[i]+'_'+tt[j]+':00:00.nc').variables['U_GDS0_ISBL'][17,:,:]
				u650[k,:,:] = Dataset('data/2018-10-'+dd[i]+'_'+tt[j]+':00:00.nc').variables['U_GDS0_ISBL'][19,:,:]
				v550[k,:,:] = Dataset('data/2018-10-'+dd[i]+'_'+tt[j]+':00:00.nc').variables['V_GDS0_ISBL'][17,:,:]
				v650[k,:,:] = Dataset('data/2018-10-'+dd[i]+'_'+tt[j]+':00:00.nc').variables['V_GDS0_ISBL'][19,:,:]
                                k = k+1


gph550 = gph550/9.806
gph650 = gph650/9.806

VWS = np.zeros((11,721,1440))

for i in range(11):
	for j in range(len(lat)):
        	for k in range(len(lon)):
                	VWS[i,j,k] = abs((((u550[i,j,k]-u650[i,j,k])**(2.)+(v550[i,j,k]-v650[i,j,k])**(2.))**(1./2.))/(gph550[i,j,k]-gph650[i,j,k]))

VWS = np.around(VWS, decimals = 4)
VWS = VWS*(1000)

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
        clevs2= np.arange(0,22,2)
	VWS[ VWS > max(clevs2) ] = max(clevs2)
	m.plot(x1, y1, color='r', marker='*')
	cmap = mpl.cm.Blues
	cs=m.contourf(x,y,VWS[n,:,:],clevs2,cmap=cmap)
	cbar=m.colorbar(cs)
	cbar.ax.set_yticklabels(['0','4','8','12','16','>20'])
	cbar.set_label('vertical windshear $(10^{3}/s)$', fontsize = 7)
	n = n+4

plt.subplots_adjust(wspace = 0.42, hspace = 0.3)
fig = plt.savefig('./figures/Paper_VWS')

