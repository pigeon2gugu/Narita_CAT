from netCDF4 import Dataset
import numpy as np
from mpl_toolkits.basemap import Basemap, shiftgrid
import matplotlib.pyplot as plt
from sys import exit
import matplotlib as mpl

lon = Dataset('data/2018-10-30_06:00:00.nc').variables['g0_lon_2'][:]
lat = Dataset('data/2018-10-30_06:00:00.nc').variables['g0_lat_1'][:]
u = np.zeros((11,721,1440))
v = np.zeros((11,721,1440))

k = 0
dd = ['28','29','30']
tt = ['00','06','12','18']
for i  in range(3):
        for j in range(4):
                        if dd[i] == '30' and tt[j] == '18':
                                break
                        else :
                                u[k,:,:] = Dataset('data/2018-10-'+dd[i]+'_'+tt[j]+':00:00.nc').variables['U_GDS0_ISBL'][12,:,:]
                                v[k,:,:] = Dataset('data/2018-10-'+dd[i]+'_'+tt[j]+':00:00.nc').variables['V_GDS0_ISBL'][12,:,:]

                                k = k+1

wind = np.zeros((11,721,1440))

for i in range(11):
        for j in range(len(lat)):
		for k in range(len(lon)):
                	wind[i,j,k] = (((u[i,j,k])**(2.))+((v[i,j,k])**(2.)))**(1./2.)

lat1 = lat[::-1]

su, newlon1 = shiftgrid(180., u, lon, start=False)
sv, newlon1 = shiftgrid(180., v, lon, start=False)
su = su[:,::-1,:]
sv = sv[:,::-1,:]

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
	x2, y2 = m(lon[562], lat[215]) #CAT incounter
	clevs1= np.arange(0, 90+5,5)
	wind[ wind > max(clevs1) ] = max(clevs1)
	cmap = mpl.cm.coolwarm
	m.plot(x2, y2, color='r', marker='*')
	cs1=m.contourf(x,y,wind[n,:,:],clevs1, cmap = cmap)
	uproj,vproj,xx,yy = m.transform_vector(su[n,:,:],sv[n,:,:],newlon1,lat1,10,10,returnxy=True,masked=True)
	x1 , y1 = m(xx, yy)
	barbs = m.barbs(x1, y1, uproj, vproj, length=6, pivot='middle', linewidth=0.5)
	cbar=m.colorbar(cs1)#, cax=cbaxes)
	cbar.ax.set_yticklabels(['0','10','20','30','40','50','60','70','80','>90'])
	cbar.set_label('Wind speed $(m/s)$', fontsize = 7)
	n = n+4

plt.subplots_adjust(wspace = 0.42, hspace = 0.30)
fig = plt.savefig('./figures/Paper_windspeed')
