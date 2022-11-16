from netCDF4 import Dataset
import numpy as np
from mpl_toolkits.basemap import Basemap, shiftgrid
import matplotlib.pyplot as plt
from sys import exit


lon = Dataset('data/2018-10-30_06:00:00.nc').variables['g0_lon_2'][:]
lat = Dataset('data/2018-10-30_06:00:00.nc').variables['g0_lat_1'][:]
gph500 = np.zeros((11,721,1440))
gph1000 = np.zeros((11,721,1440))
th = np.zeros((11,721,1440))

k = 0
dd = ['28','29','30']
tt = ['00','06','12','18']

for i  in range(3):
        for j in range(4):
                        if dd[i] == '30' and tt[j] == '18':
                                break
                        else :
                                gph500[k,:,:] = Dataset('data/2018-10-'+dd[i]+'_'+tt[j]+':00:00.nc').variables['Z_GDS0_ISBL'][16,:,:]
				gph1000[k,:,:] = Dataset('data/2018-10-'+dd[i]+'_'+tt[j]+':00:00.nc').variables['Z_GDS0_ISBL'][31,:,:]

                                k = k+1


gph500 = gph500/9.806
gph1000 = gph1000/9.806
th = gph500-gph1000

titles= ['(a) 20181028 00UTC GPH','(b) 20181028 06UTC GPH', '(c) 20181028 12UTC GPH', '(d) 20181028 18UTC GPH', '(e) 20181029 00UTC GPH', '(f) 20181029 06UTC GPH', '(g) 20181029 12UTC GPH', '(h) 20181029 18UTC GPH', '(i) 20181030 00UTC GPH', '(j) 20181030 06UTC GPH', '(k) 20181030 12UTC GPH']

fig = plt.figure(figsize=(15,10))
for i in range(11):
        plt.subplot(3, 4, i+1)
        m = Basemap(llcrnrlat=15,urcrnrlat=65,\
        llcrnrlon=100,urcrnrlon=160, resolution='i')#,resolution='c')
        plt.title(titles[i], fontsize=12, loc='left')
        m.drawcoastlines()
        m.drawcountries()
        m.drawparallels(np.arange(15,65,10),labels=[1,0,0,0])
        m.drawmeridians(np.arange(100,160,10),labels=[0,0,0,1])
        m.fillcontinents(alpha=0)
        m.drawmapboundary(fill_color='w')
        lon_2d1, lat_2d1 = np.meshgrid(lon, lat)
        x, y = m(lon_2d1, lat_2d1)
        x1, y1 = m(lon[562], lat[215]) #CAT incounter
	clevs1= np.arange(int(np.min(gph500[i,100:301,400:641])),int(np.max(gph500[i,100:301,400:641]))+60.,60.)
        clevs2= np.arange(int(np.min(th[i,100:301,400:641])),int(np.max(th[i,100:301,400:641]))+60.,60.)
	m.plot(x1, y1, color='r', marker='*')
	cs=m.contour(x,y,gph500[i,:,:],clevs1,colors='red', linewidths = 0.5, linestyles='solid')
	cs2=m.contour(x,y,th[i,:,:],clevs2,colors='blue', linewidths = 0.5,linestyles='dashed')
	plt.clabel(cs, manual = True,  colors='red')
	plt.clabel(cs2, manual = True,  colors='blue')


plt.show()
plt.subplots_adjust(wspace = 0.1, hspace = 0.3)

