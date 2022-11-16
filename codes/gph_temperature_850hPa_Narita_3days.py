from netCDF4 import Dataset
import numpy as np
from mpl_toolkits.basemap import Basemap, shiftgrid
import matplotlib.pyplot as plt
from sys import exit


lon = Dataset('data/2018-10-30_06:00:00.nc').variables['g0_lon_2'][:]
lat = Dataset('data/2018-10-30_06:00:00.nc').variables['g0_lat_1'][:]
gph850 = np.zeros((11,721,1440))
t = np.zeros((11,721,1440))


k = 0
dd = ['28','29','30']
tt = ['00','06','12','18']

# 850 hpa = 25, 700 hpa = 20

for i  in range(3):
        for j in range(4):
                        if dd[i] == '30' and tt[j] == '18':
                                break
                        else :
                                gph850[k,:,:] = Dataset('data/2018-10-'+dd[i]+'_'+tt[j]+':00:00.nc').variables['Z_GDS0_ISBL'][25,:,:]

				t[k,:,:] = Dataset('data/2018-10-'+dd[i]+'_'+tt[j]+':00:00.nc').variables['T_GDS0_ISBL'][25,:,:]
				k = k+1


gph850 = gph850/9.806


titles= ['(a) 20181028 00UTC GPH & T','(b) 20181028 06UTC GPH & T', '(c) 20181028 12UTC GPH & T', '(d) 20181028 18UTC GPH & T', '(e) 20181029 00UTC GPH & T', '(f) 20181029 06UTC GPH & T', '(g) 20181029 12UTC GPH & T', '(h) 20181029 18UTC GPH & T', '(i) 20181030 00UTC GPH & T', '(j) 20181030 06UTC GPH & T', '(k) 20181030 12UTC GPH & T']

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
	clevs1= np.arange(int(round(np.min(t[i,100:301,400:641]))),int(round(np.max(t[i,100:301,400:641])))+10.,10.)
        clevs2= np.arange(int(round(np.min(gph850[i,100:301,400:641]))),int(round(np.max(gph850[i,100:301,400:641])))+60.,60.)
	m.plot(x1, y1, color='r', marker='*')
	cs=m.contour(x,y,t[i,:,:],clevs1,colors='red', linewidths = 0.6, linestyles='dashed')
	cs2=m.contour(x,y,gph850[i,:,:],clevs2,colors='blue', linewidths = 0.6,linestyles='solid')
	plt.clabel(cs, colors='red')
	plt.clabel(cs2, colors='blue')

plt.subplots_adjust(wspace = 0.1, hspace = 0.3)

fig = plt.savefig('./figures/GPH_T_850hPa_Narita_201828_29_30')
