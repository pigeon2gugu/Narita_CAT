from netCDF4 import Dataset
import numpy as np
#from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from sys import exit

file = './Narita_lon.txt'
f= open(file)
lines=f.readlines()

lon = []

for line in lines:
        lon.append(float(line.strip()))

file = './Narita_lat.txt'
g = open(file)
lines=g.readlines()

lat = []

for line in lines:
        lat.append(float(line.strip()))


file = './Narita_wind.txt'
h= open(file)
lines=h.readlines()
DIR = []
SPD = []
for line in lines:
	DIR.append(float(line.split( )[0]))
	SPD.append(float(line.split( )[1]))

lon = np.array(lon)
lat = np.array(lat)
DIR = np.array(DIR)
SPD = np.array(SPD)
u = np.zeros((len(DIR)))
v = np.zeros((len(DIR)))
ubar = []
vbar = []
uu = np.zeros((len(DIR))-1)
vv = np.zeros((len(DIR))-1)
TKE = np.zeros(len(uu))

for i in range(len(DIR)):
	u[i] = -0.5144*(SPD[i])*(np.sin((DIR[i])*(np.pi)/180.))
	v[i] = -0.5144*(SPD[i])*(np.cos((DIR[i])*(np.pi)/180.))


CAT = np.zeros((1,2))
CAT[0,0] = 36.233566
CAT[0,1] = 140.5771 #lon=314

n=0
for i in range((len(u)/10)):
	ubar.append((u[n]+u[n+1]+u[n+2]+u[n+3]+u[n+4]+u[n+5]+u[n+6]+u[n+7]+u[n+8]+u[n+9])/10)
	vbar.append((v[n]+v[n+1]+v[n+2]+v[n+3]+v[n+4]+v[n+5]+v[n+6]+v[n+7]+v[n+8]+v[n+9])/10)
	n=n+10

ubar = np.array(ubar)
vbar = np.array(vbar)


n=0
k=0
for i in range(len(u)-1):
	if k <= 9:
		uu[i] = (u[i]-ubar[n])
		vv[i] = (v[i]-vbar[n])
		k = k+1
	if k > 9:
		k = 0
		n = n+1

		
for i in range(len(uu)):
	TKE[i] = (uu[i]**(2.)+vv[i]**(2.))/2.

exit()

plt.figure(figsize=(10,7))
m = Basemap(projection='merc',llcrnrlat=30,urcrnrlat=40,\
                        llcrnrlon=135,urcrnrlon=145,resolution='i')
m.drawcoastlines()
m.drawcountries()
m.drawparallels(np.arange(30,40,2),labels=[1,0,0,0])
m.drawmeridians(np.arange(135,145,2),labels=[0,0,0,1])
m.fillcontinents(alpha=0)
m.drawmapboundary(fill_color='w')
x, y = m(lon, lat)
for i in range(len(lon)):
        m.plot(x[i], y[i], marker='.', color='b', markersize=2)
m.plot(x[314], y[314], color='r', marker='*')
fig = plt.savefig('./figures/20181030NaritaCAT')


	
	
