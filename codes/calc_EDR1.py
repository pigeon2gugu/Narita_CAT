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

file = './Narita_time.txt'
h= open(file)
lines=h.readlines()
time = []
for line in lines:
        time.append(line.strip())

lon = np.array(lon)
lat = np.array(lat)
DIR = np.array(DIR)
SPD = np.array(SPD)
time = np.array(time)
u = np.zeros((len(DIR)))
v = np.zeros((len(DIR)))
ubar = []
vbar = []
uu = np.zeros((len(DIR))-1)
vv = np.zeros((len(DIR))-1)
TKE = np.zeros(len(DIR))
EDR = np.zeros(len(DIR))

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
	EDR[i] = ((1./311.)* (TKE[i])**(3./2.))**(1./3.)

TKE[600] = ((u[600]-ubar[-1])**(2.)+(v[600]-vbar[-1])**(2.))/2.
EDR[600] = (0.84/300. * (TKE[600])**(3./2.))**(1./3.)

plt.figure(figsize=(10,7))
plt.plot(time, TKE, color='tab:blue', linewidth=1.5, label='TKE')
plt.title('TKE Time Series [2018/10/30] ', weight='bold', size=17)
plt.ylabel('$TKE (m^2/s^2)$', color='tab:blue', weight='bold', size=14)
plt.xlabel('$Time$', weight='bold', size=14)
plt.xticks(['10:20:00', '10:21:00', '10:22:00', '10:23:00', '10:24:00' , '10:25:00','10:26:00', '10:27:00', '10:28:00', '10:29:00', '10:30:00'])
plt.gcf().autofmt_xdate()
fig = plt.savefig('./figures/TKE_timeseries')


	
	
