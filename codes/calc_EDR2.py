from netCDF4 import Dataset
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from sys import exit
from scipy import signal

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

psdu = np.zeros((10,16))
psdv = np.zeros((10,16))

for i in range(0,10):
	freqsu, psdu[i][:] = signal.welch(u[60*i:60+60*i], fs=1., window='hamming', nperseg=30)
	freqsv, psdv[i][:] = signal.welch(v[60*i:60+60*i], fs=1., window='hamming', nperseg=30)

logx = np.zeros((10,14))
logyu = np.zeros((10,14))
logyv = np.zeros((10,14))
du = np.zeros((10,14))
dv = np.zeros((10,14))
meanu = np.zeros((10))
meanv = np.zeros((10))
EDR2u = np.zeros((10))
EDR2v = np.zeros((10))

for i in range(10):
	for j in range(14):
		logx[i][j] = np.log(freqsu[j+2])
		logyu[i][j] = np.log(psdu[i][j+2])
		logyv[i][j] = np.log(psdv[i][j+2])
		du[i][j] = (logyu[i][j]+5./3.*logx[i][j])
		dv[i][j] = (logyv[i][j]+5./3.*logx[i][j])
	for k in range(13):
		du[i][k+1] = du[i][k]+du[i][k+1]
		dv[i][k+1] = dv[i][k]+dv[i][k+1]
	meanu[i] = (du[i][13])/14.
	meanv[i] = (dv[i][13])/14.
	EDR2u[i] = (np.exp(meanu[i])/0.53)**(1./2.)
	EDR2v[i] = (np.exp(meanv[i])/0.53)**(1./2.)

plt.figure(figsize=(10,7))
for i in range(10):	
	plt.loglog(freqsu,psdu[i,:], color = 'black', linewidth=2, alpha=0.5)
	plt.loglog(freqsu[2:16],np.exp(meanu[i])*freqsu[2:16]**(-5./3.),linewidth=1.7, color = 'red')
plt.title("U wind PSD [2018/10/30 QAR]", weight='bold', loc='center', size=18)
plt.xlabel('Frequency ($s^{-1}$)', size=14)
plt.ylabel('Energy', size=14)
plt.savefig('./figures/PSD_EDR_U')

'''

time = ['10:20:30', '10:21:30', '10:22:30', '10:23:30', '10:24:30' , '10:25:30','10:26:30', '10:27:30', '10:28:30', '10:29:30']
plt.figure(figsize=(10,7))
plt.plot(time, EDR2u, color='tab:blue', linewidth=1.5, label='TKE')
plt.title('EDR2_U Time Series [2018/10/30 QAR] ', weight='bold', size=18)
plt.ylabel('EDR2 $(m^{2/3}/s)$', color='tab:blue', weight='bold', size=14)
plt.xlabel('Time', weight='bold', size=14)
#plt.xticks(['10:20:00', '10:21:00', '10:22:00', '10:23:00', '10:24:00' , '10:25:00','10:26:00', '10:27:00', '10:28:00', '10:29:00'])
plt.gcf().autofmt_xdate()
plt.savefig('./figures/EDR_U_Timeseries')

'''
