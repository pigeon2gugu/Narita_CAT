from netCDF4 import Dataset
import numpy as np
#from mpl_toolkits.basemap import Basemap
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

for i in range(len(DIR)):
        u[i] = -0.5144*(SPD[i])*(np.sin((DIR[i])*(np.pi)/180.))
        v[i] = -0.5144*(SPD[i])*(np.cos((DIR[i])*(np.pi)/180.))

freq, psd = signal.welch(u, fs = 1./60.)

tlen = 60
seglen=30

ti=int(len(u)//tlen)
fi=np.zeros((ti, seglen//2+1),float)
wi=np.zeros((ti, seglen//2+1),float)

tt=np.arange(0.022, 0.21, 0.01)
def exp(x,y,c):
    z = c*x**y
    return z

c=0.53
vv=2

U = []

for i in range(0,ti):
    U.append((np.sum(u[int((i)*tlen):int((i+1)*tlen)]))/tlen)

wav=[]
EDR2_spd, meany = [], []
for i in range(0, ti) :
    freq, fi[i][:], = signal.welch(u[int((i)*tlen):int((i+1)*tlen)], fs=1., window='hamming', nperseg=seglen, noverlap=None)

    wav = freq*(2*np.pi/U[i])
    wi[i][:] = fi[i][:]*(U[i]/(2*np.pi))

    logx=[]
    logy=[]
    for j in range (0, len(wav)) :
        logx.append(np.log(wav[j]))
        logy.append(np.log(wi[i][j]))        
    bb=[]
    for k in range(len(logy)-vv):
        bb.append(logy[k+vv]+5./3.*logx[k+vv])
    meany.append(np.mean(bb))

    EDR2_spd.append((np.exp(meany[i])/c)**0.5)


print(EDR2_spd)

'''
plt.figure(figsize=(9,7.5))
for i in range(0, ti) :
    if 0.4 <= float(EDR2_spd[i]) < 0.7 :
        plt.loglog(wav, wi[i,:], c='tab:blue', linewidth=2, alpha=0.75)
        plt.loglog(wav[2:16], exp(wav[2:16], -5./3., np.exp(meany[i])), linewidth=1.7, c='tab:red', alpha=0.8)
    else:
        plt.loglog(wav, wi[i][:], c='black', alpha=0.2)
        plt.loglog(wav[2:16], exp(wav[2:16], -5./3., np.exp(meany[i])), c='orange', linewidth=1.7, alpha=0.25)
plt.title("Wind-SPD PSD [2018/10/28 0542 UTC]", weight='bold', loc='center', size=17)
plt.xlabel('Wavenumber ($rad*m^{-1}$)', size=14)
plt.ylabel('$Energy(m^{3}/s^{2})$', size=14)
#plt.xlim(10**-3,10**0)
#plt.ylim(10**-3,10**3)
'''
EDR2_spd = np.array(EDR2_spd)

time = ['10:20:30', '10:21:30', '10:22:30', '10:23:30', '10:24:30' , '10:25:30','10:26:30', '10:27:30', '10:28:30', '10:29:30']
plt.figure(figsize=(10,7))
plt.plot(time, EDR2_spd, color='tab:blue', linewidth=1.5, label='TKE')
plt.title('EDR2_U Time Series [2018/10/30 QAR] ', weight='bold', size=18)
plt.ylabel('EDR2 $(m^{2/3}/s)$', color='tab:blue', weight='bold', size=14)
plt.xlabel('Time', weight='bold', size=14)
#plt.xticks(['10:20:00', '10:21:00', '10:22:00', '10:23:00', '10:24:00' , '10:25:00','10:26:00', '10:27:00', '10:28:00', '10:29:00'])
plt.gcf().autofmt_xdate()
plt.savefig('./figures/EDR_U_Timeseries')

