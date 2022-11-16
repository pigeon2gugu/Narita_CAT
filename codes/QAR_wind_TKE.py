from netCDF4 import Dataset
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib
import matplotlib.pyplot as plt
from sys import exit

file = './Narita_time.txt'
f= open(file)
lines=f.readlines()

time = []

for line in lines:
        time.append(line.strip())

file = './Narita_wind.txt'
h= open(file)
lines=h.readlines()
DIR = []
SPD = []
for line in lines:
        DIR.append(float(line.split( )[0]))
        SPD.append(float(line.split( )[1]))

DIR = np.array(DIR)
SPD = np.array(SPD)
u = np.zeros((len(DIR)))
v = np.zeros((len(DIR)))
ubar = []
vbar = []
uu = np.zeros((len(DIR))-1)
vv = np.zeros((len(DIR))-1)
TKE = np.zeros(len(u))

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


for i in range(len(uu)):
        TKE[i] = (uu[i]**(2.)+vv[i]**(2.))/2.

TKE[600] = ((u[600]-ubar[-1])**(2.)+(v[600]-vbar[-1])**(2.))/2.

time = np.array(time)


fig, ax1 = plt.subplots()
ax1.plot(time,TKE,'k-', label='TKE')
ax1.set_xlabel('time')
ax1.set_ylabel('TKE', color='k')
ax1.tick_params('y', colors='k') 

ax2 = ax1.twinx()
ax2.spines["right"].set_position(("axes", 1.2))
ax2.plot(time, v, 'r-', label='v')
ax2.tick_params('y', colors='r')
ax2.set_ylabel('v', color='r')

ax3 = ax1.twinx()
ax3.plot(time, u, 'b-',  label='u')
ax3.set_ylabel('u', color='b')
ax3.tick_params('y', colors='b')

fig.tight_layout()
plt.xticks(['10:20:00', '10:21:00', '10:22:00', '10:23:00', '10:24:00' , '10:25:00','10:26:00', '10:27:00', '10:28:00', '10:29:00', '10:30:00'])
plt.gcf().autofmt_xdate()
fig = plt.savefig('./figures/wind and TKE')
