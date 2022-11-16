from netCDF4 import Dataset
import numpy as np
#from mpl_toolkits.basemap import Basemap
import matplotlib
import matplotlib.pyplot as plt
from sys import exit

file = './Narita_time.txt'
f= open(file)
lines=f.readlines()

time = []

for line in lines:
        time.append(line.strip())

file = './Narita_altitude.txt'
g = open(file)
lines=g.readlines()

alt = []

for line in lines:
        alt.append(float(line.strip()))


file = './Narita_verticalacceleration.txt'
h= open(file)
lines=h.readlines()
va = []

for line in lines:
        va.append(float(line.strip()))

time = np.array(time)
alt = np.array(alt)
va = np.array(va)

fig, ax1 = plt.subplots(figsize=(10,7))
plt.title('Altitude and Vertical acceleration [2018/10/30] ', weight='bold', size=17)
ax1.plot(time, alt, 'b-')
ax1.set_xlabel('$Time$', weight='bold', size=14)
ax1.set_ylabel('$Altitude($ft$)$', color='b', size=14)
ax1.tick_params('y', colors='b')

ax2 = ax1.twinx()
ax2.plot(time, va, 'r-')
ax2.set_ylabel('Vertical acceleration', color='r', size = 14)
ax2.tick_params('y', colors='r')
plt.annotate('CAT encounter at 10:25:14', xy=('10:25:14' , 0.004))

fig.tight_layout()
plt.xticks(['10:20:00', '10:21:00', '10:22:00', '10:23:00', '10:24:00' , '10:25:00','10:26:00', '10:27:00', '10:28:00', '10:29:00', '10:30:00'])
plt.gcf().autofmt_xdate()
fig = plt.savefig('./figures/altitude and vertical acceleration')
