from netCDF4 import Dataset
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from sys import exit


data1=Dataset('/home/data/ERA5/20181030/2018-10-30_10:00:00.nc')
data2=Dataset('/home/data/ERA5/20181030/2018-10-30_11:00:00.nc')

lon1 = data1.variables['g0_lon_2'][:]
lat1 = data1.variables['g0_lat_1'][:]
ggp1 = data1.variables['Z_GDS0_ISBL'][:,:,:]
uu1 = data1.variables['U_GDS0_ISBL'][:,:,:]
vv1 = data1.variables['V_GDS0_ISBL'][:,:,:]
ww1 = data1.variables['W_GDS0_ISBL'][:,:,:]
t1 = data1.variables['T_GDS0_ISBL'][:,:,:]

lon2 = data2.variables['g0_lon_2'][:]
lat2 = data2.variables['g0_lat_1'][:]
ggp2 = data2.variables['Z_GDS0_ISBL'][:,:,:]
uu2 = data2.variables['U_GDS0_ISBL'][:,:,:]
vv2 = data2.variables['V_GDS0_ISBL'][:,:,:]
ww2 = data2.variables['W_GDS0_ISBL'][:,:,:]
t2 = data2.variables['T_GDS0_ISBL'][:,:,:]



dx = np.zeros((721,1440))
dy = np.zeros((721,1440))

for i in range(len(lat1)):
        for j in range(1,len(lon1)-1):
                dx[i,j] = 2*np.pi*6371000*np.cos(lat1[i]*2*np.pi/360.)*(lon1[j+1]-lon1[j-1])/360.

for i in range(len(lat1)):
        dx[i,-1] = dx[i,-2]
        dx[i,0] = dx[i,1]

for i in range(1,len(lat1)-1):
        for j in range(len(lon1)):
                dy[i,j] = 2*np.pi*6371000*(lat1[i+1]-lat1[i-1])/360.


for j in range(len(lon1)):
        dy[-1,j] = 2*np.pi*6371000*(0.5)/360.
        dy[0,j] = 2*np.pi*6371000*(0.5)/360.


gh1 = ggp1/9.806
gh2 = ggp2/9.806

dudx = np.zeros((721,1440))
dudx2 = np.zeros((721,1440))

dvdx = np.zeros((721,1440))
dvdx2 = np.zeros((721,1440))

dudy = np.zeros((721,1440))
dudy2 = np.zeros((721,1440))

dvdy = np.zeros((721,1440))
dvdy2 = np.zeros((721,1440))

for j in range(len(lat1)):
        for k in range(1,len(lon1)-1):
                dudx[j,k] = (uu1[16,j,k+1]-uu1[16,j,k-1])/(dx[j,k])
                dudx2[j,k] = (uu2[16,j,k+1]-uu2[16,j,k-1])/(dx[j,k])
                dvdx[j,k] = (vv1[16,j,k+1]-vv1[16,j,k-1])/(dx[j,k])
                dvdx2[j,k] = (vv2[16,j,k+1]-vv2[16,j,k-1])/(dx[j,k])


for j in range(len(lat1)):
        dudx[j,-1] = (uu1[16,j,0]-uu1[16,j,-2])/(dx[j,-1])
        dudx2[j,-1] = (uu2[16,j,0]-uu2[16,j,-2])/(dx[j,-1])
        dvdx[j,-1] = (vv1[16,j,0]-vv1[16,j,-2])/(dx[j,-1])
        dvdx2[j,-1] = (vv2[16,j,0]-vv2[16,j,-2])/(dx[j,-1])
        dudx[j,0] = (uu1[16,j,-1]-uu1[16,j,1])/(dx[j,0])
        dudx2[j,0] = (uu2[16,j,-1]-uu2[16,j,1])/(dx[j,0])
        dvdx[j,0] = (vv1[16,j,-1]-vv1[16,j,1])/(dx[j,0])
        dvdx2[j,0] = (vv2[16,j,-1]-vv2[16,j,1])/(dx[j,0])

for j in range(1,len(lat1)-1):
        for k in range(len(lon1)):
                dudy[j,k] = (uu1[16,j+1,k]-uu1[16,j-1,k])/(dy[j,k])
                dudy2[j,k] = (uu2[16,j+1,k]-uu2[16,j-1,k])/(dy[j,k])
                dvdy[j,k] = (vv1[16,j+1,k]-vv1[16,j-1,k])/(dy[j,k])
                dvdy2[j,k] = (vv2[16,j+1,k]-vv2[16,j-1,k])/(dy[j,k])

for j in range(len(lon1)):
        dudy[-1,j] = dudy[-2,j]
        dudy2[-1,j] = dudy2[-2,j]
        dvdy[-1,j] = dvdy[-2,j]
        dvdy2[-1,j] = dvdy2[-2,j]
        dudy[0,j] = dudy[1,j]
        dudy2[0,j] = dudy2[1,j]
        dvdy[0,j] = dvdy[1,j]
        dvdy2[0,j] = dvdy2[1,j]

VWS = np.zeros((721,1440))
VWS2 = np.zeros((721,1440))
CVG = np.zeros((721,1440))
CVG2 = np.zeros((721,1440))
DEF = np.zeros((721,1440))
DEF2 = np.zeros((721,1440))
TI1 = np.zeros((721,1440))
TI12 = np.zeros((721,1440))
TI2 = np.zeros((721,1440))
TI22 = np.zeros((721,1440))
PT550 = np.zeros((721,1440))
PT450 = np.zeros((721,1440))
PT500 = np.zeros((721,1440))
Ns = np.zeros((721,1440))
PPT550 = np.zeros((721,1440))
PPT450 = np.zeros((721,1440))
PPT500 = np.zeros((721,1440))
Ns2 = np.zeros((721,1440))
Ri = np.zeros((721,1440))
Ri2 = np.zeros((721,1440))
vo = np.zeros((721,1440))
vo2 = np.zeros((721,1440))


for i in range(len(lat1)):
        for j in range(len(lon1)):
                VWS[i,j] = abs((((uu1[15,i,j]-uu1[17,i,j])**(2.)+(vv1[15,i,j]-vv1[17,i,j])**(2.))**(1./2.))/(gh1[15,i,j]-gh1[17,i,j]))
                VWS2[i,j] = abs((((uu2[15,i,j]-uu2[17,i,j])**(2.)+(vv2[15,i,j]-vv2[17,i,j])**(2.))**(1./2.))/(gh2[15,i,j]-gh2[17,i,j]))
                CVG[i,j] = -(dudx[i,j]+dvdy[i,j])
                CVG2[i,j] = -(dudx2[i,j]+dvdy2[i,j])
                DEF[i,j] = ((dudx[i,j]-dvdy[i,j])**(2.)+(dvdx[i,j]+dudy[i,j])**(2.))**(1./2.)
                DEF2[i,j] = ((dudx2[i,j]-dvdy2[i,j])**(2.)+(dvdx2[i,j]+dudy2[i,j])**(2.))**(1./2.)
                TI1[i,j] = (VWS[i,j])*(DEF[i,j])
                TI12[i,j] = (VWS2[i,j])*(DEF2[i,j])
                TI2[i,j] = (VWS[i,j])*(DEF[i,j]+CVG[i,j])
                TI22[i,j] = (VWS2[i,j])*(DEF2[i,j]+CVG2[i,j])
		PT550[i,j] = (t1[17,i,j])*((1013./550.)**(2./7.))
		PT450[i,j] = (t1[15,i,j])*((1013./450.)**(2./7.))
		PT500[i,j] = (t1[16,i,j])*((1013./500.)**(2./7.))
 		Ns[i,j] = ((9.806)*(PT450[i,j]-PT550[i,j]))/((PT500[i,j])*(gh1[15,i,j]-gh1[17,i,j]))	
		PPT550[i,j] = (t2[17,i,j])*((1013./550.)**(0.286))
		PPT450[i,j] = (t2[15,i,j])*((1013./450.)**(0.286))
		PPT500[i,j] = (t2[16,i,j])*((1013./500.)**(0.286))
		Ns2[i,j] = ((9.806)*(PPT450[i,j]-PPT550[i,j]))/((PPT500[i,j])*(gh2[15,i,j]-gh2[17,i,j]))
		Ri[i,j] = (Ns[i,j])/(((uu1[15,i,j]-uu1[17,i,j])/(gh1[15,i,j]-gh1[17,i,j]))**(2.))
		Ri2[i,j] = (Ns2[i,j])/(((uu2[15,i,j]-uu2[17,i,j])/(gh2[15,i,j]-gh2[17,i,j]))**(2.))
		vo[i,j] = (dvdx[i,j]-dudy[i,j])+((2.)*(2.*np.pi/86400.)*(np.sin(lat1[i]*np.pi/180.)))
		vo2[i,j] = (dvdx2[i,j]-dudy2[i,j])+((2.)*(2*np.pi/86400.)*(np.sin(lat2[i]*np.pi/180.)))

