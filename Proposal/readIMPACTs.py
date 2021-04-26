from netCDF4 import *
import matplotlib.pyplot as plt
import numpy as np
fname='/media/grecu/ExtraDrive1/IMPACTS/IMPACTS_HIWRAP_L1B_RevA_20200227T074338_to_20200227T142238.h5'

fh=Dataset(fname)

zKu=fh['Products/Ku/Combined/Data/dBZe']
zKa=fh['Products/Ka/Combined/Data/dBZe']
lat=fh['Navigation/Data/Latitude']
lon=fh['Navigation/Data/Longitude']
alt=fh['Navigation/Data/Height']
rrange=fh['/Products/Information/Range']
import datetime
dt=(datetime.datetime(2020,2,27)-datetime.datetime(1970,1,1)).days*24*3600
time=(fh['Time/Data/TimeUTC'][:]-dt)/3600.
a=np.nonzero((time-9.85)*(time-10.25)<0)
zKa_slice=zKa[a[0],:]
zKu_slice=zKu[a[0],:]
dZ=zKu_slice[:,:]-zKa_slice[:,:]
hm=alt[a[0]].mean()
nr=zKa_slice.shape[1]
nx=zKa_slice.shape[0]
rmean=rrange[:]
import scattering as sc
coeffs=[0.057072960160875125, 0.4285461505798106, -3.777856630433933]
nw=6.5+(hm-rmean[0,:])/1000.*2.5/10.0
nw[nw>9]=9
nw[nw<6.5]=6.5
swcL=[]
hsurf=[]
zsurf=[]

from numba import jit

@jit(nopython=True)
def bisection(xr,xvec):
    n=xvec.shape[0]
    if xr<=xvec[0]:
        i1,i2=0,0
        return i1,i2
    if xr>=xvec[n-1]:
        i1,i2=n-1,n-1
        return i1,i2
    i1=0
    x1=xvec[i1]
    i2=n-1
    x2=xvec[i2]
    while i2-i1>1:
        imid=int((i1+i2)/2)
        if xvec[imid]>xr:
            i2=imid
        else:
            i1=imid
    return i1,i2
xvec=np.arange(10)
i1,i2=0,9
i1,i2=bisection(4.1,xvec)
from zprofTest import *
kuka_Table['dFR']=np.array(kuka_Table['zKu'])-np.array(kuka_Table['zKa'])

dFR_Table=kuka_Table['dFR']
dm_Table=np.array(kuka_Table['dm'])
dFR=np.array(dFR)
dmRet=dFR.copy()*0
@jit(nopython=True)
def fromdfr(dFR,dmRet,dFR_table,dm_Table):
    n=dFR.shape[0]
    for i in range(n):
        i1,i2=bisection(dFR[i],dFR_table)
        dmRet[i]=dm_Table[i1]
        
fromdfr(dFR,dmRet,dFR_Table,dm_Table)
stop
for i in range(nx):
    ind=np.argmax(zKu_slice[i,500:-30])
    swc1=np.zeros((nr),float)-99
    zKa_slice[i,500+ind-20:]=np.log(0)
    a1=np.nonzero(zKa_slice[i,:]==zKa_slice[i,:])
    swc1[a1]=10**(zKa_slice[i,a1[0]]*coeffs[0]+nw[a1]*coeffs[1]+coeffs[2])
    swc1[500+ind-20:]=np.log(0)
    swcL.append(swc1)
    hsurf.append(hm-rmean[0,500+ind-20])
    zsurf.append(zKa_slice[i,500+ind-20-1])
#plt.pcolormesh(time[a[0],0],hm-rmean[0,:],zKa_slice.T,cmap='jet',vmax=30)
plt.pcolormesh(time[a[0],0],hm-rmean[0,:],dZ.T,cmap='jet',vmin=0,vmax=10)
plt.plot(time[a[0],0],hsurf)
plt.colorbar()
plt.figure()
plt.plot(zsurf)
plt.figure()
swcm=np.ma.array(np.array(swcL),mask=np.array(swcL)<0.001)
plt.pcolormesh(swcm[:,::-1].T,vmin=0,vmax=5,cmap='jet')
plt.colorbar()
