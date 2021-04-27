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

dFR_Table=np.array(kuka_Table['dFR'])
zKu_Table=np.array(kuka_Table['zKu'])
dm_Table=np.array(kuka_Table['dm'])
swc_Table=np.array(kuka_Table['swc'])
dFR=np.array(dFR)
dmRet=dFR.copy()*0
dm_h_table=np.loadtxt('dm-hgt.txt')
dm_h_table[:,0]*=0.7
@jit(nopython=True)
def fromdfr(dFR,zku,dmRet,swcRet,dFR_table,dm_Table,zKu_Table,swc_Table):
    n=dFR.shape[0]
    for i in range(n):
        i1,i2=bisection(dFR[i],dFR_table)
        dmRet[i]=dm_Table[i1]
        dn=0.1*(zku[i]-zKu_Table[i1])/10.0
        swcRet[i]=10**dn*swc_Table[i1]

@jit(nopython=True)
def from_clim_dm(zku,dmRet,swcRet,dm_Table,zKu_Table,swc_Table,dm_h_table,h):
    n=zku.shape[0]
    for i in range(n):
        i1,i2=bisection(h[i],dm_h_table[:,0])
        #print(i1,i2,h[i],dm_h_table[i1,0])
        dm_clim=dm_h_table[i1,1]*0.1
        #print(dm_clim)
        i1,i2=bisection(dm_clim,dm_Table)
        #print(dm_clim,i1,dm_Table[i1],zKu_Table[i1],zku[i1])
        dmRet[i]=dm_Table[i1]
        dn=(zku[i]-zKu_Table[i1])/10.0
        swcRet[i]=10**dn*swc_Table[i1]
        #print(dn,zku[i],zKu_Table[i1])
        #print(dm_clim,i1,dm_Table[i1],zKu_Table[i1],zku[i1],swcRet[i])
    #stop
    
#fromdfr(dFR,dmRet,dFR_Table,dm_Table)
dmRetL=[]
lwpL1=[]
lwpL2=[]
for i in range(nx):
    ind=np.argmax(zKu_slice[i,500:-30])
    swc1=np.zeros((nr),float)-99
    dm=np.zeros((nr),float)-99
    zKa_slice[i,500+ind-20:]=np.log(0)
    a1=np.nonzero(zKa_slice[i,:]==zKa_slice[i,:])
    swc1[a1]=10**(zKa_slice[i,a1[0]]*coeffs[0]+nw[a1]*coeffs[1]+coeffs[2])
    swc1[500+ind-20:]=np.log(0)
    dfr=zKu_slice[i,:500+ind-20]-zKa_slice[i,:500+ind-20]
    a1=np.nonzero(dfr==dfr)
    dmRet=dfr[a1].data.copy()*0
    zku=zKu_slice[i,:500+ind-20][a1].data
    h1=(hm-rmean[0,:500+ind-20][a1])/1000.0
    h1=h1.data-h1[-1]+0.25
    swcRet=dmRet.copy()
    swcRet_clim=dmRet.copy()
    dmRet_clim=dmRet.copy()
    fromdfr(dfr[a1].data,zku,dmRet,swcRet,dFR_Table,dm_Table,zKu_Table,swc_Table)
    from_clim_dm(zku,dmRet_clim,swcRet_clim,dm_Table,zKu_Table,swc_Table,dm_h_table,h1)
    lwpL1.append(swcRet.sum()*26.25)
    lwpL2.append(swcRet_clim.sum()*26.25)
    dm[0:500+ind-20][a1]=dmRet
    swcL.append(swc1)
    hsurf.append(hm-rmean[0,500+ind-20])
    zsurf.append(zKa_slice[i,500+ind-20-1])
    dmRetL.append(dm)

plt.rcParams.update({'font.size': 13})
#plt.pcolormesh(time[a[0],0],hm-rmean[0,:],zKa_slice.T,cmap='jet',vmax=30)
plt.subplot(211)
plt.pcolormesh(time[a[0],0],hm-rmean[0,:],zKu_slice.T,cmap='jet',vmin=0,vmax=35)
plt.title('Observed Ku and Ka reflectivity\n 27 Februarie 2020')
#plt.xlabel('Time')
plt.ylabel('Height (km)')
plt.ylim(0,7500)
cbar=plt.colorbar()
cbar.ax.set_title('dBZ')
plt.subplot(212)
plt.pcolormesh(time[a[0],0],hm-rmean[0,:],zKa_slice.T,cmap='jet',vmin=0,vmax=30)
plt.ylabel('Height (km)')
plt.ylim(0,7500)
#plt.plot(time[a[0],0],hsurf)
cbar=plt.colorbar()
cbar.ax.set_title('dBZ')
plt.xlabel('Time (UTC)')
plt.tight_layout()
plt.savefig('observedZKuKa.png')
plt.figure()
plt.plot(zsurf)
plt.figure()
swcm=np.ma.array(np.array(dmRetL),mask=np.array(dmRetL)<0.001)
plt.pcolormesh(time[a[0],0],hm-rmean[0,:],10*swcm[:,:].T,vmin=0,vmax=4,cmap='jet')
plt.ylim(0,7500)
plt.title('Retrieved Dm\n 27 February 2020')
plt.xlabel('Time')
plt.ylabel('Height (km)')
c=plt.colorbar()
c.ax.set_title('mm')
plt.tight_layout()
plt.savefig('retrievedDm.png')
plt.figure()
lwpL1=np.array(lwpL1)/1e3
lwpL2=np.array(lwpL2)/1e3
plt.plot(time[a[0],0],lwpL1)
plt.plot(time[a[0],0],lwpL2)
plt.xlim(9.85,10.25)
plt.xlabel('Time (UTC)')
plt.ylabel('Snow Water Path (mm)')
plt.legend(['Dual Frequency','Ku-only'])
plt.savefig('retrieved_LWP.png')
