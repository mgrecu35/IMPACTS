from netCDF4 import *
import matplotlib.pyplot as plt
import numpy as np
fname='/media/grecu/ExtraDrive1/IMPACTS/IMPACTS_HIWRAP_L1B_RevA_20200227T074338_to_20200227T142238.h5'
fname='/media/grecu/ExtraDrive1/IMPACTS/IMPACTS2020_HIWRAP_L1B_RevB_20200205.h5'

fname_merra='https://goldsmr5.gesdisc.eosdis.nasa.gov/opendap/MERRA2/M2I3NPASM.5.12.4/2020/02/MERRA2_400.inst3_3d_asm_Np.20200227.nc4'
import json
import pickle
iread=0
if iread==1:
    fhm=Dataset(fname_merra)
    latx,lonx=[43.25151,-76.40099]
    lonm=fhm['lon'][:]
    latm=fhm['lat'][:]
    dlat=0.5
    dlon=0.625
    iy=int((latx+90)/dlat)
    ix=int((lonx+180)/dlon)
    t1d=(fhm['T'][:,:,iy,ix].data)
    qv1d=(fhm['QV'][:,:,iy,ix].data)
    h1d=(fhm['H'][:,:,iy,ix].data)
    h_sfc=(fhm['PHIS'][:,iy,ix].data/9.81)
    p_sfc=(fhm['PS'][:,iy,ix].data)
    p_lev=(fhm['lev'][:].data)

    d={"t_1d":t1d,"qv_1d":qv1d,"h_1d":h1d,"h_sfc":h_sfc,\
       "p_sfc":p_sfc,"p_lev":p_lev}
    pickle.dump(d,open("merraEnv.pklz","wb"))
else:
    dmerra=pickle.load(open("merraEnv.pklz","rb"))
    rho=dmerra['p_lev'][:]/287/dmerra['t_1d'][:,:]*1e2
freqs=[89,166,183.3+3,183.3+7]
import combAlg as sdsu
kext_atmL=[]

for freq in freqs:
    i0=3
    iret=1
    kext1=[]
    h1L=[]
    for k in range(2,28):
        absair,abswv = sdsu.gasabsr98(freq,dmerra['t_1d'][3,k],1.0*dmerra['qv_1d'][3,k]*rho[3,k],\
                                     dmerra['p_lev'][k]*1e2,iret)
        kext1.append(absair+abswv)
        h1L.append(dmerra['h_1d'][3,k])
    kext_atmL.append(np.interp(np.arange(80)*250+375,h1L,kext1))
    t1dMerra=np.interp(np.arange(81)*250+250,h1L,dmerra['t_1d'][3,2:28])

dWRF=pickle.load(open("envData.pklz","rb"))
kextWRF=dWRF['kext']
tWRF=dWRF['T']
zWRF=dWRF['z']
kext_atmL2=[]
h1m=(np.arange(80)*250+375)/1000.
zWRFm=0.5*(zWRF[1:]+zWRF[:-1])
for kext in kextWRF:
    kexti=np.interp(h1m,zWRFm,kext)
    kext_atmL2.append(kexti)
t1dMerra=np.interp(np.arange(81)*.250+.250,zWRFm,tWRF)

kext_atmL2=np.array(kext_atmL2)
#stop
fhH=Dataset(fname)
from scatteringRet import *
from dualFreqRetr import *
zKu=fhH['Products/Ku/Combined/Data/dBZe']
zKa=fhH['Products/Ka/Combined/Data/dBZe']
lat=fhH['Navigation/Data/Latitude']
lon=fhH['Navigation/Data/Longitude']
alt=fhH['Navigation/Data/Height']
rrange=fhH['/Products/Information/Range']
import datetime
dz=1.0
dt=(datetime.datetime(2020,2,5)-datetime.datetime(1970,1,1)).days*24*3600
time=(fhH['Time/Data/TimeUTC'][:]-dt)/3600.
a=np.nonzero((time-21.7)*(time-22.1)<0)
zKa_slice=zKa[a[0],:]+dz
zKu_slice=zKu[a[0],:]+dz
dZ=zKu_slice[:,:]-zKa_slice[:,:]
hm=alt[a[0]].mean()
nr=zKa_slice.shape[1]
nx=zKa_slice.shape[0]
lon_a=lon[a[0]]
lat_a=lat[a[0]]
rmean=rrange[:]
a1=np.nonzero(zKu_slice[900:1900,300:(192*3-35)]>25)
for i1,j1 in zip(a1[0],a1[1]):
    zKu_slice[i1+900,j1+300]=25

#print(zKu_slice[900:1900,300:(192*3-35)][a1])
#stop

plt.subplot(111)
plt.pcolormesh(time[a[0]],hm-rmean[:192*3],(zKu_slice[:,:192*3]-zKa_slice[:,:192*3]).T,cmap='jet',vmin=0,vmax=10)
plt.xlim(time[a[0][1050]],time[a[0][1900]])
plt.title('Observed Ku and Ka reflectivity\n 05 February 2020')
plt.ylim(0,8000)
plt.colorbar()
plt.figure()
plt.subplot(111)
plt.pcolormesh(time[a[0]],hm-rmean[:192*3],(zKu_slice[:,:192*3]).T,cmap='jet',vmin=0,vmax=30)
plt.xlim(time[a[0][1050]],time[a[0][1900]])
plt.title('Observed Ku and Ka reflectivity\n 05 February 2020')
plt.ylim(0,8000)
plt.colorbar()
from getFlightLine import info1,info2
plt.figure(figsize=(8,12))
plt.subplot(211)
plt.pcolormesh(time[a[0]],hm-rmean[300:(192*3-35)],(zKu_slice[:,300:(192*3-35)]).T,cmap='jet',vmin=0,vmax=30)
plt.xlim(time[a[0][1050]],time[a[0][1900]])
plt.title('Observed Ku  reflectivity\n 05 February 2020')
plt.ylim(0,8000)
plt.subplot(212)
plt.pcolormesh(time[a[0]],hm-rmean[300:(192*3-35)],(zKa_slice[:,300:(192*3-35)]).T,cmap='jet',vmin=0,vmax=30)
plt.xlim(time[a[0][1050]],time[a[0][1900]])
plt.title('Observed Ka  reflectivity')
plt.ylim(0,8000)
plt.colorbar()

plt.figure()
plt.plot(lon[a[0][1050:1900]],lat[a[0][1050:1900]])
info1=np.array(info1)
plt.plot(info1[:,1],info1[:,2])
info2=np.array(info2)


#plt.figure()
#plt.plot((zKu_slice[1700,:192*3]-zKa_slice[1700,:192*3]),hm-rmean[0:192*3])

import pickle

d=pickle.load(open("knnModelsIMPACTS_RC.pklz","rb"))
knn_dm=d["knn_dm"]
knn_iwc=d["knn_iwc"]
knn_zw=d["knn_zw"]
dm_ret=zKa_slice[900:1900,300:(192*3-35)].copy()*0
iwc_ret=zKa_slice[900:1900,300:(192*3-35)].copy()*0
zw_ret=zKa_slice[900:1900,300:(192*3-35)].copy()*0
zKu_slice[900:1900,300:(192*3-35)][zKu_slice[900:1900,300:(192*3-35)]>29]=29
dfr_knn=(zKu_slice[900:1900,300:(192*3-35)]-zKa_slice[900:1900,300:(192*3-35)])
zku_knn=zKu_slice[900:1900,300:(192*3-35)]
a1=np.nonzero(zku_knn==zku_knn)
b1=np.nonzero(dfr_knn[a1]==dfr_knn[a1])
#dfr_knn[dfr_knn>1]=1
X=np.array([(zku_knn[a1][b1]-12)/16.,(dfr_knn[a1][b1]-2)/2.5]).T
iwc_p=knn_iwc.predict(X)
zw_p=knn_zw.predict(X)
it=0
for i,j in zip(a1[0][b1],a1[1][b1]):
    iwc_ret[i,j]=iwc_p[it]
    zw_ret[i,j]=zw_p[it]
    it+=1

collocZ=[]
tL=[]
for i in range(1050,1900):
    ii=np.argmin((lon[a[0][i]]-info1[:,1])**2+(lat[a[0][i]]-info1[:,2])**2)
    ik=np.argmin(abs(info1[ii,3]-(alt[a[0][i]]-rmean[420:460])))
    #print(zKu[a[0][i],ik+420],zKa[a[0][i],ik+420],info2[ii,0],info2[ii][1])
    collocZ.append([zKu[a[0][i],ik+420],zKa[a[0][i],ik+420],info2[ii,4],info2[ii][5],\
                    info2[ii,-1],info2[ii,-2],iwc_ret[150+i-1050,420+ik-300]])
    tL.append([time[a[0][i]],info1[ii,0]])

tL=np.array(tL)
collocZ=np.array(collocZ)
plt.figure()
plt.subplot(211)
plt.plot(collocZ[:,0])
plt.plot(collocZ[:,2])
plt.subplot(212)
plt.plot(collocZ[:,1])
plt.plot(collocZ[:,3])

from scipy.ndimage import gaussian_filter

iwc_ret=gaussian_filter(iwc_ret,sigma=1)
plt.figure()
plt.pcolormesh(time[a][900:1900],hm-rmean[300:(192*3-35)],iwc_ret.T,cmap='jet',vmin=0,vmax=1.2)
plt.title('Retrieved IWC\n 05 February 2020')
plt.ylabel('Height(m)')
plt.xlabel('Time')
cbar=plt.colorbar()
cbar.ax.set_title('g/m^3')
plt.savefig('retIWC_05Feb2020_2.png')

fnameCRS='/media/grecu/ExtraDrive1/IMPACTS/IMPACTS_CRS_L1B_RevA_20200205T192324_to_20200206T004925.h5'
fCRS=Dataset(fnameCRS)
zW=fCRS['Products/Data/dBZe']
latCRS=fCRS['Navigation/Data/Latitude']
lonCRS=fCRS['Navigation/Data/Longitude']
altCRS=fCRS['Navigation/Data/Height']
rrangeCRS=fCRS['/Products/Information/Range']
timeCRS=(fCRS['Time/Data/TimeUTC'][:]-dt)/3600.
aC=np.nonzero((timeCRS[:,0]-21.7)*(timeCRS[:,0]-22.1)<0)
hmCRS=altCRS[a[0]].mean()
rmeanCRS=rrangeCRS[:]
zW_slice=zW[aC[0],:]

hint=750+np.arange(54)*125
zwm1=np.zeros((54),float)
zwm2=np.zeros((54),float)
count=np.zeros((54),float)
for i1 in a[0][900:1900]:
    ind=np.argmin(abs(timeCRS[aC[0]]-time[i1]))
    #print(timeCRS[aC[0][ind]],time[i1])
    zw1=np.interp(hint,(hmCRS-rmeanCRS[0,:192*3][::-1]),zW_slice[ind,:192*3][::-1])
    zwret1=np.interp(hint,(hm-rmean[300:(192*3-35)])[::-1],zw_ret[i1-a[0][900],::-1])
    for k in range(54):
        if zw1[k]>-5 and zw1[k] and zwret1[k]>-5 and zwret1[k]<30:
            count[k]+=1
            zwm1[k]+=zw1[k]
            zwm2[k]+=zwret1[k]

ac=np.nonzero(count>0)
zwm1[ac]/=count[ac]
zwm2[ac]/=count[ac]
plt.figure()
plt.plot(zwm1,hint)
plt.plot(zwm2,hint)

zw_ret=gaussian_filter(zw_ret,sigma=1)
plt.figure(figsize=(8,12))
plt.subplot(211)
plt.pcolormesh(timeCRS[aC[0],0],hmCRS-rmeanCRS[0,:192*3],(dz+zW_slice[:,:192*3]).T,cmap='jet',vmin=0,vmax=25)
plt.xlim(time[a[0][900]],time[a[0][1900]])
plt.title('Observed W reflectivity')
plt.ylim(0,8000)
plt.colorbar()
plt.subplot(212)
plt.pcolormesh(time[a][900:1900],hm-rmean[300:(192*3-35)],(zw_ret).T,cmap='jet',vmin=0,vmax=25)
plt.title('Retrieved Zw\n 05 February 2020')
plt.ylabel('Height(m)')
plt.xlabel('Time')
cbar=plt.colorbar()
cbar.ax.set_title('g/m^3')
plt.savefig('retZw_05Feb2020_2.png')

plt.figure()
plt.plot(time[a[0][1050:1900]],collocZ[:,-1])
plt.plot(time[a[0][1050:1900]],collocZ[:,-3])
plt.plot(time[a[0][1050:1900]],collocZ[:,-2])
plt.ylabel('IWC (g/m^3)')
plt.xlabel('Time')
plt.legend(['Radar','P3 Merged'])
plt.savefig('iwcRetComp_20200205.png')
