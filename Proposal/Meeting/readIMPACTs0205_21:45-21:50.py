from netCDF4 import *
import matplotlib.pyplot as plt
import numpy as np
fname='/media/grecu/ExtraDrive1/IMPACTS/IMPACTS_HIWRAP_L1B_RevA_20200227T074338_to_20200227T142238.h5'
fname='/media/grecu/ExtraDrive1/IMPACTS/IMPACTS2020_HIWRAP_L1B_RevB_20200205.h5'

fname_merra='https://goldsmr5.gesdisc.eosdis.nasa.gov/opendap/MERRA2/M2I3NPASM.5.12.4/2020/02/MERRA2_400.inst3_3d_asm_Np.20200205.nc4'
import json
import pickle
iread=0
if iread==1:
    fhm=Dataset(fname_merra)
    latx,lonx=[43.25151,-76.40099]
    latx,lonx=[41.3883,-89.40024]
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
    pickle.dump(d,open("merraEnv_20200205.pklz","wb"))
else:
    dmerra=pickle.load(open("merraEnv_20200205.pklz","rb"))
    t1d=dmerra["t_1d"]
    h1d=dmerra["h_1d"]
    rho=dmerra['p_lev'][:]/287/dmerra['t_1d'][:,:]*1e2
    h_sfc=dmerra['h_sfc']
freqs=[89,166,183.3+3,183.3+7]
import combAlg as sdsu
kext_atmL=[]

for freq in [95.0]:#freqs:
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

#dWRF=pickle.load(open("envData.pklz","rb"))
#kextWRF=dWRF['kext']
#tWRF=dWRF['T']
#zWRF=dWRF['z']
#kext_atmL2=[]
#h1m=(np.arange(80)*250+375)/1000.
#zWRFm=0.5*(zWRF[1:]+zWRF[:-1])
#for kext in kextWRF:
#    kexti=np.interp(h1m,zWRFm,kext)
#    kext_atmL2.append(kexti)
#t1dMerra=np.interp(np.arange(81)*.250+.250,zWRFm,tWRF)

#kext_atmL2=np.array(kext_atmL2)
#stop
fhH=Dataset(fname)
#from scatteringRet import *
#from dualFreqRetr import *
zKu=fhH['Products/Ku/Combined/Data/dBZe'][:]
zKa=fhH['Products/Ka/Combined/Data/dBZe'][:]
lat=fhH['Navigation/Data/Latitude']
lon=fhH['Navigation/Data/Longitude']
alt=fhH['Navigation/Data/Height'][:]
rrange=fhH['/Products/Information/Range']
import datetime
dz=1.0
dt=(datetime.datetime(2020,2,5)-datetime.datetime(1970,1,1)).days*24*3600
time=(fhH['Time/Data/TimeUTC'][:]-dt)/3600.

ap3=np.nonzero((time-21.95)*(time-22.99)<0)
rmean=rrange[:]

t_env=np.interp((alt.mean()-rmean[:550])[::-1],h1d[7,1:25],t1d[7,1:25])[::-1]
kext_env=np.interp((alt.mean()-rmean[:576])[::-1],h1L,kext1)[::-1]
dzw_env=kext_env.cumsum()*26.25/1000*4.343*2
#stop
#load packages
from pickle import load
import tensorflow as tf

#load scalers 
scaler_X = load(open('./Chase_et_al_2021_NN-master/'+'scaler_X_V2.pkl', 'rb'))
scaler_y = load(open('./Chase_et_al_2021_NN-master/'+'scaler_y_V2.pkl', 'rb'))
model = tf.keras.models.load_model('./Chase_et_al_2021_NN-master/' + 'NN_6by8.h5',custom_objects=None,compile=True)

iwc_model = tf.keras.models.load_model('iwc_model.h5')

for i1 in zip(ap3[0]):
    ak=np.nonzero((alt[i1]-rmean-3000)*(alt[i1]-rmean-4000)<0)
    for k in ak[0]:
        if zKu[i1,k]>24:
            #print(zKu[i1,k])
            zKu[i1,k]=24


a=np.nonzero((time-21.7)*(time-22.2)<0)

zKa_slice=zKa[a[0],:]+dz
zKu_slice=zKu[a[0],:]+dz
dZ=zKu_slice[:,:]-zKa_slice[:,:]
hm=alt[a[0]].mean()
nr=zKa_slice.shape[1]
nx=zKa_slice.shape[0]
lon_a=lon[a[0]]
lat_a=lat[a[0]]

a1=np.nonzero(zKu_slice[900:1900,300:(192*3-35)]>30)
for i1,j1 in zip(a1[0],a1[1]):
    zKu_slice[i1+900,j1+300]=30


plt.subplot(111)
plt.pcolormesh(time[a[0]],hm-rmean[:192*3],(zKu_slice[:,:192*3]-zKa_slice[:,:192*3]).T,cmap='jet',vmin=0,vmax=10)
plt.xlim(time[a[0][1050]],time[a[0][3500]])
plt.title('Observed DFR\n 05 February 2020')
plt.ylim(0,8000)
plt.colorbar()

from getFlightLine import getinfo
t1,t2=21.7,22.15
info1,info2=getinfo(t1,t2)

plt.figure(figsize=(8,12))
plt.subplot(211)
plt.pcolormesh(time[a[0]],hm-rmean[300:(192*3-35)],(zKu_slice[:,300:(192*3-35)]).T,cmap='jet',vmin=0,vmax=30)
plt.xlim(time[a[0][1050]],time[a[0][3500]])
plt.title('Observed Ku  reflectivity\n 05 February 2020')
plt.ylim(0,8000)
plt.colorbar()
plt.subplot(212)
plt.pcolormesh(time[a[0]],hm-rmean[300:(192*3-35)],(zKa_slice[:,300:(192*3-35)]).T,cmap='jet',vmin=0,vmax=30)
plt.xlim(time[a[0][1050]],time[a[0][3500]])
plt.title('Observed Ka  reflectivity')
plt.ylim(0,8000)
plt.colorbar()

#stop



#plt.figure()
#plt.plot((zKu_slice[1700,:192*3]-zKa_slice[1700,:192*3]),hm-rmean[0:192*3])

import pickle

fnameCRS='/media/grecu/ExtraDrive1/IMPACTS/IMPACTS_CRS_L1B_RevA_20200205T192324_to_20200206T004925.h5'
fCRS=Dataset(fnameCRS)
zW=fCRS['Products/Data/dBZe'][:,:]+dz
latCRS=fCRS['Navigation/Data/Latitude']
lonCRS=fCRS['Navigation/Data/Longitude']
altCRS=fCRS['Navigation/Data/Height']
rrangeCRS=fCRS['/Products/Information/Range']
timeCRS=(fCRS['Time/Data/TimeUTC'][:]-dt)/3600.
aC=np.nonzero((timeCRS[:,0]-21.7)*(timeCRS[:,0]-22.2)<0)
hmCRS=altCRS[a[0]].mean()
rmeanCRS=rrangeCRS[:]
zW_slice=zW[aC[0],:]

hint=750+np.arange(54)*125
zwm1=np.zeros((54),float)
zwm2=np.zeros((54),float)
count=np.zeros((54),float)
zW_L=[]
for i1 in a[0][200:3500]:
    ind=np.argmin(abs(timeCRS[aC[0]]-time[i1]))
    zw1=np.interp((alt[i1]-rmean[:192*3])[::-1],\
                  (altCRS[aC[0][ind]]-rmeanCRS[0,:192*3][::-1]),zW_slice[ind,:192*3][::-1])
    zW_L.append(zw1[::-1]+dzw_env[:192*3])

dictZ={"zKu":zKu_slice[200:3500],"zKa":zKa_slice[200:3500],\
       "zW":np.array(zW_L),"time":time[a[0][200:3500]],\
       "rrange":rmean,"alt":alt[a[0][200:3500]],\
       "lon":lon[a[0][200:3500]],"lat":lat[a[0][200:3500]]}
pickle.dump(dictZ,open("z_obs_Feb05_22:00-22:06","wb"))
stop
zW_L=np.array(zW_L)
rs=300
iwc2d_1=np.zeros((3300,450-rs),float)
iwc2d_2=np.zeros((3300,450-rs),float)
iwc2d_3=np.zeros((3300,450-rs),float)
zPCs=np.zeros((3300,450-rs,2),float)
rho2d=np.zeros((3300,450-rs),float)-99
attw2d=np.zeros((3300,450-rs),float)
dfrkukaL=[]
dfrkawL=[]
d=pickle.load(open("knnModelsIMPACTS_2_DDA_sorted.pklz","rb"))
d3=pickle.load(open("knnModelsIMPACTS_3_DDA_sorted.pklz","rb"))
d4=pickle.load(open("knnModelsIMPACTS_4_DDA_sorted.pklz","rb"))
d3dwr=pickle.load(open("knnModelsIMPACTS_3_DDA_DWR_sorted.pklz","rb"))
knn_dm=d["knn_dm"]
knn_iwc=d["knn_iwc"]
knn_zw=d["knn_zw"]
knn_attw=d["knn_attw"]
knn_attw2=d3["knn_attw"]

plt.figure()
plt.plot(lon[a[0][200:3500]],lat[a[0][200:3500]]+0.01)
info1=np.array(info1)
plt.plot(info1[:,1],info1[:,2])
info2=np.array(info2)

zku1d=[]
zka1d=[]
zw1d=[]
#{"zpca":zpca,"scalerZ":scalerZ}
dpca=pickle.load(open("zpca.pklz","rb"))
zpca=dpca["zpca"]
scalerZ=dpca["scalerZ"]
for i in range(rs,450,1):
    zkuh1=zKu_slice[200:3500,i]
    zkah1=zKa_slice[200:3500,i]
    zwh1=zW_L[0:,i]
    x2=np.array([(zkuh1-12)/16.,(zkuh1-zkah1-2)/2.5]).T
    x3=np.array([(zkuh1-12)/16.,(zkuh1-zkah1-2)/2.5,(zwh1-2)/7.]).T
    a1=np.nonzero(zkuh1==zkuh1)
    a2=np.nonzero(zkah1[a1]==zkah1[a1])
    a3=np.nonzero(zwh1[a1][a2]==zwh1[a1][a2])
    iwc_1=knn_iwc.predict(x2[a1][a2][a3])
    iwc_2=d3["knn_iwc"].predict(x3[a1][a2][a3])
    attw1=d3["knn_attw"].predict(x2[a1][a2])

    attw2d[a1[0][a2],i-rs]=attw1
    piaW=attw2d.sum(axis=-1)*2*26.25/1000.
    piaW[a1[0][a2]]-=attw2d[a1[0][a2],i-rs]*26.25/100.
    zwh1=zW_L[0:,i]+piaW
    x2=np.array([(zkuh1-12)/16.,(zkuh1-zkah1-2)/2.5]).T
    x3=np.array([(zkuh1-12)/16.,(zkuh1-zkah1-2)/2.5,(zwh1-2)/7.]).T
    normTc=(zkuh1*0+t_env[i]- 273.15+11)/8.0
    x31=np.array([(zkuh1-12)/16.,(zkuh1-zkah1-2)/1.5,(zkah1-zwh1-2.5)/1.5,\
                  normTc]).T

    x4=np.array([(zkuh1-12)/16.,(zkuh1-zkah1-2)/2.5,(zwh1-2)/7.,normTc]).T

   
    
    #scale itr properly
    
    a1=np.nonzero(zkuh1==zkuh1)
    a2=np.nonzero(zkah1[a1]==zkah1[a1])
    a3=np.nonzero(zwh1[a1][a2]==zwh1[a1][a2])
    zku1d.extend(zkuh1[a1][a2][a3])
    zka1d.extend(zkah1[a1][a2][a3])
    zw1d.extend(zwh1[a1][a2][a3])
    zs=np.array([zkuh1[a1][a2][a3],zkah1[a1][a2][a3],zwh1[a1][a2][a3]]).T
    zs_sc=scalerZ.transform(zs)
    zPCs1=zpca.transform(zs_sc)
    zPCs[a1[0][a2][a3],i-rs,0]=zPCs1[:,0]
    zPCs[a1[0][a2][a3],i-rs,1]=zPCs1[:,1]
    iwc_1=knn_iwc.predict(x2[a1][a2][a3])
    iwc_2=d3dwr["knn_iwc"].predict(x3[a1][a2][a3])
    #iwc_2=d3dwr["knn_iwc"].predict(x31[a1][a2][a3])
    #iwc_2=iwc_model.predict(x31[a1][a2][a3])[:,0]
    #Yice=(np.log(zSims[:,-1])+1.89)/1.94
    iwc_2=(0.5*iwc_2+0.4)
    #iwc_2=d4["knn_iwc"].predict(x4[a1][a2][a3])
    #iwc_2=d["regr_iwc"].predict(x3[a1][a2][a3])
    rho1=np.log10(d4["knn_rho"].predict(x4[a1][a2][a3]))
    x11=(zkuh1-scaler_X.mean_[0])/scaler_X.scale_[0] #Ku
    x12=((zkuh1-zkah1)- scaler_X.mean_[1])/scaler_X.scale_[1] #DFR Ku - Ka
    x13=(zkuh1*0+t_env[i]- 273.15-scaler_X.mean_[2])/scaler_X.scale_[2] #T
    #
    x_nn=np.array([x11[a1][a2][a3],x12[a1][a2][a3],x13[a1][a2][a3]]).T
    #conduct the retrieval 
    yhat = model.predict(x_nn)
    yhat = scaler_y.inverse_transform(yhat)
    Nw=10**yhat[:,0]
    Dm=10**yhat[:,1]/1000.
    iwc_nn=(Nw*(Dm)**4*1000*np.pi)/4**(4)*1e3
    #stop
    iwc2d_1[a1[0][a2][a3],i-rs]=iwc_1
    iwc2d_2[a1[0][a2][a3],i-rs]=iwc_2
    iwc2d_3[a1[0][a2][a3],i-rs]=iwc_nn
    rho2d[a1[0][a2][a3],i-rs]=rho1
    a4=np.nonzero(zkuh1[a1][a2][a3]>5)
    dwr1=zkuh1[a1][a2][a3][a4]-zkah1[a1][a2][a3][a4]
    dwr1[dwr1<0]=0
    
    dwr2=zkah1[a1][a2][a3][a4]-zwh1[a1][a2][a3][a4]
    dwr2[dwr2<0]=0
    a5=np.nonzero((time[a][a1][a2][a3][a4]-22)*\
                  (time[a][a1][a2][a3][a4]-22.06)<0)
    if i>420:
        dfrkukaL.extend(dwr1[a5])
        dfrkawL.extend(dwr2[a5])
    print('Means=',iwc_1.mean(),iwc_2.mean(),iwc_nn.mean())
    print('CC=',np.corrcoef(iwc_1,iwc_2)[0,1], np.corrcoef(iwc_2,iwc_nn)[0,1],np.corrcoef(zkah1[a1][a2][a3],zwh1[a1][a2][a3])[0,1])

#zs_sc[:,0]*zpca.components_[1,0]+zs_sc[:,1]*zpca.components_[1,1]+zs_sc[:,2]*zpca.components_[1,2]

ipca=0
if ipca==1:
    from sklearn import preprocessing
    scalerZ  = preprocessing.StandardScaler()
    xz=np.array([zku1d,zka1d,zw1d]).T
    scalerZ.fit(xz)
    xz_s=scalerZ.transfor(xz)
    from sklearn.decomposition import pca
    zpca = pca.PCA()
    zpca.fit(xz_s)
    PCs=zpca.transform(xz_s)
    pickle.dump({"zpca":zpca,"scalerZ":scalerZ},open("zpca.pklz","wb"))
    stop
    
i1=550-200
i2=2600-200
collocZ=[]
tL=[]
for i in range(i1,i2):
   h1=alt[a[0][i+200]]-rmean[rs:450]
   ii=np.argmin((lon[a[0][i+200]]-info1[:,1])**2+(lat[a[0][i+200]]-info1[:,2])**2)
   #print(ii)
   ik=np.argmin(abs(info1[ii,3]-h1))
   #print(ii,ik,(lon[a[0][i+200]]-info1[ii,1])**2+(lat[a[0][i+200]]-info1[ii,2])**2)
   collocZ.append([zKu[a[0][i+200],ik+rs],zKa[a[0][i+200],ik+rs],\
                   info2[ii,4],info2[ii][5],rho2d[i,ik],zPCs[i,ik,0],\
                   zPCs[i,ik,1],\
                   iwc2d_3[i,ik],iwc2d_1[i,ik],\
                   info2[ii,-1],info2[ii,-2],iwc2d_2[i,ik]])
   tL.append([time[a[0][i+200]],info1[ii,0],time[a[0][i+200]]-info1[ii,0]])

collocZ=np.array(collocZ)
tL=np.array(tL)

import matplotlib
plt.figure()
plt.plot(tL[:,0],collocZ[:,-3])
plt.plot(tL[:,0],collocZ[:,-1])
plt.plot(tL[:,0],collocZ[:,-5])
plt.legend(['NCAR Probes','TFRet','NN'])
plt.ylabel('IWC (mm/h)')
plt.xlabel('Time')
plt.title('5 February 2020')
plt.savefig('IWC_20200205.png')


plt.figure(figsize=(8,12))
fig1=plt.subplot(311)
plt.pcolormesh(time[a[0][200:3500]],alt.mean()-rmean[rs:450],iwc2d_2.T,cmap='jet',norm=matplotlib.colors.LogNorm(),vmin=0.01,vmax=1)
#plt.ylim(450,350)
plt.title('Triple Frequency')
fig1.xaxis.set_visible(False)
c=plt.colorbar()
c.ax.set_title('g/m^3')
fig2=plt.subplot(312)
plt.pcolormesh(time[a[0][200:3500]],alt.mean()-rmean[rs:450],iwc2d_1.T,cmap='jet',norm=matplotlib.colors.LogNorm(),vmin=0.01,vmax=1)
plt.title('Dual Frequency')
fig2.xaxis.set_visible(False)
plt.ylabel('Height (m)')
plt.colorbar()
plt.subplot(313)
plt.pcolormesh(time[a[0][200:3500]],alt.mean()-rmean[rs:450],iwc2d_3.T,cmap='jet',norm=matplotlib.colors.LogNorm(),vmin=0.01,vmax=1)
#plt.ylim(450,350)
plt.title('Chase et al. 2021 NN')
plt.colorbar()
plt.savefig('IWC2D_20200205.png')
stop

zW_L=np.array(zW_L)
dm_ret=zKa_slice[900:1900,300:(192*3-35)].copy()*0
iwc_ret=zKa_slice[900:1900,300:(192*3-35)].copy()*0
zw_ret=zKa_slice[900:1900,300:(192*3-35)].copy()*0
attw_ret=zKa_slice[900:1900,300:(192*3-35)].copy()*0

zKu_slice[900:1900,300:(192*3-35)][zKu_slice[900:1900,300:(192*3-35)]>29]=29
dfr_knn=(zKu_slice[900:1900,300:(192*3-35)]-zKa_slice[900:1900,300:(192*3-35)])
zku_knn=zKu_slice[900:1900,300:(192*3-35)]
zw_knn=zW_L[:,300:(192*3-35)]
a1=np.nonzero(zku_knn==zku_knn)
b1=np.nonzero(dfr_knn[a1]==dfr_knn[a1])
c1=np.nonzero(zw_knn[a1][b1]==zw_knn[a1][b1])
#dfr_knn[dfr_knn>1]=11
X2=np.array([(zku_knn[a1][b1]-12)/16.,(dfr_knn[a1][b1]-2)/2.5]).T#, (zw_knn[a1][b1][c1]-2)/7.]).T
X3=np.array([(zku_knn[a1][b1][c1]-12)/16.,(dfr_knn[a1][b1][c1]-2)/2.5, (zw_knn[a1][b1][c1]-2)/7.]).T
iwc_p=knn_iwc.predict(X2)
knn_iwc3=d3["knn_iwc"]
iwc_p3=knn_iwc3.predict(X3)
zw_p=knn_zw.predict(X2)
attw_p=knn_attw.predict(X2)
it=0
#attw_ret[attw_ret!=attw_ret]=0

for i,j in zip(a1[0][b1],a1[1][b1]):
    iwc_ret[i,j]=iwc_p[it]
    zw_ret[i,j]=zw_p[it]
    attw_ret[i,j]=attw_p[it]
    it+=1
it=0
for i,j in zip(a1[0][b1][c1],a1[1][b1][c1]):
    iwc_ret[i,j]=iwc_p3[it]
    it+=1
piaW=np.zeros((900),float)

attw_ret[attw_ret!=attw_ret]=0
intAtt=attw_ret.copy()
for i in range(900):
    piaW[i]=attw_ret[i,:].sum()*26.25*2/1e3
    intAtt[i,:]=attw_ret[i,:].cumsum()*26.25*2/1e3
zw_knn+=intAtt

X3=np.array([(zku_knn[a1][b1][c1]-12)/16.,(dfr_knn[a1][b1][c1]-2)/2.5, (zw_knn[a1][b1][c1]-2)/7.]).T
d3=pickle.load(open("knnModelsIMPACTS_3_DDA_sorted.pklz","rb"))
d4=pickle.load(open("knnModelsIMPACTS_4_DDA_sorted.pklz","rb"))
knn_iwc3=d3["knn_iwc"]
iwc_p3=knn_iwc3.predict(X3)

it=0
for i,j in zip(a1[0][b1][c1],a1[1][b1][c1]):
    iwc_ret[i,j]=iwc_p3[it]
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
plt.pcolormesh(time[a][900:1900],hm-rmean[300:(192*3-35)],iwc_ret.T,cmap='jet',vmin=0,vmax=1.0)
plt.title('Retrieved IWC\n 05 February 2020')
plt.ylabel('Height(m)')
plt.xlabel('Time')
cbar=plt.colorbar()
cbar.ax.set_title('g/m^3')
plt.savefig('retIWC_05Feb2020_2.png')



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
plt.pcolormesh(timeCRS[aC[0],0],hmCRS-rmeanCRS[0,:192*3],(zW_slice[:,:192*3]).T,cmap='jet',vmin=0,vmax=25)
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
