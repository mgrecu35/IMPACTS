import numpy as np
from netCDF4 import Dataset
fh=Dataset("Chase_et_al_2021_NN-master/Unrimed_simulation_wholespecturm_train_V2.nc","r")

import matplotlib.pyplot as plt
zKu=fh["Z"][:]
zKa=fh["Z2"][:]
zW=fh["Z3"][:]
Dm=fh["Dm"][:]
Nw=fh["Nw"][:]
iwc=fh["IWC"][:]

import numpy as np

#X2,iwcS,dBZs,dmS,Nt,dBase=readDBabseLiu()

from readIMPACTS_PSDs import *
kextL=np.array(kextL)
#from calcZfromPSDS import *
import json

#plt.scatter(zKu-10*np.log10(Nw),Dm)
#plt.scatter(zKu-10*np.log10(Nw),np.log10(iwc/Nw))



plt.figure()
plt.scatter(zKu-zKa,Dm)
#plt.plot(zKu_rg-zKa_rg,0.01*np.array(dm_rg),color='red')
plt.scatter(zSims[:,0]-zSims[:,1],zSims[:,-6]*1e-3)
plt.figure()
import matplotlib
plt.subplot(121)
plt.hist2d(zSims[:,0]-zSims[:,1],zSims[:,-6]*1e-3,bins=(np.arange(24)*0.5,\
                                                        np.arange(40)*0.000075),\
           norm=matplotlib.colors.LogNorm(),vmin=0.1,cmap='jet')

plt.subplot(122)
plt.hist2d(zKu-zKa,Dm,bins=(np.arange(24)*0.5,\
                                                        np.arange(40)*0.000075),\
           norm=matplotlib.colors.LogNorm(),vmin=0.1,cmap='jet')
#stop
#stop
#plt.figure()
#plt.scatter(zKa-10*np.log10(Nw),Dm)
#plt.plot(zKa_rg-10*np.log10(10**6.6),0.01*np.array(dm_rg),color='red')

from sklearn.neighbors import KNeighborsRegressor

nne=15
neigh = KNeighborsRegressor(n_neighbors=nne,weights='distance')

neigh_iwc = KNeighborsRegressor(n_neighbors=nne,weights='distance')
neigh_zw = KNeighborsRegressor(n_neighbors=nne,weights='distance')
neigh_attw = KNeighborsRegressor(n_neighbors=nne,weights='distance')

neigh2 = KNeighborsRegressor(n_neighbors=nne,weights='distance')
fh=Dataset("Chase_et_al_2021_NN-master/Unrimed_simulation_wholespecturm_train_V2.nc","r")

iopt=2
if iopt==1:
    n=zKu.shape[0]
    X=np.zeros((n,3),float)
    X[:,0]=(zKu-12)/16.
    X[:,1]=(zKu-zKa-2)/2.5
    X[:,2]=fh["T_env"][:]/10.
    X[:,2]=(zW-2)/7.0
    Y=Dm
    Yice=iwc
    Yzw=zW
    Yattw=zW*0
else:
    n=zSims.shape[0]
    X=np.zeros((n,3),float)
    X[:,0]=zSims[:,0]
    X[:,1]=(zSims[:,0]-zSims[:,1])
    X[:,2]=(zSims[:,2])
    x_m=X.mean(axis=0)
    x_s=X.std(axis=0)
    Y=zSims[:,-6]
    Yice=zSims[:,-1]
    Yzw=zSims[:,2]
    Yattw=kextL[:,-1]
    X[:,0]=(X[:,0]-12)/16.
    X[:,1]=(X[:,1]-2)/2.5
    X[:,2]=(X[:,2]-2)/7.0
    

r=np.random.random(n)
a=np.nonzero(r<0.5)
b=np.nonzero(r>0.5)
X_train=X[a[0],:]
X_test=X[b[0],:]
y_train=Y[a[0]]
y_test=Y[b[0]]
yice_train=Yice[a[0]]
yice_test=Yice[b[0]]
yzw_train=Yzw[a[0]]
yzw_test=Yzw[b[0]]
yattw_train=Yattw[a[0]]
yattw_test=Yattw[b[0]]
neigh.fit(X_train[:,0:2], y_train)
neigh_iwc.fit(X_train[:,0:2], yice_train)
neigh_zw.fit(X_train[:,0:2], yzw_train)
neigh_attw.fit(X_train[:,0:2], yattw_train)
yp=neigh.predict(X_test[:,0:2])
yp_iwc=neigh_iwc.predict(X_test[:,0:2])
yp_zw=neigh_zw.predict(X_test[:,0:2])
yp_attw=neigh_attw.predict(X_test[:,0:2])

rms_1=(((yp-y_test)**2).mean())**.5

import pickle
pickle.dump({"knn_dm":neigh,"knn_iwc":neigh_iwc, "knn_zw":neigh_zw,\
             "knn_attw":neigh_attw},open("knnModelsIMPACTS_2_DDA.pklz","wb"))
plt.close('all')
plt.figure()
fig=plt.subplot(111)
fig.set_aspect('equal')
dwr_1=zSims[:,0]-zSims[:,1]
dwr_2=zSims[:,1]-zSims[:,2]

import matplotlib
plt.hist2d(dwr_2,dwr_1,bins=np.arange(30)*0.5,norm=matplotlib.colors.LogNorm(),cmap='jet')
plt.colorbar()
plt.figure()
fig=plt.subplot(111)
fig.set_aspect('equal')
dwr_1=zKu-zKa
dwr_2=zKa-zW
plt.hist2d(dwr_2,dwr_1,bins=np.arange(30)*0.5,norm=matplotlib.colors.LogNorm(),cmap='jet')
plt.colorbar()
