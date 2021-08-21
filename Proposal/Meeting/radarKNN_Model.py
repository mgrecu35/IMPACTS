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


#from readIMPACTS_PSDs import *
import pickle

import matplotlib
d=pickle.load(open("simulatedPropSortedTemp_2.pklz","rb"))
zSims=d["zSims"]
kextL=d["kextL"]
rhoL=d["rhoL"]
tempC=np.array(d["tempC"])
kextL=np.array(kextL)

from sklearn.neighbors import KNeighborsRegressor

nne=30
neigh = KNeighborsRegressor(n_neighbors=nne,weights='distance')

neigh_iwc = KNeighborsRegressor(n_neighbors=nne,weights='distance')
neigh_zw = KNeighborsRegressor(n_neighbors=nne,weights='distance')
neigh_attw = KNeighborsRegressor(n_neighbors=nne,weights='distance')
neigh_rho = KNeighborsRegressor(n_neighbors=nne,weights='distance')

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
    X=np.zeros((n,4),float)
    X[:,0]=zSims[:,0]
    X[:,1]=(zSims[:,0]-zSims[:,1])
    X[:,2]=(zSims[:,2])
    X[:,3]=(tempC[:])
    x_m=X.mean(axis=0)
    x_s=X.std(axis=0)
    Y=zSims[:,-6]
    Yice=(zSims[:,-1])
    Yzw=zSims[:,2]
    Yattw=kextL[:,-1]
    Yrho=np.array(rhoL)
    X[:,0]=(X[:,0]-12)/16.
    X[:,1]=(X[:,1]-2)/1.5
    X[:,2]=(X[:,2]-2)/7
    X[:,3]=(X[:,3]+11)/8.0


dpca=pickle.load(open("zpca.pklz","rb"))
zpca=dpca["zpca"]
scalerZ=dpca["scalerZ"]
zs_sc=scalerZ.transform(zSims[:,0:3])
zPCs1=zpca.transform(zs_sc)

#stop
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
yrho_train=Yrho[a[0]]
yrho_test=Yrho[b[0]]
neigh.fit(X_train[:,0:3], y_train)
neigh_iwc.fit(X_train[:,0:3], yice_train)
neigh_zw.fit(X_train[:,0:3], yzw_train)
neigh_attw.fit(X_train[:,0:2], yattw_train)
neigh_rho.fit(X_train[:,0:3], yrho_train)
yp=neigh.predict(X_test[:,0:3])
yp_iwc=neigh_iwc.predict(X_test[:,0:3])
yp_zw=neigh_zw.predict(X_test[:,0:3])
yp_attw=neigh_attw.predict(X_test[:,0:2])
yp_rho=neigh_rho.predict(X_test[:,0:3])

rms_1=(((yp-y_test)**2).mean())**.5

import pickle
from knnRC import knnRC

from sklearn.ensemble import RandomForestRegressor

neigh_iwcRC=knnRC(zKu,zKa,zW,Dm,iwc,nne)
regr_iwc = RandomForestRegressor(n_estimators=20,max_depth=5, random_state=0)
regr_iwc.fit(X_train[:,0:3],yice_train)
yp_iwc2=regr_iwc.predict(X_test[:,0:3])
print(np.corrcoef(yp_iwc2,yice_test))
      
pickle.dump({"knn_dm":neigh,"knn_iwc":neigh_iwc, "knn_zw":neigh_zw,\
             "knn_attw":neigh_attw,"regr_iwc":regr_iwc,\
             "knn_rho":neigh_rho},\
            open("knnModelsIMPACTS_3_DDA_2_sorted.pklz","wb"))
plt.close('all')
import matplotlib

plt.figure()
fig=plt.subplot(111)
fig.set_aspect('equal')
dwr_1=zSims[:,0]-zSims[:,1]
dwr_2=zSims[:,1]-zSims[:,2]
plt.hist2d(dwr_2,dwr_1,bins=np.arange(30)*0.5,norm=matplotlib.colors.LogNorm(),cmap='jet')

plt.colorbar()
plt.figure()
fig=plt.subplot(111)
fig.set_aspect('equal')
dwr_1=zKu-zKa
dwr_2=zKa-zW
plt.hist2d(dwr_2,dwr_1,bins=np.arange(30)*0.5,norm=matplotlib.colors.LogNorm(),cmap='jet')
plt.colorbar()

stop
import tensorflow as tf
import tensorflow.keras as keras
K = keras.backend

def iwc_model(n1,n2):
    input1 = keras.layers.Input(shape=[n1])
    z = keras.layers.Dense(40, activation="relu")(input1)
    z = keras.layers.Dropout(0.1)(z)
    z = keras.layers.Dense(40, activation="relu")(z)
    z = keras.layers.Dropout(0.1)(z)
    output = keras.layers.Dense(n2)(z)
    model = keras.models.Model(
        inputs=[input1], outputs=[output])
    return model

model=iwc_model(4,1)

model.compile(optimizer=tf.keras.optimizers.Adam(),  \
              loss='mse',\
              metrics=[tf.keras.metrics.MeanSquaredError()])

history = model.fit(X_train[:,:4], yice_train, batch_size=32,epochs=50,
                    validation_data=(X_test[:,:4], yice_test))

model.save('iwc_model.h5')

