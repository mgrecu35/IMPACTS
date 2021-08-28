from sklearn import preprocessing
from netCDF4 import Dataset

fh=Dataset('IMPACTS_PSD_Database.nc')
#fh=Dataset("IMPACTS_PSDs.nc")
import numpy as np
PSDs=fh["Nc"][:]
Z=fh["Z_sim"][:]
iwc=fh["iwc"][:]
dZ=fh["dZ"][:]
Nc=fh["Nc"][:]
Dm=fh['Dm'][:]
Nw=4**4/np.pi*iwc/(0.1*Dm)**4/1e6
rho=fh['rho'][:]
X=np.array([iwc,Dm,np.log10(Nw),np.log10(rho)]).T
Y=Z.copy()
Y[:,1]=Z[:,0]-Z[:,1]
Y[:,2]=Z[:,1]-Z[:,2]
from sklearn import preprocessing
scalerX  = preprocessing.StandardScaler()
scalerY  = preprocessing.StandardScaler()
scalerX.fit(X)
scalerY.fit(Z)
Xs=scalerX.transform(X)
Ys=scalerY.transform(Z)
n=Xs.shape[0]

r=np.random.random(n)
a=np.nonzero(r<0.5)
b=np.nonzero(r>0.5)
X_train=Xs[a[0],:]
X_test=Xs[b[0],:]
y_train=Ys[a[0],:]
y_test=Ys[b[0],:]

import tensorflow as tf
import tensorflow.keras as keras
Keras = keras.backend


def iwc_model(n1,n2):
    input1 = keras.layers.Input(shape=[n1])
    z = keras.layers.Dense(40, activation="relu")(input1)
    z = keras.layers.Dropout(0.1)(z,training=True)
    z = keras.layers.Dense(40, activation="relu")(z)
    z = keras.layers.Dropout(0.1)(z,training=True)
    output = keras.layers.Dense(n2)(z)
    model = keras.models.Model(
        inputs=[input1], outputs=[output])
    return model
itrain=1
if itrain==1:
    model=iwc_model(4,3)
    model.compile(optimizer=tf.keras.optimizers.Adam(),  \
                  loss='mse',\
                  metrics=[tf.keras.metrics.MeanSquaredError()])
    
    history = model.fit(X_train[:,:4], y_train, batch_size=32,epochs=50,
                        validation_data=(X_test[:,:4], y_test))

    model.save("DFRmodel.h5")
else:
    model=tf.keras.models.load_model("DFRmodel.h5")
    
import pickle
pickle.dump({"InputScaler":scalerX,"OutputScaler":scalerY},open("zDFR_OutputScaler.pklz","wb"))


#A=inv(K.T*inv(Se)*K+inv(Sa))*(K.T*inv(Se)*K)

#stop
#Nw=4^4/pi/rhow*lwc/
from sklearn.cluster import  MiniBatchKMeans
random_state=0
import pickle
zmin=20

az=np.nonzero((Z[:,0]-zmin)*(Z[:,0]-(zmin+1))<0)
xinput=Xs[az[0],:]
yn=model.predict(xinput)
K=[]
for i in range(4):
    xinput_i=Xs[az[0],:].copy()
    xinput_i[:,i]+=0.1
    yn_i=model.predict(xinput_i)
    K.append((yn_i-yn)/0.1)

K=np.array(K).mean(axis=1)
Se=np.zeros((3,3),float)
Se[0,0]=1/16.
Se[1,1]=2/15
Se[2,2]=2/14.

KSK=np.dot(np.dot(K,np.linalg.inv(Se)),K.T)
Sa=np.cov(Xs[:,:].T)
A1=np.linalg.pinv(np.linalg.inv(Sa)+KSK)
A=np.dot(A1,KSK)
#Z_train=Z[a[0],:]
#Z_test=Z[b[0],:]
#iwc_train=iwc[a]
#iwc_test=iwc[b]
#train=1
#if train==1:
#    kmeans = MiniBatchKMeans(n_clusters=60, batch_size=1000,random_state=random_state).fit(Z_train)
#    pickle.dump(kmeans,open("impacts_PSD_kmeans.pklz","wb"))
#else:
#    kmeans=pickle.load(open("impacts_PSD_kmeans.pklz","rb")

#labels_=kmeans.predict(Z_test)
#logdn=np.log(1.1)
#dZdN_all=[]
