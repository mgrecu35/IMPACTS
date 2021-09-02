from sklearn import preprocessing
from netCDF4 import Dataset

fh=Dataset('IMPACTS_PSD_Database.nc')
#fh=Dataset("IMPACTS_PSDs.nc")
import numpy as np
PSDs=fh["Nc"][:]
Z=fh["Z_sim"][:]
iwc=fh["iwc"][:]
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

#stop
import tensorflow as tf
import tensorflow.keras as keras
Keras = keras.backend


def iwc_model(n1,n2):
    input1 = keras.layers.Input(shape=[n1])
    z = keras.layers.Dense(40, activation="relu")(input1)
    z = keras.layers.Dropout(0.2)(z,training=True)
    z = keras.layers.Dense(40, activation="relu")(z)
    z = keras.layers.Dropout(0.2)(z,training=True)
    output = keras.layers.Dense(n2)(z)
    model = keras.models.Model(
        inputs=[input1], outputs=[output])
    return model

def iwc_modelM(n1,n2):
    input1 = keras.layers.Input(shape=[n1])
    z = keras.layers.Dense(8, activation="relu")(input1)
    z = keras.layers.Dropout(0.1) (z)
    z = keras.layers.Dense(16, activation="relu")(z)
    z = keras.layers.Dropout(0.1)(z)
    z = keras.layers.Dense(16, activation="relu")(z)
    z = keras.layers.Dropout(0.1)(z)
    #z = keras.layers.Dense(6, activation="relu")(z)
    z = keras.layers.Dense(n2)(z)
    #output = K.abs(z)
    output=z
    model = keras.models.Model(
        inputs=[input1], outputs=[output])
    return model

itrain=1
if itrain==1:
    model=iwc_modelM(3,4)
    model.compile(optimizer=tf.keras.optimizers.Adam(),  \
                  loss='mae',\
                  metrics=[tf.keras.metrics.MeanSquaredError()])
    
    history = model.fit(y_train[:,:], X_train, batch_size=32,epochs=30,
                        validation_data=(y_test[:,:], X_test))

    model.save("DFR_InvModel.h5")
else:
    model=tf.keras.models.load_model("DFR_InvModel.h5")
    
import pickle
pickle.dump({"InputScaler":scalerX,"OutputScaler":scalerY},open("iwcScalers.pklz","wb"))

stop
#A=inv(K.T*inv(Se)*K+inv(Sa))*(K.T*inv(Se)*K)

#stop
#Nw=4^4/pi/rhow*lwc/
from sklearn.cluster import  MiniBatchKMeans
random_state=0
import pickle
zmin=20

az=np.nonzero((Z[:,0]-zmin)*(Z[:,0]-(zmin+1))<0)
xinput=Ys[az[0],:]
#yn=model.predict(Y_test)

xEns=[]
for i in range(30):
    print(i)
    yn=model.predict(xinput)
    xEns.append(yn)
