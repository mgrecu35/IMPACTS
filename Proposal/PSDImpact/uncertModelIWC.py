from sklearn import preprocessing
from netCDF4 import Dataset


import tensorflow as tf
import tensorflow
import tensorflow.keras as keras
from tensorflow.keras.layers import Lambda

K = keras.backend


def iwc_modelS(n1,n2):
    input1 = keras.layers.Input(shape=[n1])
    z = keras.layers.Dense(8, activation="relu")(input1)
    z = keras.layers.Dropout(0.1) (z)
    z = keras.layers.Dense(16, activation="relu")(z)
    z = keras.layers.Dropout(0.1)(z)
    z = keras.layers.Dense(16, activation="relu")(z)
    z = keras.layers.Dropout(0.1)(z)
    #z = keras.layers.Dense(6, activation="relu")(z)
    z = keras.layers.Dense(n2)(z)
    output = K.abs(z)
    #output=z
    model = keras.models.Model(
        inputs=[input1], outputs=[output])
    return model

def iwc_model(n1,n2):
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



import numpy as np
from netCDF4 import Dataset

fh=Dataset("IMPACTS_PSD_Database.nc")

import numpy as np
PSDs=fh["Nc"][:]
Z=fh["Z_sim"][:]
Z=Z+np.random.randn(Z.shape[0],3)
iwc=fh["iwc"][:]
Nc=fh["Nc"][:]
Dm=fh['Dm'][:]
Nw=4**4/np.pi*iwc/(0.1*Dm)**4/1e6
rho=fh['rho'][:]
X=np.array([iwc,Dm,np.log10(Nw),np.log10(rho)]).T
X=X[:,[0,1,2,3]]
Y=Z.copy()
Y[:,1]=Z[:,0]-Z[:,1]
Y[:,2]=Z[:,1]-Z[:,2]
from sklearn import preprocessing
scalerX  = preprocessing.StandardScaler()
scalerY  = preprocessing.StandardScaler()
scalerX.fit(X)
scalerY.fit(Y)
Xs=scalerX.transform(X)
Ys=scalerY.transform(Y)
n=Xs.shape[0]

model_y=iwc_model(3,4)
model_sy=iwc_modelS(3,4)


def custom_loss(y_true, y_pred):
    #eps = tf.Variable(1.e-9)
    y_pred= y_pred+1e-2
    y_div=Lambda(lambda inputs: inputs[0] / inputs[1])([y_true, y_pred])
    y_div2= K.square(y_div)
    y_pred2= K.square(y_pred)
    loss = K.log(y_pred2)+y_div2
    return loss


def custom_loss_y(y_true, y_pred):
    #eps = tf.Variable(1.e-9)
    y_pred= y_pred+1e-1
    y_div=Lambda(lambda inputs: inputs[0] / inputs[1])([y_true, y_pred])
    y_div2= K.square(y_div)
    y_pred2= K.square(y_pred)
    loss = (K.log(y_pred2)+y_div2)
    loss = tf.math.reduce_mean(loss,axis=-1)

    #loss = loss.mean(axis=-1)
    return loss



model_y.compile(optimizer=tf.keras.optimizers.Adam(),  \
              loss='mae',\
              metrics=[tf.keras.metrics.MeanSquaredError()])


model_sy.compile(optimizer=tf.keras.optimizers.Adam(),  \
              loss=custom_loss_y,\
              metrics=[tf.keras.metrics.MeanSquaredError()])


import numpy as np

r=np.random.random(n)
a=np.nonzero(r<0.5)
b=np.nonzero(r>0.5)
X_train=Xs[a[0],:]
X_test=Xs[b[0],:]
y_train=Ys[a[0],:]
y_test=Ys[b[0],:]

#weight=y_t*0+1

#for it in range(20):
itrain=1
import pickle
if itrain==1:
    history_y = model_y.fit(y_train[:,:], X_train, batch_size=32,epochs=50,
                            validation_data=(y_test, X_test))
    pickle.dump({"InputScaler":scalerX,"OutputScaler":scalerY,"varX":\
                 ["iwc","Dm","logNw","rho"],"varY":["zKu","DRF1","DFR2"]},\
                open("iwcScalers.pklz","wb"))

    model_y.save("DFR_IWC_Inv_Model.h5")
else:
    model_y=tf.keras.models.load_model("DFR_IWC_Inv_Model.h5")

itrain=1
yp=model_y.predict(y_train)
yp_t=model_y.predict(y_test)
#yps=model_sy.predict(y_train)
#loss=custom_loss_y(X_train-yp,yps)

#stop
if itrain==1:
    history_sy = model_sy.fit(y_train[:,:], X_train-yp, batch_size=32,epochs=50,
                          validation_data=(y_test, yp_t-X_test))
    model_sy.save("DFR_IWC_Inv_Model_Sigma.h5")

yp_s=model_sy.predict(y_test)
    
#weight=1./(yp_s+1e-3)
    
import matplotlib.pyplot as plt


