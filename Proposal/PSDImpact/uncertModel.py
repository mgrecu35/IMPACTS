from sklearn import preprocessing
from netCDF4 import Dataset


import tensorflow as tf
import tensorflow
import tensorflow.keras as keras
from tensorflow.keras.layers import Lambda

K = keras.backend


def iwc_model(n1,n2):
    input1 = keras.layers.Input(shape=[n1])
    z = keras.layers.Dense(8, activation="relu")(input1)
    z = keras.layers.Dropout(0.01) (z)
    z = keras.layers.Dense(8, activation="relu")(z)
    z = keras.layers.Dropout(0.01)(z)
    z = keras.layers.Dense(8, activation="relu")(z)
    z = keras.layers.Dropout(0.01)(z)
    #z = keras.layers.Dense(6, activation="relu")(z)
    z = keras.layers.Dense(n2)(z)
    output = K.abs(z)
    model = keras.models.Model(
        inputs=[input1], outputs=[output])
    return model


model_y=iwc_model(1,2)
model_sy=iwc_model(1,2)
import numpy as np
x=np.random.random(15000)[:,np.newaxis]



y=np.zeros((x.shape[0],2),float)
y[:,0]=0.5*x[:,0]+0.5*(abs(x[:,0])+0.01)*np.random.randn(15000)
y[:,1]=0.5*x[:,0]+0.3*(abs(x[:,0])+0.01)*np.random.randn(15000)

#y=np.array(,\
#            0.5*x+0.3*(abs(x[:,0])+0.01)*np.random.randn(15000)]).T

x_t=np.random.random(15000)[:,np.newaxis]

y_t=np.zeros((x.shape[0],2),float)
y_t[:,0]=0.5*x_t[:,0]+0.5*(abs(x_t[:,0])+0.01)*np.random.randn(15000)
y_t[:,1]=0.5*x_t[:,0]+0.3*(abs(x_t[:,0])+0.01)*np.random.randn(15000)


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
    y_pred= y_pred+1e-2
    y_div=Lambda(lambda inputs: inputs[0] / inputs[1])([y_true, y_pred])
    y_div2= K.square(y_div)
    y_pred2= K.square(y_pred)
    loss = K.log(y_pred2)+y_div2
    loss = tf.math.reduce_mean(loss,axis=-1)
    return loss


model_y.compile(optimizer=tf.keras.optimizers.Adam(),  \
              loss='mae',\
              metrics=[tf.keras.metrics.MeanSquaredError()])


model_sy.compile(optimizer=tf.keras.optimizers.Adam(),  \
              loss=custom_loss_y,\
              metrics=[tf.keras.metrics.MeanSquaredError()])




weight=y_t[:,0]*0+1

y_error=y.copy()
y_error[:,0]-=0.5*x[:,0]
y_error[:,1]-=0.5*x[:,0]
for it in range(20):
    history_y = model_y.fit(x_t[:,:], y_t, sample_weight=weight,batch_size=32,epochs=5,
                        validation_data=(x, y))
    yp=model_y.predict(x_t)
    
    history_sy = model_sy.fit(x_t[:,:], y_t-yp, batch_size=32,epochs=5,
                              validation_data=(x,y_error))
    yp_s=model_sy.predict(x_t)
    
    #weight=1./(yp_s+1e-3)
    
import matplotlib.pyplot as plt


