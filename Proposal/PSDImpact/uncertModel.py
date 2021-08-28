from sklearn import preprocessing
from netCDF4 import Dataset


import tensorflow as tf
import tensorflow
import tensorflow.keras as keras
from tensorflow.keras.layers import Lambda

K = keras.backend


def iwc_model(n1,n2):
    input1 = keras.layers.Input(shape=[n1])
    z = keras.layers.Dense(6, activation="relu")(input1)
    z = keras.layers.Dropout(0.1) (z)
    z = keras.layers.Dense(6, activation="relu")(z)
    z = keras.layers.Dropout(0.1)(z)
    #z = keras.layers.Dense(6, activation="relu")(z)
    z = keras.layers.Dense(1)(z)
    output = K.abs(z)
    model = keras.models.Model(
        inputs=[input1], outputs=[output])
    return model


model=iwc_model(1,1)
import numpy as np
x=np.random.random(15000)[:,np.newaxis]
y=1*(abs(x)+0.01)*np.random.randn(15000)[:,np.newaxis]

x_t=np.random.random(15000)[:,np.newaxis]
y_t=1*(abs(x_t)+0.01)*np.random.randn(15000)[:,np.newaxis]


def custom_loss(y_true, y_pred):
    #eps = tf.Variable(1.e-9)
    y_pred= y_pred+1e-2
    y_div=Lambda(lambda inputs: inputs[0] / inputs[1])([y_true, y_pred])
    y_div2= K.square(y_div)
    y_pred2= K.square(y_pred)
    loss = K.log(y_pred2)+y_div2
    return loss



model.compile(optimizer=tf.keras.optimizers.Adam(),  \
              loss=custom_loss,\
              metrics=[tf.keras.metrics.MeanSquaredError()])

history = model.fit(x_t[:,:], y_t, batch_size=32,epochs=100,
                        validation_data=(x, y))




