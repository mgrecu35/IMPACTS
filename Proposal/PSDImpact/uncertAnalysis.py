import pickle
import tensorflow as tf
import tensorflow.keras as keras
from tensorflow.keras.layers import Lambda

K = keras.backend
d=pickle.load(open("iwcScalers.pklz","rb"))

def custom_loss_y(y_true, y_pred):
    #eps = tf.Variable(1.e-9)
    y_pred= y_pred+1e-2
    y_div=Lambda(lambda inputs: inputs[0] / inputs[1])([y_true, y_pred])
    y_div2= K.square(y_div)
    y_pred2= K.square(y_pred)
    loss = (K.log(y_pred2)+y_div2)
    loss = tf.math.reduce_mean(loss,axis=-1)

    #loss = loss.mean(axis=-1)
    return loss
#{"InputScaler":scalerX,"OutputScaler":scalerY,"varX":\
#["iwc","Dm","logNw","rho"],"varY":["zKu","DRF1","DFR2"
scalerX=d["InputScaler"]
scalerZ=d["OutputScaler"]

model_y=tf.keras.models.load_model("DFR_IWC_Inv_Model.h5")
model_sy=tf.keras.models.load_model("DFR_IWC_Inv_Model_Sigma.h5",compile=False)

from netCDF4 import Dataset

fh=Dataset("IMPACTS_PSD_Database.nc")

import numpy as np
PSDs=fh["Nc"][:]
Z=fh["Z_sim"][:]
iwc=fh["iwc"][:]
Nc=fh["Nc"][:]
Dm=fh['Dm'][:]
Nw=4**4/np.pi*iwc/(0.1*Dm)**4/1e6
rho=fh['rho'][:]

dfr1=Z[:,0]-Z[:,1]
dfr2=Z[:,1]-Z[:,2]
Y=np.array([Z[:,0],dfr1,dfr2]).T
X=np.array([iwc,Dm,np.log10(Nw),np.log10(rho)]).T
Xs=scalerX.transform(X)
#scalerZ.fit(Y)
Ys=scalerZ.transform(Y)
sigmaIWC=np.zeros((60,60),float)
from sklearn.neighbors import RadiusNeighborsRegressor, NearestNeighbors
neigh = RadiusNeighborsRegressor(radius=0.05,weights='distance')
neigh2 = NearestNeighbors(n_neighbors=5, radius=0.1)
neigh.fit(Ys, Xs)
neigh2.fit(Ys, Xs)
stdMap=np.zeros((25,30,30,4),float)-999.
meanMap=np.zeros((25,30,30,4),float)-999.
for i in range(10,35):
    z1=i-0.5
    z2=(i+0.5)
    a=np.nonzero((Z[:,0]-z1)*(Z[:,0]-z2)<0)
    dfr1_min=dfr1[a].min()
    dfr1_max=dfr1[a].max()
    dfr2_min=dfr2[a].min()
    dfr2_max=dfr2[a].max()
    print(z1,dfr1_min,dfr1_max,dfr2_min,dfr2_max)
    dfr1i=np.linspace(dfr1_min,dfr1_max,50)
    dfr2i=np.linspace(dfr2_min,dfr2_max,50)
    xsL=[]
    ijL=[]
    for j in np.arange(0,15,0.5):
        for k in np.arange(0,15,0.5):
            b=np.nonzero((dfr1[a]-(j-0.5))*(dfr1[a]-(j+0.5))<0)
            c=np.nonzero((dfr2[a][b]-(k-0.5))*(dfr2[a][b]-(k+0.5))<0)
            if len(c[0])>0:
                ijL.append((int(2*j),int(2*k)))
                x=np.array([Z[:,0][a][b][c],dfr1[a][b][c],dfr2[a][b][c]]).T
                xs=scalerZ.transform(x)
                xsL.append(xs)
    

    stdL=[]
    for j,xs in enumerate(xsL):
        dist,ind=neigh.radius_neighbors(xs)
        dist2,ind2=neigh2.kneighbors(xs)
        indT=[]
        for ind1 in ind2:
            indT.extend(list(ind1))
        #stop
        cov1=np.cov(Xs[indT,:].T)
        std1=Xs[indT,:].std(axis=0)
        mean1=Xs[indT,:].mean(axis=0)
        #stdL.append(std1)
    
        stdMap[i-10,ijL[j][0],ijL[j][1],:]=std1*scalerX.scale_
        meanMap[i-10,ijL[j][0],ijL[j][1],:]=scalerX.inverse_transform(mean1)
        #if mean1[0]<1e-5:
        #    stop
    #stop
#plt.pcolormesh(np.arange(30)*0.5,np.arange(30)*0.5,meanMap,norm=matplotlib.colors.LogNorm(),cmap='jet')
import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits.axes_grid1 import make_axes_locatable


import xarray as xr
stdMapX=xr.DataArray(stdMap)
meanMapX=xr.DataArray(meanMap)

d=xr.Dataset({"stdMap":stdMapX,"meanMap":meanMapX})
d.to_netcdf("iceRet&Uncert.nc")
