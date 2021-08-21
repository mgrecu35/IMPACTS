import numpy as np
from sklearn.neighbors import KNeighborsRegressor


def knnRC(zKu,zKa,zW,Dm,iwc,nne):
    n=zKu.shape[0]
    X=np.zeros((n,3),float)
    X[:,0]=(zKu-12)/16.
    X[:,1]=(zKu-zKa-2)/2.5
    X[:,2]=(zW-2)/7.0
    Y=Dm
    Yice=iwc
    Yzw=zW
    Yattw=zW*0
    neigh = KNeighborsRegressor(n_neighbors=nne,weights='distance')
    neigh_iwc = KNeighborsRegressor(n_neighbors=nne,weights='distance')
    r=np.random.random(n)
    a=np.nonzero(r<0.5)
    b=np.nonzero(r>0.5)
    X_train=X[a[0],:]
    X_test=X[b[0],:]
    y_train=Y[a[0]]
    y_test=Y[b[0]]
    yice_train=Yice[a[0]]
    yice_test=Yice[b[0]]

    neigh.fit(X_train[:,0:2], y_train)
    neigh_iwc.fit(X_train[:,0:2], yice_train)
    yp=neigh.predict(X_test[:,0:2])
    yp_iwc=neigh_iwc.predict(X_test[:,0:2])
    print(np.corrcoef(yp_iwc,yice_test))
    return neigh_iwc
