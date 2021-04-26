import scattering as sc
import numpy as np
import matplotlib.pyplot as plt

mu=2.0
freqKu=13.8
freqKa=35.5
scatTables=sc.scatTables
swc=np.arange(800)*0.005+0.01
nws=10**(3.5+np.arange(800)*0.0)
zKu,dmKu=sc.calcSnowProp(swc,nws,scatTables,mu,freqKu)
zKa,dmKa=sc.calcSnowProp(swc,nws,scatTables,mu,freqKa)


from sklearn.linear_model import LinearRegression
#X = np.array([zL,nwsL]).T
#y = np.array(swcL)
#reg = LinearRegression().fit(X, y)
#print(reg.coef_[0],reg.coef_[1],reg.intercept_)
