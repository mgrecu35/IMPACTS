import matplotlib.pyplot as plt
import numpy as np
lines=open('base_df_DDA.csv').readlines()

vL=[]
for l in lines[1:]:
    v=np.array([float(v) for v in l.split(',')[1:]])
    vL.append(v)


vL=np.array(vL)

#plt.scatter(vL[:,1],np.log10(vL[:,2]))
vL[:,1]*=1.2

dbins=[50.0,   70.0,    90.0,   112.5,   137.5,   175.0,   225.0,   275.0,   325.0,   375.0,   437.5,   512.5,   587.5,   662.5,   750.0,   850.0,   950.0,  1100.0,  1300.0,  1500.0,  1700.0,  2000.0,  2400.0,  2800.0,  3200.0,  3600.0,  4000.0,  4400.0,  4800.0,  5500.0,  6500.0,  7500.0,  8500.0,  9500.0, 11000.0, 13000.0, 15000.0, 17000.0, 19000.0, 22500.0,  27500.0]

dbins=[50.0, 70.0, 90.0, 112.5, 137.5, 175.0, 225.0, 275.0, 325.0, 375.0, 437.5, 512.5, 587.5, 662.5, 750.0, 850.0, 950.0, 1100.0, 1300.0, 1500.0, 1700.0, 2000.0, 2400.0, 2800.0, 3200.0, 3600.0, 4000.0, 4400.0, 4800.0, 5500.0, 6500.0, 7500.0, 8500.0, 9500.0, 11000.0, 13000.0, 15000.0, 17000.0, 19000.0, 22500.0, 27500.0]
dint=[40.0, 60.0, 80.0, 100.0, 125.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0, 475.0, 550.0, 625.0, 700.0, 800.0, 900.0, 1000.0, 1200.0, 1400.0, 1600.0, 1800.0, 2200.0, 2600.0, 3000.0, 3400.0, 3800.0, 4200.0, 4600.0, 5000.0, 6000.0, 7000.0, 8000.0, 9000.0, 10000.0, 12000.0, 14000.0, 16000.0, 18000.0, 20000.0, 25000.0, 30000.0]
print(len(dint))

dbins=np.array(dbins)
dint=0.5*(dbins[0:-1]+dbins[1:])
dintL=[]
dintL.append(2*dbins[0]-dint[0])
for d in dint:
    dintL.append(d)
dintL.append(2*dbins[-1]-dint[-1])
a=np.zeros((42,42),float)
b=np.zeros((42),float)
a[0,0]=1
b[0]=40
for i in range(41):
    a[i+1,i]=0.5
    a[i+1,i+1]=0.5
    b[i+1]=dbins[i]
dint=np.dot(np.linalg.pinv(a),b)
ds=dint[1:]-dint[0:-1]

massB=0.0061*(1e-4*dbins)**2.05
dbins=[   30. ,    50. ,    70. ,    90. ,   112.5,   137.5,   175. ,
         225. ,   275. ,   325. ,   375. ,   437.5,   512.5,   587.5,
         662.5,   750. ,   850. ,   950. ,  1100. ,  1300. ,  1500. ,
        1700. ,  2000. ,  2400. ,  2800. ,  3200. ,  3600. ,  4000. ,
        4400. ,  4800. ,  5500. ,  6500. ,  7500. ,  8500. ,  9500. ,
       11000. , 13000. , 15000. , 17000. , 19000. , 22500. , 27500. ]
dint=[2.00e+01, 4.00e+01, 6.00e+01, 8.00e+01, 1.00e+02, 1.25e+02,
       1.50e+02, 2.00e+02, 2.50e+02, 3.00e+02, 3.50e+02, 4.00e+02,
       4.75e+02, 5.50e+02, 6.25e+02, 7.00e+02, 8.00e+02, 9.00e+02,
       1.00e+03, 1.20e+03, 1.40e+03, 1.60e+03, 1.80e+03, 2.20e+03,
       2.60e+03, 3.00e+03, 3.40e+03, 3.80e+03, 4.20e+03, 4.60e+03,
       5.00e+03, 6.00e+03, 7.00e+03, 8.00e+03, 9.00e+03, 1.00e+04,
       1.20e+04, 1.40e+04, 1.60e+04, 1.80e+04, 2.00e+04, 2.50e+04,
       3.00e+04]

dbins=[50.0, 70.0, 90.0, 112.5, 137.5, 175.0, 225.0, 275.0, 325.0, 375.0, 437.5, 512.5, 587.5, 662.5, 750.0, 850.0, 950.0, 1100.0, 1300.0, 1500.0, 1700.0, 2000.0, 2400.0, 2800.0, 3200.0, 3600.0, 4000.0, 4400.0, 4800.0, 5500.0, 6500.0, 7500.0, 8500.0, 9500.0, 11000.0, 13000.0, 15000.0, 17000.0, 19000.0, 22500.0, 27500.0]
dint=[40.0, 60.0, 80.0, 100.0, 125.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0, 475.0, 550.0, 625.0, 700.0, 800.0, 900.0, 1000.0, 1200.0, 1400.0, 1600.0, 1800.0, 2200.0, 2600.0, 3000.0, 3400.0, 3800.0, 4200.0, 4600.0, 5000.0, 6000.0, 7000.0, 8000.0, 9000.0, 10000.0, 12000.0, 14000.0, 16000.0, 18000.0, 20000.0, 25000.0, 30000.0]
print(len(dint))

dbins=np.array(dbins)
dint=np.array(dint)
ds=dint[1:]-dint[0:-1]

mass=0.0061*(1e-4*dbins)**2.05
massInt=0.0061*(1e-4*dint)**2.05

#from scattering import *

mass=0.0042*(1e-4*dbins)**2.04
massInt=0.0042*(1e-4*dint)**2.04

sigma_ku_av=np.zeros((41,45),float)
sigma_ka_av=np.zeros((41,45),float)
sigma_w_av=np.zeros((41,45),float)
kext_ku_av=np.zeros((41,45),float)
kext_ka_av=np.zeros((41,45),float)
kext_w_av=np.zeros((41,45),float)
mass_av=np.zeros((41,45),float)
vol_av=np.zeros((41,45),float)

import random
from bhmief import bhmie
import pytmatrix.refractive
wl=[pytmatrix.refractive.wl_Ku,pytmatrix.refractive.wl_Ka,\
    pytmatrix.refractive.wl_W]

import pytmatrix.refractive as refr
rho=0.1
nang=22
def getsigma_mie(refr,wl1,rho,mass1):
    mi=refr.mi(wl1,rho)
    ds=1e3*(6*mass1/(1e3*rho*np.pi))**(1/3)
    x=ds*np.pi/wl1
    s1,s2,qext,qsca,qback,gsca=bhmie(x,mi,nang)
    return (qback*np.pi/4*(1e-3*ds)**2),(qext*np.pi/4*(1e-3*ds)**2)

coeffs=np.polyfit(np.log(vL[:,0]),np.log(vL[:,1]),1)

residA=np.log(vL[:,1])-coeffs[0]*np.log(vL[:,0])-coeffs[1]
ir=np.random.random(45)
ir=np.arange(45)/46.
import pickle
for i in range(41):
    #a=np.nonzero((vL[:,1]*1e3-massInt[i])*(vL[:,1]*1e3-massInt[i+1])<0)
    a=np.nonzero((vL[:,0]*1e6-dint[i])*(vL[:,0]*1e6-dint[i+1])<0)
    mbound=0.002938030918684807*(1e-4*dint[i+1])**1.9
    b=np.nonzero(vL[:,1][a]*1e3>1.0*mbound)
    print(len(a[0]),len(b[0]))
    #print((vL[a[0],1]*1e3).std()/(vL[a[0],1]*1e3).mean(),(vL[a[0],1]*1e3).mean())
    rho1L=[]
    indk=np.argsort(residA[a])
    #indk=range(len(a[0]))
    n2=int(len(a[0])/100.0)
    indk=a[0][indk]
    #indk=a[0]
    ir=np.arange(45)/45.
    print(n2)
    if n2<1:
        n2=1
    for k in range(45):
        #ik=a[0][indk[int(ir[k]*len(a[0]))]]
        ik=int(ir[k]*len(a[0])-n2)
        if ik<0:
            ik=0
        sigma_ku_av[i,k]=vL[indk[ik:ik+n2],2].mean()
        sigma_ka_av[i,k]=vL[indk[ik:ik+n2],3].mean()
        sigma_w_av[i,k]=vL[indk[ik:ik+n2],4].mean()
        mass_av[i,k]=vL[indk[ik:ik+n2],1].mean()
        vol_av[i,k]=(vL[indk[ik:ik+n2],0]**3).mean()
        mass1=mass_av[i,k]
        dmax1=(vL[indk[ik:ik+n2],0]).mean()
        rho1=mass1/(np.pi*dmax1**3/6)
        #rho1L.append(rho1)
        rho1=rho1*1e-3
        rho1=max(rho1,0.01)
        sigmaKu_mie,kextKu_mie=getsigma_mie(refr,wl[0],rho1,mass1)
        sigmaKa_mie,kextKa_mie=getsigma_mie(refr,wl[1],rho1,mass1)
        sigmaW_mie,kextW_mie=getsigma_mie(refr,wl[2],rho1,mass1)
        kext_ku_av[i,k]=kextKu_mie
        kext_ka_av[i,k]=kextKa_mie
        kext_w_av[i,k]=kextW_mie
        rn=np.random.random()
        rn=0.85+0.15*rn
 
        #print(ik)
    #ind=np.argsort(mass_av[i,:])
    #sigma_ku_av[i,:]=sigma_ku_av[i,ind]
    #sigma_ka_av[i,:]=sigma_ka_av[i,ind]
    #sigma_w_av[i,:]=sigma_w_av[i,ind]
    #mass_av[i,:]=mass_av[i,ind]
    #print((vL[a[0],1]*1e3).std()/(vL[a[0],1]*1e3).mean(),(vL[a[0],1]*1e3).mean(),\
    #      mass_sorted[20]/vL[a[0],1].mean(),mass_sorted[380]/vL[a[0],1].mean())
    #stop

