import numpy as np
from netCDF4 import Dataset

fh=Dataset("scattering/ice-self-similar-aggregates_035.50-GHz_scat.nc")
fhKu=Dataset("scattering/ice-self-similar-aggregates_013.80-GHz_scat.nc")
fhW=Dataset("scattering/ice-self-similar-aggregates_094.00-GHz_scat.nc")
dmax_rh=fh["maximum_dimension"][:]
mass_rh=fh["mass"][:,:]
bscat_rh=fh["bscat"][:,:]*4*np.pi
bscatKu_rh=fhKu["bscat"][:,:]*4*np.pi
bscatW_rh=fhW["bscat"][:,:]*4*np.pi

#stop
fL=open('PSDs/ess238-sup-0002-supinfo.tex').readlines()
typeL=[]
dataL=[]
for f1 in fL:
    ls=f1.split()
    typeL.append(ls[0])
    dataL.append([float(v) for v in ls[1:]])
#0-elwp,1-mass,2-Dmax,3-radius of gyr, 4-Axis ratio, 6-Rimed fract, 
#6-sigma_X, 7-ext_X, 8-sigma_Ku, 9-ext_Ku, 10-sigma_Ka, 11-ext_Ka,

dint=np.array([2.00e+01, 4.00e+01, 6.00e+01, 8.00e+01, 1.00e+02, 1.25e+02,
               1.50e+02, 2.00e+02, 2.50e+02, 3.00e+02, 3.50e+02, 4.00e+02,
               4.75e+02, 5.50e+02, 6.25e+02, 7.00e+02, 8.00e+02, 9.00e+02,
               1.00e+03, 1.20e+03, 1.40e+03, 1.60e+03, 1.80e+03, 2.20e+03,
               2.60e+03, 3.00e+03, 3.40e+03, 3.80e+03, 4.20e+03, 4.60e+03,
               5.00e+03, 6.00e+03, 7.00e+03, 8.00e+03, 9.00e+03, 1.00e+04,
               1.20e+04, 1.40e+04, 1.60e+04, 1.80e+04, 2.00e+04, 2.50e+04,
               3.00e+04])
dbins=0.5*(dint[1:]+dint[:-1])      
massB=0.0061*(1e-4*dbins)**2.05
vols=np.pi*(dbins*1e-6)**3/6.
rhos=massB*1e-3/vols
fract=fh['fraction'][:]
ibin=0
sigmaRH=np.zeros((42,3),float)
massRH=np.zeros((42),float)
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
    return (qback*np.pi/4*(1e-3*ds)**2)

f=0.2
for mass1, rho1 in zip(massB,rhos):
    fract1=rho1/930.0
    a1=np.nonzero((dint[ibin+1]-dmax_rh*1e6)*(dint[ibin]-dmax_rh*1e6)<0)
    mass_int=mass_rh[:,a1[0]].mean(axis=-1)*1e3
    ind=np.argmin(abs(mass1-mass_int))
    if ind<47:
        ind+=0
    print(dbins[ibin],ind,fract[ind], mass1*1e-3/mass_rh[ind,a1[0]].mean(axis=-1))
    massRH[ibin]=mass_rh[ind,a1[0]].mean(axis=-1)
    sigmaKa_mie=getsigma_mie(refr,wl[1],rho,mass_rh[ind,a1[0]].mean())
    sigmaKu_mie=getsigma_mie(refr,wl[0],rho,mass_rh[ind,a1[0]].mean())
    sigmaW_mie=getsigma_mie(refr,wl[2],rho,mass_rh[ind,a1[0]].mean())
    print(sigmaKa_mie,bscat_rh[-2,ind,a1[0]].mean())

    sigmaRH[ibin,0]=f*bscatKu_rh[-2,ind,a1[0]].mean(axis=-1)+(1-f)*sigmaKu_mie
    sigmaRH[ibin,1]=f*bscat_rh[-2,ind,a1[0]].mean(axis=-1)+(1-f)*sigmaKa_mie
    sigmaRH[ibin,2]=f*bscatW_rh[-2,ind,a1[0]].mean(axis=-1)+(1-f)*sigmaW_mie
    ibin+=1
import pickle
pickle.dump({"mass":massRH,"sigma":sigmaRH},open("sigmaRH_ob.pklz","wb"))
stop

dataL=np.array(dataL)

import matplotlib.pyplot as plt

plt.scatter(dataL[:,1],dataL[:,2])
typeL=np.array(typeL)
a=np.nonzero(typeL=='A')
for i in range(2,3):
    b=np.nonzero(abs(dataL[a[0],0]-i*0.5)<1e-2)
    if len(b[0])>0:
        print(i*0.5,len(b[0]))
        print(dataL[a[0][b],6].mean(axis=0),i)
        coeffs=np.polyfit(np.log10(dataL[a[0][b],2]*100),np.log10(1e3*dataL[a[0][b],1]),1)
        print('Dmax=%5.f'%(dataL[a[0][b],2].max()*1e6),coeffs[0],10**coeffs[1])
        cL=[]
        for i in range(42):
            c=np.nonzero((dataL[a[0][b],2]*1e6-dint[i])*\
                (dataL[a[0][b],2]*1e6-dint[i+1])<0)
            cL.append(len(c[0]))
        print(cL)

mass_l=dataL[:,1]
dmax_l=dataL[:,2]
bscat_l=dataL[:,-4]
plt.scatter(mass_l[a[0][b]],dmax_l[a[0][b]])

plt.figure()
for i in range(24):
    plt.scatter(np.log(mass_rh[i,:]),np.log(bscat_rh[-2,i,:]),s=2)
#plt.scatter(np.log(mass_l[a[0][b]]),np.log(bscat_l[a[0][b]]),s=1)
plt.ylim(-30,-5)
plt.xlim(-25,5)

stop
lines=open('Chase_et_al_2021_NN-master/base_df_DDA.csv').readlines()

vL=[]
for l in lines[1:]:
    v=np.array([float(v) for v in l.split(',')[1:]])
    vL.append(v)

vL=np.array(vL)
#grep bhmie ~/JGRStudy/*py
from bhmief import bhmie
import pytmatrix.refractive
wl=[pytmatrix.refractive.wl_Ku,pytmatrix.refractive.wl_Ka,\
    pytmatrix.refractive.wl_W]

import pytmatrix.refractive as refr
#mw1= pytmatrix.refractive.m_w_10C[wl[0]]
rho=0.2
mi=refr.mi(wl[1],rho)
nang=10
sigmaka_mieL=[]
massL_ab=[]
for i in a[0][b]:
    ds=1e3*(6*mass_l[i]/(1e3*rho*np.pi))**(1/3)
    x=ds*np.pi/wl[1]
    #print(ds,x)
    s1,s2,qext,qsca,qback,gsca=bhmie(x,mi,nang)
    sigmaka_mieL.append(qback*np.pi/4*(1e-3*ds)**2)
    massL_ab.append(mass_l[i])

from sklearn.neighbors import NearestNeighbors
X=np.zeros((25*100,2),float)
sigma_Ka=np.zeros((25*100),float)
for i in range(25):
    X[i*100:i*100+100,0]=np.log(mass_rh[i,:])
    X[i*100:i*100+100,1]=np.log(dmax_rh[:])
    sigma_Ka[i*100:i*100+100]=np.log(bscat_rh[-1,i,:])
    
nbrs = NearestNeighbors(n_neighbors=2, algorithm='ball_tree').fit(X)
distances, indices = nbrs.kneighbors(np.log(dataL[:,1:3]))
#plt.scatter(np.log(vL[:,1]),np.log(vL[:,-2]))
#plt.scatter(np.log(np.array(massL_ab)),np.log(np.array(sigmaka_mieL)))
