import numpy as np
from netCDF4 import Dataset

fh=Dataset("scattering/ice-self-similar-aggregates_needles_035.50-GHz_scat.nc")
fhKu=Dataset("scattering/ice-self-similar-aggregates_needles_013.80-GHz_scat.nc")
fhW=Dataset("scattering/ice-self-similar-aggregates_needles_094.00-GHz_scat.nc")
dmax_RG=fh["maximum_dimension"][:]
mass_RG=fh["mass"][:,:]
bscatKa_RG=fh["bscat"][:,:]*4*np.pi
bscatKu_RG=fhKu["bscat"][:,:]*4*np.pi
bscatW_RG=fhW["bscat"][:,:]*4*np.pi
kextKa_RG=fh["ext"][:,:]
kextKu_RG=fhKu["ext"][:,:]
kextW_RG=fhW["ext"][:,:]
#print(fh)
#stop


dint=np.array([2.00e+01, 4.00e+01, 6.00e+01, 8.00e+01, 1.00e+02, 1.25e+02,
               1.50e+02, 2.00e+02, 2.50e+02, 3.00e+02, 3.50e+02, 4.00e+02,
               4.75e+02, 5.50e+02, 6.25e+02, 7.00e+02, 8.00e+02, 9.00e+02,
               1.00e+03, 1.20e+03, 1.40e+03, 1.60e+03, 1.80e+03, 2.20e+03,
               2.60e+03, 3.00e+03, 3.40e+03, 3.80e+03, 4.20e+03, 4.60e+03,
               5.00e+03, 6.00e+03, 7.00e+03, 8.00e+03, 9.00e+03, 1.00e+04,
               1.20e+04, 1.40e+04, 1.60e+04, 1.80e+04, 2.00e+04, 2.50e+04,
               3.00e+04])

#stop
fract=fh['fraction'][:]
ibin=0
sigmaKu_RGT=np.zeros((42,95),float)
sigmaKa_RGT=np.zeros((42,95),float)
sigmaW_RGT=np.zeros((42,95),float)
kextKu_RGT=np.zeros((42,95),float)
kextKa_RGT=np.zeros((42,95),float)
kextW_RGT=np.zeros((42,95),float)
mass_RGT=np.zeros((42,95),float)
vol_av=np.zeros((42,95),float)
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
    return (qback*np.pi/4*(1e-3*ds)**2), (qext*np.pi/4*(1e-3*ds)**2)


for i in range(42):
    a1=np.nonzero((dint[i+1]-dmax_RG*1e6)*(dint[i]-dmax_RG*1e6)<0)
    print(len(a1[0]))
    for j in range(95):
        mass_RGT[i,j]=mass_RG[j,a1[0]].mean(axis=-1)
        sigmaKu_RGT[i,j]=bscatKu_RG[-2,j,a1[0]].mean(axis=-1)
        sigmaKa_RGT[i,j]=bscatKa_RG[-2,j,a1[0]].mean(axis=-1)
        sigmaW_RGT[i,j]=bscatW_RG[-2,j,a1[0]].mean(axis=-1)
        kextKu_RGT[i,j]=kextKu_RG[-2,j,a1[0]].mean(axis=-1)
        kextKa_RGT[i,j]=kextKa_RG[-2,j,a1[0]].mean(axis=-1)
        kextW_RGT[i,j]=kextW_RG[-2,j,a1[0]].mean(axis=-1)
        vol_av[i,j]=(dmax_RG[a1[0]]**3).mean()
        mass1=mass_RGT[i,j]
        dmax1=dmax_RG[a1[0]].mean()
        rho1=mass1/(np.pi*dmax1**3/6)
        rho1=rho1*1e-3
        rho1=max(rho1,0.01)
        sigmaKu_mie,kextKu_mie=getsigma_mie(refr,wl[0],rho1,mass1)
        sigmaKa_mie,kextKa_mie=getsigma_mie(refr,wl[1],rho1,mass1)
        sigmaW_mie,kextW_mie=getsigma_mie(refr,wl[2],rho1,mass1)
        #kext_ku_av[i,:]=kextKu_mie
        #kext_ka_av[i,:]=kextKa_mie
        #kext_w_av[i,:]=kextW_mie
        rn=np.random.random()
        rn=0.5
        sigmaKu_RGT[i,j]=rn*sigmaKu_RGT[i,j]+(1-rn)*sigmaKu_mie
        sigmaKa_RGT[i,j]=rn*sigmaKa_RGT[i,j]+(1-rn)*sigmaKa_mie
        sigmaW_RGT[i,j]=rn*sigmaW_RGT[i,j]+(1-rn)*sigmaW_mie
    #sigmaKu_mie=getsigma_mie(refr,wl[0],rho,mass_rh[ind,a1[0]].mean())
    #sigmaW_mie=getsigma_mie(refr,wl[2],rho,mass_rh[ind,a1[0]].mean())
    #print(sigmaKa_mie,bscat_rh[-2,ind,a1[0]].mean())

    #sigmaRH[ibin,0]=f*bscatKu_rh[-2,ind,a1[0]].mean(axis=-1)+(1-f)*sigmaKu_mie
    #sigmaRH[ibin,1]=f*bscat_rh[-2,ind,a1[0]].mean(axis=-1)+(1-f)*sigmaKa_mie
    #sigmaRH[ibin,2]=f*bscatW_rh[-2,ind,a1[0]].mean(axis=-1)+(1-f)*sigmaW_mie
 
