from numpy import *
fileList=['../PSDs/IMPACTS_MergedHorizontal-P3_20200118_sizedistributions_v01.nc',
          '../PSDs/IMPACTS_MergedHorizontal-P3_20200207_sizedistributions_v01.nc',
          '../PSDs/IMPACTS_MergedHorizontal-P3_20200125_sizedistributions_v01.nc',
          '../PSDs/IMPACTS_MergedHorizontal-P3_20200213_sizedistributions_v01.nc',
          '../PSDs/IMPACTS_MergedHorizontal-P3_20200201_sizedistributions_v01.nc',
          '../PSDs/IMPACTS_MergedHorizontal-P3_20200218_sizedistributions_v01.nc',
          '../PSDs/IMPACTS_MergedHorizontal-P3_20200205_sizedistributions_v01.nc',
          '../PSDs/IMPACTS_MergedHorizontal-P3_20200220_sizedistributions_v01.nc',
          '../PSDs/IMPACTS_MergedVertical-P3_20200207_sizedistributions_v01.nc',
          '../PSDs/IMPACTS_MergedVertical-P3_20200220_sizedistributions_v01.nc',
          '../PSDs/IMPACTS_MergedVertical-P3_20200213_sizedistributions_v01.nc',
          '../PSDs/IMPACTS_MergedVertical-P3_20200218_sizedistributions_v01.nc']
from netCDF4 import Dataset
fh=Dataset(fileList[0])

dbins=fh['CONCENTRATION'].bin_midpoints
dint=fh['CONCENTRATION'].bin_endpoints
#stop
dbins=array(dbins)
dint=array(dint)
ds=dint[1:]-dint[0:-1]

mass=0.0061*(1e-4*dbins)**2.05
massInt=0.0061*(1e-4*dint)**2.05

mass=0.0042*(1e-4*dbins)**2.04
massInt=0.0042*(1e-4*dint)**2.04
import pickle


dbinsM=(mass*6/pi/0.92)**(1./3.)*10.
dbinsI=(mass*6/pi)**(1./3.)*10.

import pytmatrix.refractive
wl=[pytmatrix.refractive.wl_Ku,pytmatrix.refractive.wl_Ka,pytmatrix.refractive.wl_W]

from readDDArho import *

mass=0.0042*(1e-4*dbins)**2.04
sback95=(wl[2]**4/pi**5)*sigma_w_av*1e6
sback37=(wl[1]**4/pi**5)*sigma_ka_av*1e6
sback13=(wl[0]**4/pi**5)*sigma_ku_av*1e6
#ext95=0.5*(10**luo_Tables[0][1]+10**luo_Tables2[0][1])
#ext37=0.5*(10**luo_Tables[1][1]+10**luo_Tables2[1][1])
#ext13=0.5*(10**luo_Tables[2][1]+10**luo_Tables2[2][1])

#sca95=0.5*(10**luo_Tables[0][2]+10**luo_Tables2[0][2])
#sca37=0.5*(10**luo_Tables[1][2]+10**luo_Tables2[1][2])
#sca13=0.5*(10**luo_Tables[2][2]+10**luo_Tables2[2][2])

#gsca95=0.5*(luo_Tables[0][3]+luo_Tables2[0][3])
#gsca37=05.*(luo_Tables[1][3]+luo_Tables2[1][3])
#gsca13=0.5*(luo_Tables[2][3]+luo_Tables2[2][3])

#luo_Tables=pickle.load(open('luo_Graupel_Tables.pklz','rb'))
#luo_Tables2=pickle.load(open('luo_Graupel_Tables.pklz','rb'))

mi1= pytmatrix.refractive.mi(wl[0],0.92)
mi2= pytmatrix.refractive.mi(wl[1],0.92)
mi3= pytmatrix.refractive.mi(wl[2],0.92)
m=array([mi1,mi2,mi3])
K2r=abs((m**2-1)/(m**2+2))**2
K2s=((m**2-1)/(m**2+2))**2
#print K2r
K2=0.93

dbins1=0.5*dbinsI*1e3+0.5*dbins
xW1=4*3.142/wl[0]*(0.30*dbins1*1e-3)
xf1=(1+0.159*xW1**2)/(1+(0.159+1./3)*xW1**2+0.164*xW1**4);
xW2=4*3.142/wl[1]*(0.30*dbins1*1e-3)
xf2=(1+0.159*xW2**2)/(1+(0.159+1./3)*xW2**2+0.164*xW2**4);
xW3=4*3.142/wl[2]*(0.30*dbins1*1e-3)

xf3=(1+0.159*xW3**2)/(1+(0.159+1./3)*xW3**2+0.164*xW3**4);
#print sca37/ext37
#print gsca37


zSims=[]
datL=[]
import random
choiceL=range(45)
kextL=[]
vol_av=vol_av*np.pi/6
rhoL=[]
iwcL=[]
tempL=[]
ncL=[]
tempCL=[]
zL=[]
dmL=[]

from numba import jit
@jit(nopython=True)
def calcProp(Nbins,ds,mass_av,vol_av,dbinsM,K2r,xf1,xf2,xf3,K2,ik,gradZ,graddm,gradiwc,gradrho,dNbins):
    dm=sum(Nbins*ds*mass_av[:,ik]*dbinsM)/sum(Nbins*ds*mass_av[:,ik])
    iwc_ps=sum(Nbins*ds*mass_av[:,ik])*1e-3
    rho=sum(Nbins*ds*mass_av[:,ik])/sum(Nbins*ds*vol_av[:,ik])
    Z1=log10(sum(Nbins*ds*1e-6*(dbinsI)**6*xf1)*K2r[0]/K2)*10.
    Z2=log10(sum(Nbins*ds*1e-6*(dbinsI)**6*xf2)*K2r[1]/K2)*10.
    Z3=log10(sum(Nbins*ds*1e-6*(dbinsI)**6*xf3)*K2r[2]/K2)*10.
    r1=0.5*np.random.randn()
    r2=0.5*np.random.randn()
    r3=0.5*np.random.randn()
    Z1_dda=log10(sum(Nbins*1e-6*ds*sback13[:,ik]/K2))*10.
    Z3_dda=log10(sum(Nbins*1e-6*ds*sback95[:,ik]/K2))*10.
    Z2_dda=log10(sum(Nbins*1e-6*ds*sback37[:,ik]/K2))*10.
    
    gradZ[42,0]=log10(sum(Nbins*1e-6*ds*sback13[:,ik+1]/K2))*10.
    gradZ[42,2]=log10(sum(Nbins*1e-6*ds*sback95[:,ik+1]/K2))*10.
    gradZ[42,1]=log10(sum(Nbins*1e-6*ds*sback37[:,ik+1]/K2))*10.

    gradiwc[42]=sum((Nbins+dNbins)*ds*mass_av[:,ik+1])*1e-3
    graddm[42]=sum((Nbins+dNbins)*ds*mass_av[:,ik+1]*dbinsM)/sum((Nbins+dNbins)*ds*mass_av[:,ik+1])
    gradrho[42]=sum((Nbins+dNbins)*ds*mass_av[:,ik+1])/sum((Nbins+dNbins)*ds*vol_av[:,ik+1])
        
    for i in range(-42):
        dNbins*=0
        dNbins[i]=Nbins[i]*0.1
        gradZ[i,0]=log10(sum((Nbins+dNbins)*1e-6*ds*sback13[:,ik]/K2))*10.
        gradZ[i,1]=log10(sum((Nbins+dNbins)*1e-6*ds*sback37[:,ik]/K2))*10.
        gradZ[i,2]=log10(sum((Nbins+dNbins)*1e-6*ds*sback95[:,ik]/K2))*10.
        gradiwc[i]=sum((Nbins+dNbins)*ds*mass_av[:,ik])*1e-3
        graddm[i]=sum((Nbins+dNbins)*ds*mass_av[:,ik]*dbinsM)/sum((Nbins+dNbins)*ds*mass_av[:,ik])
        gradrho[i]=sum((Nbins+dNbins)*ds*mass_av[:,ik])/sum((Nbins+dNbins)*ds*vol_av[:,ik])
    kext35=4.34*sum(Nbins*ds*1e-6*kext_ka_av[:,ik])*1e3
    kext95=4.34*sum(Nbins*ds*1e-6*kext_w_av[:,ik])*1e3
    kext13=4.34*sum(Nbins*ds*1e-6*kext_ku_av[:,ik])*1e3
    return Z1,Z2,Z3,Z1_dda,Z2_dda,Z3_dda,rho,dm,iwc_ps,kext13,kext35,kext95,gradZ,gradiwc,graddm,gradrho

gradZL=[]
graddmL=[]
gradiwcL=[]
gradrhoL=[]
for fname in fileList[:]:
    fh=Dataset(fname)
    date=fname.split("_")[2]
    fnameT="../PSDs/temp_%s.nc"%date
    fhT=Dataset(fnameT)
    tempC=fhT["tempC"][:]
    #if '0205' not in fname:
    #    continue
    c=fh['CONCENTRATION'][:,:].T
    iwc_a=fh['IWC'][:]
    iL=0
    print(fname,c.shape[0])
    for i in range(c.shape[0]):
        Nbins=c[i,:].data
        if Nbins.min()<0 or Nbins.sum()<1e-3:
            continue
        if tempC[i]<-998:
            continue
        if tempC[i]>0:
            continue
        #iwc_psL=[]
        #for ik in range(45):
        #    iwc_ps=sum(Nbins*ds*mass_av[:,ik])*1e-3
        #    iwc_psL.append(iwc_ps)
        #iwcL1=[iwc_a[i]]
        #iwcL1.extend(iwc_psL)
        #iwcL.append(iwcL1)
        
   
        for it in range(2):
            
            #ik=np.argmin(abs(iwc_a[i]-np.array(iwc_ps)))
            #a=np.nonzero((iwc_a[i]-2*np.array(iwc_ps))*(iwc_a[i]-0.1*np.array(iwc_ps)))
            #if len(a[0])>0:
            ik=random.choice(range(44))
            #ik=6
            #Nbins=array(dat1[3:])
                
            gradZ=np.zeros((43,3),float)
            gradiwc=np.zeros((43),float)
            graddm=np.zeros((43),float)
            gradrho=np.zeros((43),float)
            dNbins=np.zeros((42),float)
            Z1,Z2,Z3,Z1_dda,Z2_dda,Z3_dda,rho,\
                dm,iwc_ps,kext13,kext35,kext95,\
                gradZ,gradiwc,graddm,gradrho=\
                    calcProp(Nbins,ds,mass_av,vol_av,dbinsM,K2r,xf1,xf2,xf3,K2,ik,\
                             gradZ,graddm,gradiwc,gradrho,dNbins)

           
            if iwc_a[i]>0:
                ncL.append(Nbins)
                tempCL.append(tempC[i])
                zL.append([Z1_dda,Z2_dda,Z3_dda])
                kextL.append([kext13,kext35,kext95])
                dmL.append(dm)
                rhoL.append(rho)
                iwcL.append(iwc_ps)
                tempL.append(tempC[i])
                zSims.append([Z1_dda,Z2_dda,Z3_dda,dm,Z1,Z2,Z3,iwc_a[i],iwc_ps])
                gradZL.append(gradZ)
                graddmL.append(graddm)
                gradiwcL.append(gradiwc)
                gradrhoL.append(gradrho)
                iL+=1
            #if Z1_dda>10:
                #print(Z1_dda,Z1_dda1)
            #    stop
            if abs(iwc_ps/(iwc_a[i]+1e-5)-1)>1.50:
                continue

            #kextL.append([kext13,kext35,kext95])

            #rhoL.append(rho)
            #
    print(iL)
            #stop

import numpy as np
zSims=np.array(zSims)
import xarray as xr
ncLx=xr.DataArray(ncL)
iwcLx=xr.DataArray(iwcL)
tempLx=xr.DataArray(tempCL)
zLx=xr.DataArray(zL,dims=['dim_0','n3'])
kextLx=xr.DataArray(kextL,dims=['dim_0','n3'])
rhoLx=xr.DataArray(rhoL)
dmLx=xr.DataArray(dmL)
iwcLx=xr.DataArray(iwcL)
gradZLx=xr.DataArray(gradZL,dims=['dim_0','n43','n3'])
gradiwcLx=xr.DataArray(gradiwcL,dims=['dim_0','n43'])
graddmLx=xr.DataArray(graddmL,dims=['dim_0','n43'])
gradrhoLx=xr.DataArray(gradrhoL,dims=['dim_0','n43'])

ds=xr.Dataset({"Nc":ncLx,"tempC":tempLx,"Z_sim":zLx,"kext":kextLx,"rho":rhoLx,"Dm":dmLx,"iwc":iwcLx,"diwc":gradiwcLx,\
              "ddm":graddmLx,"drho":gradrhoLx,"dZ":gradZLx})

ds=xr.Dataset({"Nc":ncLx,"tempC":tempLx,"Z_sim":zLx,"kext":kextLx,"rho":rhoLx,"Dm":dmLx,\
               "iwc":iwcLx})

comp = dict(zlib=True, complevel=5)
encoding = {var: comp for var in ds.variables}
#ds.to_netcdf(filename, encoding=encoding)

ds.to_netcdf("IMPACTS_PSD_Database.nc",encoding=encoding)
    
def readAH(fname):
    lines=open(fname,'r').readlines()
    datL=[]
    datL2=[]
    for l in lines[:]:
        #print l
        try:
            dat1=[float(v) for v in l.split()]
        except:
            dat1=[]
            pass
        if len(dat1)==44:
            if dat1[2]>0.0001 and dat1[2]<9.999e3:
                Nbins=array(dat1[3:])
                dm=sum(Nbins*ds*mass*dbinsM)/sum(Nbins*ds*mass)
                Z1=log10(sum(Nbins*ds*1e-6*(dbinsI)**6*xf1)*K2r[0]/K2)*10.
                Z2=log10(sum(Nbins*ds*1e-6*(dbinsI)**6*xf2)*K2r[1]/K2)*10.
                Z3=log10(sum(Nbins*ds*1e-6*(dbinsI)**6*xf3)*K2r[2]/K2)*10.
                Z3=log10(sum(Nbins*ds*1e-6*sback95/K2))*10.
                Z1=log10(sum(Nbins*ds*1e-6*sback13/K2))*10.
                Z2=log10(sum(Nbins*ds*1e-6*sback37/K2))*10.
                kext35=4.34*sum(Nbins*ds*1e-6*ext37)*1e-3
                kext95=4.34*sum(Nbins*ds*1e-6*ext95)*1e-3
                kext13=4.34*sum(Nbins*ds*1e-6*ext13)*1e-3
                ksca13=4.34*sum(Nbins*ds*1e-6*sca13)*1e-3
                ksca35=4.34*sum(Nbins*ds*1e-6*sca37)*1e-3
                ksca95=4.34*sum(Nbins*ds*1e-6*sca95)*1e-3
                _gsca13=4.34*sum(Nbins*ds*1e-6*sca13*gsca13)*1e-3
                _gsca35=4.34*sum(Nbins*ds*1e-6*sca37*gsca37)*1e-3
                _gsca95=4.34*sum(Nbins*ds*1e-6*sca95*gsca95)*1e-3
                #if Z1>250:
                #    print Z3, kext95, kext35
                Nw=4**4/pi*dat1[2]/(0.1*dm)**4/1e6
                #if abs(dm-1)<1e-3:
                    #print sum(Nbins*ds*mass*dbinsM)
                    #print sum(Nbins*ds*mass), dat1[2]
                dat2=[kext13,ksca13/kext13,ksca35/kext35,ksca95/kext95,\
                      _gsca13/ksca13,_gsca35/ksca35,_gsca95/ksca95]
                #print dat2
            else:
                Z1=-99
                Z2=-99
                Z3=-99
                kext35=0
                kext95=0
                kext13=0
                dm=0
                Nw=-99
                dat2=[kext13,0.,0.,0.,0.,0.,0.]
                #print [Nw,kext35,kext95,Z1,Z2,Z3,dm]
            dat1.extend([Nw,kext35,kext95,Z1,Z2,Z3,dm])
            
            datL.append(dat1)
            datL2.append(dat2)
            #print dat1
    #stop
    return array(datL),array(datL2)

def readAH2(fname):
    lines=open(fname,'r').readlines()
    datL=[]
    nbinsL=[]
    for l in lines[:]:
        #print l
        try:
            dat1=[float(v) for v in l.split()]
        except:
            dat1=[]
            pass
        if len(dat1)==44:
            if dat1[2]>0.0001 and dat1[2]<9.999e3:
                Nbins=array(dat1[3:])
                dm=sum(Nbins*ds*mass*dbinsM)/sum(Nbins*ds*mass)
                Z1=log10(sum(Nbins*ds*1e-6*(dbinsI)**6*xf1)*K2r[0]/K2)*10.
                Z2=log10(sum(Nbins*ds*1e-6*(dbinsI)**6*xf2)*K2r[1]/K2)*10.
                Z3=log10(sum(Nbins*ds*1e-6*(dbinsI)**6*xf3)*K2r[2]/K2)*10.
                Z3=log10(sum(Nbins*ds*1e-6*sback95/K2))*10.
                Z1=log10(sum(Nbins*ds*1e-6*sback13/K2))*10.
                Z2=log10(sum(Nbins*ds*1e-6*sback37/K2))*10.
                kext35=4.34*sum(Nbins*ds*1e-6*ext37)*1e-3
                kext95=4.34*sum(Nbins*ds*1e-6*ext95)*1e-3
                #if Z1>250:
                #    print Z3, kext95, kext35
                Nw=4**4/pi*dat1[2]/(0.1*dm)**4/1e6
                #if abs(dm-1)<1e-3:
                    #print sum(Nbins*ds*mass*dbinsM)
                    #print sum(Nbins*ds*mass), dat1[2]
            else:
                Z1=-99
                Z2=-99
                Z3=-99
                kext35=0
                kext95=0
                dm=0
                Nw=-99
                Nbins=array(dat1[3:])*0-99
                #print [Nw,kext35,kext95,Z1,Z2,Z3,dm]
            dat1.extend([Nw,kext35,kext95,Z1,Z2,Z3,dm])

            datL.append(dat1)
            nbinsL.append(Nbins)
            #print dat1
    return array(datL), array(nbinsL), sback95, sback13, sback37, K2, ds

#datAH=readAH()
import pickle
pickle.dump({"zSims":zSims,"kextL":kextL,"rhoL":rhoL,"iwcL":iwcL,\
             "tempC":tempL},open("simulatedPropTempSorted.pklz","wb"))

#pickle.dump({"zSims":zSims,"kextL":kextL,"rhoL":rhoL,"iwcL":iwcL,\
#             "tempC":tempL},open("simulatedPropSortedTemp_2.pklz","wb"))
