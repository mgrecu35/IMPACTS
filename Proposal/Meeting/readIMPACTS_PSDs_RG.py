from numpy import *
fileList=['PSDs/IMPACTS_MergedHorizontal-P3_20200118_sizedistributions_v01.nc',
          'PSDs/IMPACTS_MergedHorizontal-P3_20200207_sizedistributions_v01.nc',
          'PSDs/IMPACTS_MergedHorizontal-P3_20200125_sizedistributions_v01.nc',
          'PSDs/IMPACTS_MergedHorizontal-P3_20200213_sizedistributions_v01.nc',
          'PSDs/IMPACTS_MergedHorizontal-P3_20200201_sizedistributions_v01.nc',
          'PSDs/IMPACTS_MergedHorizontal-P3_20200218_sizedistributions_v01.nc',
          'PSDs/IMPACTS_MergedHorizontal-P3_20200205_sizedistributions_v01.nc',
          'PSDs/IMPACTS_MergedHorizontal-P3_20200220_sizedistributions_v01.nc',
          'PSDs/IMPACTS_MergedVertical-P3_20200207_sizedistributions_v01.nc',
          'PSDs/IMPACTS_MergedVertical-P3_20200220_sizedistributions_v01.nc',
          'PSDs/IMPACTS_MergedVertical-P3_20200213_sizedistributions_v01.nc',
          'PSDs/IMPACTS_MergedVertical-P3_20200224_sizedistributions_v01.nc',
          'PSDs/IMPACTS_MergedVertical-P3_20200218_sizedistributions_v01.nc']
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

#
from read_RG_properties import mass_RGT, sigmaKu_RGT, sigmaKa_RGT, sigmaW_RGT,\
    kextKu_RGT,kextKa_RGT,kextW_RGT, vol_av
#from readDDA import sigma_ku_av, sigma_ka_av, sigma_w_av, mass_av
mass=0.0042*(1e-4*dbins)**2.04
sback95=(wl[2]**4/pi**5)*sigmaW_RGT*1e6
sback37=(wl[1]**4/pi**5)*sigmaKa_RGT*1e6
sback13=(wl[0]**4/pi**5)*sigmaKu_RGT*1e6
mass_av1=mass_RGT
#sback95[:,50:65]=(wl[2]**4/pi**5)*sigma_w_av*1e6
#sback13[:,50:65]=(wl[0]**4/pi**5)*sigma_ku_av*1e6
#sback37[:,50:65]=(wl[1]**4/pi**5)*sigma_ka_av*1e6
#mass_av1[:,50:65]=mass_av
mass_av=mass_av1

import numpy as np
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
choiceL=range(60)
kextL=[]

kext_ku_av=kextKu_RGT
kext_ka_av=kextKa_RGT
kext_w_av=kextW_RGT
rhoL=[]
for fname in fileList[:]:
    fh=Dataset(fname)
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
        #Nbins/=np.exp(-2*np.arange(42)/42)

        for it in range(1):
            iL+=1
            ik=random.choice(choiceL)
            #ik=6
            #Nbins=array(dat1[3:])
            dm=sum(Nbins*ds*mass_av[:,ik]*dbinsM)/sum(Nbins*ds*mass_av[:,ik])
            iwc_ps=sum(Nbins*ds*mass_av[:,ik])*1e-3
            rho=sum(Nbins*ds*mass_av[:,ik])/sum(Nbins*ds*vol_av[:,ik])
            Z1=log10(sum(Nbins*ds*1e-6*(dbinsI)**6*xf1)*K2r[0]/K2)*10.+np.random.randn()*1.0
            Z2=log10(sum(Nbins*ds*1e-6*(dbinsI)**6*xf2)*K2r[1]/K2)*10.+np.random.randn()*1.0
            Z3=log10(sum(Nbins*ds*1e-6*(dbinsI)**6*xf3)*K2r[2]/K2)*10.+np.random.randn()*1.0
            Z1_dda=log10(sum(Nbins*1e-6*ds*sback13[:,ik]/K2))*10.+np.random.randn()*1.0
            Z3_dda=log10(sum(Nbins*1e-6*ds*sback95[:,ik]/K2))*10.+np.random.randn()*1.0
            Z2_dda=log10(sum(Nbins*1e-6*ds*sback37[:,ik]/K2))*10.+np.random.randn()*1.0
            kext35=4.34*sum(Nbins*ds*1e-6*kext_ka_av[:,ik])*1e3
            kext95=4.34*sum(Nbins*ds*1e-6*kext_w_av[:,ik])*1e3
            kext13=4.34*sum(Nbins*ds*1e-6*kext_ku_av[:,ik])*1e3
            kextL.append([kext13,kext35,kext95])
            zSims.append([Z1_dda,Z2_dda,Z3_dda,dm,Z1,Z2,Z3,iwc_a[i],iwc_ps])
            rhoL.append(rho)
    print(iL)
            #stop

import numpy as np
zSims=np.array(zSims)
kextL=np.array(kextL)
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
pickle.dump({"zSims":zSims,"kextL":kextL,"rhoL":rhoL},open("simulatedPropRGNeedles.pklz","wb"))
