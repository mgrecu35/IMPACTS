from numpy import *
fileList=['PSDs/IMPACTS_MergedHorizontal-P3_20200125_sizedistributions_v01.nc']
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
zL=[]
timeL=[]
t1=21.20
t2=22.35
t1=21.75
t2=22.1
t1=21.25
t2=21.35
t1,t2=21+10/60.,21.3#19/60.
latL=[]
lonL=[]
iwcL2=[]

zEns=[]
iwcEns=[]
dmL=[]
dmNCAR=[]
iwcNCAR=[]
rhoNCAR=[]
rhoEns=[]

NbinL=[]
for fname in fileList[0:1]:
    fh=Dataset(fname)
    date=fname.split("_")[2]
    fnameT="PSDs/temp_%s.nc"%date
    fhT=Dataset(fnameT)
    tempC=fhT["tempC"][:]
    time=fh["time"][:]/3600.
    #if '0205' not in fname:
    #    continue
    c=fh['CONCENTRATION'][:,:].T
    iwc_a=fh['IWC'][:]
    iL=0
    a=np.nonzero((time-t1)*(time-t2)<0)
    galt=fh["GALT"][:]
    latp=fh["LAT"][:]
    lonp=fh["LON"][:]
    print(fname,c.shape[0])
    mass_avP=mass_av[:,:].mean(axis=-1)
    vol_avP=vol_av[:,:].mean(axis=-1)
    for i in a[0]:
        Nbins=c[i,:].data
        #Nbins[0:10]/=100.
        if Nbins.min()<0 or Nbins.sum()<1e-3:
            continue
        if tempC[i]<-998:
            continue
        if tempC[i]>0:
            continue
        iwc_psL=[]
        ik=random.choice(range(45))
        #ik=6
        zL.append(galt[i])
        #Nbins=array(dat1[3:])
        dmAvg=sum(Nbins*ds*mass_avP*dbinsM)/sum(Nbins*ds*mass_avP)
        iwcAvg=sum(Nbins*ds*mass_avP)*1e-3
        rhoAvg=sum(Nbins*ds*mass_avP)/sum(Nbins*ds*vol_avP)
        Z1=log10(sum(Nbins*ds*1e-6*(dbinsI)**6*xf1)*K2r[0]/K2)*10.
        Z2=log10(sum(Nbins*ds*1e-6*(dbinsI)**6*xf2)*K2r[1]/K2)*10.
        Z3=log10(sum(Nbins*ds*1e-6*(dbinsI)**6*xf3)*K2r[2]/K2)*10.
        Z1_dda=log10(sum(Nbins*1e-6*ds*sback13[:,ik]/K2))*10.+0.5*np.random.randn()
        Z3_dda=log10(sum(Nbins*1e-6*ds*sback95[:,ik]/K2))*10.+0.5*np.random.randn()
        Z2_dda=log10(sum(Nbins*1e-6*ds*sback37[:,ik]/K2))*10.+0.5*np.random.randn()
        Z1_ddaL=[]
        Z2_ddaL=[]
        Z3_ddaL=[]
        dm1=[]
        rho1=[]
        for it in range(45):
            ik=it
            Z1_dda=log10(sum(Nbins*1e-6*ds*sback13[:,ik]/K2))*10.+0.5*np.random.randn()
            Z3_dda=log10(sum(Nbins*1e-6*ds*sback95[:,ik]/K2))*10.+0.5*np.random.randn()
            Z2_dda=log10(sum(Nbins*1e-6*ds*sback37[:,ik]/K2))*10.+0.5*np.random.randn()
            Z1_ddaL.append(Z1_dda)
            Z2_ddaL.append(Z2_dda)
            Z3_ddaL.append(Z3_dda)
            iwc_ps=sum(Nbins*ds*mass_av[:,ik])*1e-3
            iwc_psL.append(iwc_ps)
            dm=sum(Nbins*ds*mass_av[:,ik]*dbinsM)/sum(Nbins*ds*mass_av[:,ik])
            dm1.append(dm)
            rho=sum(Nbins*ds*mass_av[:,ik])/sum(Nbins*ds*vol_av[:,ik])
            rho1.append(rho)
        kext35=4.34*sum(Nbins*ds*1e-6*kext_ka_av[:,ik])*1e3
        kext95=4.34*sum(Nbins*ds*1e-6*kext_w_av[:,ik])*1e3
        kext13=4.34*sum(Nbins*ds*1e-6*kext_ku_av[:,ik])*1e3
        #if abs(iwc_ps/(iwc_a[i]+1e-5)-1)>0.20:
        #    continue
        iL+=1
        kextL.append([kext13,kext35,kext95])
        zSims.append([np.mean(Z1_ddaL),np.mean(Z2_ddaL),\
                      np.mean(Z3_ddaL),dm,Z1,Z2,Z3,iwc_a[i],iwc_ps])
        if zSims[-1][0]>25:
            NbinL.append(Nbins)
        zEns.append([Z1_ddaL,Z2_ddaL,Z3_ddaL])
        iwcEns.append(iwc_psL)
        rhoEns.append(rho1)
        rhoNCAR.append(rhoAvg)
        tempL.append(tempC[i])
        timeL.append(time[i])
        lonL.append(lonp[i])
        latL.append(latp[i])
        iwcL.append(np.mean(iwc_ps))
        iwcL2.append(iwc_a[i])
        dmL.append(dm1)
        dmNCAR.append(dmAvg)
        iwcNCAR.append(iwcAvg)
    print(iL)
            #stop

import numpy as np
zSims=np.array(zSims)
         
#stop

import pickle
pickle.dump({"zSims":zSims,"kextL":kextL,"rhoL":rhoL,"iwcL":iwcL,\
             "tempC":tempL},open("simulatedProp0205_22:00-22:06.pklz","wb"))


#dictZ={"zKu":zKu_slice[200:3500],"zKa":zKa_slice[200:3500],\
#       "zW":np.array(zW_L),"time":time[a[0][200:3500]],\
#       "rrange":rmean,"alt":alt[a[0][200:3500]]}
dictZ=pickle.load(open("z_obs_Jan25_21:10-21:19","rb"))

zKu=dictZ["zKu"]
zKa=dictZ["zKa"]
zW=dictZ["zW"]
alt=dictZ["alt"]
rrange=dictZ["rrange"]
timeE=dictZ["time"]
lon=dictZ["lon"]
lat=dictZ["lat"]
ae=np.nonzero((timeE-t1)*(timeE-t2)<0)

h1=alt[ae[0][0]]-rrange
ind=np.argmin(abs(h1-zL[0]))
plt.figure(figsize=(8,12))
plt.subplot(311)
zkuL=[]
zkaL=[]
zwL=[]
rhoEns=np.array(rhoEns)
rhoCombRet=[]
zKuRet=[]
zKaRet=[]
zWRet=[]
iwcRet=[]
i=0
zEns=np.array(zEns)
iwcEns=np.array(iwcEns)
x3L=[]
dmL=np.array(dmL)
dmCombRet=[]
wku=1.
wka=1.
indL=[]
for alt1,lon1,lat1 in zip(zL,lonL,latL):
    ind1=np.argmin((lon1-lon)**2+(lat1-lat)**2)
    ind2=np.argmin(abs(alt[ind1]-rrange-alt1))
    zKa[ind1,ind2]-=1
    indL.append(ind1)
    zkuL.append(zKu[ind1,ind2])
    zkaL.append(zKa[ind1,ind2])
    zwL.append(zW[ind1,ind2])
    wku=1.
    wka=1.
    if timeL[i]>21.277 and timeL[i]<21.287:
        wku=0
        wka=0
    if zKu[ind1,ind2]==zKu[ind1,ind2] and \
       zKa[ind1,ind2]==zKa[ind1,ind2] and \
       zW[ind1,ind2]==zW[ind1,ind2]:
        indSol=np.argsort(wku*(zkuL[-1]-zEns[i,0,:])**2+wka*(zkaL[-1]-zEns[i,1,:])**2+\
                          (zwL[-1]-zEns[i,2,:])**2)
        zKuRet.append(zEns[i,0,indSol[0:3]].mean())
        zKaRet.append(zEns[i,1,indSol[0:3]].mean())
        zWRet.append(zEns[i,2,indSol[0:3]].mean())
        iwcRet.append(iwcEns[i,indSol[0:3]].mean())
        dmCombRet.append(dmL[i,indSol[0:3]].mean())
        rhoCombRet.append(rhoEns[i,indSol[0:3]].mean())
        if timeL[i]>21.277 and timeL[i]<21.287:
            zkuL[-1]=zEns[i,0,indSol[0:3]].mean()
            zkaL[-1]=zEns[i,1,indSol[0:3]].mean()
           
        #zkuL[-1]=zEns[i,0,indSol[0:1]].mean()
        #zkaL[-1]=zEns[i,1,indSol[0:1]].mean()
        #zwL[-1]=zEns[i,2,indSol[0:1]].mean()
        x3L.append([(zkuL[-1]-12)/16.,(zkuL[-1]-zkaL[-1]-2)/2.5,(zwL[-1]-2)/7.])
    else:
        zKuRet.append(np.nan)
        zKaRet.append(np.nan)
        zWRet.append(np.nan)
        iwcRet.append(np.nan)
        dmCombRet.append(np.nan)
        rhoCombRet.append(np.nan)
        x3L.append([(-12-12)/16.,(0-2)/2.5,(-12-2)/7.])
    i+=1

    
d31=pickle.load(open("knnModelsIMPACTS_3_RG_Needles.pklz","rb"))
d3=pickle.load(open("knnModelsIMPACTS_3_DDA_sorted.pklz","rb"))
x3L=np.array(x3L)
iwc_2=d3["knn_iwc"].predict(x3L)
dm_2=d3["knn_dm"].predict(x3L)
rho_2=d3["knn_rho"].predict(x3L)
plt.plot(timeL,zkuL)
plt.plot(timeL,zkaL)
plt.plot(timeL,zwL)
plt.ylim(0,35)
plt.legend(['Ku','Ka','W'])
plt.ylabel('observed Z')
plt.xlim(t1,t2)

plt.subplot(312)
#plt.plot(timeL,zSims[:,0])
#plt.plot(timeL,zSims[:,1])
#plt.plot(timeL,zSims[:,2])
plt.plot(timeL,zKuRet[:])
plt.plot(timeL,zKaRet[:])
plt.plot(timeL,zWRet[:])
plt.legend(['Ku','Ka','W'])
plt.ylabel('simulated Z')
plt.ylim(0,35)
plt.xlim(t1,t2)
plt.subplot(313)
#plt.plot(timeL,iwcL)
plt.plot(timeL,iwcL2)
plt.plot(timeL,iwcRet)
plt.plot(timeL,iwc_2)
plt.legend(['NCAR','Combined','Retrieved'])

plt.ylabel('IWC (g/m^3)')
plt.xlabel('Time')

plt.xlim(t1,t2)
plt.savefig('combinedRadarProbes_IWC_0125_21:10-21:19.png')

plt.figure(figsize=(8,8))
fig1=plt.subplot(211)
plt.plot(timeL,dmNCAR)
plt.plot(timeL,dmCombRet)
plt.plot(timeL,dm_2)
plt.xlim(t1,t2)
plt.legend(['NCAR','Combined','Retrieved'])
plt.ylabel('Dm (mm)')
fig2=plt.subplot(212)
plt.semilogy(timeL,rhoNCAR)
plt.semilogy(timeL,rhoCombRet)
plt.semilogy(timeL,rho_2)
plt.ylabel("Bulk density (kg/m^3)")
plt.xlabel('Time')
plt.xlim(t1,t2)
plt.savefig('combinedRadarProbes_RhoDm_0125_21:10-21:19.png')
