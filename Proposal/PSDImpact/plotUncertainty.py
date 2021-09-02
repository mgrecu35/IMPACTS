#d=xr.Dataset({"stdMap":stdMapX,"meanMap":meanMapX})
from netCDF4 import Dataset
import numpy as np
fh=Dataset("iceRet&Uncert.nc")

stdMap=fh["stdMap"][:,:]
meanMap=fh["meanMap"][:,:]
fh.close()
izL=[5]
import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits.axes_grid1 import make_axes_locatable

meanMap[:,:,:,2]+=5
meanMapm=np.ma.array(meanMap,mask=meanMap<-99)
stdMapm=np.ma.array(stdMap,mask=meanMap<-99)
izL=[5,10,15]
cbar=['g/m^3','mm','','kg/m^3']
tit=['IWC[g/m$^3$]','D$_m$[mm]','Nw[log10(m$^{-3}$mm$^{-1}$)]','log($\\rho$) [log(kg/m$^3$)]']
for i in izL:
    plt.figure(figsize=(8,8))
    plt.suptitle('Mean estimates, Z(Ku)=%6.2f dBZ'%(i+10))
    for k in range(4):
        fig=plt.subplot(2,2,(k+1))
        ax = plt.gca()
        ax.set_aspect('equal')
        if k==0:
            im=plt.pcolormesh(np.arange(30)*0.5,np.arange(30)*0.5,meanMapm[i,:,:,k],\
                              cmap='jet',norm=matplotlib.colors.LogNorm())
        else:
            im=plt.pcolormesh(np.arange(30)*0.5,np.arange(30)*0.5,meanMapm[i,:,:,k],\
                              cmap='jet')
        plt.xlim(0,10)
        plt.ylim(0,10)
        if k==0 or k==2:
            plt.ylabel('DFR(Ku,Ka)')
        if k==2 or k==3:
            plt.xlabel('DFR(Ka,W)')
        plt.title(tit[k])
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        
        c=plt.colorbar(im, cax=cax)
        #c.ax.set_title(cbar[k])
      
    plt.tight_layout()
    plt.savefig('retrMeanZ%5.2f.png'%(i+10))

for i in izL:
    plt.figure(figsize=(8,8))
    plt.suptitle('Uncertainty estimates, Z(Ku)=%6.2f dBZ'%(i+10))
    for k in range(4):
        fig=plt.subplot(2,2,(k+1))
        ax = plt.gca()
        ax.set_aspect('equal')
        if k==0:
            im=plt.pcolormesh(np.arange(30)*0.5,np.arange(30)*0.5,stdMapm[i,:,:,k],\
                              cmap='jet',norm=matplotlib.colors.LogNorm())
        else:
            im=plt.pcolormesh(np.arange(30)*0.5,np.arange(30)*0.5,stdMapm[i,:,:,k],\
                       cmap='jet')
        plt.xlim(0,10)
        plt.ylim(0,10)
        if k==0 or k==2:
            plt.ylabel('DFR(Ku,Ka)')
        if k==2 or k==3:
            plt.xlabel('DFR(Ka,W)')
        plt.title(tit[k])
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        
        c=plt.colorbar(im, cax=cax)
        #c.ax.set_title(cbar[k])
      
    plt.tight_layout()
    plt.savefig('retrUncertaintiesZ%5.2f.png'%(i+10))
