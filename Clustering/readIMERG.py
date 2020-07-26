import glob
from netCDF4 import Dataset
fs=glob.glob("/media/grecu/ExtraDrive1/IMERG/precipAfrica*nc")
fs=sorted(fs)
fh=Dataset(fs[0])

precip=fh['precip'][0,:,:]
from numpy import *


i=0
precipL=[]
precipL1d=[]
lon=-179.75+arange(720)*0.5
lat=-89.75+arange(360)*0.5

for f in fs[:9]:
    fh=Dataset(f,'r')
    precip=fh['precip'][:,320:400,179:220]
    for i in range(precip.shape[0]):
        precipL.append(precip[i,:,:])
        precipL1d.append(precip[i,:,:].flatten())


import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits.basemap import Basemap
m = Basemap(projection='mill',llcrnrlat=-0,urcrnrlat=20,\
            llcrnrlon=-20,urcrnrlon=20,resolution='i')

lon2d,lat2d=meshgrid(lon[320:401],lat[179:221])
x,y=m(lon2d,lat2d)
m.drawcoastlines()
m.drawcountries()
m.pcolormesh(x,y,array(precipL).mean(axis=0).T,cmap='jet')


precipL1d=array(precipL1d)
precipL1dm=precipL1d.mean(axis=0)

for k in range(precipL1dm.shape[0]):
    precipL1d[:,k]-=precipL1dm[k]

from sklearn.cluster import MiniBatchKMeans, KMeans
nc1=5
nc2=4

k_means = KMeans(init='k-means++', n_clusters=nc1*nc2,\
                random_state=0)

k_means.fit(precipL1d)

precipL=array(precipL)
it=0


fig=plt.figure(figsize=(12,8))
plt.suptitle('Precipitation Anomaly Classes')

for i in range(nc1):
    for j in range(nc2):
        ax=plt.subplot(nc1,nc2,it+1)
        m.drawcoastlines()
        a=nonzero(k_means.labels_==it)
        pm=precipL[a[0],:,:].mean(axis=0)-precipL1dm.reshape(80,41)
        im=plt.pcolormesh(x,y,(pm/48).T,vmax=5,vmin=-5,cmap='RdBu_r')
        plt.title('%2.2i'%(it+1))
        it+=1
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.875, 0.15, 0.025, 0.7])
c=fig.colorbar(im, cax=cbar_ax)
c.ax.set_title('mm/day')
plt.savefig('precipAnomalyClasses.png')

it=0
Freq=[]
nt=precipL1d.shape[0]
for i in range(nc1):
    for j in range(nc2):
        a=nonzero(k_means.labels_==it)
        Freq.append(len(a[0])/nt)
        it+=1

matplotlib.rcParams.update({'font.size': 13})


fig, ax = plt.subplots()
ax.bar((arange(20)+1).astype(int),Freq)
ax.set_xticks([1,5,9,13,17])
plt.xlabel('Class',fontsize=13)
plt.ylabel('Relative Frequency',fontsize=13)
plt.tight_layout()
plt.savefig('classFrequency.png')


transMat=zeros((20,20),float)
for k in range(nt-1):
    i0=k_means.labels_[k]
    j0=k_means.labels_[k+1]
    transMat[i0,j0]+=1

for k in range(20):
    transMat[k,:]/=sum(transMat[k,:])

plt.figure(figsize=(6,6))

transMatm=ma.array(transMat,mask=transMat<0.00001)

fig=plt.figure()
ax = plt.subplot(111)
ax.set_aspect('equal')
plt.pcolormesh(arange(0,21)+0.5,arange(0,21)+0.5,transMatm,cmap='jet')
plt.colorbar()
plt.savefig('transitionMatrix.png')


