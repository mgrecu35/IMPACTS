from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap


from pyresample import geometry, utils, image, kd_tree
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.colors as col

fh2=Dataset("slp.2001.nc",'r')
lon=fh2.variables['lon'][:]
lat=fh2.variables['lat'][:]
slp=fh2.variables['slp'][:,:,:]
lon[lon>180]-=360.
from numpy import *

lon2=-179.75+arange(720)*0.5
lat2=-89.75+arange(360)*0.5

area_id = 'ease_nh'
description = 'Arctic EASE grid'
proj_id = 'ease_nh'
x_size = 50
y_size = 50

m2 = Basemap(projection='npstere',boundinglat=30,lon_0=0,resolution='l')

#stop
area_extent = (-5326849.0625,-5326849.0625,5326849.0625,5326849.0625)
#area_extent = (m2.xmin ,m2.ymin,m2.xmax,m2.ymax)
shape=(50,50)
proj_dict = {'a': '6371228.0', 'units': 'm', 'lon_0': '0', \
             'proj': 'npstere', 'lat_0': '90'}
#proj_dict={'proj': 'stere', 'a': 6370997.0, 'units': 'm', 'lon_0': 0.0, 'lat_ts': 90.0, 'lat_0': 90.0, 'x_0': 3414207.0022622114, 'y_0': 3414207.0022622123}
center=(0,0)
radius=5000000
area_extent = (-5326849.0625,-5326849.0625,5326849.0625,5326849.0625)
proj_dict = {'a': '6371228.0', 'units': 'm', 'lon_0': '0', \
             'proj': 'laea', 'lat_0': '90'}
area_def = geometry.AreaDefinition(area_id, description, proj_id, \
                                   proj_dict, x_size, y_size, area_extent)

#area_def = geometry.AreaDefinition.from_circle(area_id, \
#                                               proj_dict, center, radius,
#                                              proj_id='npstere',\
#                                              shape=shape,
#                                               area_extent=area_extent)
print(area_def)
import numpy as np


lons,lats=np.meshgrid(lon,lat)
lons2,lats2=np.meshgrid(lon2,lat2)
grid_def = geometry.SwathDefinition(lons=lons, lats=lats)
grid_def2 = geometry.SwathDefinition(lons=lons2, lats=lats2)
wf = lambda r: 1 - r/350000.0
wf2 = lambda r: 1 - r/60000.0
resultL=[]
timeL=[]
import datetime
m = Basemap(projection='npstere',boundinglat=30,lon_0=0,resolution='l')
rtime=datetime.datetime(1800,1,1)

import glob
fs=glob.glob("slp.2*.nc")
fs=sorted(fs)
fsr=glob.glob("prate.sfc.g*2*.nc")
fsr=sorted(fsr)
fsr=glob.glob("/media/grecu/ExtraDrive1/IMERG/precipA*.nc")
fsr=sorted(fsr)

prateL=[]
for f1,f2 in zip(fs[:18],fsr[:18]):
    fh2=Dataset(f1,'r')
    slp=fh2.variables['slp'][:,:,:]
    fh3=Dataset(f2,'r')
    prate=fh3.variables['precip'][:,:,:]
    prate[:,:,:270]=-999
    t=fh2.variables['time'][:]
    a=range(59)
    a=list(a)
    nt=slp.shape[0]
    a.extend(range(nt-31,nt))
    
    for k,k1 in zip(a,range(60)):
        result  = kd_tree.resample_custom(grid_def, slp[k,:,:],
                                      area_def, \
                                      radius_of_influence=250000, \
                                      weight_funcs=wf)
        prate2 = kd_tree.resample_custom(grid_def2, prate[k1,:,:].T,
                                         area_def, \
                                         radius_of_influence=50000, \
                                         weight_funcs=wf)
        prate2=ma.array(prate2,mask=prate2<0)
        resultL.append(result.flatten())
        prateL.append(prate2)
        ctime=rtime+datetime.timedelta(hours=t[k])
        timeL.append(ctime)
        if 1==-1:
            plt.figure()
            m.drawcoastlines()
            m.drawparallels(np.arange(-80.,81.,20.))
            m.drawmeridians(np.arange(-180.,181.,20.))
            xx,yy=area_def.get_lonlats()
            x2,y2=m(xx,yy)
            cs=plt.contourf(x2,y2,prate2,linewidths=2)
            plt.colorbar()
            #stop

def condMean(precip):
    nt,nx,ny=precip.shape
    precipm=zeros((nx,ny),float)-99
    for i in range(nx):
        for j in range(ny):
            a=nonzero(precip[:,i,j]>=-1)
            if len(a[0])>0:
                precipm[i,j]=precip[a[0],i,j].mean()
    return precipm

from sklearn.cluster import KMeans
from sklearn import cluster
from sklearn.neighbors import kneighbors_graph
from numpy import *
resultL=array(resultL)
prateL=array(prateL)
resultm=resultL.mean(axis=0)
results=resultL.std(axis=0)
for k in range(resultm.shape[0]):
    resultL[:,k]-=resultm[k]
    resultL[:,k]/=results[k]
nc=35
kmeans=KMeans(n_clusters=35,random_state=0)
kmeans.fit(resultL)
#lab12_10=kmeans.labels_[1328]
#a12_10=nonzero(kmeans.labels_==lab12_10)
#print(array(timeL)[a12_10])
prateM=condMean(prateL)
anomL=[]
for lab in range(35):
    a=nonzero(kmeans.labels_==lab)
    plt.figure()
    m.drawcoastlines()
    m.drawparallels(np.arange(-80.,81.,20.))
    m.drawmeridians(np.arange(-180.,181.,20.))
    xx,yy=area_def.get_lonlats()
    x2,y2=m(xx,yy)
    #cs=plt.contourf(x2,y2,(resultL[a[0],:].mean(axis=0)*results+resultm).reshape(50,50))
    pratem=prateL[a[0],:,:].mean(axis=0)
    #pratem=ma.array(pratem,mask=pratem<0.000005)
    anom=condMean(prateL[a[0],:,:])
    anomm=ma.array(anom/2.-prateM/2,mask=anom<0)
    prateMm=ma.array(prateM/2.,mask=prateM<0)
    cs=plt.pcolormesh(x2,y2,anomm,cmap='RdBu',vmin=-20,vmax=20)
    anomL.append((anomm-prateMm).sum())
    plt.savefig('precipAnomClass%2.2i.png'%lab)
    plt.colorbar()

plt.figure()
plt.bar(arange(35),anomL)
plt.xlabel("Class")
plt.ylabel("Accumulation")
plt.title("Precipitation Anomalies")
plt.savefig('precipAnom.png')
stop
ExClass=array([10, 14, 21, 24])
ExClasses2=array([10, 14, 18, 21, 24])

X=resultL
connectivity = kneighbors_graph(
        X, n_neighbors=25, include_self=False)
    # make connectivity symmetric
connectivity = 0.5 * (connectivity + connectivity.T)
ward = cluster.AgglomerativeClustering(
    n_clusters=25, linkage='ward',
    connectivity=connectivity)
params={'n_clusters':25}   
spectral = cluster.SpectralClustering(
    n_clusters=params['n_clusters'], eigen_solver='arpack',
    affinity="nearest_neighbors")
spectral.fit(X)
for lab in range(25):
    a=nonzero(spectral.labels_==lab)
    plt.figure()
    m.drawcoastlines()
    m.drawparallels(np.arange(-80.,81.,20.))
    m.drawmeridians(np.arange(-180.,181.,20.))
    xx,yy=area_def.get_lonlats()
    x2,y2=m(xx,yy)
    cs=plt.contourf(x2,y2,(resultL[a[0],:].mean(axis=0)*results+resultm).reshape(50,50))
    plt.colorbar()

