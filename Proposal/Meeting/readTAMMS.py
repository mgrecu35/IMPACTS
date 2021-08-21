import glob
files=glob.glob("TAMMS/*ict")
import xarray as xr
for f in []:#files:
    lines=open(f).readlines()
    print(lines[47])
    timeL=[]
    latL=[]
    lonL=[]
    altL=[]
    tempCL=[]
    for l in lines[48:]:
        ls=l.split(',')
        v=[float(v) for v in ls[:7]]
        if abs(v[0]-int(v[0]))<1e-4:
            timeL.append(v[0])
            latL.append(v[1])
            lonL.append(v[2])
            altL.append(v[3])
            tempCL.append(v[6])
    fname=f.replace("ict","nc")
    timeX=xr.DataArray(timeL)
    latX=xr.DataArray(latL)
    lonX=xr.DataArray(lonL)
    altX=xr.DataArray(altL)
    tempX=xr.DataArray(tempCL)
    d={"time":timeX,"lat":latX,"lon":lonX,"alt_feet":altX,"tempC":tempX}
    dSet=xr.Dataset(d)
    dSet.to_netcdf(fname)
    #stop
    #stop
filesN=glob.glob("PSDs/*sizedist*.nc")
from netCDF4 import Dataset
import numpy as np
for f in filesN:
    fh=Dataset(f)
    timeN=fh['time'][:]
    date=f.split("_")[2]
    tammsF=glob.glob("TAMMS/*"+date+"*nc")
    fhT=Dataset(tammsF[0])
    timeT=fhT["time"][:]
    tempCT=fhT["tempC"][:]
    tempL=[]
    for t in timeN:
        a=np.nonzero(timeT==t)
        if len(a[0])==1:
            tempT=tempCT[a[0][0]]
            tempL.append(tempT)
        else:
            tempL.append(-999)
    tempL=np.array(tempL)
    fnameT="PSDs/temp_%s.nc"%date
    tempX=xr.DataArray(tempL)
    dSet=xr.Dataset({"tempC":tempX,"time":xr.DataArray(timeN)})
    dSet.to_netcdf(fnameT)
    #stop
