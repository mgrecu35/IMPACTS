import glob
import xarray as xr
files=glob.glob('/home/grecu/Citation/Cit_NCAR_Olympex/201*.comb.spectrum.1Hz')
files2=glob.glob('/home/grecu/Citation/Cit_NCAR_Olympex/olym*nav*.txt')
files=sorted(files)[:-1]
files2=sorted(files2)
npsd=0
allTemps=[]
import numpy as np
for f,f2 in zip(files,files2):
    lines=open(f).readlines()

    lsp=lines[26].split()
    dmid=[float(v1) for v1 in lsp]
    lsp=lines[10].split()
    dInt=[float(v1) for v1 in lsp]
    
    NcL=[]
    timeL=[]
    iwcL=[]
    tempL=[]
    lines2=open(f2).readlines()
    timeL2=[]
    tempL2=[]
    for l in lines2[70:]:
        lsp=l.split()
        tempL2.append(float(lsp[1]))
        timeL2.append(float(lsp[0]))
    timeL2=np.array(timeL2)
    
    for l in lines[43:]:
        lsp=l.split()
        v=[float(v1) for v1 in lsp[3:]]
        iwc=float(lsp[2])
        if iwc>0:
            ind=np.argmin(abs(float(lsp[0])-timeL2))
            #stop
            if tempL2[ind]<-5 and abs(float(lsp[0])-timeL2[ind])<1:
                timeL.append(float(lsp[0]))
                iwcL.append(iwc)
                NcL.append(v)
                tempL.append(tempL2[ind])
            
    timeL=np.array(timeL)
    allTemps.extend(tempL)
    fnameout=f.split('/')[-1]+'.nc'
    NcLx=xr.DataArray(NcL)
    timeLx=xr.DataArray(timeL)
    iwcLx=xr.DataArray(iwcL )
    tempLx=xr.DataArray(tempL)
    d=xr.Dataset({"time":timeLx,"iwc":iwcLx,"NcL":NcLx,"tempC":tempLx})
    d.to_netcdf(fnameout)
    npsd+=len(timeL)
    print(len(timeL))
    #stop
    #stop
