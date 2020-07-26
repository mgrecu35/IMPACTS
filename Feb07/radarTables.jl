using PyCall

pickle=pyimport("pickle")

push!(PyVector(pyimport("sys")["path"]), "./")

include("scatTables.jl")
using .scatTables

temp,mass,fraction,bscat,Deq,ext,scat,g,vfall,
temp_r,mass_r,bscat_r,Deq_r,ext_r,scat_r,g_r,vfall_r,
tempKu,massKu,fractionKu,bscatKu,DeqKu,extKu,scatKu,gKu,vfallKu,
tempKu_r,massKu_r,bscatKu_r,DeqKu_r,extKu_r,scatKu_r,gKu_r,vfallKu_r,
tempW,massW,fractionW,bscatW,DeqW,extW,scatW,gW,vfallW,tempW_r,massW_r,bscatW_r,DeqW_r,
extW_r,scatW_r,gW_r,vfallW_r=scatTables.init()

xr=pyimport("xarray")
function saveTables()
    nwTX=xr["DataArray"](nwT)
    pwcTX=xr["DataArray"](pwcT)
    attKaTX=xr["DataArray"](attKaT)
    attKuTX=xr["DataArray"](attKuT)
    zKaTX=xr["DataArray"](zKaT)
    zKuTX=xr["DataArray"](zKuT)
    zKaTsX=xr["DataArray"](zKaTs)
    zKuTsX=xr["DataArray"](zKuTs)
    dmTsX=xr["DataArray"](dmTs)
    dmTX=xr["DataArray"](dmT)
    pwcTsX=xr["DataArray"](pwcTs)
    attKaTsX=xr["DataArray"](attKaTs)
    attKuTsX=xr["DataArray"](attKuTs)
    attWTsX=xr["DataArray"](attWTs)
    zWTsX=xr["DataArray"](zWTs)
    nwTsX=xr["DataArray"](nwTs)
    rateX=xr["DataArray"](rateT)
    ratesX=xr["DataArray"](rateTs)
    d=Dict("nwT"=>nwTX,"pwcT"=>pwcTX,"attKaT"=>attKaTX,
           "attKuT"=>attKuTX, "dmT"=>dmTX, "zKuT"=>zKuTX, "zKaT"=>zKaTX,
           "nwTs"=>nwTsX,"pwcTs"=>pwcTsX,"attKaTs"=>attKaTsX,
           "attKuTs"=>attKuTsX, "dmTs"=>dmTsX, "zKuTs"=>zKuTsX, "zKaTs"=>zKaTsX,
           "zWTs"=>zWTsX,"attWTs"=>attWTsX, "rate"=>rateX,
           "rateS"=>ratesX)
    tables=xr["Dataset"](d)
    tables["to_netcdf"]("tables_nsv_rho02_dmTs11.nc")
end

iread=0
if iread==1
    ns=9
    scatTables.getDmNwSF(tempKu,massKu,fractionKu,bscatKu,DeqKu,extKu,scatKu,gKu,vfallKu,
    tempKu_r,massKu_r,bscatKu_r,DeqKu_r,extKu_r,scatKu_r,gKu_r,vfallKu_r,
    temp,mass,fraction,bscat,Deq,ext,scat,g,vfall,
    temp_r,mass_r,bscat_r,Deq_r,ext_r,scat_r,g_r,vfall_r,ns,
    tempW,massW,fractionW,bscatW,DeqW,extW,scatW,gW,vfallW,
    tempW_r,massW_r,bscatW_r,DeqW_r,extW_r,scatW_r,gW_r,vfallW_r)
    scatTables.getDmNwR(tempKu,massKu,fractionKu,bscatKu,DeqKu,extKu,scatKu,gKu,vfallKu,
    tempKu_r,massKu_r,bscatKu_r,DeqKu_r,extKu_r,scatKu_r,gKu_r,vfallKu_r,
    temp,mass,fraction,bscat,Deq,ext,scat,g,vfall,
    temp_r,mass_r,bscat_r,Deq_r,ext_r,scat_r,g_r,vfall_r,
    tempW,massW,fractionW,bscatW,DeqW,extW,scatW,gW,vfallW,
    tempW_r,massW_r,bscatW_r,DeqW_r,extW_r,scatW_r,gW_r,vfallW_r)
    saveTables()
    exit(1)
end  
n4=pyimport("netCDF4")
fh=n4["Dataset"]("tables_nsv_rho02_dmTs11.nc")
zKuT=fh["variables"]["get"]("zKuT")["_get"]([0],[240],[1])
zKaT=fh["variables"]["get"]("zKaT")["_get"]([0],[240],[1])
attKuT=fh["variables"]["get"]("attKuT")["_get"]([0],[240],[1])
attKaT=fh["variables"]["get"]("attKaT")["_get"]([0],[240],[1])
dmT=fh["variables"]["get"]("dmT")["_get"]([0],[240],[1])
nwT=fh["variables"]["get"]("nwT")["_get"]([0],[240],[1])
pwcT=fh["variables"]["get"]("pwcT")["_get"]([0],[240],[1])
zKuTs=fh["variables"]["get"]("zKuTs")["_get"]([0],[240],[1])
zKaTs=fh["variables"]["get"]("zKaTs")["_get"]([0],[240],[1])
attKuTs=fh["variables"]["get"]("attKuTs")["_get"]([0],[240],[1])
attKaTs=fh["variables"]["get"]("attKaTs")["_get"]([0],[240],[1])
dmTs=fh["variables"]["get"]("dmTs")["_get"]([0],[240],[1])
nwTs=fh["variables"]["get"]("nwTs")["_get"]([0],[240],[1])
pwcTs=fh["variables"]["get"]("pwcTs")["_get"]([0],[240],[1])
attWTs=fh["variables"]["get"]("attWTs")["_get"]([0],[240],[1])
zWTs=fh["variables"]["get"]("zWTs")["_get"]([0],[240],[1])
rateTs=fh["variables"]["get"]("rateS")["_get"]([0],[240],[1])
rateT=fh["variables"]["get"]("rate")["_get"]([0],[240],[1])
