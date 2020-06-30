include("scatTables.jl")
using .scatTables
using PyCall
pickle=pyimport("pickle")
fh=pybuiltin("open")("IMPACTS0201.pklz","rb")
impactsData=pickle["load"](fh)
netCDF=pyimport("netCDF4")
temp,mass,fraction,bscat,Deq,ext,scat,g,vfall,
temp_r,mass_r,bscat_r,Deq_r,ext_r,scat_r,g_r,vfall_r,
tempKu,massKu,fractionKu,bscatKu,DeqKu,extKu,scatKu,gKu,vfallKu,
tempKu_r,massKu_r,bscatKu_r,DeqKu_r,extKu_r,scatKu_r,gKu_r,vfallKu_r,
tempW,massW,fractionW,bscatW,DeqW,extW,scatW,gW,vfallW,tempW_r,massW_r,bscatW_r,DeqW_r,
extW_r,scatW_r,gW_r,vfallW_r=scatTables.init()

for v in [tempKu_r,massKu_r,bscatKu_r,DeqKu_r,extKu_r,scatKu_r,gKu_r,vfallKu_r]
    println(typeof(v))
end
xr=pyimport("xarray")
function saveTables()
     #using .scatTables
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
     d=Dict("nwT"=>nwTX,"pwcT"=>pwcTX,"attKaT"=>attKaTX, "attKuT"=>attKuTX, "dmT"=>dmTX, "zKuT"=>zKuTX, "zKaT"=>zKaTX,
     "nwTs"=>nwTsX,"pwcTs"=>pwcTsX,"attKaTs"=>attKaTsX, "attKuTs"=>attKuTsX, "dmTs"=>dmTsX, "zKuTs"=>zKuTsX, "zKaTs"=>zKaTsX,
     "zWTs"=>zWTsX,"attWTs"=>attWTsX)
     tables=xr["Dataset"](d)
     tables["to_netcdf"]("tables_nsv.nc")
 end

iread=1
if iread==0
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
else
    fh=netCDF["Dataset"]("tables_nsv.nc")
    print(fh)#exit()
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
    println(size(zKuT))
end
include("psdInt_2.jl")
function fromZKuTs(zKuC,dn)
    if zKuC-10*dn<-7
        dn=(zKuC+7)/(10)
    end
    if zKuC-10*dn>43.75
        dn=(zKuC-43.75)/(10)
    end
    zKuC1=zKuC-10*dn
    n1,n2=bisection(zKuTs,zKuC1)
    dZ=zKuTs[n2]-zKuTs[n1]+1e-5
    f=(zKuC1-zKuTs[n1])/dZ
    attKu=10^dn*((1-f)*attKuTs[n1]+f*attKuTs[n2])
    attKa=10^dn*((1-f)*attKaTs[n1]+f*attKaTs[n2])
    zKa=(1-f)*zKaTs[n1]+f*zKaTs[n2]+10*dn
    pwcRet=10^dn*((1-f)*pwcTs[n1]+f*pwcTs[n2])
    return attKu,attKa,zKa,pwcRet
end
function fromZKuT(zKuC,dn)
    if zKuC-10*dn<0
        dn=(zKuC+0)/(10)
    end
    if zKuC-10*dn>50.75
        dn=(zKuC-50.75)/(10)
    end
    zKuC1=zKuC-10*dn
    n1,n2=bisection(zKuT,zKuC1)
    dZ=zKuTs[n2]-zKuTs[n1]+1e-5
    f=(zKuC1-zKuT[n1])/dZ
    attKu=10^dn*((1-f)*attKuT[n1]+f*attKuT[n2])
    attKa=10^dn*((1-f)*attKaT[n1]+f*attKaT[n2])
    zKa=(1-f)*zKaT[n1]+f*zKaT[n2]+10*dn
    pwcRet=10^dn*((1-f)*pwcT[n1]+f*pwcT[n2])
    return attKu,attKa,zKa, pwcRet
end
function attCorrect(zKu,dr,ibbL,dn)
    piaKu=0
    attKu1=0.0
    attKu2=0
    nz=size(zKu)[1]
    zKuC=zeros(nz).+zKu
    zeta1d=zeros(nz)

    q=0.2*log(10)
    beta=0.72
    eps=1

    for it=1:3
    zeta=0.0
    global kaCorr=zeros(nz)
    piaKa=0.0
    piaKu=0.0
    for k=1:ibbL[2]-20
        if zKu[k]==zKu[k] && zKu[k]>0
            attKu,attKa,zKa,pwcS=fromZKuTs(zKuC[k],dn[k])
            piaKu=piaKu+2*attKu*dr/1e3
            piaKa=piaKa+2*attKa*dr/1e3
            kaCorr[k]=kaCorr[k]+piaKa
            zeta=zeta+attKu/10.0^(0.1*zKuC[k]*beta)*10.0^(0.1*zKu[k]*beta)*dr/1e3
            zeta1d[k]=zeta
        end
        if k==ibbL[2]-20
            attKu1=attKu
        end
    end
    k=ibbL[2]+20
    attKu2,attKa2,zKa2=fromZKuT(zKuC[k],dn[k])
    for k=ibbL[2]-19:ibbL[2]+19
        zeta=zeta+0.5*(attKu1+attKu2)/10.0^(0.1*zKuC[k]*beta)*10.0^(0.1*zKu[k]*beta)*dr/1e3
        zeta1d[k]=zeta
    end
    for k=ibbL[2]+20:399
        if zKu[k]==zKu[k] && zKu[k]>0
            attKu,attKa,zKa,pwcR=fromZKuT(zKuC[k],dn[k])
            piaKu=piaKu+2*attKu*dr/1e3
            zeta=zeta+attKu/10.0^(0.1*zKuC[k]*beta)*10.0^(0.1*zKu[k]*beta)*dr/1e3
            zeta1d[k]=zeta
        end
        if k==ibbL[2]+20
            piaKu=piaKu+2*0.5*(attKu1+attKu)*38*dr/1e3
        end
    end
    piaMax=25
    if q*beta*zeta1d[nz]>0.99999
        eps=(1-10^(-0.1*piaMax*beta))/(q*beta*zeta1d[nz])
        dn=dn*eps^(1.0/(1-beta))
    end
    #println("$(q*beta*zeta1d[nz])")

    for k=1:nz
        if zKu[k]==zKu[k] && zKu[k]>0
            zKuC[k]=zKu[k]-10/beta*log10(1-eps*q*beta*zeta1d[k])
        end
    end
    global piaHB=-10/beta*log10(1-eps*q*beta*zeta1d[nz])
    end
    return zKuC,piaHB,piaKu,dn,kaCorr
end
function prof1d(zKuC,zKa,dn,dr)
    piaKa=0.0
    nz=size(zKuC)[1]
    pwcRet=zeros(nz)
    attKa1=0.0
    zKaSim=zeros(nz).-99.9
    for k=1:ibbL[2]-20
        if zKu[k]==zKu[k] && zKu[k]>0
            attKu,attKa,zKa1,pwcS=fromZKuTs(zKuC[k],dn[k])
            pwcRet[k]=pwcS
            piaKa=piaKa+2*attKa*dr/1e3
            zKaSim[k]=zKa1-piaKa-attKa*dr/1e3
        else
            pwcRet[k]=NaN
        end
        if k==ibbL[2]-20
            global snowL=pwcRet[k]
            global zKaL=zKa1-piaKa
            global zKaL_obs=zKa[k]
            attKa1=attKa
        end
    end
    for k=ibbL[2]-19:ibbL[2]+19
        pwcRet[k]=NaN
    end

    for k=ibbL[2]+20:399
        if zKu[k]==zKu[k] && zKu[k]>0
            attKu,attKa,zKa1,pwcR=fromZKuT(zKuC[k],dn[k])
            zKaSim[k]=zKa1-piaKa-attKa*dr/1e3
            piaKa=piaKa+2*attKa*dr/1e3
            pwcRet[k]=pwcR
        else
            pwcRet[k]=NaN
        end
        if k==ibbL[2]+20
            piaKa=piaKa+2*0.5*(attKa1+attKa)*38*dr/1e3
            zKaSim[k]=zKaSim[k]-2*0.5*(attKa1+attKa)*38*dr/1e3
            global rainL=pwcR
        end
        if k==399
            global zKaLsfc=zKa1-piaKa
            global zKaLsfc_obs=zKa[k]
        end
    end
    return pwcRet,piaKa,rainL,snowL,zKaL,zKaLsfc,zKaL_obs,zKaLsfc_obs, zKaSim
end
function retrieveN(zKu,dr,ibbL,n1,n2,n3,n4)
    spl1=Spline1D([itop,ibbL[2]-20,ibbL[2]+20,399],[n1,n2,n3,n4])
    dn=zeros(399)
    dn[itop:399]=spl1(itop:399)
    zKuC,piaHB,piaKu,dn=attCorrect(zKu,dr,ibbL,dn)
    pwcRet1d,piaKa,rainSfc,snowBB,zKa1,zKa1sfc,zKa1_obs,zKa1sfc_obs,zKaSim=
    prof1d(zKuC,zKa,dn,dr)
    return pwcRet1d,piaKa,rainSfc,snowBB,zKa1,zKa1sfc,zKa1_obs,zKa1sfc_obs,zKaSim
end
include("enKF.jl")
np=pyimport("numpy")
ma=pyimport("numpy.ma")
function retrieve(impactsData)
    dmRet=zeros(300,399)
    dmRetZ=zeros(300,399)
    pwcRet=zeros(300,399)
    zKaObs=zeros(300,399)
    zKuObs=zeros(300,399)
    zW=zeros(2,300,399)
    zW.=NaN
    dmRet.=NaN
    dmRetZ.=NaN
    pwcRet.=NaN
    d1=[]
    d2=[]
    snowL=zeros(300)
    piaWL=zeros(300)
    zw1=[]
    zw2=[]
    for i=1:1:300
        global zKu,zKa,ibbL
        zKu=impactsData[1][i,1:1:end]
        zKu.=zKu.+0.5
        zKa=impactsData[2][i,1:1:end]
        ibbL=impactsData[3][i,:]
        zWObs=impactsData[6][i,1:1:end]
        ibbL2=Int(trunc(ibbL[2]/2))
        zObs=zKu[ibbL2:end]
        dn2d=zeros(399)
        k=1
        for k=1:399-2
            if zKu[k]>0 && zKu[k+1]>0
                global itop=k
                break
            end
        end
        dZs=zKuTs-zKaTs
        dn=zeros(399)
        dr=26.25
        zKuC,piaHB,piaKu,dn,kaCorr=attCorrect(zKu,dr,ibbL,dn)
        zKa1=zKa+kaCorr
        zKa1.=zKa1.-1
        for k=itop:ibbL[2]-20
            dZ=zKu[k]-zKa1[k]
            if dZ>0 && dZ==dZ
                n1,n2=bisection(dZs,dZ)
                dmRet[i,k]=dmTs[n1]
                #dn1=0.5*(zKuC[k]-zKuTs[n1])/10.
                #pwcRet[i,k]=pwcTs[n1]*10^dn1
                n1,n2=bisection(zKuTs,zKuC[k])
                dmRetZ[i,k]=dmTs[n1]
                dmM=0.5*dmRet[i,k]+0.5*dmRetZ[i,k]
                n1,n2=bisection(dmTs,dmM)
                dn1=0.95*(zKuC[k]-zKuTs[n1])/10.
                dn[k]=dn1
                pwcRet[i,k]=pwcTs[n1]*10^dn1
                #push!(d1,dmM)
                #push!(d2,dmRetZ[i,k])
                if k==ibbL[2]-20
                    snowL[i]=pwcRet[i,k]
                end
            else
                if zKu[k]>0
                    n1,n2=bisection(zKuTs,zKuC[k])
                    pwcRet[i,k]=pwcTs[n1]
                end
            end
        end
        zKuC,piaHB,piaKu,dn,kaCorr=attCorrect(zKu,dr,ibbL,dn)
        zKa1=zKa+kaCorr
        zKa1.=zKa1.-1
        piaW=0.
        for k=itop:ibbL[2]-20
            dZ=zKu[k]-zKa1[k]
            if dZ>0 && dZ==dZ
                n1,n2=bisection(dZs,dZ)
                dmRet[i,k]=dmTs[n1]
                #dn1=0.5*(zKuC[k]-zKuTs[n1])/10.
                #pwcRet[i,k]=pwcTs[n1]*10^dn1
                n1,n2=bisection(zKuTs,zKuC[k])
                dmRetZ[i,k]=dmTs[n1]
                dmM=0.85*dmRet[i,k]+0.15*dmRetZ[i,k]
                n1,n2=bisection(dmTs,dmM)
                dn1=0.95*(zKuC[k]-zKuTs[n1])/10.
                dn[k]=dn1
                pwcRet[i,k]=pwcTs[n1]*10^dn1
                attW=attWTs[n1]*10^dn1
                piaW=piaW+attW*2*dr/1e3
                zW[1,i,k]=zWTs[n1]+10*dn1-piaW
                zW[2,i,k]=zWObs[k]
                push!(d1,dmM)
                push!(d2,dmRetZ[i,k])

                if k==ibbL[2]-20
                    snowL[i]=pwcRet[i,k]
                    if zWObs[k]==zWObs[k] && zWObs[k]>0
                        push!(zw1,zWObs[k])
                        push!(zw2,zW[1,i,k])
                    end
                end
            else
                if zKu[k]>0
                    n1,n2=bisection(zKuTs,zKuC[k])
                    pwcRet[i,k]=pwcTs[n1]
                end
            end
        end
        piaWL[i]=piaW
        zKuObs[i,:]=zKu
        zKaObs[i,:]=zKa
    end
    return  zKaObs, zKuObs, dmRet, dmRetZ,d1,d2,pwcRet, snowL,zW,piaWL,zw1,zw2
end

zKaObs,zKuObs,dmRet,dmRetZ,d1,d2,pwcRet,snowL,zW,piaWL,zw1,zw2=retrieve(impactsData)

matplot=pyimport("matplotlib.colors")

plt=pyimport("matplotlib.pyplot")
h=impactsData[end-3].-impactsData[end-2]
h.=h/1000.0
plt["figure"](figsize=(8,6))
norm=matplot["LogNorm"]()
#plt["pcolormesh"](1:300,h,pwcRetTm[:,1:1:end],cmap="jet",vmin=0.1,vmax=2,norm=norm)
zKuObsT=copy(transpose(zKuObs))
plt["pcolormesh"](1:300,h,zKuObsT,cmap="jet",vmin=0.1,vmax=40)
c=plt["colorbar"]()
c["ax"]["set_title"]("g/m^3")
plt["xlabel"]("Profile #")
plt["ylabel"]("Height(km)")
plt["savefig"]("impacts0201.png")

plt["figure"](figsize=(8,6))
dmTr=copy(transpose(dmRet))
plt["pcolormesh"](1:300,h,dmTr,cmap="jet",vmin=0.1,vmax=3)
c=plt["colorbar"]()
c["ax"]["set_title"]("mm")
plt["xlabel"]("Profile #")
plt["ylabel"]("Height(km)")

plt["figure"](figsize=(8,6))
dmTZr=copy(transpose((pwcRet)))
plt["pcolormesh"](1:300,h,dmTZr,cmap="jet",vmin=0.1,vmax=2,norm=norm)
c=plt["colorbar"]()
c["ax"]["set_title"]("mm")
plt["xlabel"]("Profile #")
plt["ylabel"]("Height(km)")

plt["figure"](figsize=(8,6))
zWt=copy(transpose((zW[1,:,:])))
plt["pcolormesh"](1:300,h,zWt,cmap="jet",vmin=-10.0,vmax=20)
c=plt["colorbar"]()
c["ax"]["set_title"]("dBZ")
plt["xlabel"]("Profile #")
plt["ylabel"]("Height(km)")

plt["figure"](figsize=(8,6))
zWt=copy(transpose((zW[2,:,:])))
plt["pcolormesh"](1:300,h,zWt,cmap="jet",vmin=-10.0,vmax=20)
c=plt["colorbar"]()
c["ax"]["set_title"]("dBZ")
plt["xlabel"]("Profile #")
plt["ylabel"]("Height(km)")

plt["figure"](figsize=(8,6))
dZ=copy(transpose(zKuObs-zKaObs))
a=findall((1.5.+zKuObs.-zKaObs).<0)
for c in a
    i1=c[1]
    j1=c[2]
    dZ[j1,i1]=NaN
end
plt["pcolormesh"](1:300,h,dZ,cmap="jet",vmin=0,vmax=10)
c=plt["colorbar"]()
