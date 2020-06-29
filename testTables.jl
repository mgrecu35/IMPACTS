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
tempKu_r,massKu_r,bscatKu_r,DeqKu_r,extKu_r,scatKu_r,gKu_r,vfallKu_r=scatTables.init()
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
     nwTsX=xr["DataArray"](nwTs)
     d=Dict("nwT"=>nwTX,"pwcT"=>pwcTX,"attKaT"=>attKaTX, "attKuT"=>attKuTX, "dmT"=>dmTX, "zKuT"=>zKuTX, "zKaT"=>zKaTX,
     "nwTs"=>nwTsX,"pwcTs"=>pwcTsX,"attKaTs"=>attKaTsX, "attKuTs"=>attKuTsX, "dmTs"=>dmTsX, "zKuTs"=>zKuTsX, "zKaTs"=>zKaTsX)
     tables=xr["Dataset"](d)
     tables["to_netcdf"]("tables.nc")
 end

iread=1
if iread==0
    scatTables.getDmNwSF(tempKu,massKu,fractionKu,bscatKu,DeqKu,extKu,scatKu,gKu,vfallKu,
    tempKu_r,massKu_r,bscatKu_r,DeqKu_r,extKu_r,scatKu_r,gKu_r,vfallKu_r,
    temp,mass,fraction,bscat,Deq,ext,scat,g,vfall,
    temp_r,mass_r,bscat_r,Deq_r,ext_r,scat_r,g_r,vfall_r)

    scatTables.getDmNwR(tempKu,massKu,fractionKu,bscatKu,DeqKu,extKu,scatKu,gKu,vfallKu,
    tempKu_r,massKu_r,bscatKu_r,DeqKu_r,extKu_r,scatKu_r,gKu_r,vfallKu_r,
    temp,mass,fraction,bscat,Deq,ext,scat,g,vfall,
    temp_r,mass_r,bscat_r,Deq_r,ext_r,scat_r,g_r,vfall_r)
    saveTables()
    exit(1)
else
    fh=netCDF["Dataset"]("tables.nc")
    print(fh)
    zKuT=fh["variables"]["get"]("zKuT")["_get"]([0],[204],[1])
    zKaT=fh["variables"]["get"]("zKaT")["_get"]([0],[204],[1])
    attKuT=fh["variables"]["get"]("attKuT")["_get"]([0],[204],[1])
    attKaT=fh["variables"]["get"]("attKaT")["_get"]([0],[204],[1])
    dmT=fh["variables"]["get"]("dmT")["_get"]([0],[204],[1])
    nwT=fh["variables"]["get"]("nwT")["_get"]([0],[204],[1])
    pwcT=fh["variables"]["get"]("pwcT")["_get"]([0],[204],[1])
    zKuTs=fh["variables"]["get"]("zKuTs")["_get"]([0],[204],[1])
    zKaTs=fh["variables"]["get"]("zKaTs")["_get"]([0],[204],[1])
    attKuTs=fh["variables"]["get"]("attKuTs")["_get"]([0],[204],[1])
    attKaTs=fh["variables"]["get"]("attKaTs")["_get"]([0],[204],[1])
    dmTs=fh["variables"]["get"]("dmTs")["_get"]([0],[204],[1])
    nwTs=fh["variables"]["get"]("nwTs")["_get"]([0],[204],[1])
    pwcTs=fh["variables"]["get"]("pwcTs")["_get"]([0],[204],[1])
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
    zeta=0.0
    q=0.2*log(10)
    beta=0.72
    for k=1:ibbL[2]-20
        if zKu[k]==zKu[k] && zKu[k]>0
            attKu,attKa,zKa,pwcS=fromZKuTs(zKuC[k],dn[k])
            piaKu=piaKu+2*attKu*dr/1e3
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
    eps=1
    println("$(q*beta*zeta1d[nz])")

    for k=1:nz
        if zKu[k]==zKu[k] && zKu[k]>0
            zKuC[k]=zKu[k]-10/beta*log10(1-eps*q*beta*zeta1d[k])
        end
    end
    piaHB=-10/beta*log10(1-eps*q*beta*zeta1d[nz])
    return zKuC,piaHB,piaKu
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
include("enKF.jl")
np=pyimport("numpy")
ma=pyimport("numpy.ma")
function retrieve(impactsData)
    nEns=50
    pwcRet=zeros(nEns,300,399)
    zKaObs=zeros(300,399)
    zKaEnKF2d=zeros(300,399)
    pwcKF=zeros(300,399)
    pwcKuKa=zeros(300,399)
    dr=26.25
    snowL=zeros(300)
    rainL=zeros(300)
    zKaL=zeros(300)
    zKaL_obs=zeros(300)
    zKaLsfc=zeros(300)
    zKaLsfc_obs=zeros(300)
    piaKaL=zeros(300)
    dn=zeros(399)
    zKaSim2d=zeros(nEns,300,399).-99.9
    for i=1:1:300
        global zKu,zKa,ibbL
        zKu=impactsData[1][i,1:1:end]
        zKa=impactsData[2][i,1:1:end]
        ibbL=impactsData[3][i,:]
        ibbL2=Int(trunc(ibbL[2]/2))
        zObs=zKu[ibbL2:end]
        for iens=1:nEns
            dn=np["random"]["randn"](399)
            for ism=1:15
                dn=smooth(dn)
            end
            zKuC,piaHB,piaKu=attCorrect(zKu,dr,ibbL,dn)
            println("$(piaHB) $(piaKu)")
            piaKa=0
            attKa1=0.0
            attKa2=0
            pwcRet1d,piaKa,rainSfc,snowBB,zKa1,zKa1sfc,zKa1_obs,zKa1sfc_obs,zKaSim=
            prof1d(zKuC,zKa,dn,dr)
            indices=findall(zKaSim.>-7)
            indices2=findall(zKa[indices].==zKa[indices])
            nvars=length(indices2)
            if iens==1
                global xEns=zeros(nEns,nvars)
                global dnEns=zeros(nEns,nvars)
                global yEns=zeros(nEns,nvars)
                global yObs=zKa[indices][indices2]
                #println(yObs)
                #exit(1)
            end
            xEns[iens,:]=pwcRet1d[indices][indices2]
            dnEns[iens,:]=dn[indices][indices2]
            yEns[iens,:]=zKaSim[indices][indices2]

            if iens==nEns
                nm=nvars-10
                xm=enKF(xEns,yEns,yObs,nm)
                dnEnKF=enKF(dnEns,yEns,yObs,nm)

                #println(xm)
                ik=1
                global dnKF=zeros(399)
                for k in indices[indices2]
                    pwcKF[i,k]=xm[ik]
                    dnKF[k]=dnEnKF[ik]
                    ik=ik+1
                end
                global zKaSimEnKF
                pwcRet1d,piaKa,rainSfc,snowBB,zKa1,zKa1sfc,zKa1_obs,zKa1sfc_obs,zKaSimEnKF=
                prof1d(zKuC,zKa,dnKF,dr)
                zKaEnKF2d[i,:]=zKaSimEnKF
                zKaObs[i,:]=copy(zKa)
                #println(pwcKF[i,:])
                #exit(1)
            end
            piaKaL[i]=piaKa
            snowL[i]=snowBB
            rainL[i]=rainSfc
            zKaL_obs[i]=zKa1_obs
            pwcRet[iens,i,:].=pwcRet1d
            zKaSim2d[iens,i,:]=zKaSim
        end
    end
    return pwcRet,snowL,rainL,zKaL_obs,zKaL, zKaLsfc, zKaLsfc_obs,piaKaL,zKaSim2d, pwcKF, zKaObs, zKaEnKF2d
end

pwcRet,snowL,rainL,zKaL_obs,zKaL, zKaLsfc, zKaLsfc_obs,piaKaL,zKaSim,pwcKF,zKaObs,zKaEnKF=retrieve(impactsData)

matplot=pyimport("matplotlib.colors")
pwcRetT=copy(transpose(mean(pwcRet,dims=1)[1,:,:]))
pwcRetT=copy(transpose(pwcKF))
pwcRetTm=ma["masked_equal"](pwcRetT,0)
plt=pyimport("matplotlib.pyplot")
h=impactsData[end-1].-impactsData[end]
h.=h/1000.0
plt["figure"](figsize=(8,6))
norm=matplot["LogNorm"]()
plt["pcolormesh"](1:300,h,pwcRetTm[:,1:1:end],cmap="jet",vmin=0.1,vmax=2,norm=norm)
c=plt["colorbar"]()
c["ax"]["set_title"]("g/m^3")
plt["xlabel"]("Profile #")
plt["ylabel"]("Height(km)")
#plt["savefig"]("impacts0201.png")

#plt["show"]()
