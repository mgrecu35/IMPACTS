include("scatTables.jl")
using .scatTables
using PyCall
pickle=pyimport("pickle")
fh=pybuiltin("open")("IMPACTS0201.pklz","rb")
impactsData=pickle["load"](fh)

temp,mass,fraction,bscat,Deq,ext,scat,g,vfall,
temp_r,mass_r,bscat_r,Deq_r,ext_r,scat_r,g_r,vfall_r,
tempKu,massKu,fractionKu,bscatKu,DeqKu,extKu,scatKu,gKu,vfallKu,
tempKu_r,massKu_r,bscatKu_r,DeqKu_r,extKu_r,scatKu_r,gKu_r,vfallKu_r=scatTables.init()
for v in [tempKu_r,massKu_r,bscatKu_r,DeqKu_r,extKu_r,scatKu_r,gKu_r,vfallKu_r]
    println(typeof(v))
end

scatTables.getDmNwSF(tempKu,massKu,fractionKu,bscatKu,DeqKu,extKu,scatKu,gKu,vfallKu,
tempKu_r,massKu_r,bscatKu_r,DeqKu_r,extKu_r,scatKu_r,gKu_r,vfallKu_r,
temp,mass,fraction,bscat,Deq,ext,scat,g,vfall,
temp_r,mass_r,bscat_r,Deq_r,ext_r,scat_r,g_r,vfall_r)

scatTables.getDmNwR(tempKu,massKu,fractionKu,bscatKu,DeqKu,extKu,scatKu,gKu,vfallKu,
tempKu_r,massKu_r,bscatKu_r,DeqKu_r,extKu_r,scatKu_r,gKu_r,vfallKu_r,
temp,mass,fraction,bscat,Deq,ext,scat,g,vfall,
temp_r,mass_r,bscat_r,Deq_r,ext_r,scat_r,g_r,vfall_r)

include("psdInt_2.jl")

function attCorrect(zKu,dr,ibbL)
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
            n1,n2=bisection(zKuTs,zKuC[k])
            dZ=zKuTs[n2]-zKuTs[n1]+1e-5
            f=(zKu[k]-zKuTs[n1])/dZ
            attKu=(1-f)*attKuTs[n1]+f*attKuTs[n2]
            piaKu=piaKu+2*attKu*dr/1e3
            zeta=zeta+attKu/10.0^(0.1*zKuC[k]*beta)*10.0^(0.1*zKu[k]*beta)*dr/1e3
            zeta1d[k]=zeta
        end
        if k==ibbL[2]-20
            attKu1=attKu
        end
    end
    k=ibbL[2]+20
    n1,n2=bisection(zKuT,zKuC[k])
    dZ=zKuTs[n2]-zKuTs[n1]+1e-5
    f=(zKu[k]-zKuTs[n1])/dZ
    attKu2=(1-f)*attKuT[n1]+f*attKuT[n2]
    for k=ibbL[2]-19:ibbL[2]+19
        zeta=zeta+0.5*(attKu1+attKu2)/10.0^(0.1*zKuC[k]*beta)*10.0^(0.1*zKu[k]*beta)*dr/1e3
        zeta1d[k]=zeta
    end
    for k=ibbL[2]+20:399
        if zKu[k]==zKu[k] && zKu[k]>0
            n1,n2=bisection(zKuT,zKuC[k])
            dZ=zKuTs[n2]-zKuTs[n1]+1e-5
            f=(zKu[k]-zKuTs[n1])/dZ
            attKu=(1-f)*attKuT[n1]+f*attKuT[n2]
            piaKu=piaKu+2*attKu*dr/1e3
            zeta=zeta+attKu/10.0^(0.1*zKuC[k]*beta)*10.0^(0.1*zKu[k]*beta)*dr/1e3
            zeta1d[k]=zeta
        end
        if k==ibbL[2]+20
            piaKu=piaKu+2*0.5*(attKu1+attKu)*38*dr/1e3
        end
    end
    eps=1
    for k=1:nz
        if zKu[k]==zKu[k] && zKu[k]>0
            zKuC[k]=zKu[k]-10/beta*log10(1-eps*q*beta*zeta1d[k])
        end
    end
    piaHB=-10/beta*log10(1-eps*q*beta*zeta1d[nz])
    return zKuC,piaHB,piaKu
end

function retrieve(impactsData)
    pwcRet=zeros(300,399)
    dr=26.25
    snowL=zeros(300)
    rainL=zeros(300)
    zKaL=zeros(300)
    zKaL_obs=zeros(300)
    zKaLsfc=zeros(300)
    zKaLsfc_obs=zeros(300)
    piaKaL=zeros(300)
    for i=1:1:300
        global zKu,zKa,ibbL
        zKu=impactsData[1][i,1:1:end]
        zKa=impactsData[2][i,1:1:end]
        ibbL=impactsData[3][i,:]
        ibbL2=Int(trunc(ibbL[2]/2))
        zObs=zKu[ibbL2:end]
        zKuC,piaHB,piaKu=attCorrect(zKu,dr,ibbL)
        println("$(piaHB) $(piaKu)")
        piaKa=0
        attKa1=0.0
        attKa2=0
        for k=1:ibbL[2]-20
            if zKu[k]==zKu[k] && zKu[k]>0
                n1,n2=bisection(zKuTs,zKuC[k])
                #println("$(zKu[k]) $(n1) $(n2)")
                dZ=zKuTs[n2]-zKuTs[n1]+1e-5
                f=(zKu[k]-zKuTs[n1])/dZ
                pwcRet[i,k]=(1-f)*pwcTs[n1]+f*pwcTs[n2]
                attKa=(1-f)*attKaTs[n1]+f*attKaTs[n2]
                piaKa=piaKa+2*attKa*dr/1e3
            else
                pwcRet[i,k]=NaN
            end
            if k==ibbL[2]-20
                snowL[i]=pwcRet[i,k]
                zKa1=(1-f)*zKaTs[n1]+f*zKaTs[n2]
                zKaL[i]=zKa1-piaKa
                zKaL_obs[i]=zKa[k]
                attKa1=attKa
            end
        end
        for k=ibbL[2]-19:ibbL[2]+19
            pwcRet[i,k]=NaN
        end

        for k=ibbL[2]+20:399
            if zKu[k]==zKu[k] && zKu[k]>0
                n1,n2=bisection(zKuT,zKuC[k])
            #println("$(zKu[k]) $(n1) $(n2)")
                dZ=zKuTs[n2]-zKuTs[n1]+1e-5
                f=(zKu[k]-zKuTs[n1])/dZ
                pwcRet[i,k]=(1-f)*pwcT[n1]+f*pwcT[n2]
                attKa=(1-f)*attKaT[n1]+f*attKaT[n2]
                piaKa=piaKa+2*attKa*dr/1e3
            else
                pwcRet[i,k]=NaN
            end
            if k==ibbL[2]+20
                piaKa=piaKa+2*0.5*(attKa1+attKa)*40*dr/1e3
                rainL[i]=pwcRet[i,k]
            end
            if k==399
                zKa1=(1-f)*zKaT[n1]+f*zKaT[n2]
                zKaLsfc[i]=zKa1-piaKa
                zKaLsfc_obs[i]=zKa[k]
            end
        end
        piaKaL[i]=piaKa
    end
    return pwcRet,snowL,rainL,zKaL_obs,zKaL, zKaLsfc, zKaLsfc_obs,piaKaL
end
pwcRet,snowL,rainL,zKaL_obs,zKaL, zKaLsfc, zKaLsfc_obs,piaKaL=retrieve(impactsData)
np=pyimport("numpy")
ma=pyimport("numpy.ma")
pwcRetT=copy(transpose(pwcRet))
pwcRetTm=ma["masked_equal"](pwcRetT,0)
plt=pyimport("matplotlib.pyplot")
h=impactsData[end-1].-impactsData[end]
plt["pcolormesh"](1:300,h,pwcRetTm[:,1:1:end],cmap="jet",vmin=0.001)
plt["colorbar"]()
#plt["show"]()
