
#pickle=pyimport("pickle")
#fh=pybuiltin("open")("IMPACTS0201.pklz","rb")
#impactsData=pickle["load"](fh)
#pyimport("sys")




module scatTables
using PyCall
include("psdInt.jl")
include("psdInt_2.jl")

pushfirst!(PyVector(pyimport("sys")["path"]), ".")
export zKuT,dmT,nwT,pwcT,attKaT,attKuT,zKaT
export zKaTs,zKuTs,dmTs,pwcTs,attKuTs,attKaTs,nwTs
nTs=204
nT=204
nwT=zeros(nT)
pwcT=zeros(nT)
attKaT=zeros(nT)
attKuT=zeros(nT)
zKaT=zeros(nT)
zKuT=zeros(nT)
zKaTs=zeros(nTs)
zKuTs=zeros(nTs)
dmTs=zeros(nTs)
pwcTs=zeros(nTs)
attKaTs=zeros(nTs)
attKuTs=zeros(nTs)
nwTs=zeros(nTs)
for i=1:nT
    zKuT[i]=0+(i-1)*0.25
end
for i=1:nTs
    zKuTs[i]=-7+(i-1)*0.25
end
zCoeffs=[ 0.00159248, -0.01061376,  0.33276274]
dmT=zCoeffs[1].*zKuT.*zKuT+zCoeffs[2].*zKuT.+zCoeffs[3]
dmTs.=dmT
for i=1:16
    dmT[i]=dmT[17]-(17-i)*0.0025
    dmTs[i]=dmTs[17]-(17-i)*0.0025
end
ds=maximum(dmT)-minimum(dmT)
for i=1:nT
    dmT[i]=dmT[1]+0.75*(dmT[i]-dmT[1])
    dmTs[i]=dmTs[1]+0.75*(dmT[i]-dmT[1])
end
#
function init()
    #global readS
    readS=pyimport("readScattProf")
    fnameIce="ice-self-similar-aggregates_35-GHz_scat.nc"
    fnameIceKu="ice-self-similar-aggregates_13-GHz_scat.nc"
    fnameRain="liquid-water_35-GHz_scat.nc"
    fnameRainKu="liquid-water_13-GHz_scat.nc"

    temp,mass,fraction,bscat,Deq,
    ext,scat,g,vfall=readS["readScatProf"](fnameIce);
    temp_r,mass_r,bscat_r,Deq_r,
    ext_r,scat_r,g_r,vfall_r=readS["readScatProfR"](fnameRain);
    tempKu,massKu,fractionKu,bscatKu,DeqKu,
    extKu,scatKu,gKu,vfallKu=readS["readScatProf"](fnameIceKu);
    tempKu_r,massKu_r,bscatKu_r,DeqKu_r,
    extKu_r,scatKu_r,gKu_r,vfallKu_r=readS["readScatProfR"](fnameRainKu);
    return temp,mass,fraction,bscat,Deq,ext,scat,g,vfall,
    temp_r,mass_r,bscat_r,Deq_r,ext_r,scat_r,g_r,vfall_r,
    tempKu,massKu,fractionKu,bscatKu,DeqKu,extKu,scatKu,gKu,vfallKu,
    tempKu_r,massKu_r,bscatKu_r,DeqKu_r,extKu_r,scatKu_r,gKu_r,vfallKu_r
end
export getDmNwR, getDmNwS
function getDmNwR(tempKu::Array{Float32,1},massKu::Array{Float32,2},fractionKu::Array{Float32,1},
    bscatKu::Array{Float32,3},DeqKu::Array{Float32,2},extKu::Array{Float32,3},scatKu::Array{Float32,3},
    gKu::Array{Float32,3},vfallKu::Array{Float32,2},
    tempKu_r::Array{Float32,1},massKu_r::Array{Float32,1},bscatKu_r::Array{Float32,2},
    DeqKu_r::Array{Float32,1},extKu_r::Array{Float32,2},scatKu_r::Array{Float32,2},
    gKu_r::Array{Float32,2},vfallKu_r::Array{Float32,1},
    temp::Array{Float32,1},mass::Array{Float32,2},fraction::Array{Float32,1},
    bscat::Array{Float32,3},Deq::Array{Float32,2},ext::Array{Float32,3},scat::Array{Float32,3},
    g::Array{Float32,3},vfall::Array{Float32,2},
    temp_r::Array{Float32,1},mass_r::Array{Float32,1},bscat_r::Array{Float32,2},
    Deq_r::Array{Float32,1},ext_r::Array{Float32,2},scat_r::Array{Float32,2},
    g_r::Array{Float32,2},vfall_r::Array{Float32,1})
    nx=size(zKuT)[1]
    println(nx)
    mu=2.0
    dn=0.1
    freqKu=13.8
    wlKu=300/freqKu
    freqKa=35.5
    wlKa=300/freqKa
    zCoeffs=[ 0.00159248, -0.01061376,  0.33276274]

    for i=1:nx
        zObs1=zKuT[i]
        dnRet=dn
        for it=1:4
            dn=dnRet
            global dmc
            retZ=get_fZr(zObs1,bscatKu_r,scatKu_r,extKu_r,gKu_r,DeqKu_r,vfallKu_r,wlKu,dn,mu)
            rwc, Z, att,scatInt,gInt, vdop, dm=retZ
            dmc=0.1*dmT[i]
            dn1=dn*10
            retZ1=get_fZr(zObs1,bscatKu_r,scatKu_r,extKu_r,gKu_r,DeqKu_r,vfallKu_r,wlKu,dn1,mu)
            rwc1, Z1, att1,scatInt1,gInt1, vdop1, dm1=retZ1
            gradDm=(dm1-dm)/(log(dn1/dn))
            dnRet=dn*exp(1.0*(dmc-dm)*gradDm/(gradDm*gradDm+0.000001))
        end
        #dnRet=0.1
        retZ=get_fZr(zObs1,bscatKu_r,scatKu_r,extKu_r,gKu_r,DeqKu_r,vfallKu_r,wlKu,dnRet,mu)
        rwc, Z, att,scatInt,gInt, vdop, dm=retZ
        pwcT[i]=rwc
        attKuT[i]=att*4.343
        nwT[i]=log10(dnRet)
        zKa, attKa,scatIntKa,gIntKa, vdopKa,dmKa=get_Zr(rwc,bscat_r,scat_r,ext_r,g_r,Deq_r,vfall_r,wlKa,dnRet,mu)
        #println("$Z $(zKu[i]) $(dm) $(dmc)")
        zKaT[i]=zKa
        attKaT[i]=attKa*4.343
        #get_fZr(zKu[i],bscatKu_r,scatKu_r,extKu_r,gKu_r,DeqKu_r,vfallKu_r,wlKu,dn,mu)
    end
end
function getDmNwS(tempKu::Array{Float32,1},massKu::Array{Float32,2},fractionKu::Array{Float32,1},
    bscatKu::Array{Float32,3},DeqKu::Array{Float32,2},extKu::Array{Float32,3},scatKu::Array{Float32,3},
    gKu::Array{Float32,3},vfallKu::Array{Float32,2},
    tempKu_r::Array{Float32,1},massKu_r::Array{Float32,1},bscatKu_r::Array{Float32,2},
    DeqKu_r::Array{Float32,1},extKu_r::Array{Float32,2},scatKu_r::Array{Float32,2},
    gKu_r::Array{Float32,2},vfallKu_r::Array{Float32,1},
    temp::Array{Float32,1},mass::Array{Float32,2},fraction::Array{Float32,1},
    bscat::Array{Float32,3},Deq::Array{Float32,2},ext::Array{Float32,3},scat::Array{Float32,3},
    g::Array{Float32,3},vfall::Array{Float32,2},
    temp_r::Array{Float32,1},mass_r::Array{Float32,1},bscat_r::Array{Float32,2},
    Deq_r::Array{Float32,1},ext_r::Array{Float32,2},scat_r::Array{Float32,2},
    g_r::Array{Float32,2},vfall_r::Array{Float32,1})
    nx=size(zKaTs)[1]
    println(nx)
    mu=-1.0
    dn=0.1
    freqKu=13.8
    wlKu=300/freqKu
    freqKa=35.5
    wlKa=300/freqKa
    ns=11
    for i=1:nx
        zObs1=zKaTs[i]
        dnRet=dn
        rwc=0.1
        for it=1:4
            dn=dnRet
            global dmc
            retZ=get_fZs(zObs1,ns,bscat,scat,ext,g,Deq,vfall,wlKa,dn,mu)
            rwc, Z, att,scatInt,gInt, vdop, dm=retZ
            dmc=0.1*dmTs[i]
            dn1=dn*10
            retZ1=get_fZs(zObs1,ns,bscat,scat,ext,g,Deq,vfall,wlKa,dn1,mu)
            rwc1, Z1, att1,scatInt1,gInt1, vdop1, dm1=retZ1
            gradDm=(dm1-dm)/(log(dn1/dn))
            dnRet=dn*exp(1.0*(dmc-dm)*gradDm/(gradDm*gradDm+0.000001))
        end
        retZ=get_fZs(zObs1,ns,bscat,scat,ext,g,Deq,vfall,wlKa,dnRet,mu)
        rwc, Z, att,scatInt,gInt, vdop, dm=retZ
        pwcTs[i]=rwc
        attKaTs[i]=att*4.343
        nwTs[i]=log10(dnRet)
        zKaTs[i]=Z
        zKu, attKu,scatIntKu,gIntKu, vdopKu,dmKu=
        get_Zs(rwc,ns,bscatKu,scatKu,extKu,gKu,DeqKu,vfallKu,wlKu,dnRet,mu)
        #println("$Z $(zKu[i]) $(dm) $(dmc)")
        zKuTs[i]=zKu
        attKuTs[i]=attKu*4.343
    end
end
function getDmNwSF(tempKu::Array{Float32,1},massKu::Array{Float32,2},fractionKu::Array{Float32,1},
    bscatKu::Array{Float32,3},DeqKu::Array{Float32,2},extKu::Array{Float32,3},scatKu::Array{Float32,3},
    gKu::Array{Float32,3},vfallKu::Array{Float32,2},
    tempKu_r::Array{Float32,1},massKu_r::Array{Float32,1},bscatKu_r::Array{Float32,2},
    DeqKu_r::Array{Float32,1},extKu_r::Array{Float32,2},scatKu_r::Array{Float32,2},
    gKu_r::Array{Float32,2},vfallKu_r::Array{Float32,1},
    temp::Array{Float32,1},mass::Array{Float32,2},fraction::Array{Float32,1},
    bscat::Array{Float32,3},Deq::Array{Float32,2},ext::Array{Float32,3},scat::Array{Float32,3},
    g::Array{Float32,3},vfall::Array{Float32,2},
    temp_r::Array{Float32,1},mass_r::Array{Float32,1},bscat_r::Array{Float32,2},
    Deq_r::Array{Float32,1},ext_r::Array{Float32,2},scat_r::Array{Float32,2},
    g_r::Array{Float32,2},vfall_r::Array{Float32,1})
    nx=size(zKaTs)[1]
    println(nx)
    mu=2.0
    dn=0.1
    freqKu=13.8
    wlKu=300/freqKu
    freqKa=35.5
    wlKa=300/freqKa
    ns=11
    for i=1:nx
        zObs1=zKuTs[i]
        dnRet=dn
        rwc=0.1
        for it=1:4
            dn=dnRet
            global dmc
            retZ=get_fZs(zObs1,ns,bscatKu,scatKu,extKu,gKu,DeqKu,DeqKu_r,vfallKu,vfallKu_r,wlKu,dn,mu)
            rwc, Z, att,scatInt,gInt, vdop, dm=retZ
            dmc=0.1*dmTs[i]
            dn1=dn*10
            retZ1=get_fZs(zObs1,ns,bscatKu,scatKu,extKu,gKu,DeqKu,DeqKu_r,vfallKu,vfallKu_r,wlKu,dn1,mu)
            rwc1, Z1, att1,scatInt1,gInt1, vdop1, dm1=retZ1
            gradDm=(dm1-dm)/(log(dn1/dn))
            dnRet=dn*exp(1.0*(dmc-dm)*gradDm/(gradDm*gradDm+0.000001))
        end
        #dnRet=0.1
        retZ=get_fZs(zObs1,ns,bscatKu,scatKu,extKu,gKu,DeqKu,DeqKu_r,vfallKu,vfallKu_r,wlKu,dnRet,mu)
        rwc, Z, att,scatInt,gInt, vdop, dm=retZ
        pwcTs[i]=rwc
        attKuTs[i]=att*4.343
        nwTs[i]=log10(dnRet)
        zKuTs[i]=Z
        zKa, attKa,scatIntKa,gIntKa, vdopKa,dmKa=
        get_Zs(rwc,ns,bscat,scat,ext,g,Deq,Deq_r,vfall,vfall_r,wlKa,dnRet,mu)
        #println("$Z $(zKu[i]) $(dm) $(dmc)")
        zKaTs[i]=zKa
        attKaTs[i]=attKa*4.343
    end
end
end
