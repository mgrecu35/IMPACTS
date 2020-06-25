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
pwcRet=zeros(300,399)
include("psdInt_2.jl")
dr=26.25
snowL=zeros(300)
rainL=zeros(300)
for i=1:1:300
    global zKu,zKa,ibbL
    zKu=impactsData[1][i,1:1:end]
    zKa=impactsData[2][i,1:1:end]
    ibbL=impactsData[3][i,:]
    ibbL2=Int(trunc(ibbL[2]/2))
    zObs=zKu[ibbL2:end]
    for k=1:ibbL[2]-20
        if zKu[k]==zKu[k] && zKu[k]>0
            n1,n2=bisection(zKuTs,zKu[k])
            println("$(zKu[k]) $(n1) $(n2)")
            dZ=zKuTs[n2]-zKuTs[n1]+1e-5
            f=(zKu[k]-zKuTs[n1])/dZ
            pwcRet[i,k]=(1-f)*pwcTs[n1]+f*pwcTs[n2]
        else
            pwcRet[i,k]=NaN
        end

        if k==ibbL[2]-20
            snowL[i]=pwcRet[i,k]
        end
    end
    for k=ibbL[2]-19:ibbL[2]+19
        pwcRet[i,k]=NaN
    end
    for k=ibbL[2]+20:399
        if zKu[k]==zKu[k] && zKu[k]>0
            n1,n2=bisection(zKuT,zKu[k])
            #println("$(zKu[k]) $(n1) $(n2)")
            dZ=zKuTs[n2]-zKuTs[n1]+1e-5
            f=(zKu[k]-zKuTs[n1])/dZ
            pwcRet[i,k]=(1-f)*pwcT[n1]+f*pwcT[n2]
            if k==ibbL[2]-200
                snowL[i]=pwcRet[i,k]
            end
        else
            pwcRet[i,k]=NaN
        end
        if k==ibbL[2]+20
            rainL[i]=pwcRet[i,k]
        end
    end
end
np=pyimport("numpy")
ma=pyimport("numpy.ma")
pwcRetT=copy(transpose(pwcRet))
pwcRetTm=ma["masked_equal"](pwcRetT,0)
plt=pyimport("matplotlib.pyplot")
h=impactsData[end-1].-impactsData[end]
plt["pcolormesh"](1:300,h,pwcRetTm[:,1:1:end],cmap="jet",vmin=0.001)
plt["colorbar"]()
#plt["show"]()
