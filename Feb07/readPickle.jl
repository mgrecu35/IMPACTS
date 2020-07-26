using PyCall

pickle=pyimport("pickle")

f=pybuiltin("open")("dataVN_3D.pklz","rb")
dataOut=pickle["load"](f)

clutFL,reliabFlagNSL,reliabFlagMSL,pathAttenNSL,
pathAttenMSL,binSfcNSL,binSTopL,binZeroDegL,
zKuL,zKaL,cmbSfcRainL,gvSfcRainL,cmbPWC=dataOut;

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
    tables["to_netcdf"]("tables_nsv_rho06_dmTs11.nc")
end

iread=0
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
end  
n4=pyimport("netCDF4")
fh=n4["Dataset"]("tables_nsv_rho04_dmTs11.nc")
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
include("attCorrection.jl")

push!(PyVector(pyimport("sys")["path"]), "./")
readpick=pyimport("readPickleJ")
dataOut2=readpick["readpickle"]("dataVN_3D.pklz")
clutFL,reliabFlagNSL,reliabFlagMSL,pathAttenNSL,
pathAttenMSL,binSfcNSL,binSTopL,binZeroDegL,
zKuL,zKaL,cmbSfcRainL,gvSfcRainL,cmbPWCL,piaL=dataOut2


dr=0.125e3
p1=[]
p2=[]
p3=[]
p4=[]
p5=[]
w1=[]
w2=[]
ic=3
nt=size(zKuL[ic])[1]
using Statistics
pRateL=[]
for i=1:176
    push!(pRateL,[])
end
np=pyimport("numpy")
h=0.125*((176:-1:100).-100);
dmZ=np["loadtxt"]("dmZ.txt")

#dmS=Spline1D(dmZ[end:-1:1,2],dmZ[end:-1:1,1])


#np["loadtxt"]("dmZ.txt")
fractS=Spline1D([0.,1.25,3.85,6.0,6.35,10.],[1,1,1,0.1,0,0],k=1)
fractD=Spline1D([0.,1.25,3.85,6.35,10.],[1.8,1.8,1.7,1.6,0.7])
fract=fractS(h)
Dm=fractD(h)
zSim=zeros(77)
piaKu=0.
drk=0.125
for i=1:77
    global piaKu
    n1s,n2s=bisection(dmTs,Dm[i])
    n1,n2=bisection(dmT,Dm[i])
    att=(1-fract[i])*attKuTs[n1s]+fract[i]*attKuT[n1]
    piaKu=piaKu+att*drk
    zSim[i]=log10((1-fract[i])*10^(0.1*zKuTs[n1s])+(fract[i])*10^(0.1*zKuT[n1]))*10.0-piaKu
    piaKu+=att*drk
end
dm4=[]
for i=1:nt
    zKu=zKuL[ic][i,2,2,:]
    dn=zeros(176).-0.3
    inode=[binZeroDegL[ic][i][1]-8+1,binZeroDegL[ic][i][1]-2+1]
    zKuC,piaHB,piaKu,piaKa0,dn,kaCorr,
    pwcR0,rateR0,pRate,dmRet=attCorrect(zKu,dr,inode,clutFL[ic][i,2,2]+1,dn)
    dn1=dn.-0.1
    zKuC1,piaHB1,piaKu1,piaKa1,dn1,kaCorr1,
    pwcR1,rateR1,pRate,dmRet=attCorrect(zKu,dr,inode,clutFL[ic][i,2,2]+1,dn1)
    dpiadn=-(piaKa1-piaKa0)/0.1
    ddn=dpiadn*(pathAttenMSL[ic][i][1]-piaKa0)/(dpiadn*dpiadn+1.0)
    dn1=dn.+0.5*ddn
    global a=findall(dn1.<-3)
    #println(size(a))
    #println(dn1)
    if length(a)>0
        dn1[a].=-3
    end
    global a=findall(dn1.>3)
    #println(size(a))
    #println(dn1)
    if length(a)>0
        dn1[a].=3
    end

    global piaKa2d=zeros(3,3)
    global piaKu2d=zeros(3,3)
    pwcRC=0.0
    piaKaC=0.0
    rateRC=0
    piaKuC=0.0
    for i1=1:3
        for j1=1:3

            zKu=zKuL[ic][i,i1,j1,:]
            zKuC,piaHB,piaKu,piaKa,dn,
            kaCorr,pwcR,rateR,pRate,dmRet=attCorrect(zKu,dr,
                                               inode,clutFL[ic][i,i1,j1],dn1)
            piaKa2d[i1,j1]=piaKa
            piaKu2d[i1,j1]=piaKu
            if i1==2 && j1==2
                piaKaC=piaKa
                piaKuC=piaKu
                pwcRC=pwcR
                rateRC=rateR
                println(dmRet[[100,115,125,145,155,165]])
                v1=[0.125, 0.125, 0.125, 0.25, 0.25, 0.25]
                dm1=dmRet[[100,115,125,145,155,165]] .* (1 .+v1 .*randn(6))
                print(size(dm1))
                #print(dm1)
                #exit()
                if dmRet[100]>0
                    push!(dm4,dm1)
                end
                for k=1:176
                    if pRate[k]>1e-3
                        push!(pRateL[k],dmRet[k])
                    else
                        push!(pRateL[k],0.0)
                    end
                end
            end
        end
    end
    
    att=mean(10 .^(-0.1*piaKa2d))
    piaKaMean=-10*log10(att)
    att=mean(10 .^(-0.1*piaKu2d))
    piaKuMean=-10*log10(att)
    #println("$(piaKuMean/mean(piaKu2d))")
    if pwcRC>3.9 && piaKaC>175
        println(pwcRC," ",piaHB, " $(i)")
        
    end
    if pwcRC==pwcRC &&  pwcRC<5
        #println("$(piaHB) $(piaKu) $(pathAttenNSL[ic][i][1]) $(ddn)")
        #println("$(pwcR0) $(pwcR) $(piaKa0) $(piaKa)")
        push!(p1,piaKaC)
        push!(p3,piaKaMean)
        push!(p5,piaKuMean)
        push!(p4,piaL[ic][i,2])
        push!(p2,pathAttenMSL[ic][i][1])
        push!(w1,rateRC)
        push!(w2,cmbSfcRainL[ic][i][1])
    end
end
#println(piaHB)

#pRateM=zeros(176)
#for i=1:176
#    pRateM[i]=mean(pRateL[i])
#end
using PyPlot
dm4=np["array"](dm4)
dm4T=copy(transpose(dm4))
covT=np["cov"](dm4T)
using LinearAlgebra

#plot(pRateM[100:176],(176:-1:100).-100)
np=pyimport("numpy")
eVal,eVect=eigen(covT)
fP=pybuiltin("open")("dmCov.pklz","wb")
dmM=mean(dm4,dims=1)
pickle["dump"]([eVal,eVect,dmM],fP)
println([eVal,eVect,dmM])

#eVal=[0.04528470703061851
#      0.13603886360178852
#      0.20870289437033535
#      0.4869412600465882]
#eVect=[ 0.105793   0.632544  -0.360287  -0.677414
#        -0.523218  -0.609836  -0.133657  -0.580068
#        0.83457   -0.427956   0.111417  -0.328531
#        -0.136176   0.211759   0.916469  -0.310965]

#dmM=[ 0.861909  1.79766  1.8265  1.83898]

h1=(176:-1:1)*0.125
h1node=h1[[100,115,125,145,155,165]]
scipy=pyimport("scipy.ndimage") 

function simZ(Dm,fract,drk)
    nz=size(dmInt)[1]
    piaKu=0
    zSim=zeros(nz)
    rrateL=zeros(nz)
    noise=randn(nz)*1.0
    noise=scipy["gaussian_filter"](noise,5)
    for i=1:nz
        if fract[i]<0.5
            Dm[i]*=0.98
        else
            Dm[i]*=1-0.04*(1-fract[i])
        end
        n1s,n2s=bisection(dmTs,Dm[i])
        n1,n2=bisection(dmT,Dm[i])
        att=(1-fract[i])*attKuTs[n1s]+fract[i]*attKuT[n1]
        att=att*10^noise[i]
        piaKu=piaKu+att*drk
        zSim[i]=log10((1-fract[i])*10^(0.1*zKuTs[n1s])+(fract[i])*10^(0.1*zKuT[n1]))*10.0-piaKu+10*noise[i]
        rrateL[i]=10^noise[i]*((1-fract[i])*rateTs[n1s]+fract[i]*rateT[n1])
        piaKu+=att*drk
    end
    return piaKu,rrateL,zSim
end
nT=30000
zSimL=zeros(3*nT,66)
pRateL=zeros(3*nT,66)

eVal1=eVal.*2.0
fract=fractS(h1[100:165])
for i=1:nT
    global dmP=transpose(randn(1)*eVal1[6].*eVect[6,:]+randn(1)*eVal1[5].*eVect[5,:]+randn(1)*eVal1[4].*eVect[4,:]+
    randn(1)*eVal1[3].*eVect[3,:]+randn(1)*eVal1[2].*eVect[2,:]+randn(1)*eVal1[1].*eVect[1,:])+dmM
    dmS=Spline1D(h1node[4:-1:1],dmP[4:-1:1],k=2)
    global dmInt=dmS(h1[100:165]).+randn(66)*0.2
    nz=66
    noise=randn(nz)*1.0
    noise=scipy["gaussian_filter"](noise,10)
    dmInt=dmInt.+0.5*noise
    piaKu,rrate,zSim=simZ(0.85*dmInt,fract,drk)
    zSimL[i,:]=zSim
    pRateL[i,:]=rrate
end
nz=66
for i=1:nT
    global dmP=transpose(randn(1)*eVal1[6].*eVect[6,:]+randn(1)*eVal1[5].*eVect[5,:]+randn(1)*eVal1[4].*eVect[4,:]+
    randn(1)*eVal1[3].*eVect[3,:]+randn(1)*eVal1[2].*eVect[2,:]+randn(1)*eVal1[1].*eVect[1,:])+dmM
    dmS=Spline1D(h1node[4:-1:1],dmP[4:-1:1],k=2)
    global dmInt=dmS(h1[100:165]).+randn(66)*0.2
    noise=randn(nz)*1.0
    noise=scipy["gaussian_filter"](noise,10)
    dmInt=dmInt.+0.5*noise
    piaKu,rrate,zSim=simZ(0.95*dmInt,fract,drk)
    zSimL[i+nT,:]=zSim
    pRateL[i+nT,:]=rrate
end

for i=1:nT
    global dmP=transpose(randn(1)*eVal1[6].*eVect[6,:]+randn(1)*eVal1[5].*eVect[5,:]+randn(1)*eVal1[4].*eVect[4,:]+
    randn(1)*eVal1[3].*eVect[3,:]+randn(1)*eVal1[2].*eVect[2,:]+randn(1)*eVal1[1].*eVect[1,:])+dmM
    dmS=Spline1D(h1node[4:-1:1],dmP[4:-1:1],k=2)
    global dmInt=dmS(h1[100:165]).+randn(66)*0.2
    noise=randn(nz)*1.0
    noise=scipy["gaussian_filter"](noise,10)
    dmInt=dmInt.+0.5*noise
    piaKu,rrate,zSim=simZ(1.15*dmInt,fract,drk)
    zSimL[i+2*nT,:]=zSim
    pRateL[i+2*nT,:]=rrate
end


zSimLx=xr["DataArray"](zSimL)
pRateLx=xr["DataArray"](pRateL)
d=Dict("zSim"=>zSimLx,"pRate"=>pRateLx)
trainDataset=xr["Dataset"](d)
trainDataset["to_netcdf"]("simZDataset_2.nc")


#plot(zm[1,:],h1[100:165])
#plot(zm[1,:]-zs[1,:],h1[100:165])
#plot(zm[1,:]+zs[1,:],h1[100:165])

kneigh=pyimport("sklearn.neighbors")
neigh = kneigh["KNeighborsClassifier"](n_neighbors=3)

#>>> neigh.fit(X, y)
