include("psdInt_2.jl")
function fromZKuTs(zKuC,dn)
    if zKuC-10*dn<-7
        dn=(zKuC+7)/(10)
    end
    if zKuC-10*dn>48.75
        dn=(zKuC-48.75)/(10)
    end
    zKuC1=zKuC-10*dn
    n1,n2=bisection(zKuTs,zKuC1)
    dZ=zKuTs[n2]-zKuTs[n1]+1e-5
    f=(zKuC1-zKuTs[n1])/dZ
    attKu=10^dn*((1-f)*attKuTs[n1]+f*attKuTs[n2])
    attKa=10^dn*((1-f)*attKaTs[n1]+f*attKaTs[n2])
    zKa=(1-f)*zKaTs[n1]+f*zKaTs[n2]+10*dn
    pwcRet=10^dn*((1-f)*pwcTs[n1]+f*pwcTs[n2])
    rateRet=10^dn*((1-f)*rateTs[n1]+f*rateTs[n2])
    dm=dmTs[n1]
    return attKu,attKa,zKa,pwcRet,rateRet,dm
end
function fromZKuT(zKuC,dn)
    if zKuC-10*dn<0
        dn=(zKuC+0)/(10)
    end
    if zKuC-10*dn>54.75
        dn=(zKuC-54.75)/(10)
    end
    zKuC1=zKuC-10*dn
    n1,n2=bisection(zKuT,zKuC1)
    dZ=zKuTs[n2]-zKuTs[n1]+1e-5
    f=(zKuC1-zKuT[n1])/dZ
    attKu=10^dn*((1-f)*attKuT[n1]+f*attKuT[n2])
    attKa=10^dn*((1-f)*attKaT[n1]+f*attKaT[n2])
    zKa=(1-f)*zKaT[n1]+f*zKaT[n2]+10*dn
    pwcRet=10^dn*((1-f)*pwcT[n1]+f*pwcT[n2])
    rateRet=10^dn*((1-f)*rateT[n1]+f*rateT[n2])
    return attKu,attKa,zKa, pwcRet,rateRet,dmT[n1]
end

function attCorrect(zKu,dr,inode,clutF,dn)
    piaKu=0
    attKu1=0.0
    attKu2=0
    nz=size(zKu)[1]
    zKuC=zeros(nz).+zKu
    zeta1d=zeros(nz)

    q=0.2*log(10)
    beta=0.72
    eps=1
    pwcR=0.0
    rateR=0.0
    piaKa=0.0
    prate=zeros(nz)
    piaMax=0.0
    dmRet=zeros(nz)
    for it=1:5
        zeta=0.0

    global kaCorr=zeros(nz)
        piaKa=0.0
        piaKu=0.0

        for k=1:inode[1]
            attKu=0.0
            attKa=0.0
            if zKu[k]==zKu[k] && zKu[k]>0
                attKu,attKa,zKa,pwcS,rateS,dm=fromZKuTs(zKuC[k],dn[k])
                prate[k]=rateS
                piaKu=piaKu+2*attKu*dr/1e3
                piaKa=piaKa+2*attKa*dr/1e3
                kaCorr[k]=kaCorr[k]+piaKa
                zeta=zeta+attKu/10.0^(0.1*zKuC[k]*beta)*10.0^(0.1*zKu[k]*beta)*dr/1e3
                zeta1d[k]=zeta
                dmRet[k]=dm
                #println("$k, $(zKuC[k]) $(attKu)  $(zeta1d[k])")
            end
            if k==inode[1]
                attKu1=attKu
            end
        end
        
        for k=inode[1]+1:inode[2]-1
            if zKuC[k]>0
                f=(inode[2]-k)/(inode[2]-inode[1])
                if zKuC[k]-10*dn[k]<-7
                    dn[k]=(zKuC[k]+7)/(10)
                end
                if zKuC[k]-10*dn[k]>48.75
                    dn[k]=(zKuC[k]-48.75)/(10)
                end
                #println("$k $f $(size(dn))")
                n1,n2=bisection2(zKuTs,zKuT,f,zKuC[k]-10*dn[k])
                #println("$(n1) $(n2)")
                attKuS,attKaS,zKaS,pwcS,rateS=fromZKuTs(log10(f)*10+zKuC[k],dn[k])
                attKuR,attKaR,zKaR,pwcR,rateR=fromZKuT(log10(1-f)*10+zKuC[k],dn[k])
                dn1=dn[k]
                f1=0.
                attKuS=10^dn1*((1-f1)*attKuTs[n1]+f1*attKuTs[n2])
                attKaS=10^dn1*((1-f1)*attKaTs[n1]+f1*attKaTs[n2])
                attKuR=10^dn1*((1-f1)*attKuT[n1]+f1*attKuT[n2])
                attKaR=10^dn1*((1-f1)*attKaT[n1]+f1*attKaT[n2])
                attKu=f*attKuS+(1-f)*attKuR
                attKa=f*attKaS+(1-f)*attKaR
                rateS=10^dn1*rateTs[n1]
                rateR=10^dn1*rateT[n1]
                prate[k]=f*rateS+(1-f)*rateR
                zeta=zeta+attKu/10.0^(0.1*zKuC[k]*beta)*10.0^(0.1*zKu[k]*beta)*dr/1e3
                piaKu=piaKu+2*attKu*dr/1e3
                piaKa=piaKa+2*attKa*dr/1e3
            end
            zeta1d[k]=zeta
        end
        for k=inode[2]:clutF
            attKu=0.0
            attKa=0.0
            if zKu[k]==zKu[k] && zKuC[k]>0
                attKu,attKa,zKa,pwcR,rateR,dm=fromZKuT(zKuC[k],dn[k])
                piaKu=piaKu+2*attKu*dr/1e3
                piaKa=piaKa+2*attKa*dr/1e3
                zeta=zeta+attKu/10.0^(0.1*zKuC[k]*beta)*10.0^(0.1*zKu[k]*beta)*dr/1e3
                if rateR>400 && it==5
                    println("$(zKuC[k]) $(zKu[k]) $(dn[k]) $(rateR) $(piaMax)")
                    #exit(0)
                end
                dmRet[k]=dm
                prate[k]=rateR
                #zeta1d[k]=zeta
                #println("$k $(zKuC[k]) $(attKu)  $(zeta1d[k])")
            end
             zeta1d[k]=zeta
        end
        piaMax=55-zKu[clutF]
        
        if 1>0 #q*beta*zeta1d[clutF]>0.99999 
            eps=(1-10^(-0.1*piaMax*beta))/(q*beta*zeta1d[clutF])
            #println("eps =$(eps) $(log10(eps^(1.0/(1-beta))))")
            #println("$(zKuC[clutF]) $(pwcR) $(dn[clutF])")
            if eps>1.0
                eps=1.0
            end
            dn=dn.+log10(eps^(1.0/(1-beta)))
        end
        
        
        for k=1:clutF
            if zKu[k]==zKu[k] && zKu[k]>0
                zKuC[k]=zKu[k]-10/beta*log10(1-eps*q*beta*zeta1d[k])
            end
        end
        global piaHB=-10/beta*log10(1-eps*q*beta*zeta1d[clutF])
    end
    return zKuC,piaHB,piaKu,piaKa,dn,kaCorr,pwcR,rateR, prate,dmRet
end
