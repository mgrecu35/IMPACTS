using Statistics
using LinearAlgebra
function covar(x,y)
    nt,nx=size(x)[1],size(x)[2]
    nt,ny=size(y)[1],size(y)[2]
    covxy=zeros(nx,ny)
    for i=1:nx
        for j=1:ny
            xm=mean(x[:,i])
            ym=mean(y[:,j])
            covxy[i,j]=0
            for k=1:nt
                covxy[i,j]=covxy[i,j]+(x[k,i]-xm)*(y[k,j]-ym)
            end
            covxy[i,j]=covxy[i,j]/((nx-1)*(ny-1.0))
        end
    end
    return covxy
end
function enKF(xEns,yEns,yObs,nm)
    covyy=covar(yEns,yEns)
    covxy=covar(xEns,yEns)
    ny=size(covyy)[1]
    sigma=3.0
    e,v=eigen(covyy+sigma*Diagonal(ones(ny)))
    ymean=transpose(mean(yEns,dims=1))
    xmean=transpose(mean(xEns,dims=1))
    invCovyy=v[:,nm:end]*Diagonal(e[nm:end])*transpose(v[:,nm:end])
    x=xmean+covxy*invCovyy*(yObs-ymean)
    return x
end

function smooth(x)
    nx=size(x)[1]
    xs=copy(x)
    for i=2:nx-1
        xs[i]=0.1*x[i-1]+0.1*x[i+1]+0.8*x[i]
    end
    xs[nx]=0.9*x[nx]+0.1*x[nx-1]
    xs[1]=0.9*x[1]+0.1*x[2]
    return xs
end
#v*Diagonal(1 ./e)*transpose(v)=inv(covxx)
