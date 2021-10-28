##################
### Functions ####
##################

#Background Brownian motion in 2D from t0-tend with steps dt 
BM2D<-function(t0,tend,dt,sigma){
    t <- seq(t0,tend,dt)  # time
    nt<-ceiling(tend/dt)+1
    
    x <- c(rnorm(n = nt , sd = sigma*sqrt(dt)))
    y <- c(rnorm(n = nt , sd = sigma*sqrt(dt)))
    x <- cumsum(x)
    y <- cumsum(y)
    BM<-cbind(x,y)
    return(BM)
}

#Spherical Brownian motion 
SBM3Exact<-function(r,n,D,dt){
    d<-3
    tau<-2*D*dt/r^2
    
    update<-rnorm(3)
    update<-(update/sqrt(sum(update^2)))*r
    
    save<-update
    
    t<-2*D*dt/r^2
    beta<-(d-2)*t/2
    eta<- beta/(exp(beta)-1)
    mu<-2*eta/t
    var<-(2*eta/t)*(eta+beta)^2*(1+eta/(eta+beta)-2*eta)*(beta^(-2))
    
    M<-rnorm(n-1,mu,sd=sqrt(var))
    
    X<-rbeta(n-1,(d-1)/2,(d-1)/2+M)
    
    Y<-t(matrix(rnorm(2*(n-1)),nrow=2))
    Y<-Y/sqrt(rowSums(Y^2))
    
    for(i in 1:(n-1)){
        ed<-c(0,0,1)
        u<-(ed-update/r)/sqrt(sum((ed-update/r)^2))
        
        O<-diag(c(1,1,1))-2*u%*%t(u)
        
        update<-c(r*O%*% (c(2*sqrt(X[i]*(1-X[i]))%*%(t(Y[i,])),(1-2*X[i]))))
        save<-rbind(save,update)
    }
    
    return(save)
}


#Rotational diffusion of spherical
Rotate<-function(n,r,D,dt){
    d<-3
    update<-c(0,0,1)
    save<-update
    
    t<-2*D*dt/r^3
    beta<-(d-2)*t/2
    eta<- beta/(exp(beta)-1)
    mu<-2*eta/t
    var<-(2*eta/t)*(eta+beta)^2*(1+eta/(eta+beta)-2*eta)*(beta^(-2))
    
    M<-rnorm(n-1,mu,sd=sqrt(var))
    
    X<-rbeta(n-1,(d-1)/2,(d-1)/2+M)
    
    Y<-t(matrix(rnorm(2*(n-1)),nrow=2))
    Y<-Y/sqrt(rowSums(Y^2))
    
    ed<-c(0,0,1)
    u<-(ed-update/1)
    
    O<-diag(c(1,1,1))-2*u%*%t(u)
    
    update<-c(1*O%*% (c(2*sqrt(X[1]*(1-X[1]))%*%(t(Y[1,])),(1-2*X[1]))))
    save<-rbind(save,update)
    
    for(i in 2:(n-1)){
        ed<-c(0,0,1)
        u<-(ed-update/1)/sqrt(sum((ed-update/1)^2))
        
        O<-diag(c(1,1,1))-2*u%*%t(u)
        
        update<-c(1*O%*% (c(2*sqrt(X[i]*(1-X[i]))%*%(t(Y[i,])),(1-2*X[i]))))
        save<-rbind(save,update)
      }
    
    return((save))
}

#Simulate particle with E[N]=lambda, radius r, from t0 to tend with timesteps dt including lateral difussion of the motors and rotational difussion
SBM3D<-function(lambda, r,t0,tend,dt,D,D_R,threshold){
    n<-rpois(1,lambda)
    t <- seq(t0,tend,dt)  # time
    
    nt<-ceiling(tend/dt)+1
    
    if(D>0){
        test<-lapply(1:n, function(x){SBM3Exact(r,nt,D,tend/nt)})
        p.l<-list(do.call("rbind", lapply(test, function(x){x[1,]})))
        drift<-colSums(p.l[[1]])/r
        for (i in 1:ceiling(tend/dt)){
            p.l[[i+1]]<-do.call("rbind", lapply(test, function(x){x[1+i,]}))
            drift<-rbind(drift,colSums(p.l[[i+1]])/r)
        }
    }
    
    else{
        x <- data.frame(x=rnorm(n),y=rnorm(n),z=rnorm(n))
        
        p<-r*x/sqrt(rowSums(x^2))
        drift<-colSums(p)/r
        
        drift<- matrix(rep(drift,1+ceiling(tend/dt)),
                       nrow=1+ceiling(tend/dt),
                       byrow=T)    
    }
    
    if(D_R>0){
        
        r0<-sqrt(rowSums(drift^2))
        drift<-drift/r0
        rot<-Rotate(nt,r,D_R,dt)
        
        ed<-c(0,0,1)
        u<-t(ed-t(drift[2,])/1)/sqrt(sum(t(ed-t(drift[2,])/1)^2))
        O<-diag(c(1,1,1))-2*u%*%t(u)
        drift[2,]<-c(1*O%*% rot[2,])
        
        
        
        for(i in 2:(nt-1)){
            u<-t(ed-t(drift[(i+1),])/1)/sqrt(sum(t(ed-t(drift[(i+1),])/1)^2))
            O<-diag(c(1,1,1))-2*u%*%t(u)
            drift[(i+1),]<-c(1*O%*% rot[(i+1),])
        }
        
        
        
        drift<-drift*r0
    }
    
    
    if(threshold>0){
        vc<-(min(n,threshold)/n)
    }
    else{
        vc<-1
    }
    
    return(vc*drift)
}


##################
### Simulation ###
##################

##Simulate motility as a result of lateral and rotational diffusion##
dt<-0.04
t0<-0
tend<-10
lambda<-100
v.med<-1
vc<-v.med/(sqrt((2/3)*lambda))
threshold<-0

#D_L
sigma<-0.03

MSD3a<-list()
BM.MSD3a<-list()
M.MSD3a<-list()

#k_{B}T/8*pi*eta
sigmaR3a<-c(0,0.05,0.1,0.25,0.5,3)

j<-0
for (r in c(0.6,1,1.8)){
    for(sigma_R in sigmaR3a){
        #DT
        sigma.BM<-sqrt(2*(0.14/r))
        
        
        j<-j+1
        print(sigma_R)
        M.speed<- SBM3D(lambda, r,t0,tend,dt,sigma,sigma_R,threshold)
        M.drift<-cbind((vc)*dt*cumsum(M.speed[,1]),(vc)*dt*cumsum(M.speed[,2]))
        
        M.MSD.S<-c(M.drift[,1]^2+M.drift[,2]^2)
        
        BM.drift<-BM2D(t0,tend,dt,sigma.BM)
        BM.MSD.S<-c(BM.drift[,1]^2+BM.drift[,2]^2)
        
        T.drift<-BM.drift+M.drift
        MSD.S<-c(T.drift[,1]^2+T.drift[,2]^2)
        
        for (i in 1:9999){
            if(i%%100==0){
                print(i)
            }
            #
            M.speed<- SBM3D(lambda, r,t0,tend,dt,sigma,sigma_R,threshold)
            M.drift<-cbind((vc)*dt*cumsum(M.speed[,1]),(vc)*dt*cumsum(M.speed[,2]))
            
            M.MSD.S<-rbind(M.MSD.S,M.drift[,1]^2+M.drift[,2]^2)
            
            BM.drift<-BM2D(t0,tend,dt,sigma.BM)
            BM.MSD.S<-rbind(BM.MSD.S,BM.drift[,1]^2+BM.drift[,2]^2)
            
            T.drift<-BM.drift+M.drift
            MSD.S<-rbind(MSD.S,T.drift[,1]^2+T.drift[,2]^2)
            
            
        }
        
        M.MSD3a[[j]]<-M.MSD.S
        BM.MSD3a[[j]]<-BM.MSD.S
        MSD3a[[j]]<-MSD.S

    }
}



plot(seq(t0,tend,dt),c(0,head(apply(MSD3a[[1]],2,function(x){mean(x)}),length(seq(t0+dt,tend,dt)))),type="l",ylab="MSD",xlab="t",lty=1,lwd=2,col="black")
lines(seq(t0,tend,dt),c(0,head(apply(MSD3a[[2]],2,function(x){mean(x)}),length(seq(t0+dt,tend,dt)))),type="l",ylab="MSD",xlab="t",lty=2,lwd=2,col="black")
lines(seq(t0,tend,dt),c(0,head(apply(MSD3a[[3]],2,function(x){mean(x)}),length(seq(t0+dt,tend,dt)))),type="l",ylab="MSD",xlab="t",lty=3,lwd=2,col="black")
lines(seq(t0,tend,dt),c(0,head(apply(MSD3a[[4]],2,function(x){mean(x)}),length(seq(t0+dt,tend,dt)))),type="l",ylab="MSD",xlab="t",lty=4,lwd=2,col="black")
lines(seq(t0,tend,dt),c(0,head(apply(MSD3a[[5]],2,function(x){mean(na.omit(x))}),length(seq(t0+dt,tend,dt)))),type="l",ylab="MSD",xlab="t",lty=5,lwd=2,col="black")
lines(seq(t0,tend,dt),c(0,head(apply(MSD3a[[6]],2,function(x){mean(na.omit(x))}),length(seq(t0+dt,tend,dt)))),type="l",ylab="MSD",xlab="t",lty=6,lwd=2,col="black")

r<-0.6
Dr<-c(0,0.05,0.1,0.25,0.5,3)
Dstart<-Dr/r^3+(1+2/lambda)*0.03/r^2
v<-vc*sqrt(2/3)*sqrt(lambda)
for (i in Dstart){
    lines(seq(t0,tend,dt),sapply(seq(t0,tend,dt),function(x){4*x*(0.14/r)+0.5*v^2*(1/((i)))^2*(2*x*i+exp(-2*x*i)-1)}),type="l",ylab='y',col="green")
}


legend(0.5,65,c("0","0.05","0.1","0.25","0.5","3"),title=expression('D'[R]),lty=c(1,2,3,4,5,6))

plot(seq(t0,tend,dt),c(0,head(apply(MSD3a[[7]],2,function(x){mean(x)}),length(seq(t0+dt,tend,dt)))),type="l",ylab="MSD",xlab="t",lty=1,lwd=2,col="black")
lines(seq(t0,tend,dt),c(0,head(apply(MSD3a[[8]],2,function(x){mean(x)}),length(seq(t0+dt,tend,dt)))),type="l",ylab="MSD",xlab="t",lty=2,lwd=2,col="black")
lines(seq(t0,tend,dt),c(0,head(apply(MSD3a[[9]],2,function(x){mean(x)}),length(seq(t0+dt,tend,dt)))),type="l",ylab="MSD",xlab="t",lty=3,lwd=2,col="black")
lines(seq(t0,tend,dt),c(0,head(apply(MSD3a[[10]],2,function(x){mean(x)}),length(seq(t0+dt,tend,dt)))),type="l",ylab="MSD",xlab="t",lty=4,lwd=2,col="black")
lines(seq(t0,tend,dt),c(0,head(apply(MSD3a[[11]],2,function(x){mean(x)}),length(seq(t0+dt,tend,dt)))),type="l",ylab="MSD",xlab="t",lty=5,lwd=2,col="black")
lines(seq(t0,tend,dt),c(0,head(apply(MSD3a[[12]],2,function(x){mean(na.omit(x))}),length(seq(t0+dt,tend,dt)))),type="l",ylab="MSD",xlab="t",lty=6,lwd=2,col="black")

r<-1
Dr<-c(0,0.05,0.1,0.25,0.5,3)
Dstart<-Dr/r^3+(1+2/lambda)*0.03/r^2
v<-vc*sqrt(2/3)*sqrt(lambda)
for (i in Dstart){
    lines(seq(t0,tend,dt),sapply(seq(t0,tend,dt),function(x){4*x*(0.14/r)+0.5*v^2*(1/((i)))^2*(2*x*i+exp(-2*x*i)-1)}),type="l",ylab='y',col="green")
}

legend(0.5,80,c("0","0.05","0.1","0.25","0.5","3"),title=expression('D'[R]),lty=c(1,2,3,4,5,6))


plot(seq(t0,tend,dt),c(0,head(apply(MSD3a[[13]],2,function(x){mean(x)}),length(seq(t0+dt,tend,dt)))),type="l",ylab="MSD",xlab="t",lty=1,lwd=2,col="black")
lines(seq(t0,tend,dt),c(0,head(apply(MSD3a[[14]],2,function(x){mean(x)}),length(seq(t0+dt,tend,dt)))),type="l",ylab="MSD",xlab="t",lty=2,lwd=2,col="black")
lines(seq(t0,tend,dt),c(0,head(apply(MSD3a[[15]],2,function(x){mean(x)}),length(seq(t0+dt,tend,dt)))),type="l",ylab="MSD",xlab="t",lty=3,lwd=2,col="black")
lines(seq(t0,tend,dt),c(0,head(apply(MSD3a[[16]],2,function(x){mean(x)}),length(seq(t0+dt,tend,dt)))),type="l",ylab="MSD",xlab="t",lty=4,lwd=2,col="black")
lines(seq(t0,tend,dt),c(0,head(apply(MSD3a[[17]],2,function(x){mean(x)}),length(seq(t0+dt,tend,dt)))),type="l",ylab="MSD",xlab="t",lty=5,lwd=2,col="black")
lines(seq(t0,tend,dt),c(0,head(apply(MSD3a[[18]],2,function(x){mean(na.omit(x))}),length(seq(t0+dt,tend,dt)))),type="l",ylab="MSD",xlab="t",lty=6,lwd=2,col="black")

r<-1.8
Dr<-c(0,0.05,0.1,0.25,0.5,3)
Dstart<-Dr/r^3+(1+2/lambda)*0.03/r^2
v<-vc*sqrt(2/3)*sqrt(lambda)
for (i in Dstart){
    lines(seq(t0,tend,dt),sapply(seq(t0,tend,dt),function(x){4*x*(0.14/r)+0.5*v^2*(1/((i)))^2*(2*x*i+exp(-2*x*i)-1)}),type="l",ylab='y',col="green")
}

legend(0.5,90,c("0","0.05","0.1","0.25","0.5","3"),title=expression('D'[R]),lty=c(1,2,3,4,5,6))


