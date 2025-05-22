rmse=function(x,par,p){
  if(p==1)	out=var(x)+(mean(x)-par)**2 else{ 
    out=numeric()	
    for(i in 1:p) out[i]=var(x[,i])+(mean(x[,i])-par[i])**2
  }
  return(sqrt(out))
}

ARsimex2=function(out,yerr=0,delta,p=1,up=0.99,N=8,type="linear",weighted="equally",alpha=0.1,sigma=.05,criteria="fit"){
  j=1:N
  x=(j)*delta
  if(mean(yerr)!=0)
  {
    j=1:(N+1)
    x=(0:N)*delta
  }
  if(criteria=="extrapolation")
  {
    j=2:N
    x=(2:N)*delta
    out=matrix(out[j,1:p],ncol=p)
  }
  sal=numeric()
  for(i in 1:p){
    N1=length(j)
    y=out[1:N1,i]
    w=rep(1/N1,N1)
    if(weighted=="exponential")
    {
      w=alpha*(1-alpha)^(j-1)
      w=w/sum(w)
    }
    if(weighted=="gaussian")
    {
      w=dnorm(x=seq(delta,delta*N1,length=N1),mean=delta*N1/2,sd=2*delta)
      w=w/sum(w)    
    }
    if(weighted=="invexponential")
    {
      w=alpha*(1-alpha)^(j-1)
      w=sort(w)/sum(w)
    }
    data=data.frame(cbind(y,x))
    if(type=="linear")
      fit=lm(y~x,data=data,weights=w)
    if(type=="quadratic")
    {
      data$x2 <- data$x^2
      fit=lm(y~x+x2,data=data,weights=w)
    }
    if(type=="cubic")
    {
      data$x2 <- data$x^2
      data$x3 <- data$x^3
      fit=lm(y~x+x2+x3,data=data,weights=w)
    }
    if(type=="bi-quadratic")
    {
      data$x2 <- data$x^2
      data$x3 <- data$x^3
      data$x4 <- data$x^4
      fit=lm(y~x+x2+x3+x4,data=data,weights=w)
    }
    sal[i]=fit$coeff[1]
    if(mean(yerr)!=0)
    {
      test=data.frame(x=-delta,x2=(-delta)**2,x3=(-delta)**3,x4=(-delta)**4)
      sal[i]=predict(fit,new=test)
    }
    if(criteria=="extrapolation")
    {
      test=data.frame(x=delta,x2=(delta)**2,x3=(delta)**3,x4=(delta)**4)
      sal[i]=predict(fit,new=test)
    }
  }
  pol=-as.polynomial(c(sal[p:1],-1))
  roots=solve(pol)
  count=0
  for(i in 1:p) if(Mod(roots)[i]<1) count=count+1 else count=0
  if(count==p) return(list(par=sal,fit=fit$fitted,AIC=AIC(fit),r.squared=summary(fit)$r.squared)) else 
  {
    roots1=rep(up,p)
    sal=-poly.calc(roots1)[p:1]
    return(list(par=sal,fit=fit$fitted,AIC=AIC(fit),r.squared=summary(fit)$r.squared))
  }
}

Bestdelta=function(y,st,N,p=1,B=1,delta0=1){
  n=length(y)
  outF=matrix(NA,nr=1,nc=B)
  delta=delta0
  fit0=IARloglik(y=y,st=st)
  exit=-1
  while(exit<0){
    for(b in 1:B)
    {
      x=N*delta
      eps=x*rnorm(n)
      X=y+eps
      tryCatch({fit=IARloglik(y=X,st=st)}, error = function(e){print(e)})
      outF[,b]=fit$phi
    }
    fin=apply(outF,2,mean)
    if(fin<fit0$phi/2)
      delta=delta/2
    else
      exit=0
    dist=fit0$phi-fin
  }
  return(list(delta=delta,phimle=fit0$phi,dist=dist))		
}

BestSimex=function(out,delta,yerr=0,p=1,up=0.9999,criteria="AIC"){
  N=length(out)
  if(mean(yerr)!=0)
    N=N-1
  AIC=Inf
  R2=0
  j=1
  AICsimex=rep(NA,16)
  R2simex=rep(NA,16)
  MSE=Inf
  MSEsimex=rep(NA,16)
  for(type in c("linear","quadratic","cubic","bi-quadratic"))
  {
    for(weight in c("equally","gaussian","exponential","invexponential"))
    {
      if(criteria != "extrapolation")
      {
        simex=ARsimex2(out,delta=delta,yerr=yerr,p=p,up=up,N=N,type=type,weighted=weight)
        AICsimex[j]=simex$AIC
        R2simex[j]=simex$r.squared
      }
      if(criteria=="AIC" & AICsimex[j] < AIC)
      {
        AIC=AICsimex[j]
        weightsimex=weight
        typesimex=type
        parsimex=simex$par
        fitsimex=simex$fit
      }
      if(criteria=="r.squared" & R2simex[j] > R2)
      {
        R2=R2simex[j]
        weightsimex=weight
        typesimex=type
        parsimex=simex$par
        fitsimex=simex$fit
      }
      if(criteria=="extrapolation")
      {
        simex=ARsimex2(out,delta=delta,yerr=yerr,p=p,up=up,N=N,type=type,weighted=weight,criteria="extrapolation")
        MSEsimex[j]=(out[1]-simex$par)**2
        if(MSEsimex[j] < MSE)
        {
          MSE=MSEsimex[j]
          weightsimex=weight
          typesimex=type
        }
      }
      j=j+1
    }
  }
  crit=ifelse(criteria=="AIC",AIC,R2)
  if(criteria=="extrapolation")
  {
    simex=ARsimex2(out,delta=delta,yerr=yerr,p=p,up=up,N=N,type=typesimex,weighted=weightsimex)
    crit=MSE
    parsimex=simex$par
    fitsimex=simex$fit
  }
  return(list(par=parsimex,criteria=crit,type=typesimex,weight=weightsimex,fit=fitsimex,R2=R2simex,AIC=AICsimex))
}

SimexiAR=function(y,st,yerr=0,N,delta,p=1,B=1){
  out=matrix(NA,nr=N,nc=p)
  n=length(y)
  outF=matrix(NA,nr=N,nc=B)
  jseq=1:N
  if(mean(yerr)!=0)
  {
    out=matrix(NA,nr=N+1,nc=p)
    outF=matrix(NA,nr=N+1,nc=B)
    jseq=0:N
  }
  for(b in 1:B)
  {
    for(j in jseq)
    {
      x=j*delta
      eps=x*rnorm(n)
      if(mean(yerr)!=0)
        eps=x*rnorm(n,0,mean(yerr))  
      X=y+eps
      tryCatch({fit=IARloglik(y=X,st=st)}, error = function(e){print(e)})
      if(mean(yerr)==0)
        out[j,]=fit$phi
      if(mean(yerr)!=0)
        out[j+1,]=fit$phi
    }
    outF[,b]=out[,1]
  }
  fin=matrix(apply(outF,1,mean),nr=length(jseq),nc=p)
  return(fin)		
}


LSEestimation <- function(y, st,seed=NULL){
  set.seed(seed)
  y1 <- y[-1]
  y2 <- lag(y)[-1]
  delta <- diff(st)
  mylse <- function(x){
    z=y1-x^delta*y2
    z2=sum(z^2)
    return(z2)
  }
  u <- runif(1)
  out <- as.numeric(NA)
  out <- optim(u,mylse,lower = 0,upper = 1,method = 'L-BFGS-B')
  par <-out$par
  lse <-out$value
  return(list(par=par,lse=lse))
}

