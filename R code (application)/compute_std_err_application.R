compute_std_err_application <- function(bhat, listy, listreg, listZ, listZy, listsector, Theta_value, Sigma_value, Sigma_valueg, listdummy, nnuts,nreg, nsectors,ES, cov2Zy, listmu,listsigmatr)
     {

  f <- function(x, y){
    sum(x*y)  }
  
    f1 <- function(x, y){
    sum(x*t(y))}
  
    Omega_newg<-solve(Sigma_valueg)
    prova=vech(Omega_newg)
    nsigma=length(prova[prova!=0])
    listJ=list()
    jj=0
    for (g in  1:nsectors){
      for (h in  g:nsectors){
        if (Omega_newg[g,h]!=0)
{        jj=jj+1
        listJ[[jj]]=matrix(0,nsectors,nsectors)
        listJ[[jj]][g,h]=1
        listJ[[jj]][h,g]=1
}
      }}
    
  lambda3=list()
  lambda4=list()
  
  Bbetabeta=0
  Bbs=matrix(0,  nsigma,nreg)
  Bss=matrix(0,  nsigma,  nsigma)
  
  # the three lines below only for forbalanced datasets as in MC simulations but not in real data, 
  # in real data need to calcualet A ld and ApltA for each i
  
    for (i in seq(1, nnuts)){  
    mu=listmu[[i]]
    sigma=listsigmatr[[i]]

    ld=listdummy[[i]]

    Zy=listZy[[i]] 
    Z=listZ[[i]] 
    
    delta1 <- (Theta_value[listy[[i]]] - mu) / sigma
    delta2 <- (Theta_value[listy[[i]]+1] - mu) / sigma
    tmp1 <- (dnorm(delta1) - dnorm(delta2)) / (pnorm(delta2) - pnorm(delta1))
    
    # moments of truncated normals - horrace 2015 formula A3-A6 in the appendix
    EX <- mu + tmp1 * sigma
    EXX <- sigma^2-sigma*tmp1*EX
    EXXX =sigma*tmp1*(EX^2-EXX)
    EXXXX=2*sigma^4-3*(sigma*tmp1*EX)^2-sigma*(tmp1^(-1))*EXXX+(mu^2)*EXX

    lambda3=diag(EXXX/(EXX^(3/2)))
    lambda4=diag(EXXXX/(EXX^2)-3)
    EXX05=diag(EXX^(0.5))
    EXX=diag(EXX)
#   A=lapply(1:nsigma,function(x,listJ,ld,m)  ld%*%listJ[[x]]%*%t(ld)/m^2,ld=ld,listJ=listJ,m=m)
#   ApltA=lapply(1:nsigma,function(x,A)  A[[x]]+t(A[[x]]),A=A)
    
    ld=listdummy[[i]]
    
    ImZy=diag(length(Zy))-ld%*%diag(norm[[i]]^(-1))%*%t(ld)
    Bbb1=t(listreg[[i]])%*%ImZy%*%EXX%*%t(ImZy)%*%listreg[[i]]
    Bss3a=list()
    Bss3b=list()
    Bss5a=list()
    prova1=matrix(0,nsigma,nsigma)
    prova2=matrix(0,nsigma,nsigma)
    Bbs1a=matrix(0,nsigma, nreg)
    Bbs2a=matrix(0,nsigma, nreg)
    
    for (x in 1:nsigma){
    A=ld%*%diag(norm[[i]]^(-1))%*%listJ[[x]]%*%diag(norm[[i]]^(-1))%*%t(ld)
    ApltA=0.5*(A+t(A))
    
    # now start calculation of the B matrix:
    
    # this quantity is useful to compute A.37: (I-Z*Z'/m)
    
    # Bbb1 is the second bit in A.40, useful to calculate Bbetabeta
    
    # elements of Bbetasigma: (A.43)
    Bbs1a[x,]=t(t(listreg[[i]])%*%ImZy%*%EXX05%*%diag(diag(lambda3*EXX05*ApltA*EXX05))%*%array(1,length(Zy)))
    Bbs2a[x,]=t(t(listreg[[i]])%*%ImZy%*%EXX%*%ApltA%*%Zy)
    
    # elements of Bsigmasigma (A.44-A.48)
    #A.44
    Bss2a=EXX%*%ApltA
    prova1[x,x]=2*sum(Bss2a*t(Bss2a))
    #A.46-A.47
    Bss3a[[x]]=array(diag(lambda3*EXX05*ApltA*EXX05))
    Bss3b[[x]]=array(EXX05%*%ApltA%*%Zy)
     #A.44
    Bss4a= lambda4*EXX05*ApltA*EXX05
    Bss4b=EXX05*ApltA*EXX05
    prova2[x,x]= sum(Bss4a*t(Bss4b))
    #A.48
    Bss5a[[x]]=t(Zy)%*%ApltA%*%EXX05
    
    }
     
    Bss3=outerList(Bss3a,Bss3b,f) 
    prova3=2*Bss3+2*t(Bss3)
    prova4=4*outerList(Bss5a,Bss5a,f)
    
    Bbetabeta=Bbetabeta-t(listreg[[i]])%*%listreg[[i]]+Bbb1
    Bbs=Bbs-(1/2)*(Bbs1a+2*Bbs2a)
    
  #  prova1=outerList(Bss2a,Bss2b,f1)
   # prova2=outerList(Bss4a,Bss4b,f1)
    Bss=Bss+(1/4)*(prova1+prova2+prova3+prova4)
    }

  # first bit of A.44
  Bss1=lapply(1:nsigma,function(x,listJ,Sigma_valueg) Sigma_valueg%*%listJ[[x]],Sigma_valueg=Sigma_valueg,listJ=listJ)
  #collect results for Bsigmasigma
  Bss=Bss-(nnuts/2)*outerList(Bss1,Bss1,f1)
  
  B1=cbind(Bbetabeta,t(Bbs))
  B2=cbind(Bbs,Bss)
  B=rbind(B1,B2)
  
  variance=ginv(-B)
  varbhat=diag(variance)[1:nreg]
  }