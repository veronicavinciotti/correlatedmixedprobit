compute_std_err <- function(bhat, listy, listreg, listZ, listZy, listsector, Theta_value, Sigma_value, Sigma_valueg,Omega_newg, listdummy, nnuts,nreg, nsectors,ES, cov2Zy, listmu,listsigmatr,m)
   {

  f <- function(x, y){
    sum(x*y)  }
  
    f1 <- function(x, y){
    sum(x*t(y))}
  
  #nsigma=length(vech(Sigma_valueg))
  #listJ=list()
  #jj=0
  #for (g in  1:nsectors){
  #  for (h in  g:nsectors){
  #    jj=jj+1
  #    listJ[[jj]]=matrix(0,nsectors,nsectors)
  #    listJ[[jj]][g,h]=1
  #    listJ[[jj]][h,g]=1
  #  }}
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
  # in real data need to calculate A and ApltA for each i
  
  ld=listdummy[[1]]
  A=lapply(1:nsigma,function(x,listJ,ld,m)  ld%*%listJ[[x]]%*%t(ld)/m^2,ld=ld,listJ=listJ,m=m)
  ApltA=lapply(1:nsigma,function(x,A)  0.5*(A[[x]]+t(A[[x]])),A=A)

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
  
    # now start calculation of the B matrix:
    
    # this quantity is useful to compute A.37: (I-Z*Z'/m)
    ImZy=diag(length(Zy))-ld%*%t(ld)/m
    
    # Bbb1 is the second bit in A.40, useful to calculate Bbetabeta
    Bbb1=t(listreg[[i]])%*%ImZy%*%EXX%*%t(ImZy)%*%listreg[[i]]
    
    # elements of Bbetasigma: (A.43)
    Bbs1a=lapply(1:nsigma,function(x,ApltA,EXX05,Zy,lambda3,ld,listreg,ImZy) t(t(listreg[[i]])%*%ImZy%*%EXX05%*%diag(diag(lambda3*EXX05*ApltA[[x]]*EXX05))%*%array(1,length(Zy))),ApltA=ApltA,EXX05=EXX05,Zy=Zy,lambda3=lambda3,ld=ld,listreg=listreg,ImZy=ImZy)
    Bbs2a=lapply(1:nsigma,function(x,Zy,ApltA,listreg,ImZy) t(t(listreg[[i]])%*%ImZy%*%EXX%*%ApltA[[x]]%*%Zy),Zy=Zy,listreg=listreg,ApltA=ApltA,ImZy=ImZy)
    
    # elements of Bsigmasigma (A.44-A.48)
    #A.44
    Bss2a=lapply(1:nsigma,function(x,EXX,ApltA) EXX%*%ApltA[[x]], EXX=EXX,ApltA=ApltA)
    prova1=diag(unlist(lapply(1:nsigma,function(x,Bss2a) 2*sum(Bss2a[[x]]*t(Bss2a[[x]])), Bss2a=Bss2a)))
    #A.46-A.47
    Bss3a=lapply(1:nsigma,function(x,lambda3,EXX05,ApltA) array(diag(lambda3*EXX05*ApltA[[x]]*EXX05)),lambda3=lambda3,EXX05=EXX05,ApltA=ApltA)
    Bss3b=lapply(1:nsigma,function(x,Zy,ApltA,EXX05) array(EXX05%*%ApltA[[x]]%*%Zy), ApltA=ApltA,Zy=Zy,EXX05=EXX05)
    Bss3=outerList(Bss3a,Bss3b,f) 
    prova3=2*Bss3+2*t(Bss3)
    #A.44
    Bss4a=lapply(1:nsigma,function(x,EXX05,lambda4,ApltA) lambda4*EXX05*ApltA[[x]]*EXX05,ApltA=ApltA,EXX05=EXX05,lambda4=lambda4)
    Bss4b=lapply(1:nsigma,function(x,EXX05,ApltA) EXX05*ApltA[[x]]*EXX05,ApltA=ApltA,EXX05=EXX05)
    prova2=diag(unlist(lapply(1:nsigma,function(x,Bss4a,Bss4b) sum(Bss4a[[x]]*t(Bss4b[[x]])), Bss4a=Bss4a,Bss4b=Bss4b)))
    #A.48
    Bss5a=lapply(1:nsigma,function(x,listJ,Zy,ApltA,EXX05) t(Zy)%*%ApltA[[x]]%*%EXX05,Zy=Zy,ApltA=ApltA,EXX05=EXX05)
    prova4=4*outerList(Bss5a,Bss5a,f)
    
    Bbetabeta=Bbetabeta-t(listreg[[i]])%*%listreg[[i]]+Bbb1
    Bbs=Bbs-(1/2)*(Reduce('rbind', Bbs1a)+2*Reduce('rbind', Bbs2a))
    
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