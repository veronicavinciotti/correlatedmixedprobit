impute_S_application <- function(bhat, listy, listreg, listZ, listsector, Theta_value, Sigma_value, Sigma_valueg, listdummy, nnuts,nreg, Omega_valueg,loglik.calc,var.diag)
{            
  
  Z_new <- list()
  Zy_new <- list()
  listmu=list()
  cov2Zy=list()
  listsigmatr=list()
  Q.func=0
  
  ES=0
  Sigma_valueginv=Omega_valueg
  nsectors<-ncol(Sigma_valueginv)

  for (i in seq(1, nnuts)){  

    Zarray=array(listZ[[i]])
    prova=Zarray-array(listreg[[i]]%*%bhat)
    Sigma_v=Sigma_value[[i]]
    Gsec<- sort(listsector[[i]]) #vector of sectors in region i (sorted from smallest to largest..)
    idx<-cumsum(as.numeric(table(Gsec))) #one way of selecting one observation per sector
    veronica0=lapply(1:length(idx),function(x,i, idx,listdummy,Sigma_valueginv)  solve(Sigma_valueginv+t(listdummy[[i]][-idx[x],]) %*%listdummy[[i]][-idx[x],]), listdummy=listdummy,idx=idx, i=i,Sigma_valueginv=Sigma_valueginv)
    veronica1=lapply(1:length(idx),function(x,i, idx, Sigma_v,listdummy,veronica0)  array(Sigma_v[idx[x],-idx[x]]) - array(array(Sigma_v[idx[x],-idx[x]])%*%listdummy[[i]][-idx[x],]%*%veronica0[[x]]%*%t(listdummy[[i]][-idx[x],])), Sigma_v=Sigma_v, listdummy=listdummy,idx=idx, i=i,veronica0=veronica0)
    mutmp=lapply(1:length(Zarray),function(x,prova, veronica1,idx) veronica1[[which(idx>=x)[1]]] %*%array(prova[-x]),prova=prova,veronica1=veronica1,idx=idx)
    veronica2=lapply(1:length(idx),function(x,idx, Sigma_v,veronica1)  veronica1[[x]]%*%array(Sigma_v[-idx[x],idx[x]]) , Sigma_v=Sigma_v, idx=idx, veronica1=veronica1)
    
    veronica3<-vector(length=nsectors) 
    veronica3[unique(Gsec)]<-array(unlist(veronica2)) 
    veronica3[!unique(Gsec)]<-0
    sigma1=1+diag(Sigma_valueg)-veronica3
    
    sigma=unlist(lapply(1:length(Zarray),function(x,sigma1,listdummy) array(listdummy[[i]][x,])%*%sqrt(diag(sigma1))%*%array(listdummy[[i]][x,]) , sigma1=sigma1, listdummy=listdummy))
    
    listsigmatr[[i]]=sigma
    
    mu=array(listreg[[i]]%*%bhat)+unlist(mutmp)
    listmu[[i]]=mu
    
    obj <- truncate_component(Y=listy[[i]] , theta=Theta_value, mu=mu, sigma=sigma)      
    
    Z_new[[i]] <- obj$EX 
    Zy_new[[i]]= obj$EX - array(listreg[[i]]%*%bhat)
    
    diag_element <- obj$EXX + array(listreg[[i]]%*%bhat)^2-2 * Z_new[[i]]* array(listreg[[i]]%*%bhat)
    
    ld=listdummy[[i]]
 
    E2ystar=Zy_new[[i]]%*%t(Zy_new[[i]])
    diag(E2ystar)=diag_element
    cov2Zy[[i]]=E2ystar
    
    norm=colSums(listdummy[[i]])
    norm[norm==0]=1
    
    Zbar=diag(1/norm)%*%t(ld)%*%Zy_new[[i]]
    inner1=diag_element-Zy_new[[i]]^2
    inner2=lapply(1:nrow(ld),function(x,i,ld,inner1) inner1[x]*ld[x,], ld=ld, inner1=inner1 )
    inner3=Reduce('rbind', inner2)
    inner4=(diag(1/norm))%*%(t(inner3)%*%ld/nnuts)%*%(diag(1/norm))
    ES=ES+Zbar%*%t(Zbar)/nnuts+inner4
    
     if(loglik.calc){
      Q.func=Q.func-0.5*tr((diag(length(Zy_new[[i]]))-ld%*%(diag(1/norm))%*%t(ld))%*%E2ystar%*%t(diag(length(Zy_new[[i]]))-ld%*%(diag(1/norm))%*%t(ld)))
    }
           }
  
  #standardize diagonal
  if(sum(var.diag)!=ncol(ES)){
  diag(ES)[abs(diag(ES)) < 1e-10] <- 1e-10
  ES<-diag(1/sqrt(diag(ES)))%*%ES%*%diag(1/sqrt(diag(ES)))
  ES<-diag(sqrt(var.diag))%*%ES%*%diag(sqrt(var.diag))}
  
  
  if(loglik.calc){
    cov.mat<-ES
    Q.func=Q.func+nnuts/2*(determinant(Sigma_valueginv)$modulus-sum(Sigma_valueginv*cov.mat))
  }
  
  output <- list()
  output$ES <- ES
  output$Z <- Z_new
  output$Zy <-Zy_new
  output$Qfunc=Q.func
  output$listmu=listmu
  output$cov2Zy=cov2Zy
  output$listsigmatr=listsigmatr
  return(output)         
}

