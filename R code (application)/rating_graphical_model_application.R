rating_graphical_model_application<- function(listy, listreg, listsector, listdummy, listZ, lambda_value, bhat, Sigma_valueg,Theta_value, nnuts, nsectors,refit.MLE,var.diag,stand)
{
  Omega_valueg<-solve(Sigma_valueg)
  Sigma_value=lapply(1:nnuts ,function(x,Sigma_valueg,listdummy) listdummy[[x]]%*%Sigma_valueg%*%t(listdummy[[x]])+diag(nrow(listdummy[[x]])),Sigma_valueg=Sigma_valueg, listdummy=listdummy)
  
  num_iter1=array()
  diff_value1=array()
  diff_valuebeta1=array()
  
  nreg<-length(bhat)
  
  r=1
  s=1
  num_iter <- 0
  max_iter <- 700
  tol_value <-  0.015*(nsectors^2)
  tol_valueb  <- 0.01*nreg
  tol_valueu  <-0.1*nsectors
  diff_value <- 1e+10
  diff_u<- 1e+10
  num_iterbeta <- 0
  diff_valuebeta <- 1e+10
  max_iterbeta<-700
  estimates_sigma=list()
  estimates_omega=list()
  loglik.calc=TRUE
  
  norm=list()
  for (i in 1:nnuts){
    norm[[i]]=colSums(listdummy[[i]])
    norm[[i]][norm[[i]]==0]=1}
  
  Euy=lapply(1:nnuts,function(x, listZ,norm) 	t(listdummy[[x]])%*%(listZ[[x]]-array(listreg[[x]]%*%bhat))/norm[[x]], listZ=listZ, norm=norm)
  
  
	 # start  EM iterations:
  sden=ginv(t(matrixreg)%*%matrixreg)
  
	 while((num_iterbeta < max_iterbeta) & (diff_valuebeta > tol_valueb)  & (diff_u > tol_valueu)) {
		while((num_iter < max_iter) & (diff_value > tol_value)) {

		  obj_Estep <- impute_S_application(bhat, listy, listreg, listZ, listsector, Theta_value, Sigma_value, Sigma_valueg, listdummy, nnuts,nreg, Omega_valueg,loglik.calc,var.diag)
		  listZ  <- obj_Estep$Z
		  listZy=obj_Estep$Zy
		  
		  ES <- obj_Estep$ES
		  obj <- glasso(s=ES, rho=lambda_value, maxit=1000, penalize.diagonal=F)
		  Omega_newg <- (t(obj$wi) + obj$wi) / 2
		  Sigma_valueg <- (t(obj$w) + obj$w) / 2

		  diff_value <- sum(abs(Omega_newg - Omega_valueg))
		  
		  num_iter <- num_iter + 1
		  num_iter1[r]=num_iter 
		  diff_value1[r]=diff_value
		  Omega_valueg <- Omega_newg
		  r=r+1
		  Sigma_value=list()
		  for (x in 1:nnuts){ 
		    Sigma_value[[x]]=listdummy[[x]]%*%Sigma_valueg%*%t(listdummy[[x]])+diag(nrow(listdummy[[x]]))}
		}	   
	   
	   num_iter <- 0
	   diff_value <- 1e+10

	   Euy_new=lapply(1:nnuts,function(x, listZy,norm,listdummy) (t(listdummy[[x]])%*%listZy[[x]])/norm[[x]], listZy=listZy,norm=norm,listdummy=listdummy)
	   num=lapply(1:nnuts,function(x,listreg, listZ,listdummy,Euy_new) t(listreg[[x]])%*%(as.matrix(listZ[[x]])-listdummy[[x]]%*%Euy_new[[x]]), listreg=listreg, listZ=listZ,listdummy=listdummy, Euy_new=Euy_new)
	   bhat_new=sden%*%Reduce('+', num)
	   
	   diffu=lapply(1:nnuts,function(x, Euy_new, Euy) sum(abs(Euy_new[[x]] - Euy[[x]])), Euy_new=Euy_new, Euy=Euy)
	   diff_u= sum(unlist(diffu))
	   diff_valuebeta <- sum(abs(bhat_new - as.matrix(bhat)))
	   diff_valuebeta1[s+1]=diff_valuebeta 
	   s=s+1
	   bhat=bhat_new	
	   Euy=Euy_new
	   num_iterbeta <- num_iterbeta + 1
	 }
  
  Omega_newg[abs(Omega_newg)<0.0001]=0
  cov2Zy=obj_Estep$cov2Zy	
  listmu=obj_Estep$listmu
  listsigmatr=obj_Estep$listsigmatr
  
   
  if(refit.MLE == TRUE){
    #refit non-zero parameters of precision matrix by MLE. Only needed for eBIC
    adj<-1*(Omega_newg!=0)
    colnames(adj)<-rownames(adj)<-colnames(ES)<-rownames(ES)<-as.character(1:ncol(ES))
    net.ggm <- fitConGraph(adj, ES, nnuts)
    omega.mle <- solve(net.ggm$Shat)
    omega.mle[abs(omega.mle)<0.0001]=0
    sigma.mle <- net.ggm$Shat
    
    Sigma_value_mle=lapply(1:nnuts,function(x,sigma.mle,listdummy) listdummy[[x]]%*%sigma.mle%*%t(listdummy[[x]])+diag(nrow(listdummy[[x]])),sigma.mle=sigma.mle, listdummy=listdummy)
    
    obj_Estep.mle <- impute_S_application(bhat, listy, listreg, listZ, listsector, Theta_value, Sigma_value_mle, sigma.mle, listdummy, nnuts,nreg,omega.mle,loglik.calc,var.diag)
    
    listZ.mle  <- obj_Estep.mle$Z
    listZy.mle=obj_Estep.mle$Zy
    Euy_new_mle=lapply(1:nnuts,function(x, listZy.mle,norm,listdummy) (t(listdummy[[x]])%*%listZy.mle[[x]])/norm[[x]], listZy.mle=listZy.mle,norm=norm,listdummy=listdummy)
    num=lapply(1:nnuts,function(x,listreg, listZ.mle,listdummy,Euy_new_mle) t(listreg[[x]])%*%(as.matrix(listZ.mle[[x]])-listdummy[[x]]%*%Euy_new_mle[[x]]), listreg=listreg, listZ.mle=listZ.mle,listdummy=listdummy, Euy_new_mle=Euy_new_mle)
    bhat.mle=sden%*%Reduce('+', num)
    
    obj_Estep.mle <- impute_S_application(bhat.mle, listy, listreg, listZ.mle, listsector, Theta_value, Sigma_value_mle, sigma.mle, listdummy, nnuts,nreg,omega.mle,loglik.calc,var.diag)
    obj_Estep<-obj_Estep.mle
      }
  
  #Calculate log-lik (Q-function)
  loglik=obj_Estep$Qfunc 
  

  if (stand==1) {
    varbhat=compute_std_err_application(bhat, listy, listreg, listZ, listZy, listsector, Theta_value, Sigma_value, Sigma_valueg, listdummy, nnuts,nreg, nsectors,ES, cov2Zy, listmu,listsigmatr)
  }
  
  if (stand==0) {
    varbhat=rep(1,nreg)
  }
  
	output <- list()
	output$Sigma_value <- Sigma_value
	output$ES <- ES
	output$betahat=bhat_new
	output$varbhat<-varbhat
	output$loglik <- loglik #Q function
	output$Omega_newg<-Omega_newg 
	output$Sigma_newg<-Sigma_valueg 
	output$Eu <- Euy
	output$Z <- listZ
	output$Zy <- listZy
	
		return(output)
}

