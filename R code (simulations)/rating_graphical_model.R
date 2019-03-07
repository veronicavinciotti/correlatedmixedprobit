rating_graphical_model <- function(listy, listreg, listsector, listdummy, sden, listZ, lambda_value, bhat, Sigma_valueg, Theta_value, nnuts, nsectors, std_err, nreg,m,refit.MLE)
{
  #initial values
  Omega_valueg<-solve(Sigma_valueg)
  nreg<-length(bhat)
  
  # calculate Sigma_value: it is a list of sigma_Values, one big matrix for each region, 
  Sigma_value=lapply(1:nnuts,function(x,Sigma_valueg,listdummy) listdummy[[x]]%*%Sigma_valueg%*%t(listdummy[[x]])+diag(nrow(listdummy[[x]])),Sigma_valueg=Sigma_valueg, listdummy=listdummy)

  ##  now start the EM iteration:

  num_iter1=array()
  diff_value1=array()
  diff_valuebeta1=array()

  r=1
  s=1
  num_iter <- 0
  max_iter <- 700
  tol_value <-  0.01*(nsectors^2)
  tol_valueb  <- 0.05*nreg
  tol_valueu  <-0.01*nsectors
  diff_value <- 1e+10
  diff_u<- 1e+10
  num_iterbeta <- 0
  diff_valuebeta <- 1e+10
  max_iterbeta<-700
  loglik.calc<-TRUE #this is only needed for eBIC

	 Euy=lapply(1:nnuts,function(x, listZ) 	t(listdummy[[x]])%*%(listZ[[x]]-array(listreg[[x]]%*%bhat))/m, listZ=listZ) #1st approximation

	 while ((num_iterbeta < max_iterbeta) & (diff_valuebeta > tol_valueb) & (diff_u > tol_valueu) ) {
    while ((num_iter < max_iter) & (diff_value > tol_value)) {
		  
		  ## E-step: impute Z
		  obj_Estep <- impute_S(bhat, listy, listreg, listZ, listsector, Theta_value, Sigma_value, Sigma_valueg, listdummy, nnuts,nreg,m,loglik.calc)
		  
		  listZ  <- obj_Estep$Z
		  listZy <- obj_Estep$Zy
		  ES <- obj_Estep$ES
		  
		  ES <- obj_Estep$ES-diag(nsectors)/(m^2)
		  
		  if(lambda_value!=0){
		  obj <- glasso(s=ES, rho=lambda_value, maxit=1000, penalize.diagonal=F)
		  Sigma_valueg <- (t(obj$w) + obj$w) / 2
      Omega_newg <- (t(obj$wi) + obj$wi) / 2
		  }
		  else
		  {
		    Sigma_valueg<-ES
		    Omega_newg<-solve(ES)
		  }
		  diff_value <- sum(abs(Omega_newg - Omega_valueg))
		
		  num_iter <- num_iter + 1
		  num_iter1[r]=num_iter 
		  diff_value1[r]=diff_value
		  Omega_valueg <- Omega_newg
		  r=r+1
		  
		  Sigma_value=lapply(1:nnuts,function(x,Sigma_valueg,listdummy) listdummy[[x]]%*%Sigma_valueg%*%t(listdummy[[x]])+diag(nrow(listdummy[[x]])),Sigma_valueg=Sigma_valueg, listdummy=listdummy)
		}
	   
	   num_iter <- 0
	   diff_value <- 1e+10
	   
	# compute beta
		  Euy_new=lapply(1:nnuts,function(x, listdummy,listZy) 	t(listdummy[[x]])%*%listZy[[x]]/m, listdummy=listdummy,listZy=listZy)
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

	 if(refit.MLE == TRUE){
	   #refit non-zero parameters of precision matrix by MLE (needed for eBIC)
	   adj<-1*(Omega_newg!=0)
	   colnames(adj)<-rownames(adj)<-colnames(ES)<-rownames(ES)<-as.character(1:ncol(ES))
	   net.ggm <- fitConGraph(adj, ES, nnuts)
	   omega.mle <- solve(net.ggm$Shat)
	   omega.mle[abs(omega.mle)<0.0001]=0
	   sigma.mle <- net.ggm$Shat
	  
	   Sigma_value_mle=lapply(1:nnuts,function(x,sigma.mle,listdummy) listdummy[[x]]%*%sigma.mle%*%t(listdummy[[x]])+diag(nrow(listdummy[[x]])),sigma.mle=sigma.mle, listdummy=listdummy)
	   
	   obj_Estep.mle <- impute_S(bhat, listy, listreg, listZ, listsector, Theta_value, Sigma_value_mle, sigma.mle, listdummy, nnuts,nreg,m,loglik.calc)
	   
	   listZ.mle  <- obj_Estep.mle$Z
	   listZy.mle=obj_Estep.mle$Zy
	   Euy_new_mle=lapply(1:nnuts,function(x, listZy.mle) 	t(listdummy[[x]])%*%listZy.mle[[x]]/m, listZy.mle=listZy.mle)
	   num=lapply(1:nnuts,function(x,listreg, listZ.mle,listdummy,Euy_new_mle) t(listreg[[x]])%*%(as.matrix(listZ.mle[[x]])-listdummy[[x]]%*%Euy_new_mle[[x]]), listreg=listreg, listZ.mle=listZ.mle,listdummy=listdummy, Euy_new_mle=Euy_new_mle)
	   bhat.mle=sden%*%Reduce('+', num)
	   
	   obj_Estep.mle <- impute_S(bhat.mle, listy, listreg, listZ.mle, listsector, Theta_value, Sigma_value_mle, sigma.mle, listdummy, nnuts,nreg,m,loglik.calc)
	   }
	 
	#Calculate log-lik (Q-function)
	loglik=obj_Estep$Qfunc 
  
	cov2Zy=obj_Estep$cov2Zy	
  listmu=obj_Estep$listmu
  listsigmatr=obj_Estep$listsigmatr
 
  if (std_err==1) {
    varbhat=compute_std_err(bhat, listy, listreg, listZ, listZy, listsector, Theta_value, Sigma_value, Sigma_valueg, Omega_newg, listdummy, nnuts,nreg, nsectors,ES, cov2Zy, listmu,listsigmatr,m)
    }
  if (std_err==0) {
    varbhat=rep(1,nreg)}

	output <- list()
	output$Sigma_value <- Sigma_value
	output$ES <- ES
	output$num_iter <- num_iter1[r-1]
	output$loglik <- loglik #Q function
	output$betahat=bhat
	output$varbhat=varbhat
	output$Omega_newg<-Omega_newg 
	output$Sigma_newg<-Sigma_valueg 
	output$Eu <- Euy
	output$Z <- listZ
	output$avsigma_g=mean(diag(Sigma_valueg))

	return(output)
}

