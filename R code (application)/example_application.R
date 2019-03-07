############################################################################################
#
# Author:   Veronica Vinciotti
# e-mail:   veronica.vinciotti@brunel.ac.uk
# date:     07-03-2019

####This is an example to generate data from the mixed graphical probit and to make inference
####This code is more general and does not assume that one has the same number of observations for each sector/region


rm(list=ls())
library(glasso)
library(huge)
library(lme4)
library(OpenMx)
library(ggm)
library(rmngb)
source('rating_graphical_model_application.R')
source('truncate_component.R')
source('compute_std_err_application.R')
source('impute_S_application.R')

set.seed(1234)
nsectors<-10
nnuts<-200
p<-300
#200 regions, 300 companies in each region, 10 industrial sectors

m=p/nsectors
#number of companies per sector per region (simpler case where this is the same, but here we will use the general code)


#GENERATE MIXED GRAPHICAL PROBIT WITH CORRELATED RANDOM EFFECTS
  
b0=c(0,1)
modelg=huge.generator(nnuts, nsectors, graph = "random" ,   verbose = FALSE, prob=3/nsectors)
truesigmag=modelg$sigma
truetetag=solve(truesigmag)
truetetag[abs(truetetag)<0.0001]=0
  
truesigma=kronecker(truesigmag,matrix(1,m,m))+kronecker((diag(nsectors)),diag(m))
trueteta=solve(truesigma, sparse=TRUE)
trueteta[abs(trueteta)<0.0001]=0
  
wg=modelg$theta
w=kronecker((modelg$theta+diag(nsectors)),matrix(1,m,m))-diag(p)
  
xmodel=huge.generator(nnuts, p, graph = "random",   verbose = FALSE)
xtruesigma=xmodel$sigma
  
xd=chol(xtruesigma);
d=chol(truesigmag);

u=t(t(d)%*% matrix( rnorm(nsectors*nnuts, 0, 1),nsectors,nnuts));
error=matrix( rnorm(p*nnuts, 0, 1),nnuts,p);
xerror=t(t(xd)%*% matrix( rnorm(p*nnuts, 0, 1),p,nnuts));
  
listreg=list()
listystar=list()
listy=list()
listy01=list()
listsector=list()
listdummy=list()
sector1=kronecker((1:nsectors),array(1,m))
  
for (j in 1:nnuts) {
    listreg[[j]]= cbind(rep(1,p),as.array(xerror[j,]))
    listsector[[j]]<-sector1
    listdummy[[j]]=matrix(0,p,nsectors)
    for (h in 1:nsectors){
      listdummy[[j]][sector1==h,h]=1}
    listystar[[j]]=array(listreg[[j]]%*%b0)+array(listdummy[[j]]%*%u[j,])+as.array(error[j,])  
    listy[[j]]=as.numeric(listystar[[j]]>0)
    listy01[[j]]<-listy[[j]]
    listy[[j]][listy[[j]]==1]=2
    listy[[j]][listy[[j]]==0]=1
  }
matrixreg=Reduce('rbind', listreg)
  
y=unlist(listy)
ystar=unlist(listystar)
y01=y
  
y01[y==1]=0
y01[y==2]=1
  
Theta_value <- rep(0,3)
Theta_value[1]<- -Inf
Theta_value[3]<- Inf
  
listZ <- list()
for (j in seq(1, nnuts)){
    tmp2 <- truncate_component(Y=listy[[j]], theta=Theta_value)
    listZ[[j]] <- tmp2$EX
  }
  
nuts=kronecker((1:nnuts),array(1,p))
nreg=ncol(listreg[[1]])

#INITALISE bhat
ID=as.factor(c(kronecker(c(1:(nsectors*nnuts)),rep(1,m))))
probit=glmer(y01 ~ c(matrixreg[,2]) +  (1| ID), family=binomial(link="probit"))
bhat=fixef(probit)


norm=list()
for (i in 1:nnuts){
  norm[[i]]=colSums(listdummy[[i]])
  norm[[i]][norm[[i]]==0]=1}

ybar=matrix(0, nsectors, nnuts)

for (i in (1:nnuts)) {
  ybar[,i]=as.array(((t(listdummy[[i]])%*%listy01[[i]]/norm[[i]])))-as.array(((t(listdummy[[i]])%*%listreg[[i]]%*%bhat/norm[[i]])))
}


# initialise ES using pairwise complete covariance
ES= cov(t(ybar),use="pairwise.complete.obs")

sden=solve(t(matrixreg)%*%matrixreg)

#MLE
lambda_value=0
stand<-0
refit.MLE<-TRUE
obj <- glasso::glasso(s=ES, rho=lambda_value, maxit=1000, penalize.diagonal=F)
Sigma_valueg <- (t(obj$w) + obj$w) / 2 # this is simply ES
diag.var<-rep(1,ncol(Sigma_valueg)) #no standardization here
output.mle<-rating_graphical_model_application(listy, listreg, listsector, listdummy, listZ, lambda_value, bhat, Sigma_valueg,Theta_value, nnuts, nsectors,refit.MLE,diag.var,stand)

var.diag<-diag(output.mle$Sigma_newg) # extract variances and keep fixed throughout the path

#INITIALISE ES ACCORDING TO THE MLE ESTIMATE
ES= cor(t(ybar),use="pairwise.complete.obs")
ES<-diag(sqrt(var.diag))%*%ES%*%diag(sqrt(var.diag))

#CREATE SEQUENCE OF LAMBDAS
nlambda=10
lambda.min.ratio = 0.001
lambda.max = max(max(ES- diag(diag(ES)), -min(ES- diag(diag(ES))))) 
lambda.min = lambda.min.ratio * lambda.max
lambda = exp(seq(log(lambda.max), log(lambda.min), length = nlambda))

#INITALISE Sigma_valueg
obj <- glasso::glasso(s=ES, rho=lambda.max, maxit=1000, penalize.diagonal=F)
Sigma_valueg <- (t(obj$w) + obj$w) / 2

#INFER MODELS ACROSS PATH OF SOLUTIONS 
stand<-0
output.mod<-vector(length= nlambda+1,mode="list")
for  (s in 1: nlambda){
  print(s)
  lambda_value=lambda[s]
  output.mod[[s]]<-rating_graphical_model_application(listy, listreg, listsector, listdummy, listZ, lambda_value, bhat, Sigma_valueg,Theta_value, nnuts, nsectors,refit.MLE,var.diag,stand)
  bhat<-output.mod[[s]]$betahat
  listZ<-output.mod[[s]]$Z
  Sigma_valueg<-output.mod[[s]]$Sigma_newg
}

output.mod[[nlambda+1]]<-output.mle

df<-NULL
for(s in 1: (nlambda+1))
  df[s]<-sum(output.mod[[s]]$Omega_newg!=0)

idx<-which.min((df-sum(truetetag!=0))^2)
#For each model:
#1. estimated regression coefficients
output.mod[[idx]]$betahat
#2. estimated network
round(output.mod[[idx]]$Omega_newg,3)

#calculate standard errors of optimal model (time consuming)
lambda_value=lambda[idx]
stand<-1
output.opt=rating_graphical_model_application(listy, listreg, listsector, listdummy, listZ, lambda_value, bhat, Sigma_valueg,Theta_value, nnuts, nsectors,refit.MLE,var.diag,stand)
sqrt(output.opt$varbhat)

