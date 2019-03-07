############################################################################################
#
# Author:   Veronica Vinciotti
# e-mail:   veronica.vinciotti@brunel.ac.uk
# date:     07-03-2019

####This is an example to generate data from the mixed graphical probit and to make inference
#### The code is made more efficient by the fact that the same number of observations are generated for each sector/region

rm(list=ls())
library(glasso)
library(huge)
library(lme4)
library(OpenMx)
library(ggm)
library(rmngb)
source('rating_graphical_model.R')
source('truncate_component.R')
source('compute_std_err.R')
source('impute_S.R')

set.seed(1234)
nsectors<-10
nnuts<-200
p<-300
#200 regions, 300 companies in each region, 10 industrial sectors

m=p/nsectors
#number of companies per sector per region (simpler case where this is the same)

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

#INITIALISE ES (different options here)
y01bar=t(matrix(y01,p,nnuts))%*%listdummy[[1]]/m
ES=cor(y01bar)

sden=solve(t(matrixreg)%*%matrixreg)
  
#CREATE SEQUENCE OF LAMBDAS
nlambda=10
lambda.min.ratio = 0.001
lambda.max = max(max(ES- diag(diag(ES)), -min(ES- diag(diag(ES))))) +0.5
lambda.min = lambda.min.ratio * lambda.max
lambda = exp(seq(log(lambda.max), log(lambda.min), length = nlambda))

#INITALISE Sigma_valueg
obj <- glasso::glasso(s=ES, rho=lambda.max, maxit=1000, penalize.diagonal=F)
Sigma_valueg <- (t(obj$w) + obj$w) / 2

#INFER MODELS ACROSS PATH OF SOLUTIONS 
std_err<-0 # standard errors can be calculated just for the optimal model  
refit.MLE<-TRUE # only needed for calculations of eBIC
output.mod<-vector(length=nlambda,"list")
for  (s in 1:nlambda ){
  print(s)
lambda_value=lambda[s]
output.mod[[s]]=rating_graphical_model(listy, listreg, listsector, listdummy, sden, listZ, lambda_value, bhat, Sigma_valueg, Theta_value, nnuts, nsectors, std_err, nreg,m,refit.MLE)
bhat<-output.mod[[s]]$betahat
listZ<-output.mod[[s]]$Z
Sigma_valueg<-output.mod[[s]]$Sigma_newg
}
 
df<-NULL
for(s in 1:nlambda)
  df[s]<-sum(output.mod[[s]]$Omega_newg!=0)

idx<-which.min((df-sum(truetetag!=0))^2)
#For each model:
#1. estimated regression coefficients
output.mod[[idx]]$betahat
#2. estimated network
round(output.mod[[idx]]$Omega_newg,3)

#calculate standard errors of optimal model
lambda_value=lambda[idx]
std_err<-1
output.opt=rating_graphical_model(listy, listreg, listsector, listdummy, sden, listZ, lambda_value, bhat, Sigma_valueg, Theta_value, nnuts, nsectors, std_err, nreg,m,refit.MLE)
sqrt(output.opt$varbhat)
