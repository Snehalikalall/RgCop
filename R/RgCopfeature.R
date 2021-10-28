RgCopfeature <- function(data,classd,p,nf){
library(foreach)
library(doParallel)
library('copula')
 
 #data=a matrix format, Cells sholud be in row, Genes should be in coloumn in data
 #nf=Number of feature to be selected.
#p=Number of cores. 
 
cl <- makeCluster(p)
registerDoParallel(cl)
nf=nf
set.seed(1000)
datas2<-as.matrix(data)
classd<-as.matrix(classd)
n=nrow(datas2)
col=ncol(datas2)
count<-1:ncol(datas2)
start.time <- Sys.time()
# Feature by RgCop parallel 
fea=matrix(0, nrow=1,ncol =nf)
varlist=matrix(0,col,1)  # total variance list
for(i in 1:col)
{
  varlist[i,]<-var(datas2[,i])
}

# To find relevancy with class
mimc1<-foreach(j=count, .combine=c,.packages='copula') %dopar%  #Find the most relevant feature
  {u=pobs(cbind(classd,datas2[,j]))
  res<-mean(C.n(u, cbind(classd,datas2[,j])))
  }
mimc1<-as.matrix(mimc1)


# To find redundancy between nonselected feature and 
#other selected features. Redundancy
# Different tuning parameter gamma
coeff=0.3
# bivariate copula reg
fea[1,1]<-which.max(mimc1)
for (m in 2:nf)
{
  feas<-fea[1,1:(m-1)]
  parl<-foreach(j=count, .combine=c,.packages='copula') %dopar%
    { u1=pobs(cbind(datas2[,feas],datas2[,j]))
    res<-mean(C.n(u1, cbind(datas2[,feas],datas2[,j])))
    }
  red<-as.matrix(parl) 
  result<-(mimc1-red) + coeff*(abs(mimc1*varlist))
  result[fea]=Inf
  fea[1,m]<-which.min(result)
}

stopCluster(cl)
registerDoSEQ()

rgcopdata=data[,fea] # Feature reduced Data with Copula based feature selection
RgCopresult<- list("Feadata" = rgcopdata, "Features" = fea)
return(RgCopresult) 
}
