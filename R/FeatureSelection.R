library('MASS')
library(foreach)
library(doParallel)


# load the preprocess data
data=as.matrix(read.table("preprocessdata.csv",header=FALSE,sep=","))
#Cell type vector
cell<-as.matrix(read.table("celltype.csv",header=FALSE,sep=",")) 
n <- nrow(data)
col<-ncol(data)
count=1:ncol(data)            
#Number of feature to be selected, default
nf=100
# Number of cores, p= 20 default
p=20
# Different tuning parameter gamma
coeff=0.009


# Regularized Copula based feature selection
  datas2=data
  ## Total variance list of each gene
  varlist=matrix(0,col,1)  
  for(i in 1:col)
  {
    varlist[i,]<-var(datas2[,i])
  }
  
  
  cl <- makeCluster(p)
  registerDoParallel(cl)
  classc=cell
  # To find relevancy with class
  u <- matrix(runif(n*2), n, 2)
  mimc1<-foreach(j=count, .combine=c,.packages='copula') %dopar%  #Find the most relevant feature
    {res<-mean(C.n(u, cbind(classc,datas2[,j])))
    }
  mimc1<-as.matrix(mimc1)
  
  fea<- matrix(0, nrow=1,ncol =nf)
  fea[1,1]=which.max(mimc1)
  for (m in 2:nf)
  {u1<- matrix(runif(n*m), n,m)
  feas<-fea[1,1:(m-1)]
  parl<-foreach(j=count, .combine=c,.packages='copula') %dopar%
    {
      res<-mean(C.n(u1, cbind(datas2[,feas],datas2[,j])))
    }
  red<-as.matrix(parl) 
  result<-(mimc1-red) + coeff*(abs(mimc1*varlist))
  result[fea]=Inf
  fea[1,m]<-which.min(result)
  }
stopCluster(cl)
registerDoSEQ()
# Feature reduced data using RgCop
Selecteddata=data[,fea]
Selectedfea=fea


#Saving feature selected data for further analysis
write.table(Selecteddata,file="DatafeaSelected.csv",sep=",")

