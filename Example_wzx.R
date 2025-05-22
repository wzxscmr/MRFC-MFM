library(openxlsx)
library(MCMCpack)
library(purrr)
library(snowfall)
library(parallel)

#adjacency matrix
county<-read.csv("countyadj.csv")
index <- grep("TX", county$X)
rownames(county)<-county[,1]
county<-county[,-1]
Texas<-county[index,index]
SX<-sort(colnames(Texas),index.return=T)
Texas<-as.matrix(Texas[SX$ix,SX$ix])
distance<-Texas

#data generation for Texas counties
lambda1 <- 0
a=1
b=1
neighbour=1
n=254
beta <- c(0.5, 0.5, 0, 0)
p=length(beta)
X_grid <- rnorm(254*p, mean = 0, sd = 1)
X_grid <- matrix(X_grid, nrow = 254, ncol = p)
MoransI.Basis<-function(X,r,A){
  n = dim(X)[1]
  PsiOper = (diag(n) - X%*%solve(t(X)%*%X)%*%t(X))%*%A%*%(diag(n) - X%*%solve(t(X)%*%X)%*%t(X))
  output2<-eigen(PsiOper)
  Psi = output2$vectors[,1:r]
  return(Psi)
}
B_grid<-as.matrix(MoransI.Basis(X_grid,5,distance))
N_grid=rep(1,254)
GAMMA=1
D<-matrix(c(1,  0.1,  0.1,  0.2,  0.2,
            0.1, 1,   0.1,  0.2,  0.2,
            0.1, 0.1, 1,    0.2,  0.2,
            0.2, 0.2, 0.2,  1,    0.2,
            0.2, 0.2, 0.2, 0.2,1), nrow=5, byrow=TRUE)
phi=diag(5)
v=6
ita<-mvrnorm(1,rep(0, 5),D)
initNCluster=3
niterations=5000
Y_grid=rep(NA,254)

#cluster configuration
d1<-read.xlsx("123.xlsx",colNames = F)
d1$cluster<-rep(3,254)
for (i in 1:254){
  if (d1$X4[i] %in% c(7,8,11)){
    d1$cluster[i]<-1
  }
  else if (d1$X4[i] %in% c(1,2,3)){
    d1$cluster[i]<-2
  }
}
d1$cluster<-factor(d1$cluster)
d1<-data.frame(d1)
colnames(d1)<-c("number","county","fips","X4","X5","cluster")

lambdac<-rep(0,254)
for (i in 1:254){
  if(d1$cluster[i]==1){
    lambdac[i]=20
  }
  else if(d1$cluster[i]==2){
    lambdac[i]=5
  }
  else if(d1$cluster[i]==3){
    lambdac[i]=0.2
  }
}
for (i in 1:254){
  Y_grid[i]<-rpois(1,N_grid[i]*lambdac[i]*exp(X_grid[i,] %*% beta + B_grid[i,] %*% ita))
}

#parallel computing
sfInit (parallel=TRUE , cpus=16)
sfLibrary(MCMCpack)
sfLibrary(purrr)
sfExport("X_grid", "Y_grid", "N_grid", "B_grid", "n", "a", "b", "phi", "v", "GAMMA","initNCluster", "niterations","lambda1","neighbour","distance","wzx_Gibbs")
resultsbn <- sfLapply(1:100, function(x) wzx_Gibbs(X_grid, Y_grid, N_grid, B_grid, n, a, b, phi, v, GAMMA,initNCluster, niterations,lambda1,  neighbour,distance))

sfStop()
save(resultsbn, file = "resultsf11.rdata")