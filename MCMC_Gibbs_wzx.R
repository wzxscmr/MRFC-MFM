##' input: X_grid (covariates value for every county, n by p matrix);
##' input: B_grid (covariates value for every county, n by q matrix);
##' input: Y_grid, N_grid (n by 1 vector);
##' input: n, number of county;
##' input: hyperparameters, a, b, phi(q by q matrix), v;
##' input: GAMMA, power of the pCRP;
##' input: initNCluster, initial number of cluster;
##' input: niterations, posterior sample size;
##' input: lambda1, Markov random field parameter;
##' input: neighbour, usually set to 1;
##' input: distance, the adjacency matrix for data;
wzx_Gibbs <- function(X_grid, Y_grid, N_grid, B_grid, n, a, b, phi, v, GAMMA, 
                      initNCluster, niterations, lambda1, neighbour,distance){
  index_cluster <- c(sample(1:initNCluster, size = initNCluster, replace = FALSE),
                     sample(1:initNCluster, size = n-initNCluster, replace = TRUE))
  lambda0 <- rgamma(initNCluster, shape = a, rate = b)
  nCluster <- initNCluster
  beta <- rep(0, ncol(X_grid))
  gam <- rep(0, ncol(X_grid))
  ita <- rep(0, ncol(B_grid))
  
  #precomputation for prespecified coefficient VN
  lambda <- 1
  gamma <- GAMMA
  NN=n ## n is the number of oberservations
  VN<-0
  tmax = n+10
  for (t in 1:tmax)
  {
    r = log(0)
    for (k in t:700)
    {
      zz = sum(log((k-t+1):k))-sum(log((k*gamma):(k*gamma+NN-1))) + dpois(k-1, lambda, log = TRUE)
      m = max(zz,r)
      r = log(exp(r-m) + exp(zz-m)) + m
    }
    VN[t] = r
  }
  
  History <- vector("list", niterations)
  ## start Gibb's sampling
  for (iter in 1:niterations){
    Lambda <- N_grid * exp(X_grid %*% beta + B_grid %*% ita)
    
    ## update index_cluster
    for (i in 1:n) {
      count_cluster <- table(index_cluster)
      ## determine if i-th grid is a singleton
      if (count_cluster[index_cluster[i]] > 1) {
        ## if not a singleton, then have nCluster + 1 choice
        count_cluster[index_cluster[i]] <- count_cluster[index_cluster[i]] - 1
        ## goes to an existing cluster: z_i = 1:nCluster
        logclusterProbs = sapply(1:nCluster, function(x) {
          clusterAssign_temp = index_cluster
          clusterAssign_temp[i] = x
          return(lambda1*sum(as.numeric(distance[i,clusterAssign_temp==x]>=neighbour))+
                   Y_grid[i]*log(lambda0[x])+log(GAMMA+count_cluster[x])-
                   Lambda[i]*lambda0[x])
        })
        ## goes to a new cluster: z_i = nCluster+1
        logclusterProbs[nCluster+1]<-log(GAMMA)+VN[nCluster+1]-VN[nCluster]+
          a * log(b) + lgamma(Y_grid[i] + a) - (Y_grid[i] + a) * log(b+Lambda[i]) - 
          lgamma(a)
        logclusterProbs <- logclusterProbs - max(logclusterProbs)
        for (z in 1:(nCluster+1)){
          if(logclusterProbs[z]>=700){
            logclusterProbs[z]=700
          }
        }
        clusterProbs <- exp(logclusterProbs)
        ## get the posterior sample for Z_i
        index_i <- sample(1:(nCluster+1), size = 1,
                          prob = clusterProbs/sum(clusterProbs))
        ## if the i-th grid really goes to a new cluster
        if (index_i > nCluster) {
          lambda0_new <- rep(0, nCluster + 1)
          lambda0_new[1:nCluster] <- lambda0
          lambda0_new[nCluster+1] <- rgamma(1, shape = a, rate = b)
          lambda0 <- lambda0_new
          
          index_cluster[i] <- index_i
          nCluster <- nCluster + 1
        } else { ## if i-th grid goes to an existing cluster
          index_cluster[i] <- index_i
        }
      } 
      else { ## if grid is a singleton, then has nCluster choices
        ## move all the cluster index that are greater than index_cluster[i] 
        ## foward by one to fill in the blank, also change the count_cluster 
        ## and lambda0 coresspondingly
        index_cluster[index_cluster > index_cluster[i]] <- 
          index_cluster[index_cluster > index_cluster[i]] - 1
        lambda0[index_cluster[i] : (nCluster-1)] <- lambda0[(index_cluster[i]+1):nCluster]
        count_cluster <- table(index_cluster)
        count_cluster[index_cluster[i]] <- count_cluster[index_cluster[i]] - 1
        ## only nCluster-1 clusters remained after removing i-th grid
        lambda0 <- lambda0[1:(nCluster-1)]
        ## goes to existing cluster : 1:(nCluster-1)
        logclusterProbs <- sapply(1:(nCluster-1), function(x) {
          clusterAssign_temp = index_cluster
          clusterAssign_temp[i] = x
          return(lambda1*sum(as.numeric(distance[i,clusterAssign_temp==x]>=neighbour))+
                   Y_grid[i]*log(lambda0[x])+log(GAMMA+count_cluster[x])-Lambda[i]*lambda0[x])
        })
        ## goes to new cluster: nCluster
        logclusterProbs[nCluster] <- log(GAMMA)+VN[nCluster+1]-VN[nCluster]+
          a * log(b) + lgamma(Y_grid[i] + a) - (Y_grid[i] + a) * log(b+Lambda[i]) - 
          lgamma(a)
        logclusterProbs <- logclusterProbs - max(logclusterProbs)
        for (p in 1:(nCluster)){
          if(logclusterProbs[p]>=700){
            logclusterProbs[p]=700
          }
        }
        # get the posterior sample for Z_i
        clusterProbs <- exp(logclusterProbs)
        index_i <- sample(1:nCluster, size = 1, 
                          prob = clusterProbs/sum(clusterProbs))
        index_cluster[i] <- index_i
        if (index_i < nCluster) {
          nCluster <- nCluster-1
        } else {
          lambda0_new <- rep(0, nCluster)
          lambda0_new[1:(nCluster-1)] <- lambda0
          lambda0_new[nCluster] <- rgamma(1, shape = a, rate = b)
          lambda0 <- lambda0_new
        }
      }
    }
    
    ## update lambda
    for (i in 1:nCluster) {
      lambda0[i] <- rgamma(1, shape = a+sum(Y_grid[index_cluster == i]), 
                           rate = b + sum(Lambda[index_cluster == i]))
    } 
    
    ## update gam
    gam <- rbinom(length(gam), size = 1, prob = (1+dnorm(beta, sd = 0.01)/dnorm(beta, sd = 10))^(-1))
    
    ## update beta, using MH
    loglikelihood <- function(beta) {
      Lambda <- N_grid * exp(X_grid %*% beta + B_grid %*% ita)
      l_grid <- rep(0, n)
      for (i in 1:n) {
        l_grid[i] <- X_grid[i,] %*% beta * Y_grid[i] - Lambda[i]*lambda0[index_cluster[i]]
      }
      return(sum(l_grid))
    }
    logprior <- function(beta, gam) {
      return(sum(dnorm(beta, mean = 0, sd = sqrt(0.0001*(1-gam) + 100*gam), log = T)))
    }
    logposterior <- function(beta, gam) {
      return(loglikelihood(beta) + logprior(beta, gam))
    }
    proposalfunction <- function(beta, sd) {
      return(rnorm(length(beta), mean = beta, sd = sd))
    }
    step_beta <- rep(0.05, length(beta))
    ar_beta <- rep(0, length(beta))
    for (i in 1:length(beta)) {
      ## adaptive step
      beta_prop <- proposalfunction(beta[i], sd = step_beta[i])
      beta_proposal <- beta
      beta_proposal[i] <- beta_prop
      probab <- exp(logposterior(beta_proposal, gam) - logposterior(beta, gam))
      if (runif(1) < probab) {
        beta <- beta_proposal
        ar_beta[i] <- 1
      }
    }
    
    ## update D
    D <- riwish(v + 1, phi + ita %*% t(ita))
    
    ## update ita, using MH
    loglikelihood <- function(ita) {
      Lambda <- N_grid * exp(X_grid %*% beta + B_grid %*% ita)
      l2_grid <- rep(0, n)
      for (i in 1:n) {
        l2_grid[i] <- B_grid[i,] %*% ita * Y_grid[i] - Lambda[i]*lambda0[index_cluster[i]]
      }
      return(sum(l2_grid))
    }
    logprior <- function(ita, D) {
      return(-t(ita) %*% solve(D) %*% ita /2)
    }
    logposterior <- function(ita, D) {
      return(loglikelihood(ita) + logprior(ita, D))
    }
    proposalfunction <- function(ita, sd) {
      return(rnorm(length(ita), mean = ita, sd = sd))
    }
    step_ita <- rep(0.05, length(ita))
    ar_ita <- rep(0, length(ita))
    for (i in 1:length(ita)) {
      ## adaptive step
      ita_prop <- proposalfunction(ita[i], sd = step_ita[i])
      ita_proposal <- ita
      ita_proposal[i] <- ita_prop
      probab <- exp(logposterior(ita_proposal, D) - logposterior(ita, D))
      if (runif(1) < probab) {
        ita <- ita_proposal
        ar_ita[i] <- 1
      }
    }
    History[[iter]] <- list(gam = gam, beta = beta, lambda0 = lambda0, ita = ita,
                            index_cluster = index_cluster, D = D)
  }
  list(Iterates = History)
}