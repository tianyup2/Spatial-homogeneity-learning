require(mnormt)
require(matlib)
require(igraph)


log.dmn <- function(y_std,phi_0,d_matrix){
  t <- dim(y_std)[1]
  n <- dim(y_std)[2]
  Lambda <- diag(n)-phi_0*d_matrix
  return(-t*n/2*log(2*pi)+t/2*log(det(Lambda))-
           0.5*(sum(diag(t(y_std)%*%y_std%*%Lambda))))
}

GWCRP_rgrs_car <- function(y, phi, d_matrix, mh_sz, h, niterations, 
                           mu_0, Lambda_0, nu_0, s_0, log_alpha, initNClusters)
{
  ## Model: y_i|phi,phi_0 ,z_i,ksi_j,sigma_j^2 \sim MVN(phi_i*ksi_{z_i},sigma_{z_i}^2) ##
  ##        covariance between different y_i's at the same time point is conditional autoregressive, independent o.w.
  ##        where the covariance structure after normalizing out sigma_j^2 is (I-phi_0*d_matrix)_{n x n} (x) I_{t x t}
  ##        (x) is the Kronecker product
  ##        z|alpha \sim GWCRP(alpha, h)
  ##        ksi_j|sigma_j^2 \sim MVN(mu_0, sigma_j^2*Lambda_0^{-1})
  ##        sigma_j^2 \sim Inv-Gamma(nu_0/2,nu_0*s_0^2/2)
  ##        phi_0 \sim Unif(1/min(eig_val(d_matrix)),1/max(eig_val(d_matrix)))
  
  ################################################################
  
  ## Input: 
  ##        y, n*t matrix ##
  ##        phi, t*K matrix ##
  ##        d_matrix, n*n graph distance matrix ##
  ##        mh_sz, stepsize of sampling phi by Metropolis-Hasting, using random walk##
  ##        h = decay coefficient for graph distance ##
  ##        niterations = the total number of iterations ##
  ##        mu_0 = K*1 mean vector of the ksi prior ##
  ##        Lambda_0 = K*K precision matrix of the ksi prior ##
  ##        nu_0 = prior shape for 1/sigma_j^2 ##
  ##        s_0 = expected prior variance for sigma_j^2 ##
  ##        alpha = the parameter in Dirichlet distribution that 
  ##                controls the relative size of clusters ##
  ##        initNClusters = the initial number of clusters ##
  
  ## Output: 
  ##         zout = clustering configuration, n by 1 vector ##
  ##         ksiout = estimate ksi, K by p matrix ##
  ##         sigmaout = estimate sigma, 1 by p matrix ##
  ##         p is the number of cluster ## 
  
  #################################################################
  
  d_graph <- graph_from_adjacency_matrix(d_matrix)
  d_matrix_sp <- shortest.paths(d_graph)
  
  n <- dim(y)[1]
  t <- dim(y)[2]
  W <- exp(-d_matrix_sp*h)
  diag(W) <- 0
  W[d_matrix_sp==1] <- 1
  W <- W * n * (n-1) / sum(W)
  
  eig <- eigen(d_matrix)
  l_b <- 1/min(eig$values)
  u_b <- 1/max(eig$values)
  
  invL0 <- solve(Lambda_0)
  phiTphi <- t(phi)%*%phi
  L0mu0 <- Lambda_0%*%mu_0
  mu0TL0mu0 <- t(mu_0)%*%L0mu0
  
  #===== first initialize clustering and beta
  clusterAssign <- c(sample(1:initNClusters, size=initNClusters, replace=FALSE),
                     sample(1:initNClusters, size=n-initNClusters, replace=TRUE))
  
  #sigma2 <- 1/rgamma(initNClusters,shape = nu_0/2,rate = nu_0*s_0^2/2)
  # the initial value is trivial for conjugate prior, since we are considering marginal likelihood:
  sigma2 <- 1/rgamma(initNClusters,shape = 1/2,rate = 1*1^2/2)
  ksi <- sapply(seq(1:initNClusters), 
                function(i) return(rmnorm(1, mean=mu_0, varcov=sigma2[i]*invL0, sqrt=NULL)))
  
  phi_0 <- runif(1,min = l_b,max = u_b)
  
  History <- vector("list", niterations)
  
  pb <- txtProgressBar(min = 0, max = niterations, style = 3)
  ## start Gibb's sampling
  for (iter in 1:niterations)
  {
    ## update z ##
    clusterSizes = table(as.factor(clusterAssign))
    nClusters = length(clusterSizes)
    
    for (i in 1:n)
    { # determine whether ith component is a singleton 
      if (clusterSizes[clusterAssign[i]] > 1){
        # if not a singleton, have nClusters + 1 choices
        clusterSizes[clusterAssign[i]] = clusterSizes[clusterAssign[i]] - 1
        # the probs for choosing exsiting clusters, inference doesn't change for a single individual:
        
        clusterProbs = sapply(1:nClusters, function(x) {
          log(sum(W[i, clusterAssign == x])) + 
            sum(dnorm(y[i,], phi%*%ksi[,x,drop=FALSE], sqrt(sigma2[x]),log = TRUE))
        })
        
        # the prob for choosing a new cluster, inference doesn't change for a single individual:
        
        
        Lambda_t <- phiTphi+Lambda_0
        mu_t <- inv(Lambda_t)%*%(t(phi)%*%t(y[i,,drop=FALSE])+Lambda_0%*%mu_0)
        a_t <- (t+nu_0)/2
        b_t <- 0.5*(t(y[i,])%*%y[i,]+mu0TL0mu0+nu_0*s_0^2-t(mu_t)%*%Lambda_t%*%mu_t)
        
        clusterProbs[nClusters+1] <- log_alpha + 
          ((-t/2)*log(2*pi)+0.5*(log(det(Lambda_0))-log(det(Lambda_t)))+nu_0/2*log(0.5*nu_0*s_0^2)-
             a_t*log(b_t)+lgamma(a_t)-lgamma(nu_0/2))
        
        clusterProbs<- exp(clusterProbs-max(clusterProbs))
        
        # choose the cluster number for ith observation
        cluster.i <- sample(1:(nClusters+1), size = 1, prob = clusterProbs)
        clusterAssign[i] <- cluster.i
        
        if (cluster.i > nClusters) {
          # it will be updated later. Any value is okay here.
          ksi <- cbind(ksi, rmnorm(1, mean = mu_0, invL0))
          sigma2 <- c(sigma2,1/rgamma(1,shape = nu_0/2,rate = (nu_0*s_0^2)/2))
        }
      } else {
        # if singleton, have nClusters choices
        clusterAssign[clusterAssign > clusterAssign[i]] <- 
          clusterAssign[clusterAssign > clusterAssign[i]] - 1
        ksi <- ksi[, -clusterAssign[i],drop=FALSE]
        sigma2 <- sigma2[-clusterAssign[i]]
        # the probs for choosing exsiting clusters, inference doesn't change for a single individual:
        clusterProbs = sapply(1:(nClusters-1), function(x) {
          log(sum(W[i, clusterAssign == x])) + 
            sum(dnorm(y[i,], phi%*%ksi[,x,drop=FALSE], sqrt(sigma2[x]),log = TRUE))
        })
        
        Lambda_t <- phiTphi+Lambda_0
        mu_t <- inv(Lambda_t)%*%(t(phi)%*%t(y[i,,drop=FALSE])+Lambda_0%*%mu_0)
        a_t <- (t+nu_0)/2
        b_t <- 0.5*(t(y[i,])%*%y[i,]+mu0TL0mu0+nu_0*s_0^2-t(mu_t)%*%Lambda_t%*%mu_t)
        
        clusterProbs[nClusters] <- log_alpha + 
          ((-t/2)*log(2*pi)+0.5*(log(det(Lambda_0))-log(det(Lambda_t)))+nu_0/2*log(0.5*nu_0*s_0^2)-
             a_t*log(b_t)+lgamma(a_t)-lgamma(nu_0/2))
        
        clusterProbs <- exp(clusterProbs-max(clusterProbs))
        
        # choose the cluster number for ith observation
        cluster.i <- sample(1:nClusters, size = 1, prob = clusterProbs)
        clusterAssign[i] <- cluster.i
        
        if (cluster.i == nClusters) {
          ksi <- cbind(ksi, rmnorm(1, mean = mu_0, invL0))
          sigma2 <- c(sigma2,1/rgamma(1,shape = nu_0/2,rate = (nu_0*s_0^2)/2))
        }
      }
      
      clusterSizes <- table(as.factor(clusterAssign))
      nClusters <- length(clusterSizes)
    }
    
    ## update ksi and sigma2 ##
    for (r in 1:nClusters){
      w_j <- d_matrix[which(clusterAssign == r),which(clusterAssign == r),drop=FALSE]
      r_sz <- clusterSizes[r]
      Lambda_r <- sum(diag(r_sz)-phi_0*w_j)*phiTphi+Lambda_0
      mu_r <- inv(Lambda_r)%*%(t(phi)%*%t(y[which(clusterAssign == r),,drop=FALSE])%*%(diag(r_sz)-phi_0*w_j)%*%
                                 matrix(1,nrow = r_sz,ncol = 1)+
                                 Lambda_0%*%mu_0)
      
      a_r <- (t*clusterSizes[r]+nu_0)/2
      b_r <- 0.5*(mu0TL0mu0+nu_0*s_0^2/2+
                    sum(diag(y[which(clusterAssign == r),,drop=FALSE]%*%
                               t(y[which(clusterAssign == r),,drop=FALSE])%*%
                               (diag(r_sz)-phi_0*w_j)))-
                    t(mu_r)%*%Lambda_r%*%mu_r)
      sigma2[r] <- 1/rgamma(1,shape = a_r,rate = b_r)
      ksi[,r] = rmnorm(1, mean = mu_r, varcov = sigma2[r]*inv(Lambda_r))
    }
    
    # standardize y according to their cluster parameter
    y_std <- t(y)
    for(j in 1:n){
      y_std[,j] <- (y_std[,j]-phi%*%ksi[,clusterAssign[j]])/sqrt(sigma2[clusterAssign[j]])
    }
    
    # update phi_0:
    phi_p <- min(max(rnorm(1,mean = phi_0,sd = mh_sz),l_b+1e-4),u_b-1e-4)
    log.acc <- log.dmn(y_std = y_std,phi_0 = phi_p,d_matrix)-log.dmn(y_std = y_std,phi_0 = phi_0,d_matrix)
    if(log(runif(1)) < log.acc){
      phi_0 <- phi_p
    }
    
    History[[iter]] <- list(zout = clusterAssign, ksiout = ksi,sigma2out=sigma2,phi_0out=phi_0)
    if(iter%%1e2==0){
      cat(" iteration:", iter,"\n",clusterAssign,"\n")
    }
    
    
    #cat(" iteration:", iter,"\n")
    setTxtProgressBar(pb, iter)
  }# for loop over iterations
  
  list(Iterates = History)
}

getDahl <- function(MFMfit, burn)
{
  ################################################################
  
  ## Input: MFMfit = the result from GWCRP ##
  ##        burn = the number of burn-in iterations ##
  
  ## Output:
  ##         zout = estimated clustering configuration, a n by 1 vector##
  ##         Qout = estimated probability matrix, a k by k matrix ##
  
  #################################################################
  iters <- MFMfit$Iterates[-(1:burn)]
  n <- length(iters[[1]][[1]])
  niters <- length(iters)
  membershipMatrices <- lapply(iters, function(x){
    clusterAssign <- x$zout
    outer(clusterAssign, clusterAssign, FUN = "==")
  })
  membershipAverage <- Reduce("+", membershipMatrices)/niters
  SqError <- sapply(membershipMatrices, function(x, av) sum((x - av)^2),
                    av = membershipAverage)
  DahlIndex <- which.min(SqError)
  DahlAns <- iters[[DahlIndex]]
  attr(DahlAns, "iterIndex") <- burn + DahlIndex
  attr(DahlAns, "burnin") <- burn
  DahlAns
}

LPML <- function(y,phi,d_matrix,MFMfit,nburn){
  pst_param <- MFMfit$Iterates[-(1:nburn)]
  niter <- length(pst_param)
  k <- dim(phi)[2]
  n <- dim(y)[1]
  t <- dim(y)[2]
  
  y0 <- t(y)
  log_CPO_inv <- c()
  for(i in 1:n){
    log_b <- c()
    for(j in 1:niter){
      ksi_j <- pst_param[[j]]$ksiout[,pst_param[[j]]$zout[i]]
      sigma_j <- sqrt(pst_param[[j]]$sigma2out[pst_param[[j]]$zout[i]])
      phi_j <- pst_param[[j]]$phi_0out
      log_b[j] <- sum(dnorm(y0[,i],phi%*%ksi_j,sigma_j,log = TRUE))
    }
    log_CPO_inv[i] <- -min(log_b)+log(sum(exp(min(log_b)-log_b)))-log(niter)
  }
  lpml <- -sum(log_CPO_inv)
  return(lpml)
}