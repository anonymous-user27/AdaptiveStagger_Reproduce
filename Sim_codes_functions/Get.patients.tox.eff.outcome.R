library(truncdist)

##Need to get this patient's outcome 
##generate the tox and eff outcomes
get.pts.tox.eff.outcome <- function(d, npts, mu_E, mu_T, Sigma, DLT_window, efficacy_window, dist, alpha=c(0.5,0.5), rho.time=0){
  Z <- mvrnorm(npts, c(mu_T[d], mu_E[d]), Sigma)
  
  if(npts==1){
    T.out <-  as.numeric(Z[1] >=0)
    E.out <-  as.numeric(Z[2] >=0)
  }else{
    T.out <-  as.numeric(Z[,1] >=0)
    E.out <-  as.numeric(Z[,2] >=0)
  }
 
  d.pT <- pnorm( mu_T[d])
  d.pE <- pnorm(mu_E[d])
  ##Generate time to efficacy or toxicity
  if(dist=="Uniform"){
    T.time.random <- runif(npts, 0, DLT_window)
    E.time.random <- runif(npts,0, efficacy_window)
    

    # Define the correlation matrix
    cor_matrix <- matrix(c(1, rho.time, rho.time, 1), nrow = 2)
    
    # Perform Cholesky decomposition
    L <- chol(cor_matrix)
    
    
    # Generate independent uniform variables
    uniforms <- matrix(c(T.time.random/DLT_window, E.time.random/efficacy_window), ncol = 2)
    
    # Transform uniform to standard normal (using probit transformation)
    normals <- qnorm(uniforms)
    
    # Apply the Cholesky decomposition to induce correlation
    correlated_normals <- normals %*% L
    
    # Transform back to uniform distribution
    correlated_uniforms <- pnorm(correlated_normals)
    
    T.time.random.cor <- correlated_uniforms[,1]*DLT_window
    E.time.random.cor <- correlated_uniforms[,2]*efficacy_window
    
      
    T.time=ifelse(T.out==1, T.time.random.cor,DLT_window+0.001)
    E.time=ifelse(E.out==1, E.time.random.cor,efficacy_window+0.001)
    
  }
  ## solve parameters for Weibull given pi=1-S(T) and phalft=1-S(T/2)
  if(dist=="Weibull"){
    #Generate time to DLT
    # T.shape=log(log(1-d.pT)/log(1-d.pT/2))/log(2)
    T.pihalft=d.pT*alpha[1]
    T.shape=log(log(1-d.pT)/log(1-T.pihalft))/log(2)
    T.scale=DLT_window/exp(log(-log(1-d.pT))/T.shape)
    #  print(paste("T.scale=",T.scale, "and T.shape=",T.shape))
    T.time=ifelse(T.out==1,rtrunc(npts, "weibull", a=0, b = DLT_window, shape=T.shape,scale=T.scale),DLT_window+0.001);
    
    #   #Generate time to tumor response
    E.pihalft=d.pE*alpha[2]
    E.shape=log(log(1-d.pE)/log(1-E.pihalft))/log(2)
    E.scale=efficacy_window/exp(log(-log(1-d.pE))/E.shape)
    #print(paste("E.scale=",E.scale, "and E.shape=",E.shape))
    E.time=ifelse(E.out==1,rtrunc(npts, "weibull", a=0, b = efficacy_window, shape=E.shape,scale=E.scale),efficacy_window+0.001);
  }
  list(tox=T.out,t.tox=T.time, eff=E.out,t.eff=E.time)
}

