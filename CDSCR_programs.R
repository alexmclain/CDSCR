library(mvtnorm)
library(MASS)
library(statmod)
library(survival)


#### Wrapper to the optimization of the likelihood ####
CDSCR <- function(x.data,y.data,delta.data,knots.t=NULL,knots.y=NULL,knots.z=NULL,K1=NULL,K2=NULL,K3=NULL,B=20,theta=NULL,opt_meth="Nelder-Mead"){
  
  ##################################################################################################
  #Description: This program will run a current duration semi-competeing risk model with a  
  #             piecewise constant hazard function 
  #Usage:
  # CDSCR(T,TTFT,fail_ind,knots_T,knots_Y,knots_Z,K1,K2,K3,B,theta,opt_method)
  #
  #Arguments:
  # x.data:     Total current duration, the time from the begining of attempt to sampling.
  # y.data:     Time to the intermittent event (TTFT in the paper). If no intermittent event set =0.
  # delta.data: Indicator that the intermittent event occured.
  # knots.t:    Values of the knots for the Piecewise model for the distribution of T. Set to NULL to 
  #             use data driven knots (K1 must be set in this case).
  # knots.y:    Values of the knots for the Piecewise model for the distribution of Y. Set to NULL to 
  #             use data driven knots (K2 must be set in this case).
  # knots.z:    Values of the knots for the Piecewise model for the distribution of Z. Set to NULL to 
  #             use data driven knots (K3 must be set in this case).
  # K1:         Number of knots to use for T (ignored if knots_T is given).
  # K2:         Number of knots to use for Y (ignored if knots.y is given).
  # K3:         Number of knots to use for Z (ignored if knots.z is given).
  # B:          The number of quadrature nodes to use.  Default is 20.
  # theta:      Parameter of length 2 that gives the starting values for Corr(T,Y) and Corr(T,Z), 
  #             respectively (optional).
  # opt_meth:   the optimization method to use with 'optim'.  Default is "Nelder-Mead" other options 
  #             are "BFGS", "SANN".
  #Value:
  # The function will return a list with the following:
  # like_opt:   The results from optim.
  # knots:      A list containing the knots used for T, Y, and Z in the optimization.
  # st_vals:    A vector containing the starting values.
  
  # Cleaning the data
  y.data[delta.data==0] <- 0
  z.data <- x.data-y.data
  z.data[delta.data==0] <- 0
  
  ## Selecting knot locations (if not given). 
  if(is.null(knots.t)){
    knots.t <- quantile(x.data,probs = c(1:(K1-1))/K1)
    knots.t <- knots.t[knots.t>0]
    }
  if(is.null(knots.y)){
    knots.y <- quantile(y.data[y.data>0],probs = c(1:(K2-1))/K2)
    knots.y <- knots.y[knots.y>0]
    }
  if(is.null(knots.z)){
    knots.z <- quantile(z.data[z.data>0],probs = c(1:(K3-1))/K3)
    knots.z <- knots.z[knots.z>0]
  }
  
  ## Getting the nodes and weights for Legendre quadrature
  out <- gauss.quad(B,"legendre")
  nodes <- out$nodes
  weights <- out$weights
  
  ## Getting starting values for parameters.
  naive_est <- hazard_guess(x.data[z.data==0],knots.t) 
  naive.z   <- hazard_guess(z.data[z.data>0],knots.z)
  naive.z[naive.z<exp(-3)] <- exp(-3)
  naive.y   <- hazard_guess_ecdf(ecdf(y.data[z.data>0]),y.data[z.data>0],knots.y)
  naive.y   <- exp(log(naive.y)/0.85)
  naive.y[naive.y<exp(-3)] <- exp(-3)
  naive.yx  <- hazard_guess_ecdf(ecdf(y.data[z.data>0]),y.data[z.data>0],knots.t)
  naive.yx   <- exp(log(naive.yx)/0.85)
  naive_est_x <- naive_est - naive.yx
  naive_est_x[length(naive_est_x)] <- naive_est_x[length(naive_est_x)-1]
  naive_est_x[naive_est_x<exp(-3) | is.nan(naive_est_x) | is.na(naive_est_x)] <- exp(-3)
  
  ## Starting values correlation (if given)
  corrxy.T <- corrxz.T <- 1/4
  if(!is.null(theta)){
    corrxy.T <- max(c(abs(theta[1])/(1-abs(theta[1])),1e-4))
    corrxz.T <- max(c(abs(theta[2])/(1-abs(theta[2])),1e-4))
  }
  
  ## Starting values
  par<-c(log(c(naive_est_x,naive.y,naive.z,corrxy.T,corrxz.T)))
  
  ## Running optimization
  ## Optimization for grouped data, which uses full adaptive quadrature for the denominator (used in the data analysis).
  if(any(x.data==0) | all(is.integer(x.data))){
    cat("Grouped version running. \n")
    like_opt <-stats::optim(par,like_CDSCR_2corr_grp,x.data=x.data,y.data=y.data,z.data=z.data,delta.data=delta.data,knots.x=knots.t,knots.y=knots.y,knots.z=knots.z,nodes=nodes,weights=weights,hessian = FALSE,control = list(maxit=10000),method = opt_meth)
  }
  ## Optimization for continuous data, which uses a mixed quadrature approach for the denominator (used in the simulations).
  if(!(any(x.data==0) | all(is.integer(x.data)))){
    cat("Continuous version running. \n")
    like_opt <-stats::optim(par,like_CDSCR_2corr,x.data=x.data,y.data=y.data,z.data=z.data,delta.data=delta.data,knots.x=knots.t,knots.y=knots.y,knots.z=knots.z,nodes=nodes,weights=weights,hessian = FALSE,control = list(maxit=10000),method = opt_meth)
  }
  
  res <- list(like_opt = like_opt,knots.t=knots.t,knots.y=knots.y,knots.z=knots.z,st_vals=par)
  return(res)
}


haz_func <- function(t,alpha,lambda){
  t_0 <- c(0,t[-length(t)])
  cum_t <- (t/lambda)^alpha
  cum_t0 <- (t_0/lambda)^alpha
  return((cum_t - cum_t0)/(t-t_0))
}



#### Functions for generating data ####

CDSCR_Dep_PC <- function(n,knots.x,alpha.x,knots.y,alpha.y,knots.z,alpha.z,theta,tau){
  T        <- numeric(0)
  fail_ind <- numeric(0)
  TTF      <- numeric(0)
  mn_est   <- min(c(1/alpha.x,1/alpha.y))
  Sigma <- diag(3)
  Sigma[1,2] <- Sigma[2,1] <- theta[1]
  Sigma[1,3] <- Sigma[3,1] <- theta[3]
  Sigma[2,3] <- Sigma[3,2] <- theta[2]
  
  for(i in 1:n){
    S_t   <- 0
    S_vec <- NULL
    T_vec<- NULL
    t_len <- ceiling(tau/(mn_est*0.9))
    t_TTF <- NULL
    while(S_t < tau){
      rawvars <- mvrnorm(n=t_len, mu=c(0,0,0), Sigma=Sigma)
      u       <- pnorm(rawvars)
      r_TTP   <- PC_quan(u[,1],knots.x,alpha.x)
      r_TTF   <- PC_quan(u[,2],knots.y,alpha.y)
      r_TTP_F <- PC_quan(u[,3],knots.z,alpha.z)
      r_out   <- apply(cbind(r_TTF, r_TTP),1,min)
      r_ind   <- 1*I(r_TTP<r_TTF) + 2*I(r_TTP>r_TTF) 
      r_out   <- r_out + r_TTP_F*I(r_ind==2)
      t_l     <- cumsum(r_out)
      S_vec   <- c(S_vec, S_t + t_l)
      T_vec   <- c(T_vec, r_out)
      S_t     <- S_t + max(t_l)
      t_TTF   <- c(t_TTF,r_TTF)
    }
    A_t <- tau - S_vec[S_vec < tau]
    CD_t<- min(A_t)
    T <- c(T, CD_t)
    lastTTF <- t_TTF[S_vec > tau]
    ttf_ind <- 1*I(CD_t > lastTTF[1])
    TTF <- c(TTF,lastTTF[1]*I(ttf_ind==1))
    fail_ind <- c(fail_ind, ttf_ind)
  }
  cbind(T,TTF,fail_ind)
}


CDSCR_Dep_Weib <- function(n,alpha,lambda,theta,tau){
  T        <- numeric(0)
  fail_ind <- numeric(0)
  TTF      <- numeric(0)
  lambda.x <- lambda[1]
  lambda.y <- lambda[2]
  lambda.z <- lambda[3]
  alpha.x <- alpha[1]
  alpha.y <- alpha[2]
  alpha.z <- alpha[3]
  mn_est   <- lambda.x*gamma(1+1/alpha.x)
  Sigma <- diag(3)
  Sigma[1,2] <- Sigma[2,1] <- theta[1]
  Sigma[1,3] <- Sigma[3,1] <- theta[3]
  Sigma[2,3] <- Sigma[3,2] <- theta[2]
  
  for(i in 1:n){
    S_t   <- 0
    S_vec <- NULL
    T_vec<- NULL
    t_len <- ceiling(tau/(mn_est*0.9))
    t_TTF <- NULL
    while(S_t < tau){
      rawvars <- mvrnorm(n=t_len, mu=c(0,0,0), Sigma=Sigma)
      u       <- pnorm(rawvars)
      r_TTP   <- qweibull(u[,1],alpha.x,lambda.x)
      r_TTF   <- qweibull(u[,2],alpha.y,lambda.y)
      r_TTP_F <- qweibull(u[,3],alpha.z,lambda.z)
      r_out   <- apply(cbind(r_TTF, r_TTP),1,min)
      r_ind   <- 1*I(r_TTP<r_TTF) + 2*I(r_TTP>r_TTF) 
      r_out   <- r_out + r_TTP_F*I(r_ind==2)
      t_l     <- cumsum(r_out)
      S_vec   <- c(S_vec, S_t + t_l)
      T_vec   <- c(T_vec, r_out)
      S_t     <- S_t + max(t_l)
      t_TTF   <- c(t_TTF,r_TTF)
    }
    A_t <- tau - S_vec[S_vec < tau]
    CD_t<- min(A_t)
    T <- c(T, CD_t)
    lastTTF <- t_TTF[S_vec > tau]
    ttf_ind <- 1*I(CD_t > lastTTF[1])
    TTF <- c(TTF,lastTTF[1]*I(ttf_ind==1))
    fail_ind <- c(fail_ind, ttf_ind)
  }
  
  cbind(T,TTF,fail_ind)
}



CDSCR_Dep_Weib_T <- function(n,alpha,lambda,theta,tau){
  T        <- numeric(0)
  fail_ind <- numeric(0)
  TTF      <- numeric(0)
  lambda.x <- lambda[1]
  lambda.y <- lambda[2]
  lambda.z <- lambda[3]
  alpha.x <- alpha[1]
  alpha.y <- alpha[2]
  alpha.z <- alpha[3]
  mn_est   <- lambda.x*gamma(1+1/alpha.x)
  Sigma <- diag(3)
  Sigma[1,2] <- Sigma[2,1] <- theta[1]
  Sigma[1,3] <- Sigma[3,1] <- theta[3]
  Sigma[2,3] <- Sigma[3,2] <- theta[2]
  
  for(i in 1:n){
    S_t   <- 0
    S_vec <- NULL
    T_vec<- NULL
    t_len <- ceiling(tau/(mn_est*0.9))
    t_TTF <- NULL
    while(S_t < tau){
      rawvars <- rmvt(n=t_len,sigma=Sigma,df=20)
      u       <- pnorm(rawvars,log.p = TRUE)
      r_TTP   <- qweibull(u[,1],alpha.x,lambda.x,log.p = TRUE)
      r_TTF   <- qweibull(u[,2],alpha.y,lambda.y,log.p = TRUE)
      r_TTP_F <- qweibull(u[,3],alpha.z,lambda.z,log.p = TRUE)
      r_out   <- apply(cbind(r_TTF, r_TTP),1,min)
      r_ind   <- 1*I(r_TTP<r_TTF) + 2*I(r_TTP>r_TTF) 
      r_out   <- r_out + r_TTP_F*I(r_ind==2)
      t_l     <- cumsum(r_out)
      S_vec   <- c(S_vec, S_t + t_l)
      T_vec   <- c(T_vec, r_out)
      S_t     <- S_t + max(t_l)
      t_TTF   <- c(t_TTF,r_TTF)
    }
    A_t <- tau - S_vec[S_vec < tau]
    CD_t<- min(A_t)
    T <- c(T, CD_t)
    lastTTF <- t_TTF[S_vec > tau]
    ttf_ind <- 1*I(CD_t > lastTTF[1])
    TTF <- c(TTF,lastTTF[1]*I(ttf_ind==1))
    fail_ind <- c(fail_ind, ttf_ind)
  }
  
  cbind(T,TTF,fail_ind)
}

#### Functions of the distribution of T(x), Y(y),Z(z)    ##############

PC_quan <- function(u,knots,alpha){
  ### Quantile function.  Can be used to generate data.
  w <- -log(1-u)
  L_int <- knots - c(0,knots[-length(knots)])
  Int <- cumsum(alpha[-length(alpha)]*L_int)
  Uts = I(w<=Int[1])*(w/alpha[1]) 
  for(k in 2:length(knots)){
    Inc <- ((w-Int[k-1])/alpha[k] +knots[k-1])
    Uts = Uts + I(w <= Int[k] & w > Int[k-1])*Inc
  }
  Inc <- ((w-Int[k])/alpha[length(alpha)] +knots[k])
  Uts = Uts + I(w > Int[length(knots)])*Inc
  return(Uts)
}


PC_quan2 <- function(z,knots,alpha){
  u <- pnorm(z)
  ### Quantile function.  Can be used to generate data.
  w <- -log(1-u)
  L_int <- knots - c(0,knots[-length(knots)])
  Int <- cumsum(alpha[-length(alpha)]*L_int)
  Uts = I(w<=Int[1])*(w/alpha[1]) 
  for(k in 2:length(knots)){
    Inc <- ((w-Int[k-1])/alpha[k] +knots[k-1])
    Uts = Uts + I(w <= Int[k] & w > Int[k-1])*Inc
  }
  Inc <- ((w-Int[k])/alpha[length(alpha)] +knots[k])
  Uts = Uts + I(w > Int[length(knots)])*Inc
  return(Uts)
}



#### The main pdf and survival functions ####
PC_surv <- function(t,knots,alpha){
  ### Survival Function
  
  ###Increments between knots
  L_int <- knots - c(0,knots[-length(knots)])
  
  ###Finding the first knot bigger then t
  kn_mat <- matrix(knots,length(t),length(knots),byrow = TRUE)
  t_kn_mat <- kn_mat-t
  t_kn_mat[t_kn_mat<0] <- max(knots)
  ### Which knot (knot # counting the 0 knot) is the last smaller then t.
  first_kn_bg <- apply(t_kn_mat,1,which.min)
  first_kn_bg[t>max(knots)] <- length(knots)+1
  
  ###Getting the alpha for the interval with t.
  last_alpha <- alpha[first_kn_bg] 
  
  ###Getting the increment beyond the last knot
  last_incr <- t - c(0,knots)[first_kn_bg] 
  
  ###Calculating the survival function at the knots.  Adding one for zero too.
  ###This uses the memoryless property of the exponential distribution.
  sur_func_incr <- 1 - pexp(L_int,alpha[-length(alpha)])
  kn_sur_vals <- c(1,cumprod(sur_func_incr))
  
  ###Finding survival value at last knot smaller then t
  surv_last_kn <- kn_sur_vals[first_kn_bg]
  
  ###Adding the last increment
  surv_val <- surv_last_kn*(1-pexp(last_incr,last_alpha))
  surv_val
}



PC_pdf <- function(t,knots,alpha){
  ### Probability density function.
  ### This function is the basically same as the survival function. It returns the  
  ### survival function multiplied by alpha(t).
  
  
  ###Increments metween knots
  L_int <- knots - c(0,knots[-length(knots)])
  
  ###Finding the first knot bigger then t
  kn_mat <- matrix(knots,length(t),length(knots),byrow = TRUE)
  t_kn_mat <- kn_mat-t
  t_kn_mat[t_kn_mat<0] <- max(knots)
  first_kn_bg <- apply(t_kn_mat,1,which.min)
  first_kn_bg[t>max(knots)] <- length(knots)+1
  
  ###Getting the alpha for the interval with t.
  last_alpha <- alpha[first_kn_bg] 
  
  ###Getting the increment beyond the last knot
  last_incr <- t - c(0,knots)[first_kn_bg] 
  
  ###Calculating the survival function at the knots.  Adding one for zero too.
  sur_func_incr <- 1 - pexp(L_int,alpha[-length(alpha)])
  kn_sur_vals <- c(1,cumprod(sur_func_incr))
  
  ###Finding survival value at last knot smaller then t
  surv_last_kn <- kn_sur_vals[first_kn_bg]
  
  ###Adding the last increment
  surv_val <- surv_last_kn*(1-pexp(last_incr,last_alpha))
  
  last_alpha*surv_val
}


#### The functions for generating data ####
posi.function.inv<-function(p,knots,alpha){
  h<--log(1-p)
  knots1<-c(0,knots)
  knots.diff<-knots1[-1]-knots1[-length(knots1)]
  knots.alpha<-knots.diff*alpha[-length(alpha)]
  alpha.knots<-cumsum(knots.alpha)
  alpha.knots<-c(0,alpha.knots)
  if (h>max(alpha.knots)){return(alpha[length(alpha)])}
  else {
    for(i in (1:(length(knots1)-1))){
      if(h<=alpha.knots[i+1] & h>alpha.knots[i]) {return(alpha[i])}
    }
  }
}


posi.function.inv2 <- function(u,knots,alpha){
  ### Quantile function.  Can be used to generate data.
  w <- -log(1-u)
  L_int <- knots - c(0,knots[-length(knots)])
  Int <- cumsum(alpha[-length(alpha)]*L_int)
  Uts = I(w<=Int[1])*alpha[1]
  for(k in 2:length(knots)){
    Uts[w <= Int[k] & w > Int[k-1]] <- alpha[k]
  }
  Uts[w>max(Int)] <- alpha[length(alpha)]
  return(Uts)
}

posi.function<-function(x,knots,alpha){
  knots1<-c(0,knots)
  if (x>max(knots1)){return(alpha[length(alpha)])}
  else {
    for(i in (1:(length(knots1)-1))){
      if(x<=knots1[i+1] & x>knots1[i]) {return(alpha[i])}
    }
  }
}



#### Various joint and conditional distribution ####

surv.min.T.Y<-function(x,r,knots.x,alpha.x,knots.y,alpha.y){
  F.T<-PC_surv(x,knots.x,alpha.x)
  F.Y<-PC_surv(x,knots.y,alpha.y)
  u.Func<- qnorm(1-F.T)
  v.Func<- qnorm(1-F.Y)
  uv_mat <- cbind(-u.Func,-v.Func)
  test <- apply(uv_mat,1,pmvnorm,upper=Inf,mean = c(0,0),corr = matrix(c(1,r,r,1),2,2))
  C_bar <- F.T + F.Y -1 + test
  C_bar[C_bar<0]<-0
  C_bar[is.na(C_bar)]<-0
  return(C_bar)
}


density.T.given.Y <- function(x,y,r,knots.x,alpha.x,knots.y,alpha.y){
  u = qnorm(1-PC_surv(x,knots.x,alpha.x))
  v = qnorm(1-PC_surv(y,knots.y,alpha.y))
  u[!is.finite(u)] <- sign(u[!is.finite(u)])*1e20
  v[!is.finite(v)] <- sign(v[!is.finite(v)])*1e16
  ans1 <- dnorm(u,r*v,sqrt(1-r^2))*PC_pdf(x,knots.x,alpha.x)/dnorm(u)
  ans1[is.na(ans1)] <- 0
  ans1
}


surv.Z.given.T <- function(z,x,corrxz,knots.x,alpha.x,knots.z,alpha.z){
  F_bar_z <- PC_surv(z,knots.z,alpha.z)
  F_bar_x <- PC_surv(x,knots.x,alpha.x)
  u = qnorm(1-F_bar_z)
  v = qnorm(1-F_bar_x)
  u[!is.finite(u)] <- sign(u[!is.finite(u)])*1e20
  v[!is.finite(v)] <- sign(v[!is.finite(v)])*1e16
  ans  <- (pnorm(u,mean = corrxz*v,sd = sqrt(1-corrxz^2),lower.tail = FALSE))
  ans
}



cdf.Y.given.T <- function(x,r,knots.x,alpha.x,knots.y,alpha.y){
  F_bar_y <- PC_surv(x,knots.y,alpha.y)
  F_bar_x <- PC_surv(x,knots.x,alpha.x)
  u = qnorm(1-F_bar_y)
  v = qnorm(1-F_bar_x)
  u[!is.finite(u)] <- sign(u[!is.finite(u)])*1e16
  v[!is.finite(v)] <- sign(v[!is.finite(v)])*1e20
  ans <- (pnorm(u,mean = r*v,sd = sqrt(1-r^2),lower.tail = TRUE))
  ans
}


expect.min.T.Y<-function(r,knots.x,alpha.x,knots.y,alpha.y){
  integrate(surv.min.T.Y,0,Inf,r,knots.x,alpha.x,knots.y,alpha.y)$value
}


expect.min.T.Y_grp<-function(r,knots.x,alpha.x,knots.y,alpha.y){
  integrate(surv.min.T.Y,1,Inf,r,knots.x,alpha.x,knots.y,alpha.y)$value +1
}


#### The functions is to calculate the inner part of the integration in the denominator ####

denom.integral.inner.legendre<-function(x,corrxy,corrxz,knots.x,alpha.x,knots.y,alpha.y,knots.z,alpha.z,nodes,weights)
{
  L <- length(x)
  n <- length(nodes)
  trans_nodes <- (nodes+1)/(1-nodes)
  dp <- 1/(1-nodes)*(1+trans_nodes)
  node_vec    <- rep(trans_nodes,L)
  x_vec <- rep(x,each=n)
  inner_mat <- matrix(surv.Z.given.T(node_vec,x_vec,corrxz,knots.x,alpha.x,knots.z,alpha.z),nrow=L,ncol = n,byrow = TRUE)
  t1 <- t(t(weights*dp)%*%t(inner_mat))*cdf.Y.given.T(x,corrxy,knots.x,alpha.x,knots.y,alpha.y)*PC_pdf(x,knots.x,alpha.x)
  t1 
}

denom.integral.inner<-function(x,corrxy,corrxz,knots.x,alpha.x,knots.y,alpha.y,knots.z,alpha.z)
{
  (sapply(x,function(x,corrxz,knots.x,alpha.x,knots.z,alpha.z){integrate(surv.Z.given.T,1,Inf,x,corrxz,knots.x,alpha.x,knots.z,alpha.z)$value},corrxz,knots.x,alpha.x,knots.z,alpha.z)+1)*cdf.Y.given.T(x,corrxy,knots.x,alpha.x,knots.y,alpha.y)*PC_pdf(x,knots.x,alpha.x)
  
}

#### The functions is to calculate the inner part of the integration in the numerator ####

nume.integral.inner.legendre<-function(temp.data,nodes,corrxy,corrxz,knots.x,alpha.x,knots.y,alpha.y,knots.z,alpha.z){
  N <- dim(temp.data)[1]
  n <- length(nodes)
  trans_nodes <- (nodes+1)/(1-nodes)
  dp <- 1/(1-nodes)*(1+trans_nodes)
  x <- rep(trans_nodes,N)+rep(temp.data[,1],each=n)
  y <- rep(temp.data[,1],each=n)
  z <- rep(temp.data[,2],each=n)
  ans<-surv.Z.given.T(z,x,corrxz,knots.x,alpha.x,knots.z,alpha.z)*density.T.given.Y(x,y,corrxy,knots.x,alpha.x,knots.y,alpha.y)*PC_pdf(y,knots.y,alpha.y)*dp
  ans_mat <- matrix(ans,nrow=N,ncol=n,byrow = TRUE)
  return(ans_mat)
}

nume.integral.inner.legendre_grp <-function(temp.data,nodes,corrxy,corrxz,knots.x,alpha.x,knots.y,alpha.y,knots.z,alpha.z){
  N <- dim(temp.data)[1]
  n <- length(nodes)
  trans_nodes <- (nodes+1)/(1-nodes)
  dp <- 1/(1-nodes)*(1+trans_nodes)
  x <- rep(trans_nodes,N)+rep(temp.data[,1],each=n)
  y <- rep(temp.data[,1],each=n)
  z <- rep(temp.data[,2],each=n)
  ans<-surv.Z.given.T(z,x,corrxz,knots.x,alpha.x,knots.z,alpha.z)*density.T.given.Y(x,y,corrxy,knots.x,alpha.x,knots.y,alpha.y)*(PC_surv(y-1,knots.y,alpha.y)-PC_surv(y,knots.y,alpha.y))*dp
  ans[y==0] <- 0
  ans[x==0] <- 0
  ans_mat <- matrix(ans,nrow=N,ncol=n,byrow = TRUE)
  return(ans_mat)
}



#### Likelihood Functions ####

# Quicker version for continuous data
like_CDSCR_2corr <-function(par,x.data,y.data,z.data,delta.data,knots.x,knots.y,knots.z,nodes,weights){
  Ppar <<- par
  par <- exp(par)
  t_alpha.x <- par[1:(length(knots.x)+1)]
  t_alpha.y <- par[(length(knots.x)+2):(length(knots.x)+length(knots.y)+2)]
  t_alpha.z <- par[(length(knots.x)+length(knots.y)+3):(length(knots.x)+length(knots.y)+length(knots.z)+3)]
  corrxy <- 1*(sign(par[length(par)-1]))
  corrxz <- 1*(sign(par[length(par)]))
  if(is.finite(par[length(par)-1])) corrxy <- -(par[length(par)-1]/(1+par[length(par)-1]))
  if(is.finite(par[length(par)])) corrxz <- (par[length(par)]/(1+par[length(par)]))
  
  
  N <- length(x.data)
  numerator <- rep(0,N)
  
  numerator1<-surv.min.T.Y(x.data[delta.data==0],corrxy,knots.x,t_alpha.x,knots.y,t_alpha.y)
  
  temp.data<-cbind(y.data[delta.data==1],z.data[delta.data==1])
  ans_mat <- nume.integral.inner.legendre(temp.data,nodes,corrxy=corrxy,corrxz=corrxz,knots.x=knots.x,alpha.x=t_alpha.x,knots.y=knots.y,alpha.y=t_alpha.y,knots.z=knots.z,alpha.z=t_alpha.z)
  numerator2 <- t(t(weights)%*%t(ans_mat))
  
  numerator[delta.data==0] <- numerator1
  numerator[delta.data==1] <- numerator2
  
  t1 <- try(ETminY <- expect.min.T.Y(r=corrxy,knots.x=knots.x,alpha.x=t_alpha.x,knots.y=knots.y,alpha.y=t_alpha.y),silent = TRUE)
  if(!is.null(attr(t1,"class"))){ETminY <- Inf}
  
  t1 <- try(EZgiven <- integrate(denom.integral.inner.legendre,0,Inf,corrxy,corrxz,knots.x,t_alpha.x,knots.y,t_alpha.y,knots.z,t_alpha.z,nodes,weights)$value,silent = TRUE)
  if(!is.null(attr(t1,"class"))){EZgiven <- Inf}
  
  denominator<- ETminY + EZgiven
  
  likelihood1<-numerator/denominator
  log.likelihood1<-log(likelihood1)
  log.likelihood1[is.infinite(log.likelihood1)]<-0
  log.likelihood1[log.likelihood1==0]<- (-1000)
  sum(-log.likelihood1)
}


# Slower version for grouped data. 
like_CDSCR_2corr_grp <-function(par,x.data,y.data,z.data,delta.data,knots.x,knots.y,knots.z,nodes,weights){
  
  Ppar <<- par
  par <- exp(par)
  nsll <- Inf
  if(all(is.finite(par))){
    t_alpha.x <- par[1:(length(knots.x)+1)]
    t_alpha.y <- par[(length(knots.x)+2):(length(knots.x)+length(knots.y)+2)]
    t_alpha.z <- par[(length(knots.x)+length(knots.y)+3):(length(knots.x)+length(knots.y)+length(knots.z)+3)]
    corrxy <- 1*(sign(par[length(par)-1]))
    corrxz <- 1*(sign(par[length(par)]))
    if(is.finite(par[length(par)-1])) corrxy <- -(par[length(par)-1]/(1+par[length(par)-1]))
    if(is.finite(par[length(par)])) corrxz <- (par[length(par)]/(1+par[length(par)]))
    
    N <- length(x.data)
    numerator <- rep(0,N)
    
    numerator1<-surv.min.T.Y(x.data[delta.data==0],corrxy,knots.x,t_alpha.x,knots.y,t_alpha.y)
    
    temp.data<-cbind(y.data[delta.data==1],z.data[delta.data==1])
    ans_mat <- nume.integral.inner.legendre_grp(temp.data,nodes,corrxy=corrxy,corrxz=corrxz,knots.x=knots.x,alpha.x=t_alpha.x,knots.y=knots.y,alpha.y=t_alpha.y,knots.z=knots.z,alpha.z=t_alpha.z) 
    numerator2 <- t(t(weights)%*%t(ans_mat))
    
    numerator[delta.data==0] <- numerator1
    numerator[delta.data==1] <- numerator2
    
    t1 <- try(ETminY <- expect.min.T.Y_grp(r=corrxy,knots.x=knots.x,alpha.x=t_alpha.x,knots.y=knots.y,alpha.y=t_alpha.y),silent = TRUE)
    if(!is.null(attr(t1,"class"))){ETminY <- Inf}
    
    t1 <- try(EZgiven <- integrate(denom.integral.inner,0,Inf,corrxy,corrxz,knots.x,t_alpha.x,knots.y,t_alpha.y,knots.z,t_alpha.z)$value,silent = TRUE)
    
    if(!is.null(attr(t1,"class"))){
      Sigma <- diag(3)
      Sigma[1,2] <- Sigma[2,1] <- corrxy
      Sigma[1,3] <- Sigma[3,1] <- corrxz
      Sigma[2,3] <- Sigma[3,2] <- corrxy*corrxz
      u <- mvrnorm(n=1e6, mu=c(0,0,0), Sigma=Sigma)
      X_vals <- ceiling(gen_func(pnorm(u[,1]),knots.x,t_alpha.x))
      Y_vals <- ceiling(gen_func(pnorm(u[,2]),knots.y,t_alpha.y))
      Z_vals <- ceiling(gen_func(pnorm(u[,3]),knots.z,t_alpha.z))
      EZgiven <- mean(Z_vals*I(Y_vals<X_vals))
      ETminY <- mean(apply(cbind(X_vals,Y_vals),1,min))
    }
    denominator <- Inf
    if(is.finite(ETminY) & is.finite(EZgiven)){denominator <- ETminY + EZgiven}
    
    likelihood1<-numerator/denominator
    log.likelihood1<-log(likelihood1)
    log.likelihood1[is.infinite(log.likelihood1)]<-0
    log.likelihood1[log.likelihood1==0]<--1000
    nsll <- sum(-log.likelihood1)
  }
  return(nsll)
}


unique_knots <- function(knots){
  knots_fl <- floor(knots)
  if(length(unique(knots_fl))< length(knots)){
    knots_last <- c(knots_fl[-1],Inf)
    knots <- knots[knots_fl != knots_last]
    knots_fl <- knots_fl[knots_fl != knots_last]
  }
  return(list(knots=knots,knots_fl=knots_fl))
}


hazard_guess <- function(surv_dat,knots){
  n <- length(sort(unique(surv_dat)))
  test.z <- find.gren_cont2(surv_dat,CEN=rep(1,length(surv_dat)), plot=FALSE)
  if(length(test.z[c(0,sort(unique(surv_dat)))< sort(unique(surv_dat))[ceiling(n*0.05)]])>0){test.z[c(0,sort(unique(surv_dat)))< sort(unique(surv_dat))[ceiling(n*0.05)]] <- min(test.z[c(0,sort(unique(surv_dat)))< sort(unique(surv_dat))[ceiling(n*0.05)]])}
  GREN <- test.z/test.z[1]
  step_est.z <- stepfun(c(0,sort(unique(surv_dat))),c(0,GREN),right=FALSE)
  mx_z <- max(surv_dat)
  z_surv <- step_est.z(c(knots,mx_z))
  z_cumhaz <- -log(z_surv)
  z_haz <- (z_cumhaz - c(0,z_cumhaz[-length(z_cumhaz)]))/c(c(knots,mx_z) - c(0,knots))
  return(z_haz)
}

hazard_guess_ecdf <- function(ecdf,surv_dat,knots){
  
  mx_z <- max(surv_dat)*1.3
  y_surv <- c(1-ecdf(knots),1/length(surv_dat))
  y_cumhaz <- -log(y_surv)
  y_haz <- (y_cumhaz - c(0,y_cumhaz[-length(y_cumhaz)]))/c(c(knots,mx_z) - c(0,knots))
  return(y_haz)
}


gen_func <- function(u,knots,alpha){
  w <- -log(1-u)
  L_int <- knots - c(0,knots[-length(knots)])
  Int <- cumsum(alpha[-length(alpha)]*L_int)
  Uts = I(w<=Int[1])*(w/alpha[1]) 
  k = 1
  if(length(knots)>1){
    for(k in 2:length(knots)){
      Inc <- ((w-Int[k-1])/alpha[k] +knots[k-1])
      Uts = Uts + I(w <= Int[k] & w > Int[k-1])*Inc
    }
  }
  Inc <- ((w-Int[k])/alpha[length(alpha)] +knots[k])
  Uts = Uts + I(w > Int[length(knots)])*Inc
  return(Uts)
}








#### Functions from McLain, A. C., Sundaram, R., Thoma, M., and Louis, G. (2014). Semiparametric modeling of grouped current duration data with preferential reporting. Statistics in medicine, 33(23), 3961-3972. ####

CD_surv <- function(T,X=0,tau=Inf,weights=1,method,knots=NULL,Asymp_SE = TRUE, opt_meth = "BFGS"){
  
  ##################################################################################################
  #Description: This program will run 'Nonparametric', 'Semiparametric' and 'Piecewise' current 
  # duration analyses. 
  #Usage:
  # CD_surv(T,X,tau,w,
  #          method=c('Nonparametric','Semiparametric','Piecewise'),
  #          knots)
  #Arguments:
  # T:        Current duration values.
  # X:        Matrix of covariates not including an intercept (not required).
  # tau:      The type I censoring value (length 1), or a vector (length n) of censoring values (not 
  #           required).
  # weights:  Survey sampling weights (not required).
  # method:   What type of analysis should be done 'Nonparametric' (no covariates), 'Semiparametric' 
  #           (with or without covariates), or 'Piecewise' (with or without covariates). See 
  #           'Details' for further information.
  # knots:    values of the knots for the Piecewise model (only required for Piecewise model).
  # Asymp_SE: if TRUE standard errors will be estimated using the inverse of the estimated hessian matrix.
  #           if FALSE standard error will not be estimated.
  # opt_meth: the optimization method to use with 'optim'.  Default is "BFGS" other options are "Nelder-Mead", "SANN".
  #Details:
  # The 'Nonparametric' is uses the so-called Grenander estimator (see Jankowski and Wellner, 2009 
  # for details).  The 'Semiparametric' estimator fits a discrete time Cox model to the unobserbed 
  # survival times, as described in McLain et al. (2014).  WARNING: THIS METHOD CAN BE TIME CONSUMING, 
  # especially when there are many unique T and X values.  When this is the case the 'Piecewise' 
  # methods might be more suitable.  The Piecewise method (also described in McLain et al, 2014) fits  
  # the same discrete Cox model, but assumes the baseline log-hazard is contast over discrete partitions  
  # of the sample space of T.
  # For analyses without covariates or weights, the nonparametric is a fast and efficient method. The
  # Piecewise model is another option if smoothness in the estimated survival function is desired, 
  # for example, if the data has heaping or digit preference (see McLain et al, 2014 for details).
  # For analysis without covariates that have weights, either the semiparametric or piecewise approaches
  # should be used.
  # If censoring is desired, the tau option should be used. The uncensored data will be used to choose 
  # the truncation value for the infinite sum.  As a result, it is important that the uncensored values 
  # are inputted into the program. 'tau' must be of length one (equal censoring across all subjects) or 
  # length n (different censoring values). 
  # The starting values for all beta coefficients are obtained using the coxph function with ties=
  # 'exact'.  The starting values of the 'alpha' coefficients are obtained using the nonparametric
  # Grenander estimator. 
  #Value:
  # The function will return a list with the following:
  # Surv_est:   A function containing the fitted survival function.
  # Coef:       A data frame with the beta coefficients, standard errors and 95% confidence intervals.
  #             Standard errors are obtained using a numerical approximation to the hessian matrix.
  #             Set to 'Null' if no covariates are used.
  # Alpha_res:  A data frame with the alpha values (see equations (2) and (3) in McLain et al. for 
  #             details), standard errors and 95% confidence intervals (for Semiparametric and 
  #             Piecewise models only).  Standard errors are obtained using the delta method with a 
  #             numerical approximation to the hessian matrix. Note that the standard errors of the 
  #             alpha values have not been shown to be reiable and bootstrapping methods really are more 
  #             appropriate. If the hessian is not positive definite NA's will be returned for the 
  #             standard error and confidence intervals (more common with the semiparametric method).
  #References:
  # - Jankowski, H. K., & Wellner, J. A. (2009). Estimation of a discrete monotone distribution. 
  #   Electronic journal of statistics, 3, 1567–1605.
  # - McLain, A. C., Sundaram, R., Thoma, M., Louis, B., & Germaine, M. (2014). Semiparametric modeling 
  #   of grouped current duration data with preferential reporting. Statistics in medicine, 33(23), 
  #   3961-3972.
  
  w <- weights
  if(!all(is.integer(T))){T <- floor(T)}
  if(method == "Nonparametric"){
    if(length(X)>1){warning("Nonparametric analysis cannot incorporate covariates. They will be ignored.")}
    if(length(w)>1){warning("Nonparametric analysis cannot incorporate sampling weights They will be ignored.")}
    CEN  <- 1*I(T<=tau)
    T[T>tau] <- tau
    NP_anal <- find.gren_cont(T,CEN, plot=FALSE)
    GREN <- NP_anal/NP_anal[1]
    step_est <- stepfun(c(0:(length(GREN)-1)),c(0,GREN),right=FALSE)
    beta_res <- NULL
    alpha_res<- NULL
  }  
  if(method == "Semiparametric"){
    n    <- length(T)
    mx_T <- max(T)
    CEN  <- 1*I(T<=tau)
    T[T>tau] <- tau
    if(length(w)==1){w <- rep(1,length(T))}
    maxT <- max(c(T[CEN==0],0))
    if(maxT >= max(T[CEN==1])){
      T_vals<- c(sort(unique(T[CEN==1])),max(T[CEN==1])+1)
    }
    if(maxT <  max(T[CEN==1])){
      T_vals<- sort(unique(T[CEN==1]))
    }
    NP_anal <- find.gren_cont(T,CEN)
    step_est_NP <- stepfun(c(0:(length(NP_anal)-1)),c(0,NP_anal),right=FALSE)
    NP_surv <- step_est_NP(T_vals)/step_est_NP(0)
    NP_CH   <- -log(NP_surv)
    st_alpha <- log(c(NP_CH[1],NP_CH[-1] - NP_CH[-length(NP_CH)]))
    len_p <- length(T_vals)
    if(length(T[T==0])>0){
      len_p <- length(T_vals)-1
      st_alpha <- st_alpha[-length(st_alpha)]
    }
    st_alpha[st_alpha<(-5)] <- -5
    cov <- TRUE
    if(length(X)>1){
      if(is.null(dim(X))){stop("X must be a matrix, even if one dimensional.")}
      s_temp <- coxph(Surv(T,CEN)~X,ties = "breslow",weights = w)
      beta <- coef(s_temp)
      unq_X <- max(apply(X,2,function(X){length(unique(X))}))
    }
    if(length(X)==1){
      beta <- 0
      unq_X <- 1
      cov <- FALSE
    }
    st_vals <- c(st_alpha,beta)
    if(Asymp_SE){
      if(unq_X <= (n/2)){
        fit.res <- optim(st_vals,SP_BR_GC_like,X=X,T=T,CEN=CEN,w=w,max_T=mx_T,method=opt_meth,hessian=TRUE,control = list(maxit=10000))
      }
      if(unq_X > (n/2)){
        fit.res <- optim(st_vals,SP_BR_GC_like_cont,X=X,T=T,CEN=CEN,w=w,max_T=mx_T,method=opt_meth,hessian=TRUE,control = list(maxit=10000))
      }
      var_est <- diag(ginv(fit.res$hessian))
      beta_var<- var_est[(len_p+1):length(st_vals)]
      log_alpha_var<- var_est[1:len_p]
      log_alpha_var[log_alpha_var<0] <- 0
      
      par_est <- fit.res$par
      alpha.h <- exp(par_est[1:len_p])
      alpha.se<- alpha.h*sqrt(log_alpha_var)
      beta.h  <- par_est[(len_p+1):length(st_vals)]
      T_seq <- 1:max(T_vals) 
      ful_alpha <- T_seq*0
      ful_alpha[T_vals[T_vals>0]] <- alpha.h
      beta_res <- NULL
      if(cov){beta_res <- data.frame(beta.hat = beta.h,beta.se = sqrt(beta_var), beta.95.CI.L = beta.h - 1.96*sqrt(beta_var), beta.95.CI.U = beta.h + 1.96*sqrt(beta_var))}
      alpha_res <- data.frame(alpha.hat = alpha.h,alpha.se = alpha.se, alpha.95.CI.L = alpha.h - 1.96*alpha.se, alpha.95.CI.U = alpha.h + 1.96*alpha.se)
      if(any(log_alpha_var==0)) alpha_res[log_alpha_var==0,2:4] <- NA
    }
    
    if(!Asymp_SE){
      if(unq_X <= (n/2)){
        fit.res <- optim(st_vals,SP_BR_GC_like,X=X,T=T,CEN=CEN,w=w,max_T=mx_T,method=opt_meth,hessian=FALSE,control = list(maxit=10000))
      }
      if(unq_X > (n/2)){
        fit.res <- optim(st_vals,SP_BR_GC_like_cont,X=X,T=T,CEN=CEN,w=w,max_T=mx_T,method=opt_meth,hessian=FALSE,control = list(maxit=10000))
      }
      
      par_est <- fit.res$par
      alpha.h <- exp(par_est[1:len_p])
      beta.h  <- par_est[(len_p+1):length(st_vals)]
      T_seq <- 1:max(T_vals) 
      ful_alpha <- T_seq*0
      ful_alpha[T_vals[T_vals>0]] <- alpha.h
      beta_res <- NULL
      if(cov){beta_res <- data.frame(beta.hat = beta.h)}
      alpha_res <- data.frame(alpha.hat = alpha.h)
    }
    
    gamma       <- cumsum(ful_alpha)
    temp_vals<- c(1,exp(-gamma[-length(gamma)]),exp(-(gamma[length(gamma)-1]+c(1:ceiling(mx_T*2))*alpha.h[length(alpha.h)])))
    f_est    <- temp_vals/sum(temp_vals)
    surv_est <- f_est/f_est[1]    
    step_est <- stepfun(c(0:length(temp_vals)),c(1,surv_est,0),right=FALSE)
  }  
  if(method == "Piecewise"){
    if(is.null(knots)){stop("Knots must be supplied for the Piecewise model.")}
    if(!all(is.integer(knots))){knots <- unique(ceiling(knots))}
    knots <- knots[knots<tau]
    if(length(w)==1){w <- rep(1,length(T))}
    n    <- length(T)
    mx_T <- max(T)
    CEN  <- 1*I(T<=tau)
    T[T>tau] <- tau
    
    NP_anal <- find.gren_cont(T,CEN)
    step_est_NP <- stepfun(c(0:(length(NP_anal)-1)),c(0,NP_anal),right=FALSE)
    NP_surv <- step_est_NP(0:max(T))/step_est_NP(0)
    NP_CH   <- -log(NP_surv)
    temp_alpha <- c(NP_CH[1],NP_CH[-1] - NP_CH[-length(NP_CH)])
    st_alpha<- NULL
    temp_vec <- c(1,knots,max(T))
    for(j in 2:length(temp_vec)){
      al_j <- temp_alpha[temp_vec[j-1]:temp_vec[j]] 
      st_alpha_j <- log(mean(al_j))
      st_alpha <- c(st_alpha,st_alpha_j)
    }
    st_alpha[st_alpha<(-5) | is.nan(st_alpha) | is.na(st_alpha)] <- -5
    cov <- TRUE
    if(length(X)>1){
      if(is.null(dim(X))){stop("X must be a matrix, even if one dimensional.")}
      s_temp <- coxph(Surv(T,CEN)~X,ties = "breslow",weights = w)
      beta <- coef(s_temp)
      unq_X <- max(apply(X,2,function(X){length(unique(X))}))
    }
    if(length(X)==1){
      beta <- 0
      unq_X <- 1
      cov <- FALSE
    }
    st_vals <- c(st_alpha,beta)
    if(Asymp_SE){
      if(unq_X <= (n/2)){
        fit.res<-optim(st_vals,SP_BR_Piece_like,w=w,X=X,T=T,CEN=CEN,knots=knots,max_T=mx_T,method=opt_meth,hessian=TRUE,control = list(maxit=10000))
      }
      if(unq_X > (n/2)){
        fit.res<-optim(st_vals,SP_BR_Piece_like_cont,w=w,X=X,T=T,CEN=CEN,knots=knots,max_T=mx_T,method=opt_meth,hessian=TRUE,control = list(maxit=10000))
      }
      var_est   <- diag(ginv(fit.res$hessian))
      par_est   <- fit.res$par
      T_vals    <- sort(unique(T))
      len_p     <- length(st_alpha)
      T_seq     <- 1:max(T_vals) 
      alpha.h   <- exp(par_est[1:len_p])
      log_alpha_var<- var_est[1:len_p]
      log_alpha_var[log_alpha_var<0] <- 0
      alpha.se  <- alpha.h*sqrt(log_alpha_var)
      beta.h    <- par_est[(len_p + 1):length(var_est)]
      beta.var  <- var_est[(len_p + 1):length(par_est)]
      beta_res <- NULL
      if(cov){beta_res <- data.frame(beta.hat = beta.h,beta.se = sqrt(beta.var), beta.95.CI.L = beta.h - 1.96*sqrt(beta.var), beta.95.CI.U = beta.h + 1.96*sqrt(beta.var))}
      alpha_res <- data.frame(alpha.hat = alpha.h,alpha.se = alpha.se, alpha.95.CI.L = alpha.h - 1.96*alpha.se, alpha.95.CI.U = alpha.h + 1.96*alpha.se)
      if(any(log_alpha_var==0)) alpha_res[log_alpha_var==0,2:4] <- NA
    }
    if(!Asymp_SE){
      if(unq_X <= (n/2)){
        fit.res<-optim(st_vals,SP_BR_Piece_like,w=w,X=X,T=T,CEN=CEN,knots=knots,max_T=mx_T,method=opt_meth,hessian=FALSE,control = list(maxit=10000))
      }
      if(unq_X > (n/2)){
        fit.res<-optim(st_vals,SP_BR_Piece_like_cont,w=w,X=X,T=T,CEN=CEN,knots=knots,max_T=mx_T,method=opt_meth,hessian=FALSE,control = list(maxit=10000))
      }
      par_est   <- fit.res$par
      T_vals    <- sort(unique(T))
      len_p     <- length(st_alpha)
      T_seq     <- 1:max(T_vals) 
      alpha.h   <- exp(par_est[1:len_p])
      beta.h    <- par_est[(len_p + 1):length(par_est)]
      beta_res <- NULL
      if(cov){beta_res <- data.frame(beta.hat = beta.h)}
      alpha_res <- data.frame(alpha.hat = alpha.h)
    }
    alph_vec  <- numeric(0)
    knots2    <- c(0,knots)
    if(max(T)>=max(knots2)){knots2 <- c(knots2,max(T))}
    if(max(T)<max(knots2)){knots2 <- c(knots2)}
    for(k in 1:len_p){
      alph_vec <- c(alph_vec,rep(alpha.h[k],(knots2[k+1]-knots2[k])))
    }
    ful_alpha <- alph_vec
    gamma    <- cumsum(ful_alpha)
    temp     <- c(1,exp(-gamma),exp(-(gamma[length(gamma)]+c(1:ceiling(mx_T*2))*alph_vec[length(alph_vec)])))
    i_const  <- sum(temp)
    f_est    <- c(temp/i_const)
    surv_est <- f_est/f_est[1]
    step_est <- stepfun(c(0:length(temp)),c(1,surv_est,0),right=FALSE)
  }  
  if(!is.null(beta_res)) rownames(beta_res) <- c(paste("Beta",1:length(beta)))
  if(!is.null(alpha_res)) rownames(alpha_res) <- c(paste("Alpha",1:length(alpha.h)))
  return(list(Surv_est = step_est,Coef = beta_res, Alpha_res = alpha_res))
}










#### Nonparametric CD functions ####
# The following can be used to find the MLE of a decreasing mass function on {0,1,2, ...}.  The input of the 
# function is the vector $\widehat p_n$ and the output is $\gren(\widehat p_n)$. If \texttt{plot=TRUE}, then a # visual representation of the LCM is also given.  See the following for details:
# Jankowski, H. K., & Wellner, J. A. (2009). Estimation of a discrete monotone distribution. Electronic
# journal of statistics, 3, 1567.

find.gren  <-  function(p, plot=FALSE){
  n      <-  length(p)
  pts    <-  cbind(c(0:n,n),c(0,cumsum(p),0))
  
  hpts   <-  chull(pts)
  hpts   <-  c(hpts, hpts[1])
  hpts   <-  sort(hpts)
  
  if(plot==TRUE){
    plot(pts, cex = 1, pch=19, ylab="", xlab="", xlim=c(0,n), ylim=c(0, 1))
    lines(pts[hpts, ])
    points(pts[hpts, ], pch=21, cex=2)
  }
  
  hpairs  <-  matrix(0,length(hpts)-1,2)
  for(i in 1:(length(hpts)-1)){
    hpairs[i,]	<-	c(hpts[i+1],hpts[i])
  }
  m       <-  length(pts[,1])
  hpairs  <-  hpairs[which(hpairs[,1]!=m),]
  hpairs  <-  hpairs[which(hpairs[,2]!=m),]
  if(length(hpairs)==2){
    hpairs<-	matrix(hpairs,1,2)
  } else 	{
    hpairs	<-	hpairs[order(hpairs[,1]),]
  }
  
  s.num   <-  pts[hpairs[,1],2]-pts[hpairs[,2],2]
  s.denom <-  pts[hpairs[,1],1]-pts[hpairs[,2],1]				
  slopes  <-  s.num/s.denom
  h.hat   <-  numeric()
  for (i in 1:length(slopes)){
    h.hat	<-	c(h.hat, rep(slopes[i], hpairs[i,1]-hpairs[i,2]))
  }
  
  return(h.hat)
}


find.gren_cont  <-  function(T,CEN, plot=FALSE){
  t_vals <- 0:max(T)
  s_temp <- survfit(Surv(T,CEN)~1)
  F_step <- stepfun(c(s_temp$time),c(0,1-s_temp$surv),right=TRUE)
  
  pts   <-  cbind(c(t_vals,max(t_vals)+1,max(t_vals)+1),c(F_step(t_vals),1,0))
  
  hpts   <-  chull(pts)
  hpts   <-  c(hpts, hpts[1])
  
  if(plot==TRUE){
    plot(pts, cex = 1, pch=19, ylab="", xlab="", xlim=c(0,36), ylim=c(0, 1))
    lines(pts[hpts, ])
    points(pts[hpts, ], pch=21, cex=2)
  }
  
  hpairs  <-  matrix(0,length(hpts)-1,2)
  for(i in 1:(length(hpts)-1)){
    hpairs[i,]	<-	c(hpts[i+1],hpts[i])
  }
  m       <-  length(pts[,1])
  hpairs  <-  hpairs[which(hpairs[,1]!=m),]
  hpairs  <-  hpairs[which(hpairs[,2]!=m),]
  if(length(hpairs)==2){
    hpairs<-	matrix(hpairs,1,2)
  } else 	{
    hpairs	<-	hpairs[order(hpairs[,1]),]
  }
  
  s.num   <-  pts[hpairs[,1],2]-pts[hpairs[,2],2]
  s.denom <-  pts[hpairs[,1],1]-pts[hpairs[,2],1]				
  slopes  <-  s.num/s.denom
  h.hat   <-  numeric()
  for (i in 1:length(slopes)){
    h.hat	<-	c(h.hat, rep(slopes[i], hpairs[i,1]-hpairs[i,2]))
  }
  
  return(h.hat)
}

find.gren_cont2  <-  function(T,CEN, plot=FALSE){
  t_vals <- c(0,sort(unique(T)))
  s_temp <- survfit(Surv(T,CEN)~1)
  F_step <- stepfun(c(s_temp$time),c(0,1-s_temp$surv),right=TRUE)
  
  pts   <-  cbind(c(t_vals,max(t_vals)+1,max(t_vals)+1),c(F_step(t_vals),1,0))
  
  hpts   <-  chull(pts)
  hpts   <-  c(hpts, hpts[1])
  
  if(plot==TRUE){
    plot(pts, cex = 1, pch=19, ylab="", xlab="", xlim=c(0,36), ylim=c(0, 1))
    lines(pts[hpts, ])
    points(pts[hpts, ], pch=21, cex=2)
  }
  
  hpairs  <-  matrix(0,length(hpts)-1,2)
  for(i in 1:(length(hpts)-1)){
    hpairs[i,]	<-	c(hpts[i+1],hpts[i])
  }
  m       <-  length(pts[,1])
  hpairs  <-  hpairs[which(hpairs[,1]!=m),]
  hpairs  <-  hpairs[which(hpairs[,2]!=m),]
  if(length(hpairs)==2){
    hpairs<-	matrix(hpairs,1,2)
  } else 	{
    hpairs	<-	hpairs[order(hpairs[,1]),]
  }
  
  s.num   <-  pts[hpairs[,1],2]-pts[hpairs[,2],2]
  s.denom <-  pts[hpairs[,1],1]-pts[hpairs[,2],1]				
  slopes  <-  s.num/s.denom
  h.hat   <-  numeric()
  for (i in 1:length(slopes)){
    h.hat	<-	c(h.hat, rep(slopes[i], hpairs[i,1]-hpairs[i,2]))
  }
  
  return(h.hat)
}


########################################Semi-parametric##################################

SP_BR_GC_like <- function(par,X,T,CEN=1,w=1,max_T=400){
  
  #Semi-parametric backward recurrent discrete cox likelihood.
  
  if(length(CEN)==1){CEN <- rep(1,length(T))}
  maxT <- max(c(T[CEN==0],0))
  if(maxT >= max(T[CEN==1])){
    T_vals<- c(sort(unique(T[CEN==1])),max(T[CEN==1])+1)
    check <- TRUE
  }
  if(maxT < max(T[CEN==1])){
    T_vals<- sort(unique(T[CEN==1]))
    check <- FALSE
  }
  
  if(length(w)==1){w <- rep(1,length(T))}
  len_p <- length(T_vals)-1
  if(length(T_vals[T_vals==0])==0){len_p <- length(T_vals)}
  T_seq <- 1:max(T_vals) 
  
  alpha <- exp(par[1:len_p])
  beta  <- par[(len_p + 1):length(par)]
  if(length(beta)==1){eta <- X*beta}
  if(length(beta) >1){eta <- X%*%beta}
  
  ful_alpha <- T_seq*0
  ful_alpha[T_vals[T_vals>0]] <- alpha
  
  gamma       <- cumsum(ful_alpha)
  GAM_T       <- T
  GAM_T[T> 0] <- gamma[T[T>0]]
  
  i_const <- rep(0,length(T))
  i_surv  <- rep(0,length(T))
  ex_param = Inf
  if(check){
    ex_param <- alpha[length(alpha)]
  }
  for(i in unique(eta)){
    temp_vals<- c(1,exp(-gamma*exp(i)),exp(-(gamma[length(gamma)]+c(1:ceiling(max_T*2))*ex_param)*exp(i)))
    t_pmf    <- temp_vals/sum(temp_vals)
    t_T_vals <- T[eta==i]+1
    cs_pmf   <- cumsum(t_pmf)
    t_surv   <- 1 - cs_pmf[t_T_vals]
    i_surv[eta==i]   <- t_surv
    i_const[eta==i]  <- sum(temp_vals)
  }
  
  
  like       <- (CEN)*exp(-GAM_T*exp(eta))/i_const + (1-CEN)*i_surv
  llike      <- like
  llike[is.nan(like)] <- 0
  llike[is.na(like)]  <- 0
  llike[llike==0] <- -90000
  llike[llike>0] <- w[llike>0]*log(like[llike>0])
  -sum(llike)
}




SP_BR_GC_like_cont <- function(par,X,T,CEN=1,w=1,max_T=400){
  
  #Semi-parametric backward recurrent discrete cox likelihood.
  
  if(length(CEN)==1){CEN <- rep(1,length(T))}
  maxT <- max(c(T[CEN==0],0))
  if(maxT >= max(T[CEN==1])){
    T_vals<- c(sort(unique(T[CEN==1])),max(T[CEN==1])+1)
    check <- TRUE
  }
  if(maxT < max(T[CEN==1])){
    T_vals<- sort(unique(T[CEN==1]))
    check <- FALSE
  }
  
  if(length(w)==1){w <- rep(1,length(T))}
  len_p <- length(T_vals)-1
  if(length(T_vals[T_vals==0])==0){len_p <- length(T_vals)}
  T_seq <- 1:max(T_vals) 
  
  alpha <- exp(par[1:len_p])
  beta  <- par[(len_p + 1):length(par)]
  if(length(beta)==1){eta <- X*beta}
  if(length(beta) >1){eta <- X%*%beta}
  
  ful_alpha <- T_seq*0
  ful_alpha[T_vals[T_vals>0]] <- alpha
  
  gamma       <- cumsum(ful_alpha)
  GAM_T       <- T
  GAM_T[T> 0] <- gamma[T[T>0]]
  
  i_const <- numeric(0)
  i_surv  <- numeric(0)
  ex_param = Inf
  if(check){
    ex_param <- alpha[length(alpha)]
  }
  for(i in 1:length(T)){
    temp_vals<- c(1,exp(-gamma*exp(eta[i])),exp(-(gamma[length(gamma)]+c(1:ceiling(max_T*2))*ex_param)*exp(eta[i])))
    t_pmf    <- temp_vals/sum(temp_vals)
    t_surv   <- 1 - sum(t_pmf[1:c(T[i]+1)])
    i_surv   <- c(i_surv,t_surv)
    i_const  <- c(i_const,sum(temp_vals))
  }
  
  
  like       <- (CEN)*exp(-GAM_T*exp(eta))/i_const + (1-CEN)*i_surv
  llike      <- like
  llike[is.nan(like)] <- 0
  llike[is.na(like)]  <- 0
  llike[llike==0] <- -90000
  llike[llike>0] <- w[llike>0]*log(like[llike>0])
  -sum(llike)
}




########################################Peicewise##################################

SP_BR_Piece_like <- function(par,X,T,CEN=1,knots,w=1,max_T=400){
  
  if(length(CEN)==1){CEN <- rep(1,length(T))}	
  if(length(w)==1){w <- rep(1,length(T))}
  len_p <- length(knots)+1
  
  alpha <- exp(par[1:len_p])
  beta <- par[(len_p + 1):length(par)]
  if(length(beta)==1){eta <- X*beta}
  if(length(beta)>1){eta <- X%*%beta}
  
  alph_vec <- numeric(0)
  knots2 <- c(0,knots)
  if(max(T)>=max(knots2)){knots2 <- c(knots2,ceiling(max(T))+4)}
  if(max(T)<max(knots2)){knots2 <- c(knots2)}
  for(k in 1:(len_p)){
    alph_vec <- c(alph_vec,rep(alpha[k],(knots2[k+1]-knots2[k])))
  }
  ful_alpha <- alph_vec
  
  gamma       <- cumsum(ful_alpha)
  GAM_T       <- T
  GAM_T[T> 0] <- gamma[T[T>0]]
  
  i_const <- rep(0,length(T))
  i_surv  <- rep(0,length(T))
  for(i in unique(eta)){
    temp_vals<- c(1,exp(-gamma*exp(i)),exp(-(gamma[length(gamma)]+c(1:ceiling(max_T*2))*alph_vec[length(alph_vec)])*exp(i)))
    t_pmf    <- temp_vals/sum(temp_vals)
    t_T_vals <- as.integer(T)[eta==i]+1
    cs_pmf   <- cumsum(t_pmf)
    t_surv   <- 1 - cs_pmf[t_T_vals]
    i_surv[eta==i]   <- t_surv
    i_const[eta==i]  <- sum(temp_vals)
  }
  
  like       <- (CEN)*exp(-GAM_T*exp(eta))/i_const + (1-CEN)*i_surv
  llike      <- like
  llike[is.nan(like)] = 0
  llike[is.na(like)] = 0
  llike[llike==0] <- -90000
  llike[llike>0] <- w[llike>0]*log(like[llike>0])
  -sum(llike)
}


SP_BR_Piece_like_cont <- function(par,X,T,CEN=1,knots,w=1,max_T=400){
  
  if(length(CEN)==1){CEN <- rep(1,length(T))}	
  if(length(w)==1){w <- rep(1,length(T))}
  len_p <- length(knots)+1
  
  alpha <- exp(par[1:len_p])
  beta <- par[(len_p + 1):length(par)]
  if(length(beta)==1){eta <- X*beta}
  if(length(beta)>1){eta <- X%*%beta}
  
  alph_vec <- numeric(0)
  knots2 <- c(0,knots)
  if(max(T)>=max(knots2)){knots2 <- c(knots2,as.integer(max(T))+1)}
  if(max(T)<max(knots2)){knots2 <- c(knots2)}
  for(k in 1:len_p){
    alph_vec <- c(alph_vec,rep(alpha[k],(knots2[k+1]-knots2[k])))
  }
  ful_alpha <- alph_vec
  
  gamma       <- cumsum(ful_alpha)
  GAM_T       <- T
  GAM_T[T> 0] <- gamma[T[T>0]]
  
  i_const  <- numeric(0)
  i_surv    <- numeric(0)
  for(i in 1:length(T)){  
    temp_vals<- c(1,exp(-gamma*exp(eta[i])),exp(-(gamma[length(gamma)] + c(1:ceiling(max_T*2))*alph_vec[length(alph_vec)])*exp(eta[i])))
    t_pmf    <- temp_vals/sum(temp_vals)
    t_surv   <- 1 - sum(t_pmf[1:c(as.integer(T[i])+1)])
    i_surv   <- c(i_surv,t_surv)
    i_const <- c(i_const,sum(temp_vals))
  }
  
  like       <- (CEN)*exp(-GAM_T*exp(eta))/i_const + (1-CEN)*i_surv
  llike      <- like
  llike[is.nan(like)] = 0
  llike[is.na(like)] = 0
  llike[llike==0] <- -90000
  llike[llike>0] <- w[llike>0]*log(like[llike>0])
  -sum(llike)
}




