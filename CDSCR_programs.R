library(mvtnorm)
library(MASS)
library(statmod)
require(survival)



#### Wrapper to the optimization of the likelihood ####

CDSCR <- function(T,TTFT,fail_ind,knots_T,knots_Y,knots_Z,K1,K2,K3,B=20,theta=NULL,pl=FALSE,opt_meth="Nelder-Mead"){
  
  ##################################################################################################
  #Description: This program will run a current duration semi-competeing risk model with a  
  #Description: This program will run a current duration semi-competeing risk model with a  
  # piecewise constant hazard function 
  #Usage:
  # CDSCR(T,TTFT,fail_ind,knots_T,knots_Y,knots_Z,K1,K2,K3,B,theta,pl,opt_method)
  #
  #Arguments:
  # T:        Total current duration, the time from the begining of attempt to sampling.
  # TTFT:     Time to the intermittent event (TTFT in the paper). If no intermittent event set =0.
  # fail_ind: Indicator that the intermittent event occured.
  # tau:      The type I censoring value (length 1), or a vector (length n) of censoring values (not 
  #           required).
  # knots_T:  Values of the knots for the Piecewise model for the distribution of T. Set to NULL to 
  #           use data driven knots (K1 must be set in this case).
  # knots_Y:  Values of the knots for the Piecewise model for the distribution of Y. Set to NULL to 
  #           use data driven knots (K2 must be set in this case).
  # knots_Z:  Values of the knots for the Piecewise model for the distribution of Z. Set to NULL to 
  #           use data driven knots (K3 must be set in this case).
  # K1:       Number of knots to use for T (ignored if knots_T is given).
  # K2:       Number of knots to use for Y (ignored if knots_Y is given).
  # K3:       Number of knots to use for Z (ignored if knots_Z is given).
  # B:        The number of quadrature nodes to use.  Default is 20.
  # theta:    Parameter of length 2 that gives the starting values for Corr(T,Y) and Corr(T,Z), 
  #           respectively (optional).
  # pl:       Logical: if true the program will show plots of survival function with the starting values.
  # opt_meth: the optimization method to use with 'optim'.  Default is "Nelder-Mead" other options 
  #           are "BFGS", "SANN".
  #Value:
  # The function will return a list with the following:
  # like_opt:   The results from optim.
  # knots:      A list containing the knots used for T, Y, and Z in the optimization.
  # st_vals:    A vector containing the starting values.
  
  ## Selecting knot locations. 
  if(is.null(knots_T)){
    if(is.null(K1)){stop("Either knots or K must be specified.")}
    knots_T <- quantile(T[fail_ind==0],probs = c(1:(K1))/(K1+1))
    knots_T <- unique(knots_T[knots_T>0])
    K1 <- length(knots_T)
  } 
  if(is.null(knots_Y)){
    if(is.null(K2)){stop("Either knots or K must be specified.")}
    knots_Y <- quantile(TTFT[fail_ind==1],probs = c(1:(K2))/(K2+1))
    knots_Y <- unique(knots_Y[knots_Y>0])
    K2 <- length(knots_Y)
  } 
  if(is.null(knots_Z)){
    if(is.null(K3)){stop("Either knots or K must be specified.")}
    knots_Z <- quantile(T[fail_ind==1]-TTFT[fail_ind==1],probs = c(1:(K3))/(K3+1))
    knots_Z <- unique(knots_Z[knots_Z>0])
    K3 <- length(knots_Z)
  }
  
  
  ## Getting the nodes and weights for Legendre quadrature
  out <- gauss.quad(B,"legendre")
  nodes <- out$nodes
  weights <- out$weights
  
  ## Getting starting values for parameters.
  naive_est <- hazard_guess(T[fail_ind==0],knots_T) 
  naive.z   <- hazard_guess(T[fail_ind==1]-TTFT[fail_ind==1],knots_Z)
  naive.z[length(naive.z)] <- 0.2
  naive.y   <- hazard_guess_ecdf(ecdf(TTFT[fail_ind==1]),TTFT[fail_ind==1],knots_Y)
  naive.y   <- exp(log(naive.y)/0.8)
  
  naive.yx  <- hazard_guess_ecdf(ecdf(TTFT[fail_ind==1]),TTFT[fail_ind==1],knots_T)
  naive.yx   <- exp(log(naive.yx)/0.8)
  naive.x <- naive_est - naive.yx
  
  ## Setting minimum value for starting values
  min_val <- exp(-5)
  naive.x[naive.x<min_val] <- min_val
  naive.z[naive.z<min_val] <- min_val
  naive.y[naive.y<min_val] <- min_val
  
  ## Starting values
  if(!is.null(theta)){
    corrxy.T <- max(c(abs(theta[1])/(1-abs(theta[1])),1e-4))
    corrxz.T <- max(c(abs(theta[3])/(1-abs(theta[3])),1e-4))
  }
  if(is.null(theta)){
    corrxy.T <- 0.1/(1-0.1)
    corrxz.T <- 0.1/(1-0.1)
  }
  if(pl){
    t_vec <- seq(0,36,1)
    t_Surv_X <- PC_surv(t_vec,knots_T,naive.x)
    t_Surv_Y <- PC_surv(t_vec,knots_Y,naive.y)
    t_Surv_Z <- PC_surv(t_vec,knots_Z,naive.z)
    
    plot(t_vec,t_Surv_X,col=1,lwd=2,lty=1,type="l",ylab="Starting Survival Probability",xlab="Time",ylim=c(0,1))
    lines(t_vec,t_Surv_Y,col=2,lwd=2,lty=1)
    lines(t_vec,t_Surv_Z,col=3,lwd=2,lty=1)
  }
  
  par<-  log(c(naive.x,naive.y,naive.z,corrxy.T,corrxz.T))
  
  ## Running optimization
  like_opt <-stats::optim(par,likelihood.function.cdscr,T=T,TTF=TTFT,fail_ind=fail_ind,knots_T=knots_T,knots_Y=knots_Y,knots_Z=knots_Z,nodes=nodes,weights=weights,hessian = FALSE,method=opt_meth,control = list(maxit=10000))
  
  res <- list(like_opt = like_opt,knots=list(knots_T=knots_T,knots_Y=knots_Y,knots_Z=knots_Z),st_vals=par)
  return(res)
}



#### Peicewise constant survival function and probability functions   

PC_surv <- function(t,knots,alpha){
  
  ###Increments metween knots
  L_int <- knots - c(0,knots[-length(knots)])
  
  ###Finding the first knot bigger then t 
  kn_mat <- matrix(knots,length(t),length(knots),byrow = TRUE)
  t_kn_mat <- kn_mat-t
  t_kn_mat[t_kn_mat<0] <- max(knots)
  first_kn_bg <- apply(t_kn_mat,1,which.min)
  first_kn_bg[t>max(knots)] <- length(knots)+1
  
  ###Getting the alpha for the interval with t
  last_alpha <- alpha[first_kn_bg] 
  
  ###Getting the increment beyond the last knot
  last_incr <- t - c(0,knots)[first_kn_bg] 
  
  ###Calculating the survival function at the knots.  Adding one for zero too.###
  sur_func_incr <- 1 - pexp(L_int,alpha[-length(alpha)])
  kn_sur_vals <- c(1,cumprod(sur_func_incr))
  
  ###Finding survival value at last knot smaller then t
  surv_last_kn <- kn_sur_vals[first_kn_bg]
  
  ###Adding the last increment###
  surv_val <- surv_last_kn*(1-pexp(last_incr,last_alpha))
  surv_val
}



PC_pdf <- function(t,knots,alpha){
  ### Increments metween knots
  L_int <- knots - c(0,knots[-length(knots)])
  
  ### Finding the first knot bigger then t 
  kn_mat <- matrix(knots,length(t),length(knots),byrow = TRUE)
  t_kn_mat <- kn_mat-t
  t_kn_mat[t_kn_mat<0] <- max(knots)
  first_kn_bg <- apply(t_kn_mat,1,which.min)
  first_kn_bg[t>max(knots)] <- length(knots)+1
  
  ### Getting the alpha for the interval with t
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


#### The survival function of min(T,Y)

surv.min.T.Y<-function(x,corrxy,knots_T,alpha_T,knots_Y,alpha_Y){
  F.T<-PC_surv(x,knots_T,alpha_T)
  F.Y<-PC_surv(x,knots_Y,alpha_Y)
  u.Func<- qnorm(1-F.T)
  v.Func<- qnorm(1-F.Y)
  uv_mat <- cbind(-u.Func,-v.Func)
  test <- apply(uv_mat,1,pmvnorm,upper=Inf,mean = c(0,0),corr = matrix(c(1,corrxy,corrxy,1),2,2))
  C_bar <- F.T + F.Y -1 + test
  C_bar[C_bar<0]<-0
  C_bar[is.na(C_bar)]<-0
  return(C_bar)
}

#### The product of the conditional survival function of T|Y and the pdf of Y

cond_surv_dens_0_1 <- function(t,knots_T,alpha_T,knots_Y,alpha_Y,fCop){
  ans <- cCopula(cbind(PC_surv(t,knots_Y,alpha_Y),PC_surv(t,knots_T,alpha_T)), fCop, indices = 2)*PC_pdf(t,knots_Y,alpha_Y)
  ans[is.na(ans)] <- 0
  ans
}


####   The density function of T|Y

density.T.given.Y <- function(x,y,corrxy,knots_T,alpha_T,knots_Y,alpha_Y){
  u = qnorm(1-PC_surv(x,knots_T,alpha_T))
  v = qnorm(1-PC_surv(y,knots_Y,alpha_Y))
  u[!is.finite(u)] <- sign(u[!is.finite(u)])*1e20
  v[!is.finite(v)] <- sign(v[!is.finite(v)])*1e16
  ans1 <- dnorm(u,corrxy*v,sqrt(1-corrxy^2))*PC_pdf(x,knots_T,alpha_T)/dnorm(u)
  ans1[is.na(ans1)] <- 0
  ans1
}

density.T.given.Y.grouped <- function(x,y,corrxy,knots_T,alpha_T,knots_Y,alpha_Y){
  joint <- density.T.Y.grouped(x,y,corrxy,knots_T,alpha_T,knots_Y,alpha_Y)
  marg <- PC_surv(y,knots_Y,alpha_Y) - PC_surv(y+1,knots_Y,alpha_Y)
  joint/marg
}



#### The joint density function of T and Y

density.T.Y <- function(x,y,corrxy,knots_T,alpha_T,knots_Y,alpha_Y){
  F_y1 <- cdf.Y.given.T_2args(y,x,corrxy,knots_T,alpha_T,knots_Y,alpha_Y)
  F_y2 <- cdf.Y.given.T_2args(y+1,x,corrxy,knots_T,alpha_T,knots_Y,alpha_Y)
  ans1 <- (F_y2 - F_y1)*PC_pdf(x,knots_T,alpha_T)
  ans1[is.na(ans1)] <- 0
  ans1
}



#### The survival function of Z|T

surv.Z.given.T <- function(z,x,corrxz,knots_T,alpha_T,knots_Z,alpha_Z){
  F_bar_z <- PC_surv(z,knots_Z,alpha_Z)
  F_bar_x <- PC_surv(x,knots_T,alpha_T)
  u = qnorm(1-F_bar_z)
  v = qnorm(1-F_bar_x)
  u[!is.finite(u)] <- sign(u[!is.finite(u)])*1e20
  v[!is.finite(v)] <- sign(v[!is.finite(v)])*1e16
  ans  <- (pnorm(u,mean = corrxz*v,sd = sqrt(1-corrxz^2),lower.tail = FALSE))
  ans
}



#### The distribution function of Y|T

cdf.Y.given.T <- function(x,corrxy,knots_T,alpha_T,knots_Y,alpha_Y){
  F_bar_y <- PC_surv(x,knots_Y,alpha_Y)
  F_bar_x <- PC_surv(x,knots_T,alpha_T)
  u = qnorm(1-F_bar_y)
  v = qnorm(1-F_bar_x)
  u[!is.finite(u)] <- sign(u[!is.finite(u)])*1e16
  v[!is.finite(v)] <- sign(v[!is.finite(v)])*1e20
  ans <- (pnorm(u,mean = corrxy*v,sd = sqrt(1-corrxy^2),lower.tail = TRUE))
  ans
}


cdf.Y.given.T_2args <- function(y,x,corrxy,knots_X,alpha_T,knots_Y,alpha_Y){
  F_bar_y <- PC_surv(y,knots_Y,alpha_Y)
  F_bar_x <- PC_surv(x,knots_X,alpha_T)
  u = qnorm(1-F_bar_y)
  v = qnorm(1-F_bar_x)
  v[!is.finite(v)] <- 1e10*sign(v[!is.finite(v)])
  ans <- (pnorm(u,mean = corrxy*v,sd = sqrt(1-corrxy^2),lower.tail = TRUE))
  ans
}

#### Expectation of min T and Y

expect.min.T.Y<-function(r,knots_T,alpha_T,knots_Y,alpha_Y){
  integrate(surv.min.T.Y,0,Inf,r,knots_T,alpha_T,knots_Y,alpha_Y)$value
}


##### The function is to calculate the inner part of the integration in the denominator

expect.z.given.t<-function(x,r,knots_T,alpha_T,knots_Z,alpha_Z){
  integrate(surv.Z.given.T,0,Inf,x,r,knots_T,alpha_T,knots_Z,alpha_Z)$value
}


denom.integral.inner<-function(x,corrxy,corrxz,knots_T,alpha_T,knots_Y,alpha_Y,knots_Z,alpha_Z)
{
  sapply(x,function(x,corrxz,knots_T,alpha_T,knots_Z,alpha_Z){integrate(surv.Z.given.T,0,Inf,x,corrxz,knots_T,alpha_T,knots_Z,alpha_Z)$value},corrxz,knots_T,alpha_T,knots_Z,alpha_Z)*cdf.Y.given.T(x,corrxy,knots_T,alpha_T,knots_Y,alpha_Y)*PC_pdf(x,knots_T,alpha_T)
}


denom.integral.inner.legendre<-function(x,corrxy,corrxz,knots_T,alpha_T,knots_Y,alpha_Y,knots_Z,alpha_Z,nodes,weights)
{
  L <- length(x)
  n <- length(nodes)
  trans_nodes <- (nodes+1)/(1-nodes)
  dp <- 1/(1-nodes)*(1+trans_nodes)
  node_vec    <- rep(trans_nodes,L)
  x_vec <- rep(x,each=n)
  inner_mat <- matrix(surv.Z.given.T(node_vec,x_vec,corrxz,knots_T,alpha_T,knots_Z,alpha_Z),nrow=L,ncol = n,byrow = TRUE)
  t1 <- t(t(weights*dp)%*%t(inner_mat))*cdf.Y.given.T(x,corrxy,knots_T,alpha_T,knots_Y,alpha_Y)*PC_pdf(x,knots_T,alpha_T)
  t1 
}


#### The function is to calculate the inner part of the integration in the numerator


nume.integral.inner.legendre<-function(temp.data,nodes,corrxy,corrxz,knots_T,alpha_T,knots_Y,alpha_Y,knots_Z,alpha_Z){
  N <- dim(temp.data)[1]
  trans_nodes <- (nodes+1)/(1-nodes)
  dp <- 1/(1-nodes)*(1+trans_nodes)
  x <- rep(trans_nodes,N)+rep(temp.data[,1],each=length(nodes))
  y <- rep(temp.data[,1],each=length(nodes))
  z <- rep(temp.data[,2],each=length(nodes))
  ans<-surv.Z.given.T(z,x,corrxz,knots_T,alpha_T,knots_Z,alpha_Z)*density.T.Y(x,y,corrxy,knots_T,alpha_T,knots_Y,alpha_Y)*dp
  ans_mat <- matrix(ans,nrow=N,ncol=length(nodes),byrow = TRUE)
  return(ans_mat)
}


#### CDSCR Likelihood ####

likelihood.function.cdscr <-function(par,T,TTF,fail_ind,knots_T,knots_Y,knots_Z,nodes,weights){
  
  Ppar <<- par
  par <- exp(par)
  nsll <- Inf
  
  ### If any parameters are not finite, function returns Inf
  if(all(is.finite(par))){
    
    
    alpha_T  <- (par[1:(length(knots_T)+1)])
    alpha_Y  <- (par[(length(knots_T)+2):(length(knots_T)+length(knots_Y)+2)])
    alpha_Z  <- (par[(length(knots_T)+length(knots_Y)+3):(length(knots_T)+length(knots_Y)+length(knots_Z)+3)])
    corrxy <- -1
    corrxz <- 1
    if(is.finite(par[length(par)-1])) corrxy <- -(par[length(par)-1]/(1+par[length(par)-1]))
    if(is.finite(par[length(par)])) corrxz <- (par[length(par)]/(1+par[length(par)]))
    
    N <- length(T)
    ind_vec <- 1:N
    
    ### Calculating the numerator 
    numerator1<-surv.min.T.Y(T,corrxy,knots_T,alpha_T,knots_Y,alpha_Y)
    
    temp.data<-cbind(TTF,c(T-TTF))
    ans_mat <- nume.integral.inner.legendre(temp.data,nodes,corrxy=corrxy,corrxz=corrxz,knots_T=knots_T,alpha_T=alpha_T,knots_Y=knots_Y,alpha_Y=alpha_Y,knots_Z=knots_Z,alpha_Z=alpha_Z)
    numerator2 <- t(t(weights)%*%t(ans_mat))
    
    numerator <- numerator1*(1-fail_ind) + fail_ind*numerator2
    
    ### Calculating the expectation of min(T,Y)
    t1 <- try(ETminY <- expect.min.T.Y(r=corrxy,knots_T=knots_T,alpha_T=alpha_T,knots_Y=knots_Y,alpha_Y=alpha_Y),silent = TRUE)
    if(!is.null(attr(t1,"class"))){ETminY <- Inf}
    
    ### Calculating the expectation of Z.
    t1 <- try(EZgiven <- integrate(denom.integral.inner.legendre,0,Inf,corrxy,corrxz,knots_T=knots_T,alpha_T=alpha_T,knots_Y=knots_Y,alpha_Y=alpha_Y,knots_Z=knots_Z,alpha_Z=alpha_Z,nodes,weights)$value,silent = TRUE)
    ### If error, calculating using Monte Carlo. Note: the optimization is unlikely to be sucessful 
    #   if Monte Carlo is required at the MLE.
    if(!is.null(attr(t1,"class"))){
      cat("Int approx \n")
      Sigma <- diag(3)
      Sigma[1,2] <- Sigma[2,1] <- corrxy
      Sigma[1,3] <- Sigma[3,1] <- corrxz
      Sigma[2,3] <- Sigma[3,2] <- corrxy*corrxz
      u <- mvrnorm(n=1e7, mu=c(0,0,0), Sigma=Sigma)
      Y_vals <- gen_func(pnorm(u[,2]),knots_Y,alpha_Y)
      Z_vals <- gen_func(pnorm(u[,3]),knots_Z,alpha_Z)
      EZgiven <- mean(Z_vals*I(Y_vals<X_vals))
    }
    denominator <- Inf
    if(is.finite(ETminY) & is.finite(EZgiven)){denominator <- ETminY + EZgiven}
    
    likelihood1<-numerator/denominator
    log.likelihood1<-log(likelihood1)
    nsll <- sum(-log.likelihood1)
  }
  return(nsll)
}



##### Functions to get starting values from nonpaametric CD estimator.

hazard_guess <- function(surv_dat,knots){
  n <- length(sort(unique(surv_dat)))
  test.z <- find.gren_cont2(surv_dat,CEN=rep(1,length(surv_dat)), plot=FALSE)
  GREN <- test.z/test.z[1]
  step_est.z <- stepfun(c(0,sort(unique(surv_dat))),c(0,GREN),right=FALSE)
  mx_z <- max(surv_dat)
  z_surv <- step_est.z(c(knots,mx_z))
  z_cumhaz <- -log(z_surv)
  z_haz <- (z_cumhaz - c(0,z_cumhaz[-length(z_cumhaz)]))/c(c(knots,mx_z) - c(0,knots))
  return(z_haz)
}

##### Functions to get starting values from the ECDF.

hazard_guess_ecdf <- function(ecdf,surv_dat,knots){
  
  mx_z <- max(surv_dat)*1.3
  y_surv <- c(1-ecdf(knots),1/length(surv_dat))
  y_cumhaz <- -log(y_surv)
  y_haz <- (y_cumhaz - c(0,y_cumhaz[-length(y_cumhaz)]))/c(c(knots,mx_z) - c(0,knots))
  return(y_haz)
}


#### Function for generating PC data.

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










# R code to find the MLE of a decreasing mass function on {0,1,2, ...}.  The input of the function  is the vector 
# $\widehat p_n$ and the output is $\gren(\widehat p_n)$.  If \texttt{plot=TRUE}, then a visual representation of 
# the LCM is also given.

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



SP_BR_Piece_like <- function(par,X,T,CEN=1,knots,w=1,max_T=400){
  
  #Semi-parametric backward recurrent discrete cox likelihood.
  #The baseline is model using a piecewise constant specification.
  
  
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


########################################Semi-parametric##################################

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



SP_BR_Piece_like_cont <- function(par,X,T,CEN=1,knots,w=1,max_T=400){
  
  #Semi-parametric backward recurrent discrete cox likelihood.
  #The baseline is model using a piecewise constant specification.
  
  
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








