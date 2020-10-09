
#### Load the programs
source("CDSCR_programs.R")

ex_data <- read.csv("CDSCR_data.csv")

#### Extracting data
# Total current duration time
T <- ex_data$T 
# Time to fertility treatment (can be set to zero if fertility treatment did not occur)
TTFT <- ex_data$TTF 
# Indicator that fertility treatment occured.
fail_ind <- ex_data$fail_ind 

#Set the number of knots, locations will be choosen within function.
K1 <- 8
K2 <- 8
K3 <- 8

#### Run Maximum Likelihood Optimization
CDSCR_opt <- CDSCR(T,TTFT,fail_ind,knots_T=NULL,knots_Y=NULL,knots_Z=NULL,K1=K1,K2=K2,K3=K2,B=20,theta=NULL,opt_meth = "BFGS")

### Extract knots
knots_T <- CDSCR_opt$knots$knots_T  
knots_Y <- CDSCR_opt$knots$knots_Y  
knots_Z <- CDSCR_opt$knots$knots_Z  

### Likelihood value and optimization details.
CDSCR_opt$like_opt$value
CDSCR_opt$like_opt$counts
CDSCR_opt$like_opt$convergence

### Extract parameter estimates
par_est <- exp(CDSCR_opt$like_opt$par)
alpha_T_est  <- (par_est[1:(length(knots_T)+1)])
alpha_Y_est  <- (par_est[(length(knots_T)+2):(length(knots_T)+length(knots_Y)+2)])
alpha_Z_est  <- (par_est[(length(knots_T)+length(knots_Y)+3):(length(knots_T)+length(knots_Y)+length(knots_Z)+3)])
corrxy_est <- -(par_est[length(par_est)-1]/(1+par_est[length(par_est)-1]))
corrxz_est <- (par_est[length(par_est)]/(1+par_est[length(par_est)]))


### Calculating the survival functions of T, Y, Z and min(T,Y)
t_vec <- seq(0,36)
Surv_T_est <- PC_surv(t_vec,knots_T,alpha_T_est)
Surv_Y_est <- PC_surv(t_vec,knots_Y,alpha_Y_est)
Surv_Z_est <- PC_surv(t_vec,knots_Z,alpha_Z_est)
Min_TY_est <- surv.min.T.Y(t_vec,corrxy_est,knots_T,alpha_T_est,knots_Y,alpha_Y_est)


### Estimating the survival using the standard (McLain et al. 2014) approach

Surv_T_naive <- CD_surv(T,method = 'Piecewise', knots = knots_T)$Surv_est

#### Plotting the results 

par(mar=c(5,5 , 2, 2) + 0.1)
plot(NULL,xlim=c(0,max(t_vec)),ylim=c(0,1),axes=FALSE, ann=FALSE)

lines(t_vec,Surv_T_est,col=1,lwd=3,lty=1)
lines(t_vec,Surv_Y_est,col=2,lwd=3,lty=1)
lines(t_vec,Surv_Z_est,col=3,lwd=3,lty=1)
lines(t_vec,Min_TY_est,col=4,lwd=3,lty=1)
lines(t_vec,Surv_T_naive(t_vec),col=1,lwd=3,lty=2)


axis(2,las=1, at=seq(0,1,0.2), cex.axis=1.5)
tsq <- seq(0.1,0.9,0.1)
axis(2,las=1, at=tsq,lab=c(rep("",length(tsq))), cex.axis=1.5)

axis(1,las=1, at=seq(0,max(t_vec),6), cex.axis=1.5)
tsq <- seq(0,max(t_vec),3)
axis(1,las=1, at=tsq,lab=c(rep("",length(tsq))), cex.axis=1.5)
title(ylab="Survival Function",cex.lab=2,mgp=c(3,2,0))
title(xlab="Time (in months)",cex.lab=2,mgp=c(3,2,0))
legend("topright",legend = c("Survival of T","Survival of Y","Survival of Z","Survival of min(T,Y)","Survival of T (naive)"),lwd=3,col=c(1,2,3,4,1),lty=c(rep(1,4),2))
box()


