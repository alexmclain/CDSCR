
#### Description
# Gives two examples of the methods.  The first example uses continuous data and the second example uses 
# grouped data. The second example is more time consuming than the first as it uses full adaptive 
# quadrature to calculate the integrals in the denominator. 
# For more information see the description in the "Programs_CDSCR.R" file.



#### Load the programs
source("CDSCR_programs.R")

#### Start of the first example ####
ex_data <- read.csv("CDSCR_data.csv")

#### Extracting data
# Total current duration time
x.data <- ex_data$T 
# Time to intermittent event 
y.data <- ex_data$TTFT 
# Indicator that intermittent event occured.
delta.data <- ex_data$fail_ind 

#Setting the number of peicewise parameters, locations will be choosen within function (knots can be given)
K1 <- 5 #For T
K2 <- 5 #For Y
K3 <- 5 #For Z


#### Run Maximum Likelihood Optimization
CDSCR_opt <- CDSCR(x.data,y.data,delta.data,K1=K1,K2=K2,K3=K3,B=20)

### Extract knots
knots.t <- CDSCR_opt$knots.t
knots.y <- CDSCR_opt$knots.y
knots.z <- CDSCR_opt$knots.z

### Likelihood value and optimization details.
CDSCR_opt$like_opt$value
CDSCR_opt$like_opt$counts
CDSCR_opt$like_opt$convergence

### Extract parameter estimates
par_est <- CDSCR_opt$like_opt$par

#Transform to corresponding values
alpha_T_est <- exp(par_est[1:(length(knots.t)+1)])
alpha_Y_est <- exp(par_est[(length(knots.t)+2):(length(c(knots.t,knots.y))+2)])
alpha_Z_est <- exp(par_est[(length(c(knots.t,knots.y))+3):(length(c(knots.t,knots.y,knots.z))+3)])
corrxy_est <-  -(exp(par_est[length(par_est)-1])/(1+exp(par_est[length(par_est)-1])))
corrxz_est <-  (exp(par_est[length(par_est)])/(1+exp(par_est[length(par_est)])))

  
### Calculating the survival functions of T, Y, Z and min(T,Y)
t_vec <- seq(0,36)
Surv_T_est <- PC_surv(t_vec,knots.t,alpha_T_est)
Surv_Y_est <- PC_surv(t_vec,knots.y,alpha_Y_est)
Surv_Z_est <- PC_surv(t_vec,knots.z,alpha_Z_est)
Min_TY_est <- surv.min.T.Y(t_vec,corrxy_est,knots.t,alpha_T_est,knots.y,alpha_Y_est)


### Estimating the survival using the standard (McLain et al. 2014) approach
Surv_T_naive <- CD_surv(x.data,method = 'Piecewise', knots = knots.t)$Surv_est

#### Plotting the results ####
plot(t_vec,Surv_T_est,lwd=3,type = "l",xlim=c(0,max(t_vec)),ylim=c(0,1),ylab="Survival Function",xlab="Time",las=1)
lines(t_vec,Surv_Y_est,col=2,lwd=3)
lines(t_vec,Surv_Z_est,col=3,lwd=3)
lines(t_vec,Min_TY_est,col=4,lwd=3)
lines(t_vec,Surv_T_naive(t_vec),col=1,lwd=3,lty=2)
legend("topright",legend = c("Survival of T","Survival of Y","Survival of Z","Survival of min(T,Y)","Survival of T (naive)"),lwd=3,col=c(1:4,1),lty=c(rep(1,4),2))

#### Table of selected values ####
t_want <- c(3,6,12,24,36)
surv_res <- data.frame(time=t_want,Est_surv_T=Surv_T_est[t_vec %in% t_want],Naive_surv_T=Surv_T_naive(t_want),Est_surv_Y=Surv_Y_est[t_vec %in% t_want],Est_surv_Z=Surv_Z_est[t_vec %in% t_want],Est_surv_minTY=Min_TY_est[t_vec %in% t_want])
surv_res











#### Example of the methods with grouped data ####
ex_data <- read.csv("CDSCR_data_grp.csv")


#### Extracting data
# Total current duration time
x.data <- ex_data$T 
# Time to intermittent event 
y.data <- ex_data$TTFT 
# Indicator that intermittent event occured.
delta.data <- ex_data$fail_ind 

#Setting the number of peicewise parameters, locations will be choosen within function (knots can be given)
K1 <- 5 #For T
K2 <- 5 #For Y
K3 <- 5 #For Z


#### Run Maximum Likelihood Optimization
CDSCR_opt <- CDSCR(x.data,y.data,delta.data,K1=K1,K2=K2,K3=K3,B=20)

### Extract knots
knots.t <- CDSCR_opt$knots.t
knots.y <- CDSCR_opt$knots.y
knots.z <- CDSCR_opt$knots.z

### Likelihood value and optimization details.
CDSCR_opt$like_opt$value
CDSCR_opt$like_opt$counts
CDSCR_opt$like_opt$convergence

### Extract parameter estimates
par_est <- CDSCR_opt$like_opt$par

#Transform to corresponding values
alpha_T_est <- exp(par_est[1:(length(knots.t)+1)])
alpha_Y_est <- exp(par_est[(length(knots.t)+2):(length(c(knots.t,knots.y))+2)])
alpha_Z_est <- exp(par_est[(length(c(knots.t,knots.y))+3):(length(c(knots.t,knots.y,knots.z))+3)])
corrxy_est <-  -(exp(par_est[length(par_est)-1])/(1+exp(par_est[length(par_est)-1])))
corrxz_est <-  (exp(par_est[length(par_est)])/(1+exp(par_est[length(par_est)])))


### Calculating the survival functions of T, Y, Z and min(T,Y)
t_vec <- seq(0,36)
Surv_T_est <- PC_surv(t_vec,knots.t,alpha_T_est)
Surv_Y_est <- PC_surv(t_vec,knots.y,alpha_Y_est)
Surv_Z_est <- PC_surv(t_vec,knots.z,alpha_Z_est)
Min_TY_est <- surv.min.T.Y(t_vec,corrxy_est,knots.t,alpha_T_est,knots.y,alpha_Y_est)


### Estimating the survival using the standard (McLain et al. 2014) approach
Surv_T_naive <- CD_surv(x.data,method = 'Piecewise', knots = knots.t)$Surv_est

#### Plotting the results ####
plot(t_vec,Surv_T_est,lwd=3,type = "l",xlim=c(0,max(t_vec)),ylim=c(0,1),ylab="Survival Function",xlab="Time",las=1)
lines(t_vec,Surv_Y_est,col=2,lwd=3)
lines(t_vec,Surv_Z_est,col=3,lwd=3)
lines(t_vec,Min_TY_est,col=4,lwd=3)
lines(t_vec,Surv_T_naive(t_vec),col=1,lwd=3,lty=2)
legend("topright",legend = c("Survival of T","Survival of Y","Survival of Z","Survival of min(T,Y)","Survival of T (naive)"),lwd=3,col=c(1:4,1),lty=c(rep(1,4),2))

#### Table of selected values ####
t_want <- c(3,6,12,24,36)
surv_res <- data.frame(time=t_want,Est_surv_T=Surv_T_est[t_vec %in% t_want],Naive_surv_T=Surv_T_naive(t_want),Est_surv_Y=Surv_Y_est[t_vec %in% t_want],Est_surv_Z=Surv_Z_est[t_vec %in% t_want],Est_surv_minTY=Min_TY_est[t_vec %in% t_want])
surv_res

