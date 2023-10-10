#========================================================#
# Quantitative ALM, Financial Econometrics & Derivatives 
# ML/DL using R, Python, Tensorflow by Sang-Heon Lee 
#
# https://kiandlee.blogspot.com
#——————————————————#
# Lam’s generalized Hamilton model in Kim (1994)
#========================================================#

graphics.off()  # clear all graphs
rm(list = ls()) # remove all files from your workspace

#=======================================
#             Functions
#=======================================

# parameter constraints
trans = function(b0){
  b1 = b0
  
  # probability
  b1[1:2] = exp(b0[1:2])/(1+exp(b0[1:2]))
  
  # variance
  b1[5] = b0[5]^2
  
  # AR(1), AR(2) coefficients
  XX6 = b0[6]/(1+abs(b0[6]))
  XX7 = b0[7]/(1+abs(b0[7]))
  b1[6] = XX6 + XX7
  b1[7] = -1*XX6*XX7
  
  return(b1)
}

# Kim filter = Hamilton + Kalman
loglike <- function(param_unc, y) {
  
  # param_unc = init_param_unc; y = y
  nobs <- length(y)
  
  param <- trans(param_unc)
  
  p11 <- param[1] # Pr[St=1/St-1=1]
  p00 <- param[2] # Pr[St=0/St-1=0]
  MU0 <- param[3]       # delta0
  MU1 <- MU0 + param[4] # delta0 + delta1
  VAR_W <- param[5]     # sigma^2
  phi1 <- param[6]      # AR(1)
  phi2 <- param[7]      # AR(2)
  
  R <- matrix(0,1,1)
  Q <- rbind(c(VAR_W, 0), c(0, 0))
  H <- t(c(1, -1))
  f <- rbind(c(phi1, phi2), c(1, 0))
  
  # initialization
  B_LL0 <- B_LL1 <- param[8:9] #c(0,0)
  P_LL1 <- P_LL0 <- 
  matrix(solve(diag(4) - kronecker(f,f))%*%
           as.vector(Q), nrow = 2)
  
  # initialization with steady state probabilities
  PROB_1 <- (1-p00)/(2-p11-p00) #Pr[S_0=1/Y_0]
  PROB_0 <- 1-PROB_1            #Pr[S_0=0/Y_0]
  
  #——————————————
  #  START ITERATION
  #——————————————
  likv <- 0.0
  for(t in 1:nobs) {
    
    #———————————
    # KALMAN FILTER
    #——————————— 
    
    # state prediction
    B_TL00 <- f %*% B_LL0  # S=0, S’=0
    B_TL01 <- f %*% B_LL0  # S=0, S’=1
    B_TL10 <- f %*% B_LL1  # S=1, S’=0
    B_TL11 <- f %*% B_LL1  # S=1, S’=1
    
    # state prediction uncertainty
    P_TL00 <- f %*% P_LL0 %*% t(f) + Q
    P_TL01 <- f %*% P_LL0 %*% t(f) + Q
    P_TL10 <- f %*% P_LL1 %*% t(f) + Q
    P_TL11 <- f %*% P_LL1 %*% t(f) + Q
    
    # forecast   
    fit00 <- H %*% B_TL00 + MU0
    fit01 <- H %*% B_TL01 + MU1
    fit10 <- H %*% B_TL10 + MU0
    fit11 <- H %*% B_TL11 + MU1
    
    # forecast error
    F_CAST00 <- y[t] - fit00
    F_CAST01 <- y[t] - fit01
    F_CAST10 <- y[t] - fit10
    F_CAST11 <- y[t] - fit11
    
    # forecast error variance
    SS00 <-  H %*% P_TL00 %*% t(H) + R
    SS01 <-  H %*% P_TL01 %*% t(H) + R
    SS10 <-  H %*% P_TL10 %*% t(H) + R
    SS11 <-  H %*% P_TL11 %*% t(H) + R
    
    # Kalman gain
    kg00 <- P_TL00 %*% t(H) %*% solve(SS00)
    kg01 <- P_TL01 %*% t(H) %*% solve(SS01)
    kg10 <- P_TL10 %*% t(H) %*% solve(SS10)
    kg11 <- P_TL11 %*% t(H) %*% solve(SS11)
    
    # updating
    B_TT00 <- B_TL00 + kg00 %*% F_CAST00
    B_TT01 <- B_TL01 + kg01 %*% F_CAST01
    B_TT10 <- B_TL10 + kg10 %*% F_CAST10
    B_TT11 <- B_TL11 + kg11 %*% F_CAST11
    
    P_TT00 <- (diag(2) - kg00 %*% H ) %*% P_TL00
    P_TT01 <- (diag(2) - kg01 %*% H ) %*% P_TL01
    P_TT10 <- (diag(2) - kg10 %*% H ) %*% P_TL10
    P_TT11 <- (diag(2) - kg11 %*% H ) %*% P_TL11
    
    #——————————— 
    # HAMILTON FILTER
    #——————————— 
    
    # Pr[St,Yt/Yt-1]
    NORM_PDF00 <- dnorm(F_CAST00,0,sqrt(SS00))
    NORM_PDF01 <- dnorm(F_CAST01,0,sqrt(SS01))
    NORM_PDF10 <- dnorm(F_CAST10,0,sqrt(SS10))
    NORM_PDF11 <- dnorm(F_CAST11,0,sqrt(SS11))
    
    PR_VL00 <- NORM_PDF00*(  p00)*PROB_0
    PR_VL01 <- NORM_PDF01*(1-p00)*PROB_0
    PR_VL10 <- NORM_PDF10*(1-p11)*PROB_1
    PR_VL11 <- NORM_PDF11*(  p11)*PROB_1
    
    # f(y_t|Y_{t-1})
    PR_VAL <- PR_VL00 + PR_VL01 + PR_VL10 + PR_VL11
    
    # Pr[St,St-1/Yt]
    PRO_00 <- as.numeric(PR_VL00/PR_VAL)
    PRO_01 <- as.numeric(PR_VL01/PR_VAL)
    PRO_10 <- as.numeric(PR_VL10/PR_VAL)
    PRO_11 <- as.numeric(PR_VL11/PR_VAL)
    
    # Pr[St/Yt]
    PROB_0 <- PRO_00 + PRO_10
    PROB_1 <- PRO_01 + PRO_11
    
    #——————————— 
    # Kim’s Collapsing 
    #——————————— 
    B_LL0 <- (PRO_00*B_TT00 + PRO_10*B_TT10)/PROB_0
    B_LL1 <- (PRO_01*B_TT01 + PRO_11*B_TT11)/PROB_1
    
    temp00 <- (B_LL0-B_TT00) %*% t(B_LL0-B_TT00)
    temp10 <- (B_LL0-B_TT10) %*% t(B_LL0-B_TT10)
    temp01 <- (B_LL1-B_TT01) %*% t(B_LL1-B_TT01)
    temp11 <- (B_LL1-B_TT11) %*% t(B_LL1-B_TT11)
    
    P_LL0 <- (PRO_00*(P_TT00 + temp00) +
                PRO_10*(P_TT10 + temp10))/PROB_0
    P_LL1 <- (PRO_01*(P_TT01 + temp01) +
                PRO_11*(P_TT11 + temp11))/PROB_1
    
    likv = likv -log(PR_VAL)
  } 
  return(likv)
}

#=======================================
#                  Run
#=======================================

# Lam’s Real GNP Data Set : #1952.3 ~ 1984.4
rgnp <- c(  1378.2, 1406.8, 1431.4, 1444.9, 1438.2, 
            1426.6, 1406.8, 1401.2, 1418  , 1438.8, 1469.6, 
            1485.7, 1505.5, 1518.7, 1515.7, 1522.6, 1523.7, 
            1540.6, 1553.3, 1552.4, 1561.5, 1537.3, 1506.1, 
            1514.2, 1550  , 1586.7, 1606.4, 1637  , 1629.5, 
            1643.4, 1671.6, 1666.8, 1668.4, 1654.1, 1671.3, 
            1692.1, 1716.3, 1754.9, 1777.9, 1796.4, 1813.1, 
            1810.1, 1834.6, 1860  , 1892.5, 1906.1, 1948.7, 
            1965.4, 1985.2, 1993.7, 2036.9, 2066.4, 2099.3, 
            2147.6, 2190.1, 2195.8, 2218.3, 2229.2, 2241.8, 
            2255.2, 2287.7, 2300.6, 2327.3, 2366.9, 2385.3, 
            2383  , 2416.5, 2419.8, 2433.2, 2423.5, 2408.6, 
            2406.5, 2435.8, 2413.8, 2478.6, 2478.4, 2491.1, 
            2491  , 2545.6, 2595.1, 2622.1, 2671.3, 2734  , 
            2741  , 2738.3, 2762.8, 2747.4, 2755.2, 2719.3, 
            2695.4, 2642.7, 2669.6, 2714.9, 2752.7, 2804.4, 
            2816.9, 2828.6, 2856.8, 2896  , 2942.7, 3001.8, 
            2994.1, 3020.5, 3115.9, 3142.6, 3181.6, 3181.7, 
            3178.7, 3207.4, 3201.3, 3233.4, 3157  , 3159.1, 
            3199.2, 3261.1, 3250.2, 3264.6, 3219  , 3170.4, 
            3179.9, 3154.5, 3159.3, 3186.6, 3258.3, 3306.4, 
            3365.1, 3444.7, 3487.1, 3507.4, 3520.4) 

# # GNP growth rate : 1952.4 ~ 1984.4
y <- diff(log(rgnp),1)*100

# initial guess
init_param_unc <- c(3.5,0.0,-0.4,1.7,0.5,1,0.7,1,1)

# optimization
mle <- optim(init_param_unc, loglike, method=c('BFGS'),
             control=list(maxit=50000,trace=2), y=y)
mle <- optim(mle$par, loglike, method=c('Nelder-Mead'),
             control=list(maxit=50000,trace=2), y=y)

# compare estimation results with Kim (1994) paper
est_param <- cbind(c(trans(mle$par), -mle$value),
                   c(0.954, 0.456, -1.457, 2.421, 0.773, 
                     1.246, -0.367, 5.224, 0.535, -176.33))

est_param[5,1] <- sqrt(est_param[5,1]) # var -> std

colnames(est_param) <- c('our est', 'Kim est')
rownames(est_param) <- c('p00','p11','delta0', 'delta1', 'sigma', 
   'phi1', 'phi2','x0','x-1','likv')

est_param
