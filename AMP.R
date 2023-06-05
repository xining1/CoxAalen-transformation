
##################################################################################################################################################    
# For the real-data application (the AMP dataset) of the main paper                                                                              #
# we consider two multiplicative covariates, Z1 and Z2, and one additive covariate, X:                                                           #
# Z1: the treatment assignment, i.e., the placebo, low-dose and high-dose                                                                        #
# Z2: the age group indicator                                                                                                                    #
# X:  the region indicator                                                                                                                       #
#                                                                                                                                                #                                                                                                                                            
# by Xi Ning, Yinghao Pan, Yanqing Sun, and Peter B.Gilbert                                                                                      #
##################################################################################################################################################


# libraries needed #

library(matrixStats)
library(MASS)
library(limSolve)


# functions needed #

trans_log = function(r){
  
    transform = list()
    if (r == 0){
        transform[[1]] = function(x){
            return(x)
        }
        transform[[2]] = function(x){
            return(1)
        }
        transform[[3]] = function(x){
            return(0)
        }
        transform[[4]] = function(x){
            return(0)
        }
    }else{
        transform[[1]] = function(x){
            return(1 / r * log(1 + r * x))
        }
        transform[[2]] = function(x){
            return(1 / (1 + r * x))
        }
        transform[[3]] = function(x){
            return(- r / (1 + r * x)^2)
        }
        transform[[4]] = function(x){
            return((2 * r^2) / (1 + r * x)^3)
        }
    }
  
    return (transform);
  
}



# the log-likelihood function 
loglik = function(beta, a1, Z, X11, TindL, TindR, G){
  
  
    Zbeta = matrix(0, n, m)                       
  
    # get the beta1 * Z1 + beta2 * Z2
    for(j in 1:p){
        Zbeta = Zbeta + beta[j] * Z[[j]]
    }
  
    EZbeta = exp(Zbeta)
  
    L11 = 0
  
    # get a1 + a2 * X1 + a3 * X2
    for (j in 1:(J+1)){
     
        L11 = L11 + matrix(rep(a1[[j]], each=n), n, m) * X11[[j]]
    
    }
  
    L1 = L11 * EZbeta
  
    # S1 and S2 for each ith subject 
    S1 = rowSums(L1 * (TindL == 1))
    S2 = rowSums(L1 * (TindR == 1))
  
    # summation of t <= T and t <= C
    SL1 = exp(-G(S1))
    SR1 = exp(-G(S2))
    dS1 = sapply(S1, dG)
  
    T_sum = rowSums(L1 * t(matrix(rep(dS1, each=m), m, n)) * t(matrix(rep(SL1, each=m), m, n)) * TindLst)
  
    C_sum = SR1 * (Delta == 0)
  

    f = sum(log(T_sum[T_sum != 0])) + sum(log(C_sum[C_sum != 0]))
  
    return(f)
  
}



# download the AMP dataset (the HIV-1 trials) and read it in R, for example #

amp <- read.csv('/Users/xining/Desktop/AMP_survival_20211102.csv',header = TRUE, row.names=NULL)


# the placebo is indexed at 0, low-dose is 1, and high-dose is 2
tx_index = (amp$tx == "T1") * 1 + (amp$tx == "T2") * 2

# four age groups
age_index = (amp$age <= 20) * 1 +  (amp$agegrp == "21 - 30") * 2 + (amp$agegrp == "31 - 40") * 3 + (amp$age >=41) * 4

# four regions
country_index = (amp$country =="United States" | amp$country=="Switzerland" ) * 1 + (amp$country=="Peru"| amp$country=="Brazil") * 2 + (amp$country=="South Africa") * 3
country_index[country_index ==0] = 4 



Time.star = amp$hiv1survday
Delta.star = amp$hiv1event
Z1.star = tx_index
Z2.star = age_index
X.star = country_index
protocol = (amp$protocol == "HVTN 704") * 1 + (amp$protocol == "HVTN 703") * 2

# the HIV data 
Data = cbind(Time.star, Delta.star, Z1.star, Z2.star, X.star, protocol)
Data1 = Data[(Time.star != 0), ]        # remove 0 in survival time


Time = Data1[,1]/7
Delta1 = Data1[,2]
Z1 = Data1[,3]
Z2 = Data1[,4]
X1 = Data1[,5]
pro.ind = Data1[,6]


tau = 85.9                          # the end of the study
U1 = pmin(Time, tau)                # the observed time
Delta = Delta1 * (Time <= tau)      # the indicator of censorship

# the unique jump times 
tstar = sort(unique(U1 * (Delta==1)))
t = tstar[-1]  
m =length(t)
n = nrow(Data1)

TindL = matrix(rep(t,each=n), n, m) <= matrix(rep(U1,m), n, m)
TindR = matrix(rep(t,each=n), n, m) <= matrix(rep(U1,m), n, m)
TindLst = ((matrix(rep(t,each=n), n, m) == matrix(rep(U1,m), n, m))) * (Delta==1)   # the index of t_{k} == T


# the covariate Z1 and Z2 in the main paper
p = 5
Z = list(p)    
Z[[1]] = t(matrix(rep((Z1==1)*1, each=m), m, n))
Z[[2]] = t(matrix(rep((Z1==2)*1, each=m), m, n))
Z[[3]] = t(matrix(rep((Z2==2)*1, each=m), m, n))
Z[[4]] = t(matrix(rep((Z2==3)*1, each=m), m, n))
Z[[5]] = t(matrix(rep((Z2==4)*1, each=m), m, n))


# the covariate X in the main paper
X = matrix(NA, n, 4)      

X[, 1] = rep(1, n)
X[, 2] = ((X1 == 2))*1 
X[, 3] = ((X1 == 3))*1 
X[, 4] = ((X1 == 4))*1 



J = 3
X11 = list(J+1)

for (j in 1:(J+1)){
  
    X11[[j]] = t(matrix(rep(X[, j], each=m), m, n))
  
}


# the transformation model used in the main paper indexed by r
r = 0

obj = trans_log(r)
G = obj[[1]]; dG = obj[[2]]; ddG = obj[[3]]; dddG = obj[[4]]



# set initial values for beta and ak

beta = c(0, 0, 0, 0, 0)


a1 = list(J+1)
a1_new = list(J+1)

a1[[1]] = rep(1/m, m)

for (j in 2:(J+1)){
  
    a1[[j]] = rep(0, m)
  
}


######### ES algorithm starts ########


maxiter = 500; error = 1; iter = 0
while(iter < maxiter & error > 0.0001){
  
    #E step a1
    Zbeta = matrix(0, n, m)
  
    for (j in 1:p){
        Zbeta = Zbeta + Z[[j]] * beta[j]
    }
  
    EZbeta = exp(Zbeta)
  
    L11 = 0
  
    for (j in 1:(J+1)){
    
        L11 = L11 + matrix(rep(a1[[j]], each=n), n, m) * X11[[j]]
    
    }
  
    L1 = L11 * EZbeta
  
    S1 = rowSums(L1 * (TindL == 1)) 
    S2 = rowSums(L1 * (TindR == 1)) 
  
    SL1 = exp(-G(S1))
    SR1 = exp(-G(S2))
  
    # saaply function function takes list or vector as input and gives output in vector or matrix
    # get dG(S1) and dG(S2)
    dS1 = sapply(S1, dG)
    dS2 = sapply(S2, dG)
  
    ddS1 = sapply(S1, ddG)
    
    # the posterior mean for exact failure time observations
    xi1 = dS1 - ddS1/dS1 
  
    # the posterior mean for censored observations
    xi2 = dS2 
  
    # the posterior mean for \xi
    xi =  xi1 * (Delta == 1) +  xi2 * (Delta == 0)
    print(xi)
  
    # the posterior mean for W, w = 1 when t_{k} = T, otherwise, W = 0.
    W = TindLst * 1         
  
    # update the parameter in the E-step
    a1pp = a1
    betapp = beta
  
    # S step
  
    innermaxiter = 500; innererror = 1; inneriter = 0
    while(inneriter < innermaxiter & innererror > 0.001){
    
    
    
        IndRst = matrix(rep(t, each=n), n, m) <= t(matrix(rep(U1, each=m), m, n))
    
    
        # replicate \xi for each t_{k}
        xi_star = t(matrix(rep(xi, each=m), m, n))
    
        Zbeta = matrix(0, n, m)
    
        for (j in 1:p){
            Zbeta = Zbeta + Z[[j]] * beta[j]
        }
    
        EZbeta = exp(Zbeta)
    
        U22 = xi_star * EZbeta 
    
        # update the parameters ak in the S-step
        a1p = a1
   
        for (k in 1:m){
      
            sum1 = matrix(0, (J+1), (J+1))
            sum2 = rep(0, (J+1))     
      
            for(i in 1:n){
        
                sum1 = sum1 + IndRst[i,k] * U22[i,k] * c(X11[[1]][i,k], X11[[2]][i,k], X11[[3]][i,k], X11[[4]][i,k]) %*% t(c(X11[[1]][i,k], X11[[2]][i,k], X11[[3]][i,k], X11[[4]][i,k]))
                sum2 = sum2 + IndRst[i,k] * W[i,k] * c(X11[[1]][i,k], X11[[2]][i,k], X11[[3]][i,k], X11[[4]][i,k])
            } 

            sum1.inv = try(solve(sum1), silent = TRUE)
      
            if ('try-error' %in% class(sum1.inv)){
                print("solve problem") 
            }else{
                sum11.inv = sum1.inv
            }
       
            a.new = sum11.inv %*% sum2
      
            a1[[1]][k] = a.new[1]
            a1[[2]][k] = a.new[2]
            a1[[3]][k] = a.new[3]
            a1[[4]][k] = a.new[4]
      
        }
    
    
    
        a11 = list(J+1)
    
        for(j in 1:(J+1)){
      
            a11[[j]] = matrix(rep(a1[[j]], each=n), n ,m) 
      
        }
    
        L11 = 0
    
        for (j in 1:(J+1)){
            L11 = L11 + a11[[j]] * X11[[j]]
        }
    
        L1 = L11 * EZbeta
    
        score = rep(0, p)
        for (j in 1:p){
            score[j] = sum(IndRst * (W - xi_star * L11 * EZbeta) * Z[[j]])
        }
    
        # Ubeta takes derivative w.r.t beta
        dscore = matrix(0, p, p)
    
        for (k in 1:p){
            for(j in 1:p){
                dscore[k, j] = sum(-xi_star * L1 * Z[[k]] * Z[[j]] * (IndRst == 1))
            }
        }
    
    
        # update parameter beta
        betap = beta
        beta = betap - Solve(dscore, score)
    
        logf = loglik(beta, a1, Z, X11, TindL, TindR, G)
        print(logf)
    
        innererrorstar = 0
        for (j in 1:(J+1)){
      
          innererrorstar = innererrorstar + sum(abs(a1[[j]] - a1p[[j]])) 
      
        }
    
        innererror = innererrorstar + sum(abs(beta - betap))
        inneriter = inneriter + 1
    
    }
  
  
    errorstar = 0
    for (j in 1:(J+1)){
    
        errorstar = errorstar + sum(abs(a1[[j]] - a1pp[[j]])) 
    
    }
  
    error = errorstar + sum(abs(beta - betapp))
    iter = iter + 1
}




######### variance estimator ########


if (iter < maxiter){
  
    a11 = list(J+1)
  
    for(j in 1:(J+1)){
    
        a11[[j]] = matrix(rep(a1[[j]], each=n), n ,m) 
    
    }
  
    L11 = 0
  
    for (j in 1:(J+1)){
        L11 = L11 + a11[[j]] * X11[[j]]
    }
  
  
    Zbeta = matrix(0, n, m)
    for (j in 1:p){
        Zbeta = Zbeta + Z[[j]] * beta[j]
    }
  
    EZbeta = exp(Zbeta)
    L1 = L11 * EZbeta
  
    # S1 and S2 for each ith subject for exact failure time and censored observations, repsectively
    S1 = rowSums(L11 * EZbeta * (TindL == 1)) 
    S2 = rowSums(L11 * EZbeta * (TindR == 1)) 
  
  
  
    # get dG(S1), dG(S2), dd(S1) and dddG(S1)
    dS1 = sapply(S1, dG)
    dS2 = sapply(S2, dG)
    ddS1 = sapply(S1, ddG)
    ddS2 = sapply(S2, ddG)
    dddS1 = sapply(S1, dddG)
  
    # the posterior mean for exact failure time observations
    xi1 = dS1 - ddS1/dS1 
  
    # the posterior mean for censored observations
    xi2 = dS2 
  
    # the posterior mean for \xi
    xi =  xi1 * (Delta == 1) +  xi2 * (Delta == 0)
    xi_star = t(matrix(rep(xi, each=m), m, n))
  
    SS = ddS1 - dddS1/dS1 - (ddS1/dS1)^2
    SS_star = t(matrix(rep(SS, each=m), m, n))
  
    # the sum of U_{i} %*% U^{\top}_{i}
    U11 = W - xi_star * L11 * EZbeta
  
    Ubeta = matrix(NA, n, p)
  
    for (j in 1:p){
        Ubeta[,j] = rowSums(U11 * Z[[j]] * (IndRst == 1))
    }
  
  
    Uscore = matrix(0, (J+1)*m+p, (J+1)*m+p)
    
    for(i in 1:n){
    
        Uscore_indiv = rep(NA, (J+1)*m+p)
        for(k in 1:m){
      
            Uscore_indiv[((J+1)*k-J):((J+1)*k)] = (IndRst[i,k] == 1) * U11[i,k] * c(X11[[1]][i,k], X11[[2]][i,k], X11[[3]][i,k], X11[[4]][i,k])
        }
    
        Uscore_indiv[((J+1)*m+1):((J+1)*m+p)] = Ubeta[i,]
        Uscore = Uscore + Uscore_indiv %*% t(Uscore_indiv)
    }

  
    # the derivative of U_obs
    dUa = matrix(0, (J+1)*m+p, (J+1)*m+p)
  
    # the derivative of E(\xi) w.r.t ak
    Uxi.ak = list(J+1)
    for (j in 1:(J+1)){
        Uxi.ak[[j]] = (IndRst == 1) * SS_star *  EZbeta * X11[[j]] * (Delta == 1) + (IndRst == 1) * ddS2 * EZbeta * X11[[j]] * (Delta == 0)
    }
  
    # the derivative of E(\xi) w.r.t \beta
    Uxi.beta = list(p)
    for (j in 1:p){
        Uxi.beta[[j]] = SS * rowSums(IndRst * L11 * EZbeta * Z[[j]]) * (Delta == 1) + ddS2 * rowSums(IndRst * L11 * EZbeta * Z[[j]]) * (Delta == 0)
    
    }
  
  
    U22 = xi_star * EZbeta
  
    # the derivative of ak w.r.t ak
    for(k in 1:m){
    
        dUa1_indiv = matrix(0, (J+1), (J+1))
        for (i in 1:n){
      
            dUa1_indiv = dUa1_indiv - IndRst[i,k] * (L11[i,k] * EZbeta[i,k] * c(X11[[1]][i,k], X11[[2]][i,k], X11[[3]][i,k], X11[[4]][i,k]) %*% t(c(Uxi.ak[[1]][i,k], Uxi.ak[[2]][i,k], Uxi.ak[[3]][i,k], Uxi.ak[[4]][i,k])) + U22[i,k] * c(X11[[1]][i,k], X11[[2]][i,k], X11[[3]][i,k], X11[[4]][i,k]) %*% t(c(X11[[1]][i,k], X11[[2]][i,k], X11[[3]][i,k], X11[[4]][i,k])))  
      
        }
    
        dUa[((J+1)*k-J):((J+1)*k), ((J+1)*k-J):((J+1)*k)] = dUa1_indiv
    
    }  
  
  
    # the derivative of Uak w.r.t aj, where j > k.
    for (k in 1:(m-1)){
        for (j in (k+1):m){
      
            dUa1_indiv = matrix(0, (J+1), (J+1))
      
            for (i in 1:n){
        
                dUa1_indiv = dUa1_indiv - IndRst[i,k] * (L11[i,k] * EZbeta[i,k] * c(X11[[1]][i,k], X11[[2]][i,k], X11[[3]][i,k], X11[[4]][i,k]) %*% t(c(Uxi.ak[[1]][i,j], Uxi.ak[[2]][i,j], Uxi.ak[[3]][i,j], Uxi.ak[[4]][i,j])))
        
            }
      
            dUa[((J+1)*k-J):((J+1)*k), ((J+1)*j-J):((J+1)*j)] = dUa1_indiv
        }
    }
  
    #the derivative of ak w.r.t aj, where j < k.
    for (k in 2:m){
        for (j in 1:(k-1)){
      
            dUa1_indiv = matrix(0, (J+1), (J+1))
      
            for (i in 1:n){
        
                dUa1_indiv = dUa1_indiv - IndRst[i,k] * (L11[i,k] * EZbeta[i,k] * c(X11[[1]][i,k], X11[[2]][i,k], X11[[3]][i,k], X11[[4]][i,k]) %*% t(c(Uxi.ak[[1]][i,j], Uxi.ak[[2]][i,j], Uxi.ak[[3]][i,j], Uxi.ak[[4]][i,j])))
        
            }
      
            dUa[((J+1)*k-J):((J+1)*k), ((J+1)*j-J):((J+1)*j)] = dUa1_indiv
        }
    }
  
  
  
    # take derivative of ak w.r.t \beta, for k = 1, 2, ...m
    U33 = xi_star * L11 * EZbeta
  
    for(k in 1:m){
    
        dUa2_indiv = matrix(0, (J+1), p)
        for (i in 1:n){
      
            dUa2_indiv = dUa2_indiv - IndRst[i,k] * (L11[i,k] * EZbeta[i,k] * c(X11[[1]][i,k], X11[[2]][i,k], X11[[3]][i,k], X11[[4]][i,k]) %*% t(c(Uxi.beta[[1]][i], Uxi.beta[[2]][i], Uxi.beta[[3]][i], Uxi.beta[[4]][i], Uxi.beta[[5]][i])) + U33[i,k] * c(X11[[1]][i,k], X11[[2]][i,k], X11[[3]][i,k], X11[[4]][i,k]) %*% t(c(Z[[1]][i,k], Z[[2]][i,k], Z[[3]][i,k], Z[[4]][i,k], Z[[5]][i,k])))  
      
        }
    
        dUa[((J+1)*k-J):((J+1)*k), ((J+1)*m+1):((J+1)*m+p)] = dUa2_indiv
    
    }  
  
  
    # the derivative of \beta w.r.t ak, for k = 1, 2, ...m
  
  
    for(k in 1:m){
    
        dbeta.base = matrix(0, p, (J+1))
        dbeta.a1 = matrix(0, p, (J+1))
    
        for (i in 1:n){
      
            for (j in 1:m){
        
                dbeta.base = dbeta.base - IndRst[i,j] * L11[i,j] * EZbeta[i,j] * c(Z[[1]][i,j], Z[[2]][i,j], Z[[3]][i,j], Z[[4]][i,j], Z[[5]][i,j]) %*% t(c(Uxi.ak[[1]][i,k], Uxi.ak[[2]][i,k], Uxi.ak[[3]][i,k], Uxi.ak[[4]][i,k]))
        
            }
      
            dbeta.a1 = dbeta.a1 - IndRst[i,k] * U22[i,k] * c(Z[[1]][i,k], Z[[2]][i,k], Z[[3]][i,k], Z[[4]][i,k], Z[[5]][i,k]) %*% t(c(X11[[1]][i,k], X11[[2]][i,k], X11[[3]][i,k], X11[[4]][i,k]))
        }
    
        dUa[ ((J+1)*m+1):((J+1)*m+p), ((J+1)*k-J):((J+1)*k)] = dbeta.base + dbeta.a1
    
    }  
  
  
  
  
    dUxi3 = list(p)
  
    for (j in 1:p){
        dUxi3[[j]] = t(matrix(rep(Uxi.beta[[j]], each=m), m, n))
    }
  
    dUbeta = matrix(NA, p, p)
    for (k in 1:p){
        for(j in 1:p){
      
            # U33 = xi_star * L11 * EZbeta
            dUbeta[k, j] = -sum((L11 * EZbeta * Z[[k]] * dUxi3[[j]] + U33 * Z[[k]] * Z[[j]]) * (IndRst == 1))
        }
    }
  
    dUa[((J+1)*m+1):((J+1)*m+p), ((J+1)*m+1):((J+1)*m+p)] = dUbeta
  
  
    x.inv = try(solve(dUa), silent = TRUE)
  
    if ('try-error' %in% class(x.inv)){
        print("variance problem") 
    }else{
        dUa.inv = x.inv
    }
  
    cov = dUa.inv %*% Uscore %*% t(dUa.inv)
  
    se = sqrt(c(cov[((J+1)*m+p-4), ((J+1)*m+p-4)], cov[((J+1)*m+p-3), ((J+1)*m+p-3)],cov[((J+1)*m+p-2), ((J+1)*m+p-2)], cov[((J+1)*m+p-1), ((J+1)*m+p-1)], cov[((J+1)*m+p), ((J+1)*m+p)]))
  
}






