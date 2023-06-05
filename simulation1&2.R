

##################################################################################################################################################    
# For Scenarios 1 and 2 of the main paper,                                                                                                       #
# we consider two multiplicative covariates, Z1(t) and Z2, and one additive covariate, X2 (time-independent):                                    #
# Z1(t) = B1 * I(t \leq V) + B2 * I(t > V), where B1 and B2 are independent Bernoulli(0.5) and V ~ Unif(0, 3)                                    #
# Z2 follows a uniform distribution on the interval [0, 1].                                                                                      #
# X2 follows either a Bernoulli distribution with parameter 0.4 (Scenario 1) or a uniform distribution (Scenario 2) on the interval [0, 1].      #
#                                                                                                                                                #
# The simulation study uses an ES algorithm, which involves two steps:                                                                           #
# E-step: Calculate the posterior means of $\xi_i$ and $W_{ik}$ of the main paper.                                                               #
# S-step: Update the jump sizes $a_1$ and $a_2$ using their explicit forms and update $\beta$ using the Newton-Raphson method.                   #                
# Finally, the variance is estimated using a sandwich estimator.                                                                                 #
#                                                                                                                                                #      
# by Xi Ning, Yinghao Pan, Yanqing Sun, and Peter B.Gilbert                                                                                      #
##################################################################################################################################################




# libraries needed #
library(matrixStats)
library(MASS)
library(limSolve)
library(extremefit)



# functions needed #

# logarithmic transformation 
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

# used to generate the time-dependent covariates 
Z1func = function(B1, B2, V) {
	
    Z1 = function(s) {
	    ifelse(s <= V, B1, B2)
    }
	
    return (Z1)
	
}

# hazard functions
dA1 = function(s){

    1/(s + 4)

}

dA2 = function(s){

    0.1

}

# cumulative hazard functions
A1_true = function(s){

    log(1 + s/4)

}

A2_true = function(s){

    0.1 * s

}


# the log-likelihood function 
loglik = function(beta, a1, a2, Z, X, TindL, TindR, G){
  
    Zbeta = matrix(0, n, m)                       
 
    for(j in 1:p){
        Zbeta = Zbeta + beta[j] * Z[[j]]
    }
 
    EZbeta = exp(Zbeta)
  
    L11 = matrix(rep(a1,each=n), n, m) + matrix(rep(a2, each=n), n, m) * X
    L1 = L11 * EZbeta
  
    # S1 and S2 for each ith subject
    S1 = rowSums(L1 * (TindL == 1))
    S2 = rowSums(L1 * (TindR == 1))
  
    SL1 = exp(-G(S1))
    SR1 = exp(-G(S2))
    dS1 = sapply(S1, dG)
        
    T_sum = rowSums(L1 * t(matrix(rep(dS1, each=m), m, n)) * t(matrix(rep(SL1, each=m), m, n)) * TindLst)
        
    C_sum = SR1 * (Delta == 0)
        
    f = sum(log(T_sum[T_sum != 0])) + sum(log(C_sum[C_sum != 0]))
        
  
    return(f)
  
}



# simulation starts #

set.seed(10)


ns = 1   # replicate number 

convInd = matrix(NA, ns, 1)  # indicator for convergence

case = matrix(1, ns, 1)    # indicator for Scenario 1 

tau = 1  # the duration of the study 


t_new = seq(0, tau, by=0.01)
n_new = length(t_new)

beta.est = matrix(0, ns, 6 + 12 * n_new)  # records of the parameters


for (l in 1:ns)
{
  
    n = 200   # sample size
 
    #regression parameters
    beta0 = c(0.5, -0.5)
    za = qnorm(0.975)
    

    B1 = rbinom(n, size = 1, prob = 0.5)
    B2 = rbinom(n, size = 1, prob = 0.5)
    V = runif(n, 0, 3)

    Z1 = vector("list", n)                    # time-dependent
    Z2 = runif(n, 0, 1)                       # time-independent
    X1 = rbinom(n, size =1, prob = 0.4)       # Scenario 1
    # X1 = runif(n, 0, 1)                     # Scenario 2

    r = 0                                     # transformation function
    obj = trans_log(r)
    G = obj[[1]]; dG = obj[[2]]; ddG = obj[[3]]; dddG = obj[[4]] 


    ######################## generate failure times for the Cox-Aalen transformation model ################################################
    T = rep(NA, n)
    
    for (i in 1:n) {
	
        # Z1(s) for this subject #
	    Z1.temp = Z1func(B1[i], B2[i], V[i])
	
	    Z1[[i]] = Z1.temp
 
 	    integrand = function(s){
  	 	    exp(beta0[1] * Z1.temp(s) + beta0[2] * Z2[i]) * (dA1(s) + X1[i] * dA2(s))
 	    }
 	 
 	    integral = function(t){
 	 	    integrate(integrand, lower = 0, upper = t)$value
 	    }
 	
	    # survival function
 	    S = function(t){
  		    exp(-G(integral(t)))
  	    }
  	
 	    U = runif(n = 1, min = 0, max = 1)
 	
 	    InvS = function(U, lower = 0, upper = 10000) {
 	 	    uniroot((function (t) S(t) - U), lower = lower, upper = upper)$root
 	    }

 	    ## note if U is very small, for instance 0.00000001, then we will get an error message since there is no solution between 0 and 10000 ## 	
 	    T.value = try(InvS(U), silent = TRUE)
 	    if ('try-error' %in% class(T.value)) {
 	 	    T[i] <- as.numeric(10000)
 	    } else {
 		    T[i] <- T.value
 	    }
 	
    }


    #############################################################################################################################################

    ID = rep(1:n)
    U1_star = rexp(n, 0.5)
    U1 = pmin(U1_star, tau)
    Delta = ((T-U1) <= 0) * 1 
    
    before_tau = sum((U1_star < tau) * (Delta == 0))
    at_tau = sum((U1_star > tau) * (Delta == 0))
   
    
    Tstar = T * (Delta==1) + U1 * (Delta == 0)   
    data = cbind(ID, Tstar, Delta)
    
    
    # unique jump points
    tstar = sort(unique(T * (Delta==1)))
    t = tstar[-1]
    m = length(t)
    
    # indicator function of t<=T for observed failure times T
    TindL = (matrix(rep(t,each=n), n, m) <= matrix(rep(T, m), n, m)) 
    
    # indicator function of t<=T for censored observations
    TindR = (matrix(rep(t,each=n), n, m) <= matrix(rep(U1, m), n, m)) 
    
    # the index of (t_{k} == T) for exact observed failure times T
    TindLst = matrix(rep(t,each=n), n, m) == matrix(rep(T,m), n, m)
    

    # time-dependent covariate at jump points
    Z11 = matrix(0, n, m)
    for(i in 1:n){
  
        for(j in 1:m){
            Z11[i,j]=Z1[[i]](t[j])
        }
  
    }

    p = length(beta0)
    Z = list(p)
    Z[[1]] = Z11
    Z[[2]] = t(matrix(rep(Z2, each=m), m, n))
    X = t(matrix(rep(X1, each=m), m, n))

    # initial values
    a1 = rep(1/m, m)
    a2 = rep(0, m)
    beta = c(0, 0)


    maxiter = 500; error = 1; iter = 0
    while(iter < maxiter & error > 0.0001){
	
        # E step 
        Zbeta = matrix(0, n, m)
  
        for (j in 1:p){
            Zbeta = Zbeta + Z[[j]] * beta[j]
        }
  
        EZbeta = exp(Zbeta)
        L11 = matrix(rep(a1,each=n), n, m) + matrix(rep(a2, each=n), n, m) * X
        L1 = L11 * EZbeta
  
        S1 = rowSums(L1 * (TindL == 1))
        S2 = rowSums(L1 * (TindR == 1)) 
        
        SL1 = exp(-G(S1))
        SR1 = exp(-G(S2))
        
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
        
        a1pp = a1
        a2pp = a2
        betapp = beta
        
        maxiter = 500; innererror = 1; inneriter = 0
        while(inneriter < maxiter & innererror > 0.001){

            # S step
            
            Rst = Tstar  # R star
            
            #Indicator function of (t <= Rst) 
            IndRst = matrix(rep(t, each=n), n, m) <= t(matrix(rep(Rst, each=m), m, n))
        
        
            xi_star = t(matrix(rep(xi, each=m), m, n))
            
            Zbeta = matrix(0, n, m)
            
            for (j in 1:p){
                Zbeta = Zbeta + Z[[j]] * beta[j]
            }
            
            EZbeta = exp(Zbeta)
        
            D00 = colSums(IndRst * W)
            D01 = colSums(IndRst * W * X)
            D11 = colSums(IndRst * xi_star * EZbeta)
            D22 = colSums(IndRst * xi_star * EZbeta * X)
            D33 = colSums(IndRst * xi_star * EZbeta * X^2)
        
      
            a1p = a1
            a2p = a2

            a2 = (D00 * D22 - D01 * D11) / (D22^2 - D33 * D11)
            a1 = (D00 - a2 * D22) / D11
          
            a11 = matrix(rep(a1, each=n), n ,m) 
            a22 = matrix(rep(a2, each=n), n ,m) 
            
            L11 = (a11 + a22 * X)
            L1 = L11 * EZbeta
            
            score = rep(0, p)
            for (j in 1:p){
                score[j] = sum(IndRst * (W - xi_star * L1) * Z[[j]])
            }
            
            # Ubeta takes derivative w.r.t beta
            dscore = matrix(0, p, p)
            for (k in 1:p){
                for(j in 1:p){
                    dscore[k, j] = sum(-xi_star * L1 * Z[[k]] * Z[[j]] * (IndRst == 1))
                }
            }
            
            
            #update parameter  
            betap = beta
            beta = betap - Solve(dscore, score)
            
            logf = loglik(beta, a1, a2, Z, X, TindL, TindR, G)
            print(logf)
        
            innererror = sum(abs(a1 - a1p)) + sum(abs(beta - betap)) + sum(abs(a2 - a2p))
            inneriter = inneriter + 1
        
        }
        
        
        logf = loglik(beta, a1, a2, Z, X, TindL, TindR, G)
        print(logf)
  
        error = sum(abs(a1 - a1pp)) + sum(abs(beta - betapp)) + sum(abs(a2 - a2pp))
        iter = iter + 1
        
        
    }
    
    
    # records of convergence
    if (iter < maxiter){
      
        convInd[l, 1] = 1 
      
    }else{
        convInd[l, 1] = 0 
    }
    
    
    # variance estimator
    
    if (iter < maxiter){
      
        a11 = matrix(rep(a1, each=n), n ,m) 
        a22 = matrix(rep(a2, each=n), n ,m) 
      
        Zbeta = matrix(0, n, m)
      
        for (j in 1:p){
            Zbeta = Zbeta + Z[[j]] * beta[j]
        }
      
        EZbeta = exp(Zbeta)
        L11 = a11 + a22 * X
        
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
        

        ################################
        # the sum of U_{i} %*% U^{\top}_{i}
        U11 = W - xi_star * L11 * EZbeta
      
        Ubeta = matrix(NA, n, p)
      
        for (j in 1:p){
            Ubeta[,j] = rowSums(U11 * Z[[j]] * (IndRst == 1))
        }
          
        
        Uscore = matrix(0, 2*m+p, 2*m+p)
        
        for(i in 1:n){
          
            Uscore_indiv = rep(NA, 2*m+p)
            for(k in 1:m){
              
                Uscore_indiv[(2*k-1):(2*k)] = (IndRst[i,k] == 1) * U11[i,k] * c(1, X[i,k]) 
            }
          
            Uscore_indiv[(2*m+p-1):(2*m+p)] = Ubeta[i,]
            Uscore = Uscore + Uscore_indiv %*% t(Uscore_indiv)
        }
        
        #################################
      
        # the derivative of U_obs
        dUa = matrix(0, 2*m+p, 2*m+p)
        
        # the derivative of E(\xi) w.r.t ak
        X11 = matrix(1, n, m)
        Uxi1 = (IndRst == 1) * SS_star *  EZbeta * X11 * (Delta == 1) + (IndRst == 1) * ddS2 * EZbeta * X11 * (Delta == 0)
        Uxi2 = (IndRst == 1) * SS_star *  EZbeta * X * (Delta == 1) + (IndRst == 1) * ddS2 * EZbeta * X * (Delta == 0)
        
        # the derivative of E(\xi) w.r.t \beta
        Uxi3 = SS * rowSums(IndRst * L11 * EZbeta * Z[[1]]) * (Delta == 1) + ddS2 * rowSums(IndRst * L11 * EZbeta * Z[[1]]) * (Delta == 0)
        Uxi4 = SS * rowSums(IndRst * L11 * EZbeta * Z[[2]]) * (Delta == 1) + ddS2 * rowSums(IndRst * L11 * EZbeta * Z[[2]]) * (Delta == 0)
        
        
        U22 = xi_star * EZbeta
        
        # the derivative of ak w.r.t ak
        for(k in 1:m){
          
            dUa1_indiv = matrix(0, 2, 2)
            for (i in 1:n){
              
                dUa1_indiv = dUa1_indiv - IndRst[i,k] * (L11[i,k] * EZbeta[i,k] * c(1, X[i,k]) %*% t(c(Uxi1[i,k], Uxi2[i,k]))  + U22[i,k] * c(1, X[i,k]) %*% t(c(1, X[i,k])))  
              
            }
            
            dUa[(2*k-1):(2*k), (2*k-1):(2*k)] = dUa1_indiv
      
        }  
        
        
        # the derivative of ak w.r.t aj, where j > k.
        for (k in 1:(m-1)){
            for (j in (k+1):m){
            
                dUa1_indiv = matrix(0, 2, 2)
              
                for (i in 1:n){
              
                    dUa1_indiv = dUa1_indiv - IndRst[i,k] * (L11[i,k] * EZbeta[i,k] * c(1, X[i,k]) %*% t(c(Uxi1[i,j], Uxi2[i,j])))
              
                }
                
                dUa[((2*k-1):(2*k)), (2*j-1):(2*j)] = dUa1_indiv
            }
        }
        
        #the derivative of ak w.r.t aj, where j < k.
        for (k in 2:m){
          for (j in 1:(k-1)){
            
            dUa1_indiv = matrix(0, 2, 2)
            
            for (i in 1:n){
              
                dUa1_indiv = dUa1_indiv - IndRst[i,k] * (L11[i,k] * EZbeta[i,k] * c(1, X[i,k]) %*% t(c(Uxi1[i,j], Uxi2[i,j])))
              
            }
            
            dUa[(2*k-1):(2*k), (2*j-1):(2*j)] = dUa1_indiv
          }
        }
        
            
        # take derivative of ak w.r.t \beta, for k = 1, 2, ...m
        U33 = xi_star * L11 * EZbeta
          
        for(k in 1:m){
          
            dUa2_indiv = matrix(0, 2, 2)
            for (i in 1:n){
            
                dUa2_indiv = dUa2_indiv - IndRst[i,k] * (L11[i,k] * EZbeta[i,k] * c(1, X[i,k]) %*% t(c(Uxi3[i], Uxi4[i])) + U33[i,k] * c(1, X[i,k]) %*% t(c(Z[[1]][i,k], Z[[2]][i,k])))  
            
            }
          
            dUa[((2*k-1):(2*k)), (2*m+p-1):(2*m+p)] = dUa2_indiv
          
        }  
        
        
        # the derivative of \beta w.r.t ak, for k = 1, 2, ...m

        for(k in 1:m){
            
            dbeta.base = matrix(0, 2, 2)
            dbeta.a1 = matrix(0, 2, 2)
          
            for (i in 1:n){
                
                for (j in 1:m){
                  
                    dbeta.base = dbeta.base - IndRst[i,j] * L11[i,j] * EZbeta[i,j] * c(Z[[1]][i,j], Z[[2]][i,j]) %*% t(c(Uxi1[i,k], Uxi2[i,k]))
                  
                }
              
                #U22 = xi_star * EZbeta
                dbeta.a1 = dbeta.a1 - IndRst[i,k] * U22[i,k] * c(Z[[1]][i,k], Z[[2]][i,k]) %*% t(c(1, X[i,k]))  
            }
          
            dUa[ (2*m+p-1):(2*m+p), (2*k-1):(2*k)] = dbeta.base + dbeta.a1
          
        }  
        
             
        dUxi3 = list(2)
        dUxi3[[1]] = t(matrix(rep(Uxi3, each=m), m, n))
        dUxi3[[2]] = t(matrix(rep(Uxi4, each=m), m, n))
        
        dUbeta = matrix(NA, p, p)
        for (k in 1:p){
            for(j in 1:p){
                # U33 = xi_star * L11 * EZbeta
                dUbeta[k, j] = -sum((L11 * EZbeta * Z[[k]] * dUxi3[[j]] + U33 * Z[[k]] * Z[[j]]) * (IndRst == 1))
            }
        }
        
        dUa[(2*m+p-1):(2*m+p),(2*m+p-1):(2*m+p)] = dUbeta
          
            
        x.inv = try(solve(dUa), silent = TRUE)
      
        if ('try-error' %in% class(x.inv)){
            print("variance problem") 
        }else{
            dUa.inv = x.inv
        }
      
        cov = dUa.inv %*% Uscore %*% t(dUa.inv)
        
        se = sqrt(c(cov[(2*m+p-1), (2*m+p-1)], cov[(2*m+p), (2*m+p)]))
      
    }
    
  
    
    mat1 = t(matrix(rep(t_new, each = m), m, n_new))
    mat2 = matrix(rep(t, each = n_new), n_new, m)
    
    # indicator function of t_new <= t
    matInd = (mat1 > mat2)
    
    # cumulative hazards when t = seq(0, tau, by=0.01)
    A1 = rowSums(matInd * matrix(rep(a1, each = n_new), n_new, m))
    A2 = rowSums(matInd * matrix(rep(a2, each = n_new), n_new, m))
    
    var_A1 = rep(NA, n_new)
    var_A2 = rep(NA, n_new)
    

    # the largest index such that n_new <= t
    index = rowSums(matInd)
    
    for (k in 1:n_new){
      
        if(index[k] == 0){
            var_A1[k] = 0
            var_A2[k] = 0
        }else{
        
            cov_indiv = cov[1:(p*index[k]),1:(p*index[k])]
        
            var_A1[k] = sum(cov_indiv[seq(1, p*index[k], 2),seq(1, p*index[k], 2)])
            var_A2[k] = sum(cov_indiv[seq(2, p*index[k], 2),seq(2, p*index[k], 2)])
        }
    }
    
    sd_A1 = sqrt(var_A1)
    sd_A2 = sqrt(var_A2)   
    
    # records of beta, SEE, CP, A1 and A2
    beta.est[l, 1] = beta[1]
    beta.est[l, 2] = se[1]
    beta.est[l, 3] = (abs(beta[1] - beta0[1]) <= za * se[1])
    
    
    beta.est[l, 4] = beta[2]
    beta.est[l, 5] = se[2]
    beta.est[l, 6] = (abs(beta[2] - beta0[2]) <= za * se[2])
    
    beta.est[l, 7 : (n_new + 6)] = A1
    beta.est[l, (n_new + 7) : (6 + 2 * n_new)] = A2 
    
    # records of the standard deviation of A1 and A2
    beta.est[l, (6 + 2 * n_new + 1): (6 + 3 * n_new)] = sd_A1
    beta.est[l, (6 + 3 * n_new + 1): (6 + 4 * n_new)] = sd_A2
    
    # records of the CP of A1 and A2
    
    A10 = sapply(t_new, A1_true)
    A20 = sapply(t_new, A2_true)
    
    for (j in 1:n_new){
        beta.est[l, (6 + 4 * n_new + j)] = (abs(A1[j] - A10[j]) <= za * sd_A1[j])
    }
    
    for (j in 1:n_new){
        beta.est[l, (6 + 5 * n_new + j)] = (abs(A2[j] - A20[j]) <= za * sd_A2[j])
    }


    # kernel estimation for alpha
    
    h = 0.1   # bandwidth for alpha
    
    tk = matrix(NA, m, n_new)
    y.ker = matrix(NA, m, n_new)
    alpha1 = rep(NA, n_new)
    alpha2 = rep(NA, n_new)
    
    for (j in 1:n_new){
      
      tk[,j] = (t_new[j] - t) / h
      y.ker[,j] = Epa.kernel(tk[,j])
      
    }
    
    for (j in 1:n_new){
      
        alpha1[j] = sum(1/h * y.ker[,j] * a1)
        alpha2[j] = sum(1/h * y.ker[,j] * a2)
      
    }
    
    var_alpha1 = rep(0, n_new)
    var_alpha2 = rep(0, n_new)
    cov_0 = cov[1:(2*m), 1:(2*m)]
    
    for (i in 1: n_new){
        for (k in 1:m){
            for (j in 1:m){
              
                var_alpha1[i] = var_alpha1[i] + 1/ (h^2) * y.ker[k,i] * y.ker[j,i] * cov_0[(2*k-1), (2*j-1)]
                var_alpha2[i] = var_alpha2[i] + 1/ (h^2) * y.ker[k,i] * y.ker[j,i] * cov_0[2*k, 2*j]
            }
        }
    }
    
    sd_alpha1 = sqrt(var_alpha1)
    sd_alpha2 = sqrt(var_alpha2)
    
    alpha10 = sapply(t_new, dA1)
    alpha20 = sapply(t_new, dA2)
    
    # records of alpha1 and alpha2
    beta.est[l, (6 + 6 * n_new + 1): (6 + 7 * n_new)] = alpha1
    beta.est[l, (6 + 7 * n_new + 1): (6 + 8 * n_new)] = alpha2
    
    # records of the standard deviation of alpha1 and alpha2
    beta.est[l, (6 + 8 * n_new + 1): (6 + 9 * n_new)] = sd_alpha1
    beta.est[l, (6 + 9 * n_new + 1): (6 + 10 * n_new)] = sd_alpha2
    
    for (j in 1:n_new){
        beta.est[l, (6 + 10 * n_new + j)] = (abs(alpha1[j] - alpha10[j]) <= za * sd_alpha1[j])
    }
    
    for (j in 1:n_new){
        beta.est[l, (6 + 11 * n_new + j)] = (abs(alpha2[j] - alpha20[j]) <= za * sd_alpha2[j])
    }
    
    
    output = c(beta.est[l,1:6], convInd[l,1], case[l,1], beta.est[l, 7:(6 + 12 * n_new)])
  

  
}




