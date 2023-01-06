
# =============  X的分布不一样================================= 

gen.dat1 <- function(n1, n2, p){
  
  # ---------------  RCT  ---------------- # 
  X <- matrix(rnorm(n1*p), nrow = n1)      # X
  U <- rnorm(n1)                        # U, unmeasured confounder 
  Tr <- rbinom(n1, 1, prob = 0.5)  # T
  
  beta.s = rep(1, p) 
  S <- Tr + 2 * (X %*% beta.s) + U  + rnorm(n1) # S
  
  beta.y = rep(1, p)
  Y <- Tr + 3*(X%*%beta.y) + S + rnorm(n1)
  
  dat.rct <- data.frame(X, Tr, S, Y)
  names(dat.rct) <- c(paste0('X',1:p), 'Tr', 'S', 'Y')
  
  # ---------------  OBS  ---------------- # 
  X <- matrix(rnorm(n2*p, mean = 1, sd = 2), nrow = n2)   # X
  U <- rnorm(n2)                        # U, unmeasured confounder 
  
  beta.t = rep(1, p)
  prob.t = 1/(1 + exp(-X%*%beta.t - U))
  Tr <- rbinom(n2, 1, prob = prob.t )    # T
  
  beta.s = rep(1, p) 
  S <- Tr + 2 * (X %*% beta.s) + U  + rnorm(n2) # S
  
  beta.y = rep(1, p)
  prob.y = 1/(1 + exp(-Tr- 3*(X%*%beta.y) - S))
  Y <- Tr + 3*(X%*%beta.y) + S + rnorm(n2)
  
  dat.obs <- data.frame(X, Tr, S, Y)
  names(dat.obs) <- c(paste0('X',1:p), 'Tr', 'S', 'Y')
  
  return(list(dat.rct = dat.rct, dat.obs = dat.obs))
}


gen.dat2 <- function(n1, n2, p){
  
  # ---------------  RCT  ---------------- # 
  X <- matrix(rnorm(n1*p), nrow = n1)      # X
  U <- rnorm(n1)                        # U, unmeasured confounder 
  Tr <- rbinom(n1, 1, prob = 0.5)  # T
  
  beta.s = rep(1, p) 
  S <- Tr + 2 * (X^2 %*% beta.s) + U^2  + rnorm(n1) # S
  
  beta.y = rep(1, p)
  Y <- Tr + 3*(X%*%beta.y) + S + rnorm(n1)
  
  dat.rct <- data.frame(X, Tr, S, Y)
  names(dat.rct) <- c(paste0('X',1:p), 'Tr', 'S', 'Y')
  
  # ---------------  OBS  ---------------- # 
  X <- matrix(rnorm(n2*p, mean = 1, sd = 2), nrow = n2)   # X
  U <- rnorm(n2)                        # U, unmeasured confounder 
  
  beta.t = rep(1, p)
  prob.t = 1/(1 + exp(-X%*%beta.t - U))
  Tr <- rbinom(n2, 1, prob = prob.t )    # T
  
  beta.s = rep(1, p) 
  S <- Tr + 2 * (X^2 %*% beta.s) + U^2  + rnorm(n2) # S
  
  beta.y = rep(1, p)
  Y <- Tr + 3*(X%*%beta.y) + S + rnorm(n2)
  
  dat.obs <- data.frame(X, Tr, S, Y)
  names(dat.obs) <- c(paste0('X',1:p), 'Tr', 'S', 'Y')
  
  return(list(dat.rct = dat.rct, dat.obs = dat.obs))
}



# ================ X和U的分布都不一样 ============================= 

gen.dat3 <- function(n1, n2, p){
  
  # ---------------  RCT  ---------------- # 
  X <- matrix(rnorm(n1*p), nrow = n1)      # X
  U <- rnorm(n1)                        # U, unmeasured confounder 
  Tr <- rbinom(n1, 1, prob = 0.5)  # T
  
  beta.s = rep(1, p) 
  S <- Tr + 2 * (X %*% beta.s) + U  + rnorm(n1) # S
  
  beta.y = rep(1, p)
  Y <- Tr + 3*(X%*%beta.y) + S + rnorm(n1)
  
  dat.rct <- data.frame(X, Tr, S, Y)
  names(dat.rct) <- c(paste0('X',1:p), 'Tr', 'S', 'Y')
  
  # ---------------  OBS  ---------------- # 
  X <- matrix(rnorm(n2*p, mean = 1, sd = 2), nrow = n2)   # X
  # U <- runif(n2, 1,5)      # U, unmeasured confounder 
  U <- rnorm(n2, 1, 2)      # U, unmeasured confounder
  
  beta.t = rep(1, p)
  prob.t = 1/(1 + exp(-X%*%beta.t - U))
  Tr <- rbinom(n2, 1, prob = prob.t )    # T
  
  beta.s = rep(1, p) 
  S <- Tr + 2 * (X %*% beta.s) + U  + rnorm(n2) # S
  
  beta.y = rep(1, p)
  Y <- Tr + 3*(X%*%beta.y) + S + rnorm(n2)
  
  dat.obs <- data.frame(X, Tr, S, Y)
  names(dat.obs) <- c(paste0('X',1:p), 'Tr', 'S', 'Y')
  
  return(list(dat.rct = dat.rct, dat.obs = dat.obs))
}



gen.dat4 <- function(n1, n2, p){
  
  # ---------------  RCT  ---------------- # 
  X <- matrix(rnorm(n1*p), nrow = n1)      # X
  U <- rnorm(n1)                        # U, unmeasured confounder 
  Tr <- rbinom(n1, 1, prob = 0.5)  # T
  
  beta.s = rep(1, p) 
  S <- Tr + 2 * (X^2 %*% beta.s) + U^2  + rnorm(n1) # S
  
  beta.y = rep(1, p)
  Y <- Tr + 3*(X%*%beta.y) + S + rnorm(n1)
  
  dat.rct <- data.frame(X, Tr, S, Y)
  names(dat.rct) <- c(paste0('X',1:p), 'Tr', 'S', 'Y')
  
  # ---------------  OBS  ---------------- # 
  X <- matrix(rnorm(n2*p, mean = 1, sd = 2), nrow = n2)   # X
  # U <- runif(n2, 1,5)                        # U, unmeasured confounder 
  U <- rnorm(n2, 1, 2)      # U, unmeasured confounder
  
  beta.t = rep(1, p)
  prob.t = 1/(1 + exp(-X%*%beta.t - U))
  Tr <- rbinom(n2, 1, prob = prob.t )    # T
  
  beta.s = rep(1, p) 
  S <- Tr + 2 * (X^2 %*% beta.s) + U^2  + rnorm(n2) # S
  
  beta.y = rep(1, p)
  Y <- Tr + 3*(X%*%beta.y) + S + rnorm(n2)
  
  dat.obs <- data.frame(X, Tr, S, Y)
  names(dat.obs) <- c(paste0('X',1:p), 'Tr', 'S', 'Y')
  
  return(list(dat.rct = dat.rct, dat.obs = dat.obs))
}





# =======  X的分布不一样, U同时影响S与Y ====================


gen.dat5 <- function(n1, n2, p){
  # ---------------  RCT  ---------------- # 
  X <- matrix(rnorm(n1*p), nrow = n1)      # X
  U <- rnorm(n1)                        # U, unmeasured confounder 
  Tr <- rbinom(n1, 1, prob = 0.5)  # T
  
  beta.s = rep(1, p) 
  S <- Tr + 2 * (X %*% beta.s) + U  + rnorm(n1) # S
  
  beta.y = rep(1, p)
  Y <- Tr + 3*(X%*%beta.y) + S + U + rnorm(n1)
  
  dat.rct <- data.frame(X, Tr, S, Y)
  names(dat.rct) <- c(paste0('X',1:p), 'Tr', 'S', 'Y')
  
  # ---------------  OBS  ---------------- # 
  X <- matrix(rnorm(n2*p, mean = 1, sd = 2), nrow = n2)   # X
  U <- rnorm(n2)  # U, unmeasured confounder
  
  beta.t = rep(1, p)
  prob.t = 1/(1 + exp(-X%*%beta.t))
  Tr <- rbinom(n2, 1, prob = prob.t )    # T
  
  beta.s = rep(1, p) 
  S <- Tr + 2 * (X %*% beta.s) + U  + rnorm(n2) # S
  
  beta.y = rep(1, p)
  Y <- Tr + 3*(X%*%beta.y) + S + U + rnorm(n2)
  
  dat.obs <- data.frame(X, Tr, S, Y)
  names(dat.obs) <- c(paste0('X',1:p), 'Tr', 'S', 'Y')
  
  return(list(dat.rct = dat.rct, dat.obs = dat.obs))
}


gen.dat6 <- function(n1, n2, p){
  
  # ---------------  RCT  ---------------- # 
  X <- matrix(rnorm(n1*p), nrow = n1)      # X
  U <- rnorm(n1)                        # U, unmeasured confounder 
  Tr <- rbinom(n1, 1, prob = 0.5)  # T
  
  beta.s = rep(1, p) 
  S <- Tr + 2 * (X^2 %*% beta.s) + U^2  + rnorm(n1) # S
  
  beta.y = rep(1, p)
  Y <- Tr + 3*(X%*%beta.y) + S + U + rnorm(n1)
  
  dat.rct <- data.frame(X, Tr, S, Y)
  names(dat.rct) <- c(paste0('X',1:p), 'Tr', 'S', 'Y')
  
  # ---------------  OBS  ---------------- # 
  X <- matrix(rnorm(n2*p, mean = 1, sd = 2), nrow = n2)   # X
  U <- rnorm(n2)  # U, unmeasured confounder
  
  beta.t = rep(1, p)
  prob.t = 1/(1 + exp(-X%*%beta.t))
  Tr <- rbinom(n2, 1, prob = prob.t )    # T
  
  beta.s = rep(1, p) 
  S <- Tr + 2 * (X^2 %*% beta.s) + U^2  + rnorm(n2) # S
  
  beta.y = rep(1, p)
  Y <- Tr + 3*(X%*%beta.y) + S + U + rnorm(n2)
  
  dat.obs <- data.frame(X, Tr, S, Y)
  names(dat.obs) <- c(paste0('X',1:p), 'Tr', 'S', 'Y')
  
  return(list(dat.rct = dat.rct, dat.obs = dat.obs))
}


# ================== X, S都是离散变量 =========================== 

gen.dat7 <- function(n1, n2, p){
  
  # ---------------  RCT  ---------------- # 
  X <- matrix(rbinom(n1*p, 1, prob = 0.5), nrow = n1)      # X
  U <- rnorm(n1)                        # U, unmeasured confounder 
  Tr <- rbinom(n1, 1, prob = 0.5)  # T
  
  beta.s = rep(1, p) 
  prob.s = 1/(1 + exp(-Tr + 2 * (X %*% beta.s) - U))
  S <- rbinom(n1, 1, prob = prob.s) 
  
  beta.y = rep(1, p)
  Y <- Tr + 3*(X%*%beta.y) + S + rnorm(n1)
  
  dat.rct <- data.frame(X, Tr, S, Y)
  names(dat.rct) <- c(paste0('X',1:p), 'Tr', 'S', 'Y')
  
  # ---------------  OBS  ---------------- # 
  X <- matrix(rbinom(n2*p, 1, prob = 0.5), nrow = n2)   # X
  U <- rnorm(n2)    # U, unmeasured confounder 
  
  beta.t = rep(1, p)
  prob.t = 1/(1 + exp(-X%*%beta.t - U))
  Tr <- rbinom(n2, 1, prob = prob.t )    # T
  
  beta.s = rep(1, p) 
  prob.s = 1/(1 + exp(-Tr + 2 * (X %*% beta.s) - U))
  S <- rbinom(n2, 1, prob = prob.s) 
  
  beta.y = rep(1, p)
  Y <- Tr + 3*(X%*%beta.y) + S + rnorm(n2)
  
  dat.obs <- data.frame(X, Tr, S, Y)
  names(dat.obs) <- c(paste0('X',1:p), 'Tr', 'S', 'Y')
  
  return(list(dat.rct = dat.rct, dat.obs = dat.obs))
}



gen.dat8 <- function(n1, n2, p){
  
  # ---------------  RCT  ---------------- # 
  X <- matrix(rbinom(n1*p, 1, prob = 0.5), nrow = n1)      # X
  U <- rnorm(n1)                        # U, unmeasured confounder 
  Tr <- rbinom(n1, 1, prob = 0.5)  # T
  
  beta.s = rep(1, p) 
  prob.s = 1/(1 + exp(-Tr + 2 * (X^2 %*% beta.s) - U^2))
  S <- rbinom(n1, 1, prob = prob.s) 
  
  beta.y = rep(1, p)
  Y <- Tr + 3*(X%*%beta.y) + S + rnorm(n1)
  
  dat.rct <- data.frame(X, Tr, S, Y)
  names(dat.rct) <- c(paste0('X',1:p), 'Tr', 'S', 'Y')
  
  # ---------------  OBS  ---------------- # 
  X <- matrix(rbinom(n2*p, 1, prob = 0.5), nrow = n2)   # X
  U <- rnorm(n2)                        # U, unmeasured confounder 
  
  beta.t = rep(1, p)
  prob.t = 1/(1 + exp(-X%*%beta.t - U))
  Tr <- rbinom(n2, 1, prob = prob.t )    # T
  
  beta.s = rep(1, p) 
  prob.s = 1/(1 + exp(-Tr + 2 * (X^2 %*% beta.s) - U^2))
  S <- rbinom(n2, 1, prob = prob.s) 
  
  beta.y = rep(1, p)
  Y <- Tr + 3*(X%*%beta.y) + S + rnorm(n2)
  
  dat.obs <- data.frame(X, Tr, S, Y)
  names(dat.obs) <- c(paste0('X',1:p), 'Tr', 'S', 'Y')
  
  return(list(dat.rct = dat.rct, dat.obs = dat.obs))
}

