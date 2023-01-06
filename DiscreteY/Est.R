
Est.bin <- function(dat.rct, dat.obs){

  n1 <- nrow(dat.rct)
  n2 <- nrow(dat.obs)
  n <- n1 + n2
  p <- ncol(dat.obs) - 3 
  
  # ----- propensity score in RCT ----- # 
  ps.rct <- rep(1/2, n1)
  
  # ----- estimate h function: binary outcome  ------ # 
  mod.ho <- glm(Y ~. , data = dat.obs,family = binomial(link = 'logit'))
  h.obs <- unname(mod.ho$fitted.values)
  
  h.rct <- predict.glm(mod.ho, newdata = dat.rct, type = 'response')
  
  # ----- estimate ATE  ------- #
  Tr.rct <- dat.rct$Tr
  ate <- mean( h.rct * (Tr.rct/ps.rct) ) -  
    mean(  h.rct * (1 - Tr.rct)/(1- ps.rct) ) 
  
  # ------ asymptotic variance ------ # 
  X.obs <- as.matrix(cbind(1, dat.obs[, -(p+3)] ))
  X.rct <- as.matrix(cbind(1, dat.rct ))
  
  I_a <-  t(X.obs) %*% diag( h.obs*(1 - h.obs))  %*% X.obs / n2
  B <-  t(h.rct*(1-h.rct)*(Tr.rct - ps.rct)/(ps.rct*(1-ps.rct))) %*%
    X.rct / n1
  
  var.1 <- B %*% solve(I_a) %*% t(B)  / n2 
  var.2 <- var( h.rct * (Tr.rct/ps.rct) - h.rct * (1 - Tr.rct)/(1- ps.rct) ) / n1
  std <- sqrt(var.1 + var.2)
  
  # -------- return the result  ------- # 
  return( list(ate = ate, std = std) )
}


Est.bin.reduce <- function(dat.rct, dat.obs){
  
  n1 <- nrow(dat.rct)
  n2 <- nrow(dat.obs)
  n <- n1 + n2
  p <- ncol(dat.obs) - 3 
  
  # ----- propensity score in RCT: logistic regression ----- # 
  mod.ps <- glm(Tr ~ . - S, family = binomial(link = 'logit'),
                data = dat.rct)
  ps.rct <- unname(mod.ps$fitted.values)
  
  
  # ----- estimate h function: binary outcome  ------ # 
  mod.ho <- glm(Y ~. , data = dat.obs,family = binomial(link = 'logit'))
  h.obs <- unname(mod.ho$fitted.values)
  
  h.rct <- predict.glm(mod.ho, newdata = dat.rct, type = 'response')
  
  # ----- estimate ATE  ------- #
  Tr.rct <- dat.rct$Tr
  ate <- mean( h.rct * (Tr.rct/ps.rct) ) -  
    mean(  h.rct * (1 - Tr.rct)/(1- ps.rct) ) 
  
  # ------ asymptotic variance ------ # 
  X.obs <- as.matrix(cbind(1, dat.obs[, -(p+3)] ))
  X.rct <- as.matrix(cbind(1, dat.rct ))
  
  I_a <-  t(X.obs) %*% diag( h.obs*(1 - h.obs))  %*% X.obs / n2
  B <-  t(h.rct*(1-h.rct)*(Tr.rct - ps.rct)/(ps.rct*(1-ps.rct))) %*%
    X.rct / n1
  
  XX.rct <- as.matrix(cbind(1, dat.rct[, 1:p] ))
  tmp <- Tr.rct*h.rct*(1-ps.rct)/ps.rct + (1-Tr.rct)*h.rct*ps.rct/(1-ps.rct)
  D <- t( tmp ) %*% XX.rct / n1
  I_b <- t(XX.rct) %*% diag( ps.rct*(1-ps.rct) ) %*% XX.rct / n1
  
  var.1 <- B %*% solve(I_a) %*% t(B)  / n2 
  var.2 <- var( h.rct * (Tr.rct/ps.rct) - h.rct * (1 - Tr.rct)/(1- ps.rct)
                - c(D %*% solve(I_b) %*% ( t(XX.rct)  %*% diag(Tr.rct - ps.rct) )) 
                ) / n1
  
  # var.3 <- D %*% solve(I_b) %*% t(D)  / n1
  # std <- sqrt(var.1 + var.2 - var.3)
  std <- sqrt(var.1 + var.2)
  
  # -------- return the result  ------- # 
  return( list(ate = ate, std = std) )
}

