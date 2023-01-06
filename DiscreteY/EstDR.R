# binary outcome 


Est.DR <- function(dat.rct, dat.obs){
  
  n1 <- nrow(dat.rct)
  n2 <- nrow(dat.obs)
  n <- n1 + n2
  p <- ncol(dat.obs) - 3 
  
  # ----- propensity score model in both RCT and OBS ----- # 
  # P(T=1|X=x,G=1)
  ps.rct <- rep(1/2, n1)  
  ps.obs <- rep(1/2, n2)
  
  # ----- propensity score model in OBS ------- # 
  # P(T=1|X=x,S=s,G=0), only needs to caculate it in OBS data.
  # mod.ps.obs <- glm(Tr ~ . -Y, family = binomial(link = 'logit'),
  #                   data = dat.obs)
  # r.obs <- unname(mod.ps.obs$fitted.values)
  
  # ----- estimate mu: continuous outcome  ------ # 
  # mu1 = E(Y|X=x,S=s,T=1,G=1), mu0 = E(Y|X=x,S=s,T=0,G=1) 
  # h1  = E(Y|X=x,T=1,G=1),      h0 = E(Y|X=x,T=0,G=1) 
  # needs to caculate mu1 and mu0 in both RCT and OBS data.  
  # needs to caculate h1 and h0 in RCT data.
  # both of them are modelled with OBS data. 
  
  mod.mu1 <- glm(Y ~. -Tr, family = binomial(link = 'logit'),
                 data = dat.obs, subset = Tr == 1)
  mod.mu0 <- glm(Y ~. -Tr, family = binomial(link = 'logit'),
                 data = dat.obs, subset = Tr == 0)
  
  mu1.obs <- predict.glm(mod.mu1, newdata = dat.obs, type = 'response')
  mu1.rct <- predict.glm(mod.mu1, newdata = dat.rct, type = 'response')
  mu0.obs <- predict.glm(mod.mu0, newdata = dat.obs, type = 'response')
  mu0.rct <- predict.glm(mod.mu0, newdata = dat.rct, type = 'response')
  
  mod.h1 <- lm(mu1.rct ~. -Tr-S, data = dat.rct)
  mod.h0 <- lm(mu0.rct ~. -Tr-S, data = dat.rct)
  h1.rct <- predict.lm(mod.h1, newdata = dat.rct) 
  h0.rct <- predict.lm(mod.h0, newdata = dat.rct)
  # h1.rct <- ifelse( h1.rct <= 0, 0, ifelse(h1.rct >= 1, 1, h1.rct))
  # h0.rct <- ifelse( h0.rct <= 0, 0, ifelse(h0.rct >= 1, 1, h0.rct))

  # ----- sampling (weighting) score  ---------- #
  # gt(s,x) = P(G=1|X=x,S=s, T=t), only needs to calculate it in OBS data.
  dat.rct$Y <- rep(99, n1) 
  
  dat.combine <- rbind(dat.rct, dat.obs)
  dat.combine$G <- c( rep(1, n1), rep(0, n2))
  mod.g <- glm(G ~ . - Y, family = binomial(link = 'logit'),
               data = dat.combine)
  g.obs <- unname(mod.g$fitted.values)[(n1+1):(n1+n2)]
  
  
  # ----- estimate ATE  ------- #
  q <- n1 / n
  
  # efficient score in RCT data
  Tr.rct <- dat.rct$Tr
  term1 <- Tr.rct * (mu1.rct - h1.rct) / ps.rct 
  term2 <- (1-Tr.rct) * (mu0.rct - h0.rct) / (1-ps.rct)
  term3 <- h1.rct - h0.rct
  eff.score.rct <- (term1 - term2 + term3) / q
  
  # efficient score in OBS data 
  Y.obs <- dat.obs$Y
  Tr.obs <- dat.obs$Tr
  tmp1 <- g.obs / ((1-g.obs)*q) 
  tmp2 <- Tr.obs*(Y.obs - mu1.obs) / ps.obs 
  tmp3 <- (1-Tr.obs)*(Y.obs-mu0.obs) / (1- ps.obs)
  eff.score.obs <- tmp1 * ( tmp2 - tmp3 )
  
  psi <- c(eff.score.rct, eff.score.obs)
  ate <- mean(psi)
  
  # ------ asymptotic variance ------ # 
  std <- sd(psi) / sqrt(n)
  
  # -------- return the result  ------- # 
  return( list(ate = ate, std = std) )
}

