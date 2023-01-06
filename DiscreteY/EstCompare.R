
Est.Susan <- function(dat.rct, dat.obs){

  n1 <- nrow(dat.rct)
  n2 <- nrow(dat.obs)
  n <- n1 + n2
  p <- ncol(dat.obs) - 3 
  
  # ----- propensity score model in RCT ----- # 
  # P(T=1|X=x,G=1), needs to caculate it in both RCT and OBS data.
  ps.rct <- rep(1/2, n1)
  ps.obs <- rep(1/2, n2)   # extrapolate it to OBS data.  
  
  # ----- surrogate score in RCT ------- # 
  # P(T=1|X=x,S=s,G=1), only needs to caculate it in OBS data.
  mod.ss <- glm(Tr ~ ., family = binomial(link = 'logit'),
                data = dat.rct)
  ss.obs <- predict(mod.ss, newdata = dat.obs, type = 'response')
  
  
  # ----- estimate h function: binary outcome  ------ # 
  # h = E(Y|X=x,S=s,G=0), needs to caculate it in both RCT and OBS data.  
  mod.ho <- glm(Y ~. -Tr, data = dat.obs,family = binomial(link = 'logit'))
  h.obs <- unname(mod.ho$fitted.values)
  h.rct <- predict.glm(mod.ho, newdata = dat.rct, type = 'response')
  
  # ----- sampling (weighting) score  ---------- #
  # pi(s,x) = P(G=1|S=s, X=x), only needs to caculate it in OBS data.
  dat.rct$Y <- rep(99, n1)
  
  dat.combine <- rbind(dat.rct, dat.obs)
  dat.combine$G <- c( rep(1, n1), rep(0, n2))
  mod.pi <- glm(G ~ . - Y - Tr, family = binomial(link = 'logit'),
               data = dat.combine)
  pi.obs <- unname(mod.pi$fitted.values)[(n1+1):(n1+n2)]
  
  
  # ----- estimate ATE  ------- #
  q <- n1 / n
  
  # efficient score in RCT data
  Tr.rct <- dat.rct$Tr
  eff.score.rct <- h.rct * (Tr.rct/ps.rct - (1 - Tr.rct)/(1- ps.rct) ) / q
  
  # efficient score in OBS data 
  Y.obs <- dat.obs$Y
  tmp1 <- pi.obs / ((1-pi.obs)*q) 
  tmp2 <- (Y.obs - h.obs)*(pi.obs - ps.obs) / ( pi.obs*(1-pi.obs))
  eff.score.obs <- tmp1 * tmp2
  
  psi <- c(eff.score.rct, eff.score.obs)
  ate <- mean(psi)
  
  # ------ asymptotic variance ------ # 
  std <- sd(psi) / sqrt(n)
  
  # -------- return the result  ------- # 
  return( list(ate = ate, std = std) )
}
