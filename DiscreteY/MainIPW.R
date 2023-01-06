
rm(list = ls( ))
library(dplyr)

# ---------------- Functions for generating dataset ------------ # 

source('GenData.R')
GenData <- list(gen.dat1, gen.dat2, gen.dat3, gen.dat4,
                gen.dat5, gen.dat6, gen.dat7, gen.dat8)

# ---------------- True values of ATE  ------------------------- # 

set.seed(12)
true.ate = rep(NA, 8) 
for(ss in 1:8){   
  DAT <- GenData[[ss]](100000, 20, p = 2)
  dat.rct <- DAT$dat.rct    
  Tr <- dat.rct$Tr
  Y <- dat.rct$Y
  true.ate[ss] <- mean(Y[Tr == 1]) - mean(Y[Tr == 0] )
}

# ------------------  Simulation   ------------------------------#   

# Function for estimating ATE
source('Est.R')
# seed <- round(runif(1, 1, 1000))
seed = 755
set.seed(seed)

# Simulation 
for(ss in 1:8){     
  
  result <- data.frame(matrix(nrow = 3, ncol= 8))
  names(result) <- c('size.rct','size.obs', 'bias','std', 
                     'ese', 'CP95', 'ese.boot','CP95.boot')
  
  ate0 <- true.ate[ss]
  
  for(j in 1:3){
    
    n1 <- c(50, 100, 200)[j]
    n2 <- 500
    
    PARA <- rep(NA, 1000)  # used to save 1000 results of simulation 
    ESE <- rep(NA, 1000)
    COUNT95 <- rep(NA, 1000)
    ESE.boot <- rep(NA, 1000)
    COUNT95.boot <- rep(NA, 1000)
    cat('------ case ',ss, '-- n1 = ',n1,'-- n2 = ',n2, '------','\n')
    
    for(loop in 1:1000){
      # --------- generate data --------- #
      DAT <- GenData[[ss]](n1, n2, p = 2)
      dat.rct <- DAT$dat.rct     
      dat.obs <- DAT$dat.obs
      dat.rct$Y <- NULL       # mask Y for RCT data 
      
      # --------- estimate ATE and the associated std  -------- # 
      res <- Est.bin(dat.rct, dat.obs)
      
      # save the result 
      para <- res$ate
      ese <- res$std
      count95 <- (ate0 <= para + 1.96*ese) &  (ate0 >= para - 1.96*ese)
      
      PARA[loop] <- para
      ESE[loop] <- ese
      COUNT95[loop] <- count95
      
      # bootstrap
      B <- 50
      ate.tmp <- rep(NA, B)
      for(b in 1:B){
        dat.rct.tmp <- dat.rct[sample(1:n1, replace = T),]
        dat.obs.tmp <- dat.obs[sample(1:n2, replace = T),]
        res <- Est.bin(dat.rct.tmp, dat.obs.tmp)
        ate.tmp[b] <- res$ate
      }
      ese.boot <- sd(ate.tmp)
      count95.boot <- (ate0 <= para + 1.96*ese.boot) &  (ate0 >= para - 1.96*ese.boot)
      ESE.boot[loop] <- ese.boot
      COUNT95.boot[loop] <- count95.boot
      
      cat(loop, '\r')
      
    }  
    

    bias <- mean(PARA) - ate0
    std <- sd(PARA)
    ese <- mean(ESE)
    CP95 <- mean(COUNT95)
    ese.boot <- mean(ESE.boot)
    CP95.boot <- mean(COUNT95.boot)
    result[j, ] <- c(n1, n2, bias, std, ese, CP95, ese.boot, CP95.boot)
  }
  write.csv(result,
            file = paste0('result.IPW.case',ss,'-seed',seed,'.csv') )
}
