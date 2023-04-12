library(magrittr)
library(BOIN)
library(scales)
library(dfcomb)

# posterior probability of pj >= phi given data
post.prob.fn <- function(phi, y, n, alp.prior=0.1, bet.prior=0.1){
  alp <- alp.prior + y 
  bet <- bet.prior + n - y
  1 - pbeta(phi, alp, bet)
}


prob.int <- function(phi, y1, n1, y2, n2, alp.prior, bet.prior){
  alp1 <- alp.prior + y1
  alp2 <- alp.prior + y2
  bet1 <- alp.prior + n1 - y1
  bet2 <- alp.prior + n2 - y2
  # message(paste('y1,y2,n1,n2: ',y1, n1, y2, n2))
  fn.min <- function(x){
    dbeta(x, alp1, bet1)*(1-pbeta(x, alp2, bet2))
  }
  fn.max <- function(x){
    pbeta(x, alp1, bet1)*dbeta(x, alp2, bet2)
  }
  const.min <- integrate(fn.min, lower=0, upper=1)$value
  const.max <- integrate(fn.max, lower=0, upper=1)$value
  p1 <- integrate(fn.min, lower=0, upper=phi)$value/const.min
  p2 <- integrate(fn.max, lower=0, upper=phi)$value/const.max
  
  list(p1=p1, p2=p2)
}


OR.values <- function(phi, y1, n1, y2, n2, alp.prior, bet.prior, type){
  ps <- prob.int(phi, y1, n1, y2, n2, alp.prior, bet.prior)
  if (type=="L"){
    pC <- 1 - ps$p2
    pL <- 1 - ps$p1
    oddsC <- pC/(1-pC)
    oddsL <- pL/(1-pL)
    OR <- oddsC*oddsL
  }else if (type=="R"){
    pC <- 1 - ps$p1
    pR <- 1 - ps$p2
    oddsC <- pC/(1-pC)
    oddsR <- pR/(1-pR)
    OR <- (1/oddsC)/oddsR
  }else if (type=="D"){
    pC <- 1 - ps$p2
    pD <- 1 - ps$p1
    oddsC <- pC/(1-pC)
    oddsD <- pD/(1-pD)
    OR <- oddsC*oddsD
  }else if (type=="U"){
    pC <- 1 - ps$p1
    pU <- 1 - ps$p2
    oddsC <- pC/(1-pC)
    oddsU <- pU/(1-pU)
    OR <- (1/oddsC)/oddsU
  }
  return(OR)
}



All.OR.table <- function(phi, n1, n2, type, alp.prior, bet.prior){
  ret.mat <- matrix(rep(0, (n1+1)*(n2+1)), nrow=n1+1)
  for (y1cur in 0:n1){
    for (y2cur in 0:n2){
      ret.mat[y1cur+1, y2cur+1] <- OR.values(phi, y1cur, n1, y2cur, n2, alp.prior, bet.prior, type)
    }
  }
  ret.mat
}

# compute the marginal prob when lower < phiL/phiC/phiR < upper
# i.e., Pr(Y=y|lower<phi<upper)
margin.phi <- function(y, n, lower, upper){
  C <- 1/(upper-lower)
  fn <- function(phi) {
    dbinom(y, n, phi)*C
  }
  integrate(fn, lower=lower, upper=upper)$value
}

# Obtain the table of marginal distribution of (y1, y2)
# after intergrate out (phi1, phi2)
# under H0 and H1
# H0: phi1=phi, phi < phi2 < 2phi
# H1: phi2=phi, 0   < phi1 < phi
margin.ys.table <- function(n1, n2, phi, hyperthesis){
  if (hyperthesis=="H0"){
    p.y1s <- dbinom(0:n1, n1, phi)
    p.y2s <- sapply(0:n2, margin.phi, n=n2, lower=phi, upper=2*phi)
  }else if (hyperthesis=="H1"){
    p.y1s <- sapply(0:n1, margin.phi, n=n1, lower=0, upper=phi)
    p.y2s <- dbinom(0:n2, n2, phi)
  }
  p.y1s.mat <- matrix(rep(p.y1s, n2+1), nrow=n1+1)
  p.y2s.mat <- matrix(rep(p.y2s, n1+1), nrow=n1+1, byrow=TRUE)
  margin.ys <- p.y1s.mat * p.y2s.mat
  margin.ys
}

# Obtain the optimal gamma for the hypothesis test
optim.gamma.fn <- function(n1, n2, phi, type, alp.prior, bet.prior){
  # message('calculating OR table')
  OR.table <- All.OR.table(phi, n1, n2, type, alp.prior, bet.prior) 
  # message('calculation done')
  ys.table.H0 <- margin.ys.table(n1, n2, phi, "H0")
  ys.table.H1 <- margin.ys.table(n1, n2, phi, "H1")
  
  argidx <- order(OR.table)
  sort.OR.table <- OR.table[argidx]
  sort.ys.table.H0 <- ys.table.H0[argidx]
  sort.ys.table.H1 <- ys.table.H1[argidx]
  n.tol <- length(sort.OR.table)
  
  if (type=="L"){
    errs <- rep(0, n.tol-1)
    for (i in 1:(n.tol-1)){
      err1 <- sum(sort.ys.table.H0[1:i])
      err2 <- sum(sort.ys.table.H1[(i+1):n.tol])
      err <- err1 + err2
      errs[i] <- err
    }
    min.err <- min(errs)
    if (min.err > 1){
      gam <- 0
      min.err <- 1
    }else {
      minidx <- which.min(errs)
      gam <- sort.OR.table[minidx]
    }
  }else if (type=='R'){
    errs <- rep(0, n.tol-1)
    for (i in 1:(n.tol-1)){
      err1 <- sum(sort.ys.table.H1[1:i])
      err2 <- sum(sort.ys.table.H0[(i+1):n.tol])
      err <- err1 + err2
      errs[i] <- err
    }
    min.err <- min(errs)
    if (min.err > 1){
      gam <- 0
      min.err <- 1
    }else {
      minidx <- which.min(errs)
      gam <- sort.OR.table[minidx]
    }
    
  }
  list(gamma=gam, min.err=min.err)
}

make.decision.1dCFO.fn <- function(phi, cys, cns, alp.prior, bet.prior, cover.doses, diag=FALSE){
  if (cover.doses[2] == 1){
    return(1)
  }else{
    if (is.na(cys[1]) & (cover.doses[3]==1)){
      return(2)
    }else  if (is.na(cys[1]) & (!(cover.doses[3]==1))){
      gam2 <- optim.gamma.fn(cns[2], cns[3], phi, "R", alp.prior, bet.prior)$gamma
      message(paste('gam2: ', gam2))
      OR.v2 <- OR.values(phi, cys[2], cns[2], cys[3], cns[3], alp.prior, bet.prior, type="R")
      message(paste('OR.v2: ', OR.v2))
      if (OR.v2>gam2){
        return(3)
      }else{
        return(2)
      }
    }else  if (is.na(cys[3]) | (cover.doses[3]==1)){
      gam1 <- optim.gamma.fn(cns[1], cns[2], phi, "L", alp.prior, bet.prior)$gamma
      message(paste('gam1: ', gam1))
      OR.v1 <- OR.values(phi, cys[1], cns[1], cys[2], cns[2], alp.prior, bet.prior, type="L")
      message(paste('OR.v1: ', OR.v1))
      if (OR.v1>gam1){
        return(1)
      }else{
        return(2)
      }

    }else  if (!(is.na(cys[1]) | is.na(cys[3]) | cover.doses[3]==1)){
      gam1 <- optim.gamma.fn(cns[1], cns[2], phi, "L", alp.prior, bet.prior)$gamma
      gam2 <- optim.gamma.fn(cns[2], cns[3], phi, "R", alp.prior, bet.prior)$gamma
      message(paste('gam1,gam2: ',gam1,gam2))
      OR.v1 <- OR.values(phi, cys[1], cns[1], cys[2], cns[2], alp.prior, bet.prior, type="L")
      OR.v2 <- OR.values(phi, cys[2], cns[2], cys[3], cns[3], alp.prior, bet.prior, type="R")
      message(paste('OR.v1,OR.v2: ',OR.v1,OR.v2))
      v1 <- OR.v1 > gam1
      v2 <- OR.v2 > gam2
      if (v1 & !v2){
        return(1)
      }else if (!v1 & v2){
        return(3)
      }else{
        return(2)
      }
    }
  }
}



make.decision.2dCFO.fn <- function(phi, cys, cns, alp.prior, bet.prior, cover.doses, diag=FALSE){
  cidx.A <- 0
  cidx.B <- 0
  # horizontal direction
  # message(paste('cys[2,]: ', cys[2,]))
  # message(paste('cns[2,]: ', cns[2,]))
  idx.chg.A <- make.decision.1dCFO.fn(phi, cys[2,], cns[2,], alp.prior, bet.prior, cover.doses[2,]) - 2
  message(paste('decision A: ',idx.chg.A))
  # vertical direction
  # message(paste('cys[,2]: ', cys[,2]))
  # message(paste('cns[,2]: ', cns[,2]))
  # message(paste('alp, bet: ', alp.prior, bet.prior))
  # message(paste('cover.doses: ',cover.doses[,2]))
  idx.chg.B <- make.decision.1dCFO.fn(phi, cys[,2], cns[,2], alp.prior, bet.prior, cover.doses[,2]) - 2
  message(paste('decision B: ',idx.chg.B))
  
  if (idx.chg.A == 1 & idx.chg.B == 1){
    ### horizontal and vertical only
    
    OR.R <- OR.values(phi, cys[2,2], cns[2,2], cys[2,3], cns[2,3], alp.prior, bet.prior, type="R")
    OR.U <- OR.values(phi, cys[2,2], cns[2,2], cys[3,2], cns[3,2], alp.prior, bet.prior, type="R")
    
    message(paste('OR.R: ',OR.R))
    message(paste('OR.U: ',OR.U))
    
    if (OR.R == OR.U){
      rand <- rbinom(1,1,0.5)
      if(rand == 0){
        cidx.A <- 1
      } else {
        cidx.B <- 1
      }
    } else if (OR.R > OR.U){
      cidx.B <- 1
    } else {
      cidx.A <- 1
    }
    
    ### diagonal direction
    # cidx.A <- idx.chg.A
    # cidx.B <- idx.chg.B
    
  } else if (idx.chg.A == -1 & idx.chg.B == -1){
    if (is.na(cys[2,1]) & is.na(cys[1,2])){
      cidx.A <- 0
      cidx.B <- 0
    } else if (is.na(cys[2,1])){
      cidx.A <- -1
    } else if (is.na(cys[1,2])){
      cidx.B <- -1
    } else {
      OR.L <- OR.values(phi, cys[2,2], cns[2,2], cys[2,1], cns[2,1], alp.prior, bet.prior, type="L")
      OR.D <- OR.values(phi, cys[2,2], cns[2,2], cys[1,2], cns[1,2], alp.prior, bet.prior, type="L")
      
      message(paste('OR.L: ',OR.L))
      message(paste('OR.D: ',OR.D))
      
      if (OR.L == OR.D){
        rand <- rbinom(1,1,0.5)
        if(rand == 0){
          cidx.A <- -1
        } else {
          cidx.B <- -1
        }
      } else if (OR.L > OR.D){
        cidx.B <- -1
      } else {
        cidx.A <- -1
      }
    }
  } else if (idx.chg.A == 1 & idx.chg.B == -1){
    # message('rare case occurs')
    DCR <- make.decision.1dCFO.fn(phi, c(cys[1,2],cys[2,2],cys[2,3]), c(cns[1,2],cns[2,2],cns[2,3]), alp.prior, 
                                  bet.prior, c(cover.doses[1,2],cover.doses[2,2],cover.doses[2,3])) - 2
    if (DCR == 1){
      cidx.B <- 1
    } else if (DCR == -1){
      cidx.A <- -1
    }
  } else if (idx.chg.A == -1 & idx.chg.B == 1){
    # message('rare case occurs')
    LCU <- make.decision.1dCFO.fn(phi, c(cys[2,1],cys[2,2],cys[3,2]), c(cns[2,1],cns[2,2],cns[3,2]), alp.prior, 
                                  bet.prior, c(cover.doses[2,1],cover.doses[2,2],cover.doses[3,2])) - 2
    if (LCU == 1){
      cidx.A <- 1
    } else if (LCU == -1){
      cidx.B <- -1
    }
  } else if (idx.chg.A == 1 & idx.chg.B == 0){
    cidx.B <- 1
  } else if (idx.chg.A == 0 & idx.chg.B == 1){
    cidx.A <- 1
  } else if (idx.chg.A == -1 & idx.chg.B == 0){
    cidx.B <- -1
  } else if (idx.chg.A == 0 & idx.chg.B == -1){
    cidx.A <- -1
  }
  
  return (c(cidx.A, cidx.B))

}
  

overdose.fn <- function(phi, add.args=list()){
  y <- add.args$y
  n <- add.args$n
  alp.prior <- add.args$alp.prior
  bet.prior <- add.args$bet.prior
  pp <- post.prob.fn(phi, y, n, alp.prior, bet.prior)
  # message(pp)
  if ((pp >= 0.95) & (add.args$n>=3)){
    return(TRUE)
  }else{
    return(FALSE)
  }
  
}
  
# preliminary <- function(p.true, ndose.A, ndose.B, tys, tns, seed=seed){
#   set.seed(seed)
#   for (i in 1:ndose.A) {
#     tns[i,1] <- tns[i,1] + 1
#     p <- runif(1)
#     if (p<p.true[i,1]){
#       tys[i,1] <- 1
#       return(list(tys=tys, tns=tns, position=which(tys==1, arr.ind = T)))
#     }
#   }
#   
#   for (j in 2:ndose.B) {
#     tns[1,j] <- tns[1,j] + 1
#     p <- runif(1)
#     if (p<p.true[1,j]){
#       tys[1,j] <- 1
#       return(list(tys=tys, tns=tns, position=which(tys==1, arr.ind = T)))
#     }
#   }
#   
#   return(list(tys=tys, tns=tns, position=c(1,1)))
# }
  
preliminary <- function(p.true, ndose.A, ndose.B, tys, tns, seed=seed){
  set.seed(seed)
  for (i in 1:ndose.A) {
    tns[i,1] <- tns[i,1] + 1
    p <- runif(1)
    if (p<p.true[i,1]){
      tys[i,1] <- 1
      position <- which(tys==1, arr.ind = T)
      if(position[1]!=1){
        position[1] <- position[1] - 1
      }
      return(list(tys=tys, tns=tns, position=position))
    }
  }
  
  for (j in 2:ndose.B) {
    tns[1,j] <- tns[1,j] + 1
    p <- runif(1)
    if (p<p.true[1,j]){
      tys[1,j] <- 1
      position <- which(tys==1, arr.ind = T)
      position[2] <- position[2] - 1
      return(list(tys=tys, tns=tns, position=position))
    }
  }
  
  return(list(tys=tys, tns=tns, position=c(1,1)))
}
  


line.sampling <- function(phi, ndose, tys, tns, n, lower, upper){
  samples <- as.list(numeric(ndose))
  p <- replicate(ndose, 0)
  nsamples <- 0
  while (nsamples<n) {
    s <- replicate(ndose, 0)
    for (i in 1:ndose) {
      s[i] <- rbeta(1,phi+tys[i],1-phi+tns[i]-tys[i])
    }
    if (all(diff(s) >= 0)){
      for (j in 1:ndose) {
        samples[[j]] <- c(samples[[j]],s[j])
      }
      nsamples <- nsamples + 1
    } else {
      next
    }
  }
  for (k in 1:ndose) {
    p[k] <- length(samples[[k]][samples[[k]]>lower & samples[[k]]<upper])/length(samples[[k]])
  }
  return(p)
}


cross.sampling <- function(phi, ndose.A, ndose.B, tys, tns, n, lower, upper){
  p <- matrix(0,ndose.A,ndose.B)
  for (i in 1:ndose.A) {
    p[i,] <- p[i,] + line.sampling(phi, ndose.B, tys[i,], tns[i,], n, lower, upper)
  }
  for (j in 1:ndose.B) {
    p[,j] <- p[,j] + line.sampling(phi, ndose.A, tys[,j], tns[,j], n, lower, upper)
  }
  return(p)
}



# Simulation function for CFO
CFO.simu.fn <- function(phi, p.true, prelim=0, ncohort=20, cohortsize=1, init.level.A=1, init.level.B=1, add.args=list(), seed=NULL){
  # phi: Target DIL rate
  # p.true: True DIL rates under the different dose levels
  # ncohort: The number of cohorts
  # cohortsize: The sample size in each cohort
  # alp.prior, bet.prior: prior parameters
  set.seed(seed)
  earlystop <- 0
  ndose.A <- length(p.true[,1])
  ndose.B <- length(p.true[1,])
  cidx.A <- init.level.A
  cidx.B <- init.level.B
  
  tys <- matrix(0, ndose.A, ndose.B) # number of responses for different doses.
  tns <- matrix(0, ndose.A, ndose.B) # number of subject for different doses.
  tover.doses <- matrix(0, ndose.A, ndose.B) # Whether each dose is overdosed or not, 1 yes
  cohorts.left <- ncohort
  if(prelim==1){
    pre <- preliminary(p.true, ndose.A, ndose.B, tys, tns, seed)
    cidx.A <- pre$position[1]
    cidx.B <- pre$position[2]
    tys <- pre$tys
    tns <- pre$tns
    if(((ncohort*cohortsize) - sum(tns))%%3==0){
      cohorts.left <- ((ncohort*cohortsize) - sum(tns))%/%3
    } else {
      cohorts.left <- ((ncohort*cohortsize) - sum(tns))%/%3 + 1
    }
  }
  
  
  
  for (i in 1:cohorts.left){
    message(paste(i, '-th cohort:'))
    message(paste('cidx (A,B): (', cidx.A, ',', cidx.B, ')'))
    
    pc <- p.true[cidx.A, cidx.B] 
    
    # sample from current dose
    
    
    if(prelim==1){
      if(i==cohorts.left){
        cres <- rbinom((ncohort*cohortsize - sum(tns)), 1, pc)
        tys[cidx.A, cidx.B] <- tys[cidx.A, cidx.B] + sum(cres)
        tns[cidx.A, cidx.B] <- tns[cidx.A, cidx.B] + (ncohort*cohortsize - sum(tns))
      } else {
        cres <- rbinom(cohortsize, 1, pc)
        tys[cidx.A, cidx.B] <- tys[cidx.A, cidx.B] + sum(cres)
        tns[cidx.A, cidx.B] <- tns[cidx.A, cidx.B] + cohortsize
      }
    } else {
      cres <- rbinom(cohortsize, 1, pc)
      # message(paste('#DLT: ', sum(cres)))
      tys[cidx.A, cidx.B] <- tys[cidx.A, cidx.B] + sum(cres)
      tns[cidx.A, cidx.B] <- tns[cidx.A, cidx.B] + cohortsize
    }
    
    cy <- tys[cidx.A, cidx.B]
    cn <- tns[cidx.A, cidx.B]
    
    
    add.args <- c(list(y=cy, n=cn, tys=tys, tns=tns, cidx.A=cidx.A, cidx.B=cidx.B), add.args)
    
    # quick escalation
    # if (sum(tys)==0){
    #   if (cidx.A == ndose.A & cidx.B != ndose.B){
    #     cidx.B <- cidx.B +1
    #   } else if (cidx.A != ndose.A & cidx.B == ndose.B){
    #     cidx.A <- cidx.A +1
    #   } else if (cidx.A != ndose.A & cidx.B != ndose.B){
    #     rand <- rbinom(1,1,0.5)
    #     if(rand == 0){
    #       cidx.A <- cidx.A + 1
    #     } else {
    #       cidx.B <- cidx.B + 1
    #     }
    #   } else {
    #     next
    #   }
    #   next
    # }
    
    if (overdose.fn(phi, add.args)){
      tover.doses[cidx.A:ndose.A, cidx.B:ndose.B] <- 1
    }

    
    if (tover.doses[1,1] == 1){
      earlystop <- 1
      break()
    }
    
    if (cidx.A!=1 & cidx.B!=1 & cidx.A!=ndose.A & cidx.B!=ndose.B){
      # no boundary
      cys <- tys[(cidx.A-1):(cidx.A+1), (cidx.B-1):(cidx.B+1)]
      cns <- tns[(cidx.A-1):(cidx.A+1), (cidx.B-1):(cidx.B+1)]
      cover.doses <- tover.doses[(cidx.A-1):(cidx.A+1), (cidx.B-1):(cidx.B+1)]
    } else if (cidx.A==1 & cidx.B==1){
      # (1, 1)
      cys <- rbind(c(NA,NA,NA),cbind(c(NA,NA),tys[1:2,1:2]))
      cns <- rbind(c(NA,NA,NA),cbind(c(NA,NA),tns[1:2,1:2]))
      cover.doses <- rbind(c(NA,NA,NA),cbind(c(NA,NA),tover.doses[1:2,1:2]))
    } else if (cidx.A==ndose.A & cidx.B==ndose.B){
      # (nA, nB)
      cys <- rbind(cbind(tys[(cidx.A-1):cidx.A,(cidx.B-1):cidx.B],c(NA,NA)), c(NA,NA,NA))
      cns <- rbind(cbind(tns[(cidx.A-1):cidx.A,(cidx.B-1):cidx.B],c(NA,NA)), c(NA,NA,NA))
      cover.doses <- rbind(cbind(tover.doses[(cidx.A-1):cidx.A,(cidx.B-1):cidx.B],c(NA,NA)), c(NA,NA,NA))
    } else if (cidx.A==1 & cidx.B==ndose.B){
      # (1, nB) 
      cys <- rbind(c(NA,NA,NA),cbind(tys[1:2,(cidx.B-1):cidx.B],c(NA,NA)))
      cns <- rbind(c(NA,NA,NA),cbind(tns[1:2,(cidx.B-1):cidx.B],c(NA,NA)))
      cover.doses <- rbind(c(NA,NA,NA),cbind(tover.doses[1:2,(cidx.B-1):cidx.B],c(NA,NA)))
    } else if (cidx.A==ndose.A & cidx.B==1){
      # (nA, 1) 
      cys <- rbind(cbind(c(NA,NA), tys[(cidx.A-1):cidx.A,1:2]),c(NA,NA,NA))
      cns <- rbind(cbind(c(NA,NA), tns[(cidx.A-1):cidx.A,1:2]),c(NA,NA,NA))
      cover.doses <- rbind(cbind(c(NA,NA), tover.doses[(cidx.A-1):cidx.A,1:2]),c(NA,NA,NA))
    } else if (cidx.A==1 & cidx.B!=1){
      # (1, 2:(nB-1))
      cys <- rbind(c(NA,NA,NA), tys[1:2, (cidx.B-1):(cidx.B+1)])
      cns <- rbind(c(NA,NA,NA), tns[1:2, (cidx.B-1):(cidx.B+1)])
      cover.doses <- rbind(c(NA,NA,NA), tover.doses[1:2, (cidx.B-1):(cidx.B+1)])
    } else if (cidx.A!=1 & cidx.B==1){
      # (2:(nA-1), 1)
      cys <- cbind(c(NA,NA,NA), tys[(cidx.A-1):(cidx.A+1), 1:2])
      cns <- cbind(c(NA,NA,NA), tns[(cidx.A-1):(cidx.A+1), 1:2])
      cover.doses <- cbind(c(NA,NA,NA), tover.doses[(cidx.A-1):(cidx.A+1), 1:2])
    } else if (cidx.A==ndose.A & cidx.B!=ndose.B){
      # (nA, 2:(nB-1))
      cys <- rbind(tys[(ndose.A-1):ndose.A, (cidx.B-1):(cidx.B+1)], c(NA,NA,NA))
      cns <- rbind(tns[(ndose.A-1):ndose.A, (cidx.B-1):(cidx.B+1)], c(NA,NA,NA))
      cover.doses <- rbind(tover.doses[(ndose.A-1):ndose.A, (cidx.B-1):(cidx.B+1)], c(NA,NA,NA))
    } else if (cidx.A!=ndose.A & cidx.B==ndose.B){
      # (2:(nA-1), nB)
      cys <- cbind(tys[(cidx.A-1):(cidx.A+1), (cidx.B-1):cidx.B], c(NA,NA,NA))
      cns <- cbind(tns[(cidx.A-1):(cidx.A+1), (cidx.B-1):cidx.B], c(NA,NA,NA))
      cover.doses <- cbind(tover.doses[(cidx.A-1):(cidx.A+1), (cidx.B-1):cidx.B], c(NA,NA,NA))
    } else {
      message('no such case')
    }

    idx.chg <- make.decision.2dCFO.fn(phi, cys, cns, add.args$alp.prior, add.args$bet.prior, cover.doses)
    cidx.A <- cidx.A + idx.chg[1]
    cidx.B <- cidx.B + idx.chg[2]
  }
  
  est_p <- select.mtd.comb(phi, tns, tys, boundMTD=TRUE)
  if (earlystop==0){
     MTD <- select.mtd.comb(phi, tns, tys)$MTD
    # p <- cross.sampling(phi, ndose.A, ndose.B, tys, tns, n=100, lower=0.2, upper=0.4)
    # MTD <- which(p==max(p),arr.ind = T)[1,]
  }else{
    #p <- matrix(0,ndose.A,ndose.B)
    MTD <- c(99,99)
  }
  
  correct <- 0
  if(MTD[1]>ndose.A | MTD[2]>ndose.B){
    correct <- 0
  } else if (length(MTD)!=2){
    correct <- 0
  }else if (p.true[MTD[1],MTD[2]]==phi){
    correct <- 1
  }

  npercent <- 0
  for (j in 1:ndose.A) {
    for (k in 1:ndose.B) {
      if (p.true[j,k]==phi){
        npercent <- npercent + tns[j,k]
      }
    }
  }
  npercent <- percent(npercent/(ncohort*cohortsize))
  
  ptoxic <- 0
  for (j in 1:ndose.A) {
    for (k in 1:ndose.B) {
      if (p.true[j,k]>phi){
        ptoxic <- ptoxic + tns[j,k]
      }
    }
  }
  ptoxic <- percent(ptoxic/(ncohort*cohortsize))
  list(MTD=MTD, dose.ns=tns, DLT.ns=tys, p.true=p.true, target=phi, over.doses=tover.doses, correct=correct, npercent=npercent, ptoxic=ptoxic, ntox=sum(tys), est_p=est_p)
}






select.mtd.comb <- function (target, npts, ntox, cutoff.eli = 0.95, extrasafe = FALSE,
                             offset = 0.05, boundMTD=FALSE, p.tox=1.4*target,mtd.contour = FALSE)
{
  lambda_d = log((1 - target)/(1 - p.tox))/log(p.tox * (1 -target)/(target * (1 - p.tox)))
  y = ntox
  n = npts
  if (nrow(n) > ncol(n) | nrow(y) > ncol(y)) {
    stop("npts and ntox should be arranged in a way (i.e., rotated) such that for each of them, the number of rows is less than or equal to the number of columns.")
    
  }
  elimi = matrix(0, dim(n)[1], dim(n)[2])
  if (extrasafe) {
    if (n[1, 1] >= 3) {
      if (1 - pbeta(target, y[1, 1] + 1, n[1, 1] - y[1,
                                                     1] + 1) > cutoff.eli - offset) {
        elimi[, ] = 1
      }
    }
  }
  for (i in 1:dim(n)[1]) {
    for (j in 1:dim(n)[2]) {
      if (n[i, j] >= 3) {
        if (1 - pbeta(target, y[i, j] + 1, n[i, j] -
                      y[i, j] + 1) > cutoff.eli) {
          elimi[i:dim(n)[1], j] = 1
          elimi[i, j:dim(n)[2]] = 1
          break
        }
      }
    }
  }
  
  selectdose=NULL
  
  if (elimi[1] == 1) {
    selectdose = c(99, 99)
    selectdoses = matrix(selectdose, nrow = 1)
  }else {
    phat = (y + 0.05)/(n + 0.1)
    phat = round(Iso::biviso(phat, n + 0.1, warn = TRUE)[, ],2)
    # phat.out = phat
    lower.mat=qbeta(0.025,y+0.05,n-y+0.05)
    lower.mat=round(Iso::biviso(lower.mat),2)
    
    upper.mat=qbeta(0.975,y+0.05,n-y+0.05)
    upper.mat=round(Iso::biviso(upper.mat),2)
    phat.out<-matrix(paste0(format(phat,digits=1),"(",lower.mat,", ",upper.mat,")"),byrow=FALSE,nrow=dim(phat)[1])
    colnames(phat.out)=paste0("B",1:dim(n)[2])
    rownames(phat.out)=paste0("A",1:dim(n)[1])
    phat.out.noCI=round(phat,2)
    phat.out[n == 0] = "NA"
    phat[elimi == 1] = 1.1
    phat = phat * (n != 0) + (1e-05) * (matrix(rep(1:dim(n)[1],
                                                   each = dim(n)[2], len = length(n)), dim(n)[1], byrow = T) +
                                          matrix(rep(1:dim(n)[2], each = dim(n)[1], len = length(n)),
                                                 dim(n)[1]))
    
    if(boundMTD){
      if(all(phat[n!=0]>lambda_d)){
        selectdose = c(99, 99)
        selectdoses = matrix(selectdose, nrow = 1)
      }else{
        phat[phat>lambda_d]=10}}
    
    if(is.null(selectdose)){
      phat[n == 0] = 10
      selectdose = which(abs(phat - target) == min(abs(phat -
                                                         target)), arr.ind = TRUE)
      
      
      if (length(selectdose) > 2)
        selectdose = selectdose[1, ]
      aa = function(x) as.numeric(as.character(x))
      if (mtd.contour == TRUE) {
        selectdoses = cbind(row = 1:dim(n)[1], col = rep(99,
                                                         dim(n)[1]))
        for (k in dim(n)[1]:1) {
          kn = n[k, ]
          ky = y[k, ]
          kelimi = elimi[k, ]
          kphat = phat[k, ]
          if (kelimi[1] == 1 || sum(n[kelimi == 0]) ==
              0) {
            kseldose = 99
          }else {
            adm.set = (kn != 0) & (kelimi == 0)
            adm.index = which(adm.set == T)
            y.adm = ky[adm.set]
            n.adm = kn[adm.set]
            selectd = sort(abs(kphat[adm.set] - target),
                           index.return = T)$ix[1]
            kseldose = adm.index[selectd]
          }
          selectdoses[k, 2] = ifelse(is.na(kseldose), 99,
                                     kseldose)
          if (k < dim(n)[1])
            if (selectdoses[k + 1, 2] == dim(n)[2])
              selectdoses[k, 2] = dim(n)[2]
          if (k < dim(n)[1])
            if (aa(selectdoses[k + 1, 2]) == dim(n)[2] &
                aa(selectdoses[k + 1, 2]) == aa(selectdoses[k,
                                                            2]))
              selectdoses[k, 2] = 99
        }
      }else {
        selectdoses = matrix(99, nrow = 1, ncol = 2)
        selectdoses[1, ] = matrix(selectdose, nrow = 1)
      }
      
      selectdoses = matrix(selectdoses[selectdoses[, 2] !=
                                         99, ], ncol = 2)
    }
    
    colnames(selectdoses) = c("DoseA", "DoseB")
    
  }
  if (mtd.contour == FALSE) {
    if (selectdoses[1, 1] == 99 && selectdoses[1, 2] == 99) {
      cat("All tested doses are overly toxic. No MTD is selected! \n")
     out=list(target = target, MTD = 99, p_est = matrix(NA,nrow = dim(npts)[1], ncol = dim(npts)[2]))
    }
    else {
    
    out=list(target = target, MTD = selectdoses, p_est=phat.out.noCI,p_est_CI = phat.out)
    }
    
    class(out)<-"boin"
    return(out)
  }
  
  else {
    if (length(selectdoses) == 0) {
      cat("All tested doses are overly toxic. No MTD is selected! \n")
      out=list(target = target, MTD = 99, p_est = matrix(NA,nrow = dim(npts)[1], ncol = dim(npts)[2]))
    }
    else {
      
      out=list(target = target, MTD = selectdoses,  p_est=phat.out.noCI,p_est_CI = phat.out)
    }
    
    class(out)<-"boin"
    return(out)
  }
}


boin.simu.fn = function(target, p.true, ncohort, cohortsize, seed){
  res <- get.oc.comb(target, p.true, ncohort, cohortsize, ntrial=1, seed=seed)
  correct <- 0
  if(sum(res$selpercent)!=100){
    correct <- 0
  } else {
    MTD <- which(res$selpercent == 100, arr.ind = T)
    if(p.true[MTD[1],MTD[2]]==target){
      correct <- 1
    }
  }
  
  npercent <- 0
  for (j in 1:length(p.true[,1])) {
    for (k in 1:length(p.true[1,])) {
      if (p.true[j,k]==target){
        npercent <- npercent + res$npatients[j,k]
      }
    }
  }
  npercent <- percent(npercent/(ncohort*cohortsize))
  
  ptoxic <- 0
  for (j in 1:length(p.true[,1])) {
    for (k in 1:length(p.true[1,])) {
      if (p.true[j,k]>target){
        ptoxic <- ptoxic + res$npatients[j,k]
      }
    }
  }
  ptoxic <- percent(ptoxic/(ncohort*cohortsize))
  
  list(correct=correct, npercent=npercent, ntox=res$totaltox, ptoxic=ptoxic)
}

dfcomb.simu.fn = function(ndose_a1, ndose_a2, p_tox, target, target_min, target_max, prior_tox_a1, prior_tox_a2, n_cohort,
                          cohort, tite=FALSE, time_full=0, poisson_rate=0, c_e=0.85, c_d=0.45, c_stop=0.95, c_t=0.5,
                          c_over=0.25, cmin_overunder=2, cmin_mtd=3, cmin_recom=1, startup=0, alloc_rule=1, early_stop=1,
                          nburn=2000, niter=5000, seed=1){
  
  res <- CombIncrease_sim(ndose_a1=ndose_a1, ndose_a2=ndose_a2, p_tox=p_tox, target=target, target_min=target_min, target_max=target_max, prior_tox_a1=prior_tox_a1, prior_tox_a2=prior_tox_a2, n_cohort=n_cohort,
                          cohort=cohort, tite=tite, time_full=time_full, poisson_rate=poisson_rate, nsim=1, c_e=c_e, c_d=c_d, c_stop=c_stop, c_t=c_t,
                          c_over=c_over, cmin_overunder=cmin_overunder, cmin_mtd=cmin_mtd, cmin_recom=cmin_recom, startup=startup, alloc_rule=alloc_rule, early_stop=early_stop,
                          nburn=nburn, niter=niter, seed=seed)
  
  correct <- 0
  if(sum(res$rec_dose)!=100){
    correct <- 0
  } else {
    MTD <- which(res$rec_dose == 100, arr.ind = T)
    if(p_tox[MTD[1],MTD[2]]==target){
      correct <- 1
    }
  }
  npercent <- 0
  for (j in 1:ndose_a1) {
    for (k in 1:ndose_a2) {
      if (p_tox[j,k]==target){
        npercent <- npercent + res$n_pat_dose[j,k]
      }
    }
  }
  npercent <- npercent/(n_cohort*cohort)
  
  ptoxic <- 0
  for (j in 1:ndose_a1) {
    for (k in 1:ndose_a2) {
      if (p_tox[j,k]>target){
        ptoxic <- ptoxic + res$n_pat_dose[j,k]
      }
    }
  }
  ptoxic <- ptoxic/(n_cohort*cohort)
  
  list(correct=correct, npercent=npercent, ntox=sum(res$n_tox_dose), ptoxic=ptoxic)
}




######### Random scenario
rand <- function(dim1,dim2,nMTD,e){
  repeat{
    s <- sample(c(runif_gap(dim1*dim2-nMTD,0,1,e),rep(0.3,nMTD)))
    if(length(s[s<=(0.3+e) & s>=(0.3-e)])==nMTD){
      break
    }
  }
  
  m <- matrix(s,dim1,dim2)
  for (i in 1:dim1) {
    m[i,] <- sort(m[i,])
  }
  for (j in 1:dim2) {
    m[,j] <- sort(m[,j])
  }
  return(m)
}

random <- function(dim1,dim2,nMTD,e=0.05){
  i <- 0
  while (i==0) {
    p <- rand(dim1,dim2,nMTD,e)
    if(!check_monoto(0.3,p)){
      return(p)
    }
  }
}

runif_gap <- function(n,min, max, e){
  repeat{
    s <- sort(runif(n,min,max))
    mindiff <- min(diff(s))
    if (mindiff >= e){
      break
    }
  }
  return(s)
}

check_monoto <- function(phi, p){
  a <- which(p==phi, arr.ind = T)
  
  if(length(a[,1])==1){
    return(FALSE)
  }
  
  for (i in 1:length(a[1,])) {
    if(any(duplicated(a[,i]))){
      return(TRUE)
    }
  }
  return(FALSE)
}
