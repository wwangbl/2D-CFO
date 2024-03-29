source("2D_CFO_MTD.R")


target <- 0.30
ncohort <- 20
cohortsize <- 3
init.level.A <- 1
init.level.B <- 1

add.args <- list(alp.prior=target, bet.prior=1-target)

p.trues.1 <- list()
p.trues.1[[1]] <- c(0.05, 0.10, 0.15, 0.30, 0.45)
p.trues.1[[2]] <- c(0.10, 0.15, 0.30, 0.45, 0.55)
p.trues.1[[3]] <- c(0.15, 0.30, 0.45, 0.50, 0.60)
p.true.1 <- rbind(p.trues.1[[1]],p.trues.1[[2]],p.trues.1[[3]])


p.trues.2 <- list()
p.trues.2[[1]] <- c(0.15, 0.30, 0.45, 0.50, 0.60)
p.trues.2[[2]] <- c(0.30, 0.45, 0.50, 0.60, 0.75)
p.trues.2[[3]] <- c(0.45, 0.55, 0.60, 0.70, 0.80)
p.true.2 <- rbind(p.trues.2[[1]],p.trues.2[[2]],p.trues.2[[3]])


p.trues.3 <- list()
p.trues.3[[1]] <- c(0.02, 0.07, 0.10, 0.15, 0.30)
p.trues.3[[2]] <- c(0.7, 0.10, 0.15, 0.30, 0.45)
p.trues.3[[3]] <- c(0.10, 0.15, 0.30, 0.45, 0.55)
p.true.3 <- rbind(p.trues.3[[1]],p.trues.3[[2]],p.trues.3[[3]])


p.trues.4 <- list()
p.trues.4[[1]] <- c(0.30, 0.45, 0.60, 0.70, 0.80)
p.trues.4[[2]] <- c(0.45, 0.55, 0.65, 0.75, 0.85)
p.trues.4[[3]] <- c(0.50, 0.60, 0.70, 0.80, 0.90)
p.true.4 <- rbind(p.trues.4[[1]],p.trues.4[[2]],p.trues.4[[3]])





# res <- CFO.simu.fn(target, p.true.4, ncohort=ncohort, cohortsize=cohortsize, init.level.A, init.level.B, add.args=add.args)
# res

n.MTD <- 0

for (i in 1:1000){
  res <- CFO.simu.fn(target, p.true.4, ncohort=ncohort, cohortsize=cohortsize, init.level.A, init.level.B, add.args=add.args)
  if (res$MTD[1]==99 | res$MTD[2]==99){
    next
  }
  if (p.true.2[res$MTD[1],res$MTD[2]] == target){
    n.MTD <- n.MTD + 1
  }
  if (i%%100==0){
    message(i)
  }
}

n.MTD/1000
