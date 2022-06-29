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

# res <- CFO.simu.fn(target, p.true.1, ncohort=ncohort, cohortsize=cohortsize, init.level.A, init.level.B, add.args=add.args)
# res

n.MTD <- 0

for (i in 1:1000){
  res <- CFO.simu.fn(target, p.true.1, ncohort=ncohort, cohortsize=cohortsize, init.level.A, init.level.B, add.args=add.args)
  if (p.true.1[res$MTD[1],res$MTD[2]] == target){
    n.MTD <- n.MTD + 1
  }
}

n.MTD/1000
