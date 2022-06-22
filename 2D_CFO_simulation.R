source("2D_CFO_MTD.R")


target <- 0.33
ncohort <- 10
cohortsize <- 3
init.level.A <- 1
init.level.B <- 1

add.args <- list(alp.prior=target, bet.prior=1-target)
p.trues <- list()
p.trues[[1]] <- c(0.45, 0.52, 0.62, 0.70, 0.80)
p.trues[[2]] <- c(0.30, 0.40, 0.52, 0.60, 0.70)
p.trues[[3]] <- c(0.12, 0.20, 0.33, 0.40, 0.50)
p.trues[[4]] <- c(0.01, 0.06, 0.13, 0.23, 0.50)
p.trues[[5]] <- c(0.00, 0.02, 0.05, 0.10, 0.45)

p.true <- rbind(p.trues[[5]],p.trues[[4]],p.trues[[3]],p.trues[[2]],p.trues[[1]])


res <- CFO.simu.fn(target, p.true, ncohort=ncohort, cohortsize=cohortsize, init.level.A, init.level.B, add.args=add.args)
res