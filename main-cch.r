rm(list = ls())
require(survival)
require(stats)
require(splines)
require(xtable)
require(gam)
options(warn=-1)


main.path <- '~/Dropbox/UW/dissertation/paper2/biostatistics-submission/biostatistics-revision-2/code-for-biostatistics-16245-R1/code-final-data-analysis/' # mac
setwd(main.path)

source('code-helper-functions.r')

load('dataset-cch.rdata')
tau <- 6 # prediction timeframe
dc.list      <- NULL
results.list <- NULL
results.mx   <- NULL

for(si.i in seq(6, 36, 6)){
    # P = pert = 500 # ncc =T/F NCC/CCH analysis ie. ncc = T = NCC, ncc = F = CCH. 
    dc <- set.data.control(si = si.i, ti = si.i + tau, P = 500, seed = 1, ncc = F) 
    results.list <- c(results.list, main(dc, d.full = d.cch.cohort))
    dc.list <- c(dc.list, list(dc))
    # save(list = ls(), file = 'results-cch.rdata')
}

for(i in 1:length(results.list)){
    res.out <- with(results.list[[i]], cbind(sim.stats, sim.stats.se, n.case, n.ctrl, n.cens))
    dc.out  <-  with(dc.list[[i]], c(si, ti, pred.time, P))
    res.dc.out <- cbind(matrix(rep(dc.out, each = dc$n.stats), nrow = dc$n.stats, byrow = F), round(res.out, 4))
    results.mx <- rbind(results.mx, res.dc.out)
}
colnames(results.mx) <- c('s', 't', 'tau_0', 'Pert', 'Est', 'SE', 'n.case', 'n.ctrl', 'n.cens')
results.mx

save(list = ls(), file = paste('results-cch-table-p', dc$P, '-tau', tau, '.rdata', sep = ''))















