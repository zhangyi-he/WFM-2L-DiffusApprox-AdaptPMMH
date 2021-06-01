#' @title Estimating selection coefficients and testing their changes from ancient DNA data
#' @author Xiaoyang Dai, Mark Beaumont, Feng Yu, Zhangyi He

#' version 1.0
#' Single-gene phenotypes under constant natural selection and constant demographic histories
#' Time series data of genotype frequencies

#' Horse base coat colours (ASIP & MC1R) and white coat patterns (KIT13 & KIT16 & TRMP1)

# set the directory
setwd("~/Dropbox/Jeffery He/iResearch/Publications/2019/HE2021-WFM-2L-DiffusApprox-MHwGibbs-MolEcol")

#install.packages("RColorBrewer")
library("RColorBrewer")

#install.packages("viridis")
library("viridis")

#install.packages("ggplot2")
library("ggplot2")

#install.packages("plot3D")
library("plot3D")

#install.packages("emdbook")
library("emdbook")

# call R functions
source("./Code/Code v1.0/Code 1L/Code v1.0/RFUN_GEN.R")

################################################################################

#' Wutke et al. (2016) from 9320 BC (Holocene 9700 BC)
load("./Data/REAL_COL_RAW.rda")

set.seed(1)

sel_cof <- 0e+00
dom_par <- 0
pop_siz <- 1.6e+04
smp_gen <- smp_gen[-(1:3)]
smp_siz <- smp_siz[-(1:3)]
smp_cnt <- smp_gen_cnt_ASIP[, -(1:3)]
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04

system.time(sel_cof_chn <- cmprunPMMH(sel_cof, dom_par, pop_siz, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num))

save(sel_cof, dom_par, pop_siz, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, sel_cof_chn,
     file = "./Output/Output v1.0/REAL v1.0/REAL_COL_RAW_PMMH2g_ASIP.rda")

load("./Output/Output v1.0/REAL v1.0/REAL_COL_RAW_PMMH2g_ASIP.rda")

pdf(file = "./Output/Output v1.0/REAL v1.0/REAL_COL_RAW_PMMH2g_ASIP_Traceplot_SelCoeff.pdf", width = 8, height = 4)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:itn_num, sel_cof_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient")
dev.off()

# brn_num <- 1e+04
brn_num <- 0.5 * length(sel_cof_chn) + 1
sel_cof_chn <- sel_cof_chn[brn_num:length(sel_cof_chn)]

thn_num <- 5e+00
sel_cof_chn <- sel_cof_chn[(1:round(length(sel_cof_chn) / thn_num)) * thn_num]

sel_cof_est <- mean(sel_cof_chn)

sel_cof_hpd <- HPDinterval(as.mcmc(sel_cof_chn), prob = 0.95)

pdf(file = "./Output/Output v1.0/REAL v1.0/REAL_COL_RAW_PMMH2g_ASIP_Posterior_SelCoeff.pdf", width = 8, height = 6)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(sel_cof_chn, breaks = seq(min(sel_cof_chn), max(sel_cof_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the selection coefficient")
lines(density(sel_cof_chn), lwd = 2, col = 'black')
abline(v = sel_cof_est, col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2], col = 'blue', lty = 2, lwd = 2)
dev.off()

####################

#' Wutke et al. (2016) from 9320 BC (Holocene 9700 BC)
load("./Data/REAL_COL_GRP.rda")

set.seed(1)

sel_cof <- 0e+00
dom_par <- 0
pop_siz <- 1.6e+04
smp_gen
smp_siz
smp_cnt <- smp_gen_cnt_ASIP
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04

system.time(sel_cof_chn <- cmprunPMMH(sel_cof, dom_par, pop_siz, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num))

save(sel_cof, dom_par, pop_siz, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, sel_cof_chn,
     file = "./Output/Output v1.0/REAL v1.0/REAL_COL_GRP_PMMHg_ASIP.rda")

load("./Output/Output v1.0/REAL v1.0/REAL_COL_GRP_PMMHg_ASIP.rda")

pdf(file = "./Output/Output v1.0/REAL v1.0/REAL_COL_GRP_PMMHg_ASIP_Traceplot_SelCoeff.pdf", width = 8, height = 4)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:itn_num, sel_cof_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient")
dev.off()

# brn_num <- 1e+04
brn_num <- 0.5 * length(sel_cof_chn) + 1
sel_cof_chn <- sel_cof_chn[brn_num:length(sel_cof_chn)]

thn_num <- 5e+00
sel_cof_chn <- sel_cof_chn[(1:round(length(sel_cof_chn) / thn_num)) * thn_num]

sel_cof_est <- mean(sel_cof_chn)

sel_cof_hpd <- HPDinterval(as.mcmc(sel_cof_chn), prob = 0.95)

pdf(file = "./Output/Output v1.0/REAL v1.0/REAL_COL_GRP_PMMHg_ASIP_Posterior_SelCoeff.pdf", width = 8, height = 6)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(sel_cof_chn, breaks = seq(min(sel_cof_chn), max(sel_cof_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the selection coefficient")
lines(density(sel_cof_chn), lwd = 2, col = 'black')
abline(v = sel_cof_est, col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2], col = 'blue', lty = 2, lwd = 2)
dev.off()

########################################

#' Wutke et al. (2016) from 9320 BC (Holocene 9700 BC)
load("./Data/REAL_COL_RAW.rda")

set.seed(1)

sel_cof <- 0e+00
dom_par <- 0
pop_siz <- 1.6e+04
smp_gen <- smp_gen[-(1:3)]
smp_siz <- smp_siz[-(1:3)]
smp_cnt <- smp_gen_cnt_MC1R[, -(1:3)]
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04

system.time(sel_cof_chn <- cmprunPMMH(sel_cof, dom_par, pop_siz, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num))

save(sel_cof, dom_par, pop_siz, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, sel_cof_chn,
     file = "./Output/Output v1.0/REAL v1.0/REAL_COL_RAW_PMMH2g_MC1R.rda")

load("./Output/Output v1.0/REAL v1.0/REAL_COL_RAW_PMMH2g_MC1R.rda")

pdf(file = "./Output/Output v1.0/REAL v1.0/REAL_COL_RAW_PMMH2g_MC1R_Traceplot_SelCoeff.pdf", width = 8, height = 4)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:itn_num, sel_cof_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient")
dev.off()

# brn_num <- 1e+04
brn_num <- 0.5 * length(sel_cof_chn) + 1
sel_cof_chn <- sel_cof_chn[brn_num:length(sel_cof_chn)]

thn_num <- 5e+00
sel_cof_chn <- sel_cof_chn[(1:round(length(sel_cof_chn) / thn_num)) * thn_num]

sel_cof_est <- mean(sel_cof_chn)

sel_cof_hpd <- HPDinterval(as.mcmc(sel_cof_chn), prob = 0.95)

pdf(file = "./Output/Output v1.0/REAL v1.0/REAL_COL_RAW_PMMH2g_MC1R_Posterior_SelCoeff.pdf", width = 8, height = 6)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(sel_cof_chn, breaks = seq(min(sel_cof_chn), max(sel_cof_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the selection coefficient")
lines(density(sel_cof_chn), lwd = 2, col = 'black')
abline(v = sel_cof_est, col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2], col = 'blue', lty = 2, lwd = 2)
dev.off()

####################

#' Wutke et al. (2016) from 9320 BC (Holocene 9700 BC)
load("./Data/REAL_COL_GRP.rda")

set.seed(1)

sel_cof <- 0e+00
dom_par <- 0
pop_siz <- 1.6e+04
smp_gen
smp_siz
smp_cnt <- smp_gen_cnt_MC1R
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04

system.time(sel_cof_chn <- cmprunPMMH(sel_cof, dom_par, pop_siz, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num))

save(sel_cof, dom_par, pop_siz, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, sel_cof_chn,
     file = "./Output/Output v1.0/REAL v1.0/REAL_COL_GRP_PMMHg_MC1R.rda")

load("./Output/Output v1.0/REAL v1.0/REAL_COL_GRP_PMMHg_MC1R.rda")

pdf(file = "./Output/Output v1.0/REAL v1.0/REAL_COL_GRP_PMMHg_MC1R_Traceplot_SelCoeff.pdf", width = 8, height = 4)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:itn_num, sel_cof_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient")
dev.off()

# brn_num <- 1e+04
brn_num <- 0.5 * length(sel_cof_chn) + 1
sel_cof_chn <- sel_cof_chn[brn_num:length(sel_cof_chn)]

thn_num <- 5e+00
sel_cof_chn <- sel_cof_chn[(1:round(length(sel_cof_chn) / thn_num)) * thn_num]

sel_cof_est <- mean(sel_cof_chn)

sel_cof_hpd <- HPDinterval(as.mcmc(sel_cof_chn), prob = 0.95)

pdf(file = "./Output/Output v1.0/REAL v1.0/REAL_COL_GRP_PMMHg_MC1R_Posterior_SelCoeff.pdf", width = 8, height = 6)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(sel_cof_chn, breaks = seq(min(sel_cof_chn), max(sel_cof_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the selection coefficient")
lines(density(sel_cof_chn), lwd = 2, col = 'black')
abline(v = sel_cof_est, col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2], col = 'blue', lty = 2, lwd = 2)
dev.off()

########################################

#' Wutke et al. (2016) from 9320 BC (Holocene 9700 BC)
load("./Data/REAL_PTN_RAW.rda")

set.seed(1)

sel_cof <- 0e+00
dom_par <- 1
pop_siz <- 1.6e+04
smp_gen <- smp_gen[-(1:13)]
smp_siz <- smp_siz[-(1:13)]
smp_cnt <- smp_gen_cnt_KIT13[, -(1:13)]
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04

system.time(sel_cof_chn <- cmprunPMMH(sel_cof, dom_par, pop_siz, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num))

save(sel_cof, dom_par, pop_siz, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, sel_cof_chn,
     file = "./Output/Output v1.0/REAL v1.0/REAL_PTN_RAW_PMMH1g_KIT13.rda")

load("./Output/Output v1.0/REAL v1.0/REAL_PTN_RAW_PMMH1g_KIT13.rda")

pdf(file = "./Output/Output v1.0/REAL v1.0/REAL_PTN_RAW_PMMH1g_KIT13_Traceplot_SelCoeff.pdf", width = 8, height = 4)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:itn_num, sel_cof_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient")
dev.off()

# brn_num <- 1e+04
brn_num <- 0.5 * length(sel_cof_chn) + 1
sel_cof_chn <- sel_cof_chn[brn_num:length(sel_cof_chn)]

thn_num <- 5e+00
sel_cof_chn <- sel_cof_chn[(1:round(length(sel_cof_chn) / thn_num)) * thn_num]

sel_cof_est <- mean(sel_cof_chn)

sel_cof_hpd <- HPDinterval(as.mcmc(sel_cof_chn), prob = 0.95)

pdf(file = "./Output/Output v1.0/REAL v1.0/REAL_PTN_RAW_PMMH1g_KIT13_Posterior_SelCoeff.pdf", width = 8, height = 6)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(sel_cof_chn, breaks = seq(min(sel_cof_chn), max(sel_cof_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the selection coefficient")
lines(density(sel_cof_chn), lwd = 2, col = 'black')
abline(v = sel_cof_est, col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2], col = 'blue', lty = 2, lwd = 2)
dev.off()

####################

#' Wutke et al. (2016) from 9320 BC (Holocene 9700 BC)
load("./Data/REAL_PTN_GRP.rda")

set.seed(1)

sel_cof <- 0e+00
dom_par <- 1
pop_siz <- 1.6e+04
smp_gen <- smp_gen[-(1:2)]
smp_siz <- smp_siz[-(1:2)]
smp_cnt <- smp_gen_cnt_KIT13[, -(1:2)]
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04

system.time(sel_cof_chn <- cmprunPMMH(sel_cof, dom_par, pop_siz, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num))

save(sel_cof, dom_par, pop_siz, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, sel_cof_chn,
     file = "./Output/Output v1.0/REAL v1.0/REAL_PTN_GRP_PMMHg_KIT13.rda")

load("./Output/Output v1.0/REAL v1.0/REAL_PTN_GRP_PMMHg_KIT13.rda")

pdf(file = "./Output/Output v1.0/REAL v1.0/REAL_PTN_GRP_PMMHg_KIT13_Traceplot_SelCoeff.pdf", width = 8, height = 4)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:itn_num, sel_cof_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient")
dev.off()

# brn_num <- 1e+04
brn_num <- 0.5 * length(sel_cof_chn) + 1
sel_cof_chn <- sel_cof_chn[brn_num:length(sel_cof_chn)]

thn_num <- 5e+00
sel_cof_chn <- sel_cof_chn[(1:round(length(sel_cof_chn) / thn_num)) * thn_num]

sel_cof_est <- mean(sel_cof_chn)

sel_cof_hpd <- HPDinterval(as.mcmc(sel_cof_chn), prob = 0.95)

pdf(file = "./Output/Output v1.0/REAL v1.0/REAL_PTN_GRP_PMMHg_KIT13_Posterior_SelCoeff.pdf", width = 8, height = 6)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(sel_cof_chn, breaks = seq(min(sel_cof_chn), max(sel_cof_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the selection coefficient")
lines(density(sel_cof_chn), lwd = 2, col = 'black')
abline(v = sel_cof_est, col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2], col = 'blue', lty = 2, lwd = 2)
dev.off()

########################################

#' Wutke et al. (2016) from 9320 BC (Holocene 9700 BC)
load("./Data/REAL_PTN_RAW.rda")

set.seed(1)

sel_cof <- 0e+00
dom_par <- 1
pop_siz <- 1.6e+04
smp_gen <- smp_gen[-(1:13)]
smp_siz <- smp_siz[-(1:13)]
smp_cnt <- smp_gen_cnt_KIT16[, -(1:13)]
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04

system.time(sel_cof_chn <- cmprunPMMH(sel_cof, dom_par, pop_siz, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num))

save(sel_cof, dom_par, pop_siz, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, sel_cof_chn,
     file = "./Output/Output v1.0/REAL v1.0/REAL_PTN_RAW_PMMH1g_KIT16.rda")

load("./Output/Output v1.0/REAL v1.0/REAL_PTN_RAW_PMMH1g_KIT16.rda")

pdf(file = "./Output/Output v1.0/REAL v1.0/REAL_PTN_RAW_PMMH1g_KIT16_Traceplot_SelCoeff.pdf", width = 8, height = 4)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:itn_num, sel_cof_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient")
dev.off()

# brn_num <- 1e+04
brn_num <- 0.5 * length(sel_cof_chn) + 1
sel_cof_chn <- sel_cof_chn[brn_num:length(sel_cof_chn)]

thn_num <- 5e+00
sel_cof_chn <- sel_cof_chn[(1:round(length(sel_cof_chn) / thn_num)) * thn_num]

sel_cof_est <- mean(sel_cof_chn)

sel_cof_hpd <- HPDinterval(as.mcmc(sel_cof_chn), prob = 0.95)

pdf(file = "./Output/Output v1.0/REAL v1.0/REAL_PTN_RAW_PMMH1g_KIT16_Posterior_SelCoeff.pdf", width = 8, height = 6)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(sel_cof_chn, breaks = seq(min(sel_cof_chn), max(sel_cof_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the selection coefficient")
lines(density(sel_cof_chn), lwd = 2, col = 'black')
abline(v = sel_cof_est, col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2], col = 'blue', lty = 2, lwd = 2)
dev.off()

####################

#' Wutke et al. (2016) from 9320 BC (Holocene 9700 BC)
load("./Data/REAL_PTN_GRP.rda")

set.seed(1)

sel_cof <- 0e+00
dom_par <- 1
pop_siz <- 1.6e+04
smp_gen <- smp_gen[-(1:2)]
smp_siz <- smp_siz[-(1:2)]
smp_cnt <- smp_gen_cnt_KIT16[, -(1:2)]
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04

system.time(sel_cof_chn <- cmprunPMMH(sel_cof, dom_par, pop_siz, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num))

save(sel_cof, dom_par, pop_siz, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, sel_cof_chn,
     file = "./Output/Output v1.0/REAL v1.0/REAL_PTN_GRP_PMMHg_KIT16.rda")

load("./Output/Output v1.0/REAL v1.0/REAL_PTN_GRP_PMMHg_KIT16.rda")

pdf(file = "./Output/Output v1.0/REAL v1.0/REAL_PTN_GRP_PMMHg_KIT16_Traceplot_SelCoeff.pdf", width = 8, height = 4)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:itn_num, sel_cof_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient")
dev.off()

# brn_num <- 1e+04
brn_num <- 0.5 * length(sel_cof_chn) + 1
sel_cof_chn <- sel_cof_chn[brn_num:length(sel_cof_chn)]

thn_num <- 5e+00
sel_cof_chn <- sel_cof_chn[(1:round(length(sel_cof_chn) / thn_num)) * thn_num]

sel_cof_est <- mean(sel_cof_chn)

sel_cof_hpd <- HPDinterval(as.mcmc(sel_cof_chn), prob = 0.95)

pdf(file = "./Output/Output v1.0/REAL v1.0/REAL_PTN_GRP_PMMHg_KIT16_Posterior_SelCoeff.pdf", width = 8, height = 6)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(sel_cof_chn, breaks = seq(min(sel_cof_chn), max(sel_cof_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the selection coefficient")
lines(density(sel_cof_chn), lwd = 2, col = 'black')
abline(v = sel_cof_est, col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2], col = 'blue', lty = 2, lwd = 2)
dev.off()

########################################

#' Wutke et al. (2016)
load("./Data/REAL_PTN_RAW_LP.rda")

set.seed(1)

sel_cof <- 0e+00
dom_par <- 1
pop_siz <- 1.6e+04
smp_gen <- smp_gen[-(1:2)]
smp_siz <- smp_siz[-(1:2)]
smp_cnt <- smp_gen_cnt_TRPM1[, -(1:2)]
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04

system.time(sel_cof_chn <- cmprunPMMH(sel_cof, dom_par, pop_siz, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num))

save(sel_cof, dom_par, pop_siz, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, sel_cof_chn,
     file = "./Output/Output v1.0/REAL v1.0/REAL_PTN_RAW_PMMHg_TRPM1.rda")

load("./Output/Output v1.0/REAL v1.0/REAL_PTN_RAW_PMMHg_TRPM1.rda")

pdf(file = "./Output/Output v1.0/REAL v1.0/REAL_PTN_RAW_PMMHg_TRPM1_Traceplot_SelCoeff.pdf", width = 8, height = 4)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:itn_num, sel_cof_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient")
dev.off()

# brn_num <- 1e+04
brn_num <- 0.5 * length(sel_cof_chn) + 1
sel_cof_chn <- sel_cof_chn[brn_num:length(sel_cof_chn)]

thn_num <- 5e+00
sel_cof_chn <- sel_cof_chn[(1:round(length(sel_cof_chn) / thn_num)) * thn_num]

sel_cof_est <- mean(sel_cof_chn)

sel_cof_hpd <- HPDinterval(as.mcmc(sel_cof_chn), prob = 0.95)

pdf(file = "./Output/Output v1.0/REAL v1.0/REAL_PTN_RAW_PMMHg_TRPM1_Posterior_SelCoeff.pdf", width = 8, height = 6)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(sel_cof_chn, breaks = seq(min(sel_cof_chn), max(sel_cof_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the selection coefficient")
lines(density(sel_cof_chn), lwd = 2, col = 'black')
abline(v = sel_cof_est, col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2], col = 'blue', lty = 2, lwd = 2)
dev.off()

####################

#' Wutke et al. (2016)
load("./Data/REAL_PTN_GRP_LP.rda")

set.seed(1)

sel_cof <- 0e+00
dom_par <- 1
pop_siz <- 1.6e+04
smp_gen
smp_siz
smp_cnt <- smp_gen_cnt_TRPM1
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04

system.time(sel_cof_chn <- cmprunPMMH(sel_cof, dom_par, pop_siz, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num))

save(sel_cof, dom_par, pop_siz, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, sel_cof_chn,
     file = "./Output/Output v1.0/REAL v1.0/REAL_PTN_GRP_PMMHg_TRPM1.rda")

load("./Output/Output v1.0/REAL v1.0/REAL_PTN_GRP_PMMHg_TRPM1.rda")

pdf(file = "./Output/Output v1.0/REAL v1.0/REAL_PTN_GRP_PMMHg_TRPM1_Traceplot_SelCoeff.pdf", width = 8, height = 4)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:itn_num, sel_cof_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient")
dev.off()

# brn_num <- 1e+04
brn_num <- 0.5 * length(sel_cof_chn) + 1
sel_cof_chn <- sel_cof_chn[brn_num:length(sel_cof_chn)]

thn_num <- 5e+00
sel_cof_chn <- sel_cof_chn[(1:round(length(sel_cof_chn) / thn_num)) * thn_num]

sel_cof_est <- mean(sel_cof_chn)

sel_cof_hpd <- HPDinterval(as.mcmc(sel_cof_chn), prob = 0.95)

pdf(file = "./Output/Output v1.0/REAL v1.0/REAL_PTN_GRP_PMMHg_TRPM1_Posterior_SelCoeff.pdf", width = 8, height = 6)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(sel_cof_chn, breaks = seq(min(sel_cof_chn), max(sel_cof_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the selection coefficient")
lines(density(sel_cof_chn), lwd = 2, col = 'black')
abline(v = sel_cof_est, col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2], col = 'blue', lty = 2, lwd = 2)
dev.off()

################################################################################
