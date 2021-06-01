#' @title Estimating selection coefficients and testing their changes from ancient DNA data
#' @author Xiaoyang Dai, Mark Beaumont, Feng Yu, Zhangyi He

#' version 1.0
#' Two-gene phenotypes under constant natural selection and constant demographic histories
#' Horse white coat patterns (KIT13 & KIT116)

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
source("./Code/Code v1.0/Code 2L/Code v1.0/RFUN_PTN.R")

################################################################################

#' Wutke et al. (2016) from 3648 BC (Domestication 3500 BC)
load("./Data/REAL_PTN_RAW.rda")

set.seed(1)

sel_cof <- c(0e+00, 0e+00, 0e+00)
rec_rat <- 1e-08 * 4688
pop_siz <- 1.6e+04
smp_gen <- smp_gen[-(1:13)]
smp_siz <- smp_siz[-(1:13)]
smp_cnt <- smp_cnt[, -(1:13)]
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04
stp_siz <- (1:itn_num)^(-2 / 3)
apt_rto <- 4e-01

system.time(sel_cof_chn <- cmprunAdaptPMMH(sel_cof, rec_rat, pop_siz, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, stp_siz, apt_rto))

save(sel_cof, rec_rat, pop_siz, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, stp_siz, apt_rto, sel_cof_chn,
     file = "./Output/Output v1.0/REAL v1.0/REAL_PTN_RAW_PMMH1.rda")

load("./Output/Output v1.0/REAL v1.0/REAL_PTN_RAW_PMMH1.rda")

pdf(file = "./Output/Output v1.0/REAL v1.0/REAL_PTN_RAW_PMMH1_Traceplot_SelCoeff.pdf", width = 8, height = 8)
par(mfrow = c(2, 1), mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.5, cex.sub = 1.25, cex.axis = 1.25, cex.lab = 1.25)
plot(1:itn_num, sel_cof_chn[1, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient of tobiano")

plot(1:itn_num, sel_cof_chn[2, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient of sabino")
dev.off()

# brn_num <- 1e+04
brn_num <- 0.5 * dim(sel_cof_chn)[2] + 1
sel_cof_chn <- sel_cof_chn[, brn_num:dim(sel_cof_chn)[2]]

thn_num <- 5e+00
sel_cof_chn <- sel_cof_chn[, (1:round(dim(sel_cof_chn)[2] / thn_num)) * thn_num]

sel_cof_est <- rowMeans(sel_cof_chn)

sel_cof_hpd <- matrix(NA, nrow = 2, ncol = 2)
sel_cof_hpd[1, ] <- HPDinterval(as.mcmc(sel_cof_chn[1, ]), prob = 0.95)
sel_cof_hpd[2, ] <- HPDinterval(as.mcmc(sel_cof_chn[2, ]), prob = 0.95)

pdf(file = "./Output/Output v1.0/REAL v1.0/REAL_PTN_RAW_PMMH1_Posterior_SelCoeff.pdf", width = 16, height = 8)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
layout(matrix(c(1, 1, 2, 3), nrow = 2, ncol = 2))
grd_num <- 1e+03
sel_cof_pdf <- kde2d(sel_cof_chn[1, ], sel_cof_chn[2, ], n = grd_num)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
      xlab = "Selection coefficient of tobiano", ylab = "Selection coefficient of sabino",
      main = "Posterior for selection coefficient")
abline(v = sel_cof_est[1], col = 'black', lty = 2, lwd = 2)
abline(h = sel_cof_est[2], col = 'black', lty = 2, lwd = 2)
# HPDregionplot(as.mcmc(sel_cof_chn), vars = 1:2, n = grd_num, prob = 0.95, col = "blue", lwd = 2, add = TRUE)

hist(sel_cof_chn[1, ], breaks = seq(min(sel_cof_chn[1, ]), max(sel_cof_chn[1, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient of tobiano")
lines(density(sel_cof_chn[1, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[1], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 2], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_chn[2, ], breaks = seq(min(sel_cof_chn[2, ]), max(sel_cof_chn[2, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient of sabino")
lines(density(sel_cof_chn[2, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[2], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 2], col = 'blue', lty = 2, lwd = 2)
dev.off()

########################################

#' Wutke et al. (2016) from 3648 BC (Domestication 3500 BC)
load("./Data/REAL_PTN_RAW.rda")

set.seed(1)

sel_cof <- c(0e+00, 0e+00, 0e+00)
rec_rat <- 1e-08 * 4688
pop_siz <- 1.6e+04
smp_gen <- smp_gen[-(1:12)]
smp_siz <- smp_siz[-(1:12)]
smp_cnt <- smp_cnt[, -(1:12)]
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04
stp_siz <- (1:itn_num)^(-2 / 3)
apt_rto <- 4e-01

system.time(sel_cof_chn <- cmprunAdaptPMMH(sel_cof, rec_rat, pop_siz, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, stp_siz, apt_rto))

save(sel_cof, rec_rat, pop_siz, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, stp_siz, apt_rto, sel_cof_chn,
     file = "./Output/Output v1.0/REAL v1.0/REAL_PTN_RAW_PMMH2.rda")

load("./Output/Output v1.0/REAL v1.0/REAL_PTN_RAW_PMMH2.rda")

pdf(file = "./Output/Output v1.0/REAL v1.0/REAL_PTN_RAW_PMMH2_Traceplot_SelCoeff.pdf", width = 8, height = 8)
par(mfrow = c(2, 1), mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.5, cex.sub = 1.25, cex.axis = 1.25, cex.lab = 1.25)
plot(1:itn_num, sel_cof_chn[1, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient of tobiano")

plot(1:itn_num, sel_cof_chn[2, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient of sabino")
dev.off()

# brn_num <- 1e+04
brn_num <- 0.5 * dim(sel_cof_chn)[2] + 1
sel_cof_chn <- sel_cof_chn[, brn_num:dim(sel_cof_chn)[2]]

thn_num <- 5e+00
sel_cof_chn <- sel_cof_chn[, (1:round(dim(sel_cof_chn)[2] / thn_num)) * thn_num]

sel_cof_est <- rowMeans(sel_cof_chn)

sel_cof_hpd <- matrix(NA, nrow = 2, ncol = 2)
sel_cof_hpd[1, ] <- HPDinterval(as.mcmc(sel_cof_chn[1, ]), prob = 0.95)
sel_cof_hpd[2, ] <- HPDinterval(as.mcmc(sel_cof_chn[2, ]), prob = 0.95)

pdf(file = "./Output/Output v1.0/REAL v1.0/REAL_PTN_RAW_PMMH2_Posterior_SelCoeff.pdf", width = 16, height = 8)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
layout(matrix(c(1, 1, 2, 3), nrow = 2, ncol = 2))
grd_num <- 1e+03
sel_cof_pdf <- kde2d(sel_cof_chn[1, ], sel_cof_chn[2, ], n = grd_num)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
      xlab = "Selection coefficient of tobiano", ylab = "Selection coefficient of sabino",
      main = "Posterior for selection coefficient")
abline(v = sel_cof_est[1], col = 'black', lty = 2, lwd = 2)
abline(h = sel_cof_est[2], col = 'black', lty = 2, lwd = 2)
# HPDregionplot(as.mcmc(sel_cof_chn), vars = 1:2, n = grd_num, prob = 0.95, col = "blue", lwd = 2, add = TRUE)

hist(sel_cof_chn[1, ], breaks = seq(min(sel_cof_chn[1, ]), max(sel_cof_chn[1, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient of tobiano")
lines(density(sel_cof_chn[1, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[1], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 2], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_chn[2, ], breaks = seq(min(sel_cof_chn[2, ]), max(sel_cof_chn[2, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient of sabino")
lines(density(sel_cof_chn[2, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[2], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 2], col = 'blue', lty = 2, lwd = 2)
dev.off()

########################################

#' Wutke et al. (2016) from 4304 BC (Domestication 3500 BC)
load("./Data/REAL_PTN_RAW.rda")

set.seed(1)

sel_cof <- c(0e+00, 0e+00, 0e+00)
rec_rat <- 1e-08 * 4688
pop_siz <- 1.6e+04
smp_gen <- smp_gen[-(1:11)]
smp_siz <- smp_siz[-(1:11)]
smp_cnt <- smp_cnt[, -(1:11)]
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04
stp_siz <- (1:itn_num)^(-2 / 3)
apt_rto <- 4e-01

system.time(sel_cof_chn <- cmprunAdaptPMMH(sel_cof, rec_rat, pop_siz, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, stp_siz, apt_rto))

save(sel_cof, rec_rat, pop_siz, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, stp_siz, apt_rto, sel_cof_chn,
     file = "./Output/Output v1.0/REAL v1.0/REAL_PTN_RAW_PMMH3.rda")

load("./Output/Output v1.0/REAL v1.0/REAL_PTN_RAW_PMMH3.rda")

pdf(file = "./Output/Output v1.0/REAL v1.0/REAL_PTN_RAW_PMMH3_Traceplot_SelCoeff.pdf", width = 8, height = 8)
par(mfrow = c(2, 1), mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.5, cex.sub = 1.25, cex.axis = 1.25, cex.lab = 1.25)
plot(1:itn_num, sel_cof_chn[1, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient of tobiano")

plot(1:itn_num, sel_cof_chn[2, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient of sabino")
dev.off()

# brn_num <- 1e+04
brn_num <- 0.5 * dim(sel_cof_chn)[2] + 1
sel_cof_chn <- sel_cof_chn[, brn_num:dim(sel_cof_chn)[2]]

thn_num <- 5e+00
sel_cof_chn <- sel_cof_chn[, (1:round(dim(sel_cof_chn)[2] / thn_num)) * thn_num]

sel_cof_est <- rowMeans(sel_cof_chn)

sel_cof_hpd <- matrix(NA, nrow = 2, ncol = 2)
sel_cof_hpd[1, ] <- HPDinterval(as.mcmc(sel_cof_chn[1, ]), prob = 0.95)
sel_cof_hpd[2, ] <- HPDinterval(as.mcmc(sel_cof_chn[2, ]), prob = 0.95)

pdf(file = "./Output/Output v1.0/REAL v1.0/REAL_PTN_RAW_PMMH3_Posterior_SelCoeff.pdf", width = 16, height = 8)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
layout(matrix(c(1, 1, 2, 3), nrow = 2, ncol = 2))
grd_num <- 1e+03
sel_cof_pdf <- kde2d(sel_cof_chn[1, ], sel_cof_chn[2, ], n = grd_num)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
      xlab = "Selection coefficient of tobiano", ylab = "Selection coefficient of sabino",
      main = "Posterior for selection coefficient")
abline(v = sel_cof_est[1], col = 'black', lty = 2, lwd = 2)
abline(h = sel_cof_est[2], col = 'black', lty = 2, lwd = 2)
# HPDregionplot(as.mcmc(sel_cof_chn), vars = 1:2, n = grd_num, prob = 0.95, col = "blue", lwd = 2, add = TRUE)

hist(sel_cof_chn[1, ], breaks = seq(min(sel_cof_chn[1, ]), max(sel_cof_chn[1, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient of tobiano")
lines(density(sel_cof_chn[1, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[1], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 2], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_chn[2, ], breaks = seq(min(sel_cof_chn[2, ]), max(sel_cof_chn[2, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient of sabino")
lines(density(sel_cof_chn[2, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[2], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 2], col = 'blue', lty = 2, lwd = 2)
dev.off()

########################################

#' Wutke et al. (2016) (Domestication 3500 BC)
load("./Data/REAL_PTN_GRP.rda")

set.seed(1)

sel_cof <- c(0e+00, 0e+00, 0e+00)
rec_rat <- 1e-08 * 4688
pop_siz <- 1.6e+04
smp_gen <- smp_gen[-(1:2)]
smp_siz <- smp_siz[-(1:2)]
smp_cnt <- smp_cnt[, -(1:2)]
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04
stp_siz <- (1:itn_num)^(-2 / 3)
apt_rto <- 4e-01

system.time(sel_cof_chn <- cmprunAdaptPMMH(sel_cof, rec_rat, pop_siz, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, stp_siz, apt_rto))

save(sel_cof, rec_rat, pop_siz, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, stp_siz, apt_rto, sel_cof_chn,
     file = "./Output/Output v1.0/REAL v1.0/REAL_PTN_GRP_PMMH.rda")

load("./Output/Output v1.0/REAL v1.0/REAL_PTN_GRP_PMMH.rda")

pdf(file = "./Output/Output v1.0/REAL v1.0/REAL_PTN_GRP_PMMH_Traceplot_SelCoeff.pdf", width = 8, height = 8)
par(mfrow = c(2, 1), mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.5, cex.sub = 1.25, cex.axis = 1.25, cex.lab = 1.25)
plot(1:itn_num, sel_cof_chn[1, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient of tobiano")

plot(1:itn_num, sel_cof_chn[2, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient of sabino")
dev.off()

# brn_num <- 1e+04
brn_num <- 0.5 * dim(sel_cof_chn)[2] + 1
sel_cof_chn <- sel_cof_chn[, brn_num:dim(sel_cof_chn)[2]]

thn_num <- 5e+00
sel_cof_chn <- sel_cof_chn[, (1:round(dim(sel_cof_chn)[2] / thn_num)) * thn_num]

sel_cof_est <- rowMeans(sel_cof_chn)

sel_cof_hpd <- matrix(NA, nrow = 2, ncol = 2)
sel_cof_hpd[1, ] <- HPDinterval(as.mcmc(sel_cof_chn[1, ]), prob = 0.95)
sel_cof_hpd[2, ] <- HPDinterval(as.mcmc(sel_cof_chn[2, ]), prob = 0.95)

pdf(file = "./Output/Output v1.0/REAL v1.0/REAL_PTN_GRP_PMMH_Posterior_SelCoeff.pdf", width = 16, height = 8)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
layout(matrix(c(1, 1, 2, 3), nrow = 2, ncol = 2))
grd_num <- 1e+03
sel_cof_pdf <- kde2d(sel_cof_chn[1, ], sel_cof_chn[2, ], n = grd_num)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
      xlab = "Selection coefficient of tobiano", ylab = "Selection coefficient of sabino",
      main = "Posterior for selection coefficient")
abline(v = sel_cof_est[1], col = 'black', lty = 2, lwd = 2)
abline(h = sel_cof_est[2], col = 'black', lty = 2, lwd = 2)
# HPDregionplot(as.mcmc(sel_cof_chn), vars = 1:2, n = grd_num, prob = 0.95, col = "blue", lwd = 2, add = TRUE)

hist(sel_cof_chn[1, ], breaks = seq(min(sel_cof_chn[1, ]), max(sel_cof_chn[1, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient of tobiano")
lines(density(sel_cof_chn[1, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[1], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 2], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_chn[2, ], breaks = seq(min(sel_cof_chn[2, ]), max(sel_cof_chn[2, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient of sabino")
lines(density(sel_cof_chn[2, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[2], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 2], col = 'blue', lty = 2, lwd = 2)
dev.off()

################################################################################
