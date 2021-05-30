#' @title Estimating selection coefficients and testing their changes from ancient DNA data
#' @author Xiaoyang Dai, Mark Beaumont, Feng Yu, Zhangyi He

#' version 1.3
#' Two-gene phenotypes under non-constant natural selection and non-constant demographic histories conditional on genetic polymorphism
#' Horse base coat colours (ASIP & MC1R)

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
source("./Code/Code v1.0/Code 2L/Code v1.3/RFUN_COL.R")

################################################################################

#' Wutke et al. (2016) from 9320 BC (Holocene 9700 BC)
load("./Data/REAL_COL_RAW.rda")

set.seed(3)

sel_cof <- matrix(c(0e+00, 0e+00, 0e+00, 0e+00), nrow = 2, ncol = 2)
rec_rat <- 5e-01
pop_siz <- pop_siz[-(1:(smp_gen[5] - smp_gen[1]))]
ref_siz <- 1.6e+04
evt_gen <- round((-3500 - 2000) / 8)
smp_gen <- smp_gen[-(1:4)]
smp_siz <- smp_siz[-(1:4)]
smp_cnt <- smp_cnt[, -(1:4)]
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04
stp_siz <- (1:itn_num)^(-2 / 3)
apt_rto <- 4e-01

system.time(sel_cof_chn <- cmprunAdaptPMMH(sel_cof, rec_rat, pop_siz, ref_siz, evt_gen, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, stp_siz, apt_rto))

save(sel_cof, rec_rat, pop_siz, ref_siz, evt_gen, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, stp_siz, apt_rto, sel_cof_chn,
     file = "./Output/Output v1.0/REAL v1.3/REAL_COL_RAW_PMMH1.rda")

load("./Output/Output v1.0/REAL v1.3/REAL_COL_RAW_PMMH1.rda")

pdf(file = "./Output/Output v1.0/REAL v1.3/REAL_COL_RAW_PMMH1_Traceplot_SelCoeff.pdf", width = 16, height = 8)
par(mfrow = c(2, 2), mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.5, cex.sub = 1.25, cex.axis = 1.25, cex.lab = 1.25)
plot(1:itn_num, sel_cof_chn[1, 1, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient of black before domestication")

plot(1:itn_num, sel_cof_chn[2, 1, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient of chestnut before domestication")

plot(1:itn_num, sel_cof_chn[1, 2, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient of black after domestication")

plot(1:itn_num, sel_cof_chn[2, 2, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient of chestnut after domestication")
dev.off()

# brn_num <- 1e+04
brn_num <- 0.5 * dim(sel_cof_chn)[3] + 1
sel_cof_chn <- sel_cof_chn[, , brn_num:dim(sel_cof_chn)[3]]

thn_num <- 5e+00
sel_cof_chn <- sel_cof_chn[, , (1:round(dim(sel_cof_chn)[3] / thn_num)) * thn_num]

sel_cof_est <- matrix(NA, nrow = 2, ncol = 2)
sel_cof_est[, 1] <- rowMeans(sel_cof_chn[, 1, ])
sel_cof_est[, 2] <- rowMeans(sel_cof_chn[, 2, ])

sel_cof_hpd <- array(NA, dim = c(2, 2, 2))
sel_cof_hpd[1, , 1] <- HPDinterval(as.mcmc(sel_cof_chn[1, 1, ]), prob = 0.95)
sel_cof_hpd[2, , 1] <- HPDinterval(as.mcmc(sel_cof_chn[2, 1, ]), prob = 0.95)
sel_cof_hpd[1, , 2] <- HPDinterval(as.mcmc(sel_cof_chn[1, 2, ]), prob = 0.95)
sel_cof_hpd[2, , 2] <- HPDinterval(as.mcmc(sel_cof_chn[2, 2, ]), prob = 0.95)

pdf(file = "./Output/Output v1.0/REAL v1.3/REAL_COL_RAW_PMMH1_Posterior_SelCoeff.pdf", width = 16, height = 16)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
layout(matrix(c(1, 1, 4, 4, 2, 3, 5, 6), nrow = 4, ncol = 2))
grd_num <- 1e+03
sel_cof_pdf <- kde2d(sel_cof_chn[1, 1, ], sel_cof_chn[2, 1, ], n = grd_num)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
      xlab = "Selection coefficient of black", ylab = "Selection coefficient of chestnut",
      main = "Posterior for selection coefficient before domestication")
abline(v = sel_cof_est[1, 1], col = 'black', lty = 2, lwd = 2)
abline(h = sel_cof_est[2, 1], col = 'black', lty = 2, lwd = 2)
# HPDregionplot(as.mcmc(sel_cof_chn), vars = 1:2, n = grd_num, prob = 0.95, col = "blue", lwd = 2, add = TRUE)

hist(sel_cof_chn[1, 1, ], breaks = seq(min(sel_cof_chn[1, 1, ]), max(sel_cof_chn[1, 1, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient of black before domestication")
lines(density(sel_cof_chn[1, 1, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[1, 1], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 1, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 2, 1], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_chn[2, 1, ], breaks = seq(min(sel_cof_chn[2, 1, ]), max(sel_cof_chn[2, 1, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient of chestnut before domestication")
lines(density(sel_cof_chn[2, 1, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[2, 1], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 1, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 2, 1], col = 'blue', lty = 2, lwd = 2)

grd_num <- 1e+03
sel_cof_pdf <- kde2d(sel_cof_chn[1, 2, ], sel_cof_chn[2, 2, ], n = grd_num)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
      xlab = "Selection coefficient of black", ylab = "Selection coefficient of chestnut",
      main = "Posterior for selection coefficient after domestication")
abline(v = sel_cof_est[1, 2], col = 'black', lty = 2, lwd = 2)
abline(h = sel_cof_est[2, 2], col = 'black', lty = 2, lwd = 2)
# HPDregionplot(as.mcmc(sel_cof_chn), vars = 1:2, n = grd_num, prob = 0.95, col = "blue", lwd = 2, add = TRUE)

hist(sel_cof_chn[1, 2, ], breaks = seq(min(sel_cof_chn[1, 2, ]), max(sel_cof_chn[1, 2, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient of black after domestication")
lines(density(sel_cof_chn[1, 2, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[1, 2], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 1, 2], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 2, 2], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_chn[2, 2, ], breaks = seq(min(sel_cof_chn[2, 2, ]), max(sel_cof_chn[2, 2, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient of chestnut after domestication")
lines(density(sel_cof_chn[2, 2, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[2, 2], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 1, 2], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 2, 2], col = 'blue', lty = 2, lwd = 2)
dev.off()

dif_sel_chn <- sel_cof_chn[, 2, ] - sel_cof_chn[, 1, ]

dif_sel_est <- rowMeans(dif_sel_chn)

dif_sel_hpd <- matrix(NA, nrow = 2, ncol = 2)
dif_sel_hpd[1, ] <- HPDinterval(as.mcmc(dif_sel_chn[1, ]), prob = 0.95)
dif_sel_hpd[2, ] <- HPDinterval(as.mcmc(dif_sel_chn[2, ]), prob = 0.95)

pdf(file = "./Output/Output v1.0/REAL v1.3/REAL_COL_RAW_PMMH1_Posterior_SelChange.pdf", width = 16, height = 16)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
layout(matrix(c(1, 1, 3, 3, 2, 2, 4, 5), nrow = 4, ncol = 2))
grd_num <- 1e+03
sel_cof_pdf <- kde2d(sel_cof_chn[1, 1, ], sel_cof_chn[1, 2, ], n = grd_num)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
      xlab = "Selection coefficient before domestication", ylab = "Selection coefficient after domestication",
      main = "Posterior for selection coefficient of black")
abline(v = sel_cof_est[1, 1], col = 'black', lty = 2, lwd = 2)
abline(h = sel_cof_est[1, 2], col = 'black', lty = 2, lwd = 2)
# HPDregionplot(as.mcmc(sel_cof_chn), vars = 1:2, n = grd_num, prob = 0.95, col = "blue", lwd = 2, add = TRUE)

grd_num <- 1e+03
sel_cof_pdf <- kde2d(sel_cof_chn[2, 1, ], sel_cof_chn[2, 2, ], n = grd_num)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
      xlab = "Selection coefficient before domestication", ylab = "Selection coefficient after domestication",
      main = "Posterior for selection coefficient of chestnut")
abline(v = sel_cof_est[2, 1], col = 'black', lty = 2, lwd = 2)
abline(h = sel_cof_est[2, 2], col = 'black', lty = 2, lwd = 2)
# HPDregionplot(as.mcmc(sel_cof_chn), vars = 1:2, n = grd_num, prob = 0.95, col = "blue", lwd = 2, add = TRUE)

grd_num <- 1e+03
sel_cof_pdf <- kde2d(dif_sel_chn[1, ], dif_sel_chn[2, ], n = grd_num)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
      xlab = "Change in selection coefficient of black", ylab = "Change in selection coefficient of chestnut",
      main = "Posterior for change in selection coefficient before and after domestication")
abline(v = dif_sel_est[1], col = 'black', lty = 2, lwd = 2)
abline(h = dif_sel_est[2], col = 'black', lty = 2, lwd = 2)
# HPDregionplot(as.mcmc(t(dif_sel_chn)), vars = 1:2, n = grd_num, prob = 0.95, col = "blue", lwd = 2, add = TRUE)

hist(dif_sel_chn[1, ], breaks = seq(min(dif_sel_chn[1, ]), max(dif_sel_chn[1, ]), length.out = 50), freq = FALSE,
     xlab = "Change in selection coefficient",
     main = "Posterior for change in selection coefficient of black")
lines(density(dif_sel_chn[1, ]), lwd = 2, col = 'black')
abline(v = dif_sel_est[1], col = 'black', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[1, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[1, 2], col = 'blue', lty = 2, lwd = 2)

hist(dif_sel_chn[2, ], breaks = seq(min(dif_sel_chn[2, ]), max(dif_sel_chn[2, ]), length.out = 50), freq = FALSE,
     xlab = "Change in selection coefficient",
     main = "Posterior for change in selection coefficient of chestnut")
lines(density(dif_sel_chn[2, ]), lwd = 2, col = 'black')
abline(v = dif_sel_est[2], col = 'black', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[2, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[2, 2], col = 'blue', lty = 2, lwd = 2)
dev.off()

########################################

#' Wutke et al. (2016) from 12496 BC (Holocene 9700 BC)
load("./Data/REAL_COL_RAW.rda")

set.seed(2)

sel_cof <- matrix(c(0e+00, 0e+00, 0e+00, 0e+00), nrow = 2, ncol = 2)
rec_rat <- 5e-01
pop_siz <- pop_siz[-(1:(smp_gen[4] - smp_gen[1]))]
ref_siz <- 1.6e+04
evt_gen <- round((-3500 - 2000) / 8) # 3500 BCE
smp_gen <- smp_gen[-(1:3)]
smp_siz <- smp_siz[-(1:3)]
smp_cnt <- smp_cnt[, -(1:3)]
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04
stp_siz <- (1:itn_num)^(-2 / 3)
apt_rto <- 4e-01

system.time(sel_cof_chn <- cmprunAdaptPMMH(sel_cof, rec_rat, pop_siz, ref_siz, evt_gen, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, stp_siz, apt_rto))

save(sel_cof, rec_rat, pop_siz, ref_siz, evt_gen, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, stp_siz, apt_rto, sel_cof_chn,
     file = "./Output/Output v1.0/REAL v1.3/REAL_COL_RAW_PMMH2.rda")

load("./Output/Output v1.0/REAL v1.3/REAL_COL_RAW_PMMH2.rda")

pdf(file = "./Output/Output v1.0/REAL v1.3/REAL_COL_RAW_PMMH2_Traceplot_SelCoeff.pdf", width = 16, height = 8)
par(mfrow = c(2, 2), mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.5, cex.sub = 1.25, cex.axis = 1.25, cex.lab = 1.25)
plot(1:itn_num, sel_cof_chn[1, 1, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient of black before domestication")

plot(1:itn_num, sel_cof_chn[2, 1, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient of chestnut before domestication")

plot(1:itn_num, sel_cof_chn[1, 2, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient of black after domestication")

plot(1:itn_num, sel_cof_chn[2, 2, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient of chestnut after domestication")
dev.off()

# brn_num <- 1e+04
brn_num <- 0.5 * dim(sel_cof_chn)[3] + 1
sel_cof_chn <- sel_cof_chn[, , brn_num:dim(sel_cof_chn)[3]]

thn_num <- 5e+00
sel_cof_chn <- sel_cof_chn[, , (1:round(dim(sel_cof_chn)[3] / thn_num)) * thn_num]

sel_cof_est <- matrix(NA, nrow = 2, ncol = 2)
sel_cof_est[, 1] <- rowMeans(sel_cof_chn[, 1, ])
sel_cof_est[, 2] <- rowMeans(sel_cof_chn[, 2, ])

sel_cof_hpd <- array(NA, dim = c(2, 2, 2))
sel_cof_hpd[1, , 1] <- HPDinterval(as.mcmc(sel_cof_chn[1, 1, ]), prob = 0.95)
sel_cof_hpd[2, , 1] <- HPDinterval(as.mcmc(sel_cof_chn[2, 1, ]), prob = 0.95)
sel_cof_hpd[1, , 2] <- HPDinterval(as.mcmc(sel_cof_chn[1, 2, ]), prob = 0.95)
sel_cof_hpd[2, , 2] <- HPDinterval(as.mcmc(sel_cof_chn[2, 2, ]), prob = 0.95)

pdf(file = "./Output/Output v1.0/REAL v1.3/REAL_COL_RAW_PMMH2_Posterior_SelCoeff.pdf", width = 16, height = 16)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
layout(matrix(c(1, 1, 4, 4, 2, 3, 5, 6), nrow = 4, ncol = 2))
grd_num <- 1e+03
sel_cof_pdf <- kde2d(sel_cof_chn[1, 1, ], sel_cof_chn[2, 1, ], n = grd_num)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
      xlab = "Selection coefficient of black", ylab = "Selection coefficient of chestnut",
      main = "Posterior for selection coefficient before domestication")
abline(v = sel_cof_est[1, 1], col = 'black', lty = 2, lwd = 2)
abline(h = sel_cof_est[2, 1], col = 'black', lty = 2, lwd = 2)
# HPDregionplot(as.mcmc(sel_cof_chn), vars = 1:2, n = grd_num, prob = 0.95, col = "blue", lwd = 2, add = TRUE)

hist(sel_cof_chn[1, 1, ], breaks = seq(min(sel_cof_chn[1, 1, ]), max(sel_cof_chn[1, 1, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient of black before domestication")
lines(density(sel_cof_chn[1, 1, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[1, 1], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 1, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 2, 1], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_chn[2, 1, ], breaks = seq(min(sel_cof_chn[2, 1, ]), max(sel_cof_chn[2, 1, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient of chestnut before domestication")
lines(density(sel_cof_chn[2, 1, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[2, 1], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 1, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 2, 1], col = 'blue', lty = 2, lwd = 2)

grd_num <- 1e+03
sel_cof_pdf <- kde2d(sel_cof_chn[1, 2, ], sel_cof_chn[2, 2, ], n = grd_num)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
      xlab = "Selection coefficient of black", ylab = "Selection coefficient of chestnut",
      main = "Posterior for selection coefficient after domestication")
abline(v = sel_cof_est[1, 2], col = 'black', lty = 2, lwd = 2)
abline(h = sel_cof_est[2, 2], col = 'black', lty = 2, lwd = 2)
# HPDregionplot(as.mcmc(sel_cof_chn), vars = 1:2, n = grd_num, prob = 0.95, col = "blue", lwd = 2, add = TRUE)

hist(sel_cof_chn[1, 2, ], breaks = seq(min(sel_cof_chn[1, 2, ]), max(sel_cof_chn[1, 2, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient of black after domestication")
lines(density(sel_cof_chn[1, 2, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[1, 2], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 1, 2], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 2, 2], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_chn[2, 2, ], breaks = seq(min(sel_cof_chn[2, 2, ]), max(sel_cof_chn[2, 2, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient of chestnut after domestication")
lines(density(sel_cof_chn[2, 2, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[2, 2], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 1, 2], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 2, 2], col = 'blue', lty = 2, lwd = 2)
dev.off()

dif_sel_chn <- sel_cof_chn[, 2, ] - sel_cof_chn[, 1, ]

dif_sel_est <- rowMeans(dif_sel_chn)

dif_sel_hpd <- matrix(NA, nrow = 2, ncol = 2)
dif_sel_hpd[1, ] <- HPDinterval(as.mcmc(dif_sel_chn[1, ]), prob = 0.95)
dif_sel_hpd[2, ] <- HPDinterval(as.mcmc(dif_sel_chn[2, ]), prob = 0.95)

pdf(file = "./Output/Output v1.0/REAL v1.3/REAL_COL_RAW_PMMH2_Posterior_SelChange.pdf", width = 16, height = 16)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
layout(matrix(c(1, 1, 3, 3, 2, 2, 4, 5), nrow = 4, ncol = 2))
grd_num <- 1e+03
sel_cof_pdf <- kde2d(sel_cof_chn[1, 1, ], sel_cof_chn[1, 2, ], n = grd_num)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
      xlab = "Selection coefficient before domestication", ylab = "Selection coefficient after domestication",
      main = "Posterior for selection coefficient of black")
abline(v = sel_cof_est[1, 1], col = 'black', lty = 2, lwd = 2)
abline(h = sel_cof_est[1, 2], col = 'black', lty = 2, lwd = 2)
# HPDregionplot(as.mcmc(sel_cof_chn), vars = 1:2, n = grd_num, prob = 0.95, col = "blue", lwd = 2, add = TRUE)

grd_num <- 1e+03
sel_cof_pdf <- kde2d(sel_cof_chn[2, 1, ], sel_cof_chn[2, 2, ], n = grd_num)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
      xlab = "Selection coefficient before domestication", ylab = "Selection coefficient after domestication",
      main = "Posterior for selection coefficient of chestnut")
abline(v = sel_cof_est[2, 1], col = 'black', lty = 2, lwd = 2)
abline(h = sel_cof_est[2, 2], col = 'black', lty = 2, lwd = 2)
# HPDregionplot(as.mcmc(sel_cof_chn), vars = 1:2, n = grd_num, prob = 0.95, col = "blue", lwd = 2, add = TRUE)

grd_num <- 1e+03
sel_cof_pdf <- kde2d(dif_sel_chn[1, ], dif_sel_chn[2, ], n = grd_num)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
      xlab = "Change in selection coefficient of black", ylab = "Change in selection coefficient of chestnut",
      main = "Posterior for change in selection coefficient before and after domestication")
abline(v = dif_sel_est[1], col = 'black', lty = 2, lwd = 2)
abline(h = dif_sel_est[2], col = 'black', lty = 2, lwd = 2)
# HPDregionplot(as.mcmc(t(dif_sel_chn)), vars = 1:2, n = grd_num, prob = 0.95, col = "blue", lwd = 2, add = TRUE)

hist(dif_sel_chn[1, ], breaks = seq(min(dif_sel_chn[1, ]), max(dif_sel_chn[1, ]), length.out = 50), freq = FALSE,
     xlab = "Change in selection coefficient",
     main = "Posterior for change in selection coefficient of black")
lines(density(dif_sel_chn[1, ]), lwd = 2, col = 'black')
abline(v = dif_sel_est[1], col = 'black', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[1, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[1, 2], col = 'blue', lty = 2, lwd = 2)

hist(dif_sel_chn[2, ], breaks = seq(min(dif_sel_chn[2, ]), max(dif_sel_chn[2, ]), length.out = 50), freq = FALSE,
     xlab = "Change in selection coefficient",
     main = "Posterior for change in selection coefficient of chestnut")
lines(density(dif_sel_chn[2, ]), lwd = 2, col = 'black')
abline(v = dif_sel_est[2], col = 'black', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[2, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[2, 2], col = 'blue', lty = 2, lwd = 2)
dev.off()

########################################

#' Wutke et al. (2016) from 12496 BC (Holocene 9700 BC)
load("./Data/REAL_COL_RAW.rda")

set.seed(1)

sel_cof <- matrix(c(0e+00, 0e+00, 0e+00, 0e+00), nrow = 2, ncol = 2)
rec_rat <- 5e-01
pop_siz <- pop_siz[-(1:(smp_gen[4] - smp_gen[1]))]
pop_siz <- rep(1.6e+04, length.out = length(pop_siz))
ref_siz <- 1.6e+04
evt_gen <- round((-3500 - 2000) / 8) # 3500 BCE
smp_gen <- smp_gen[-(1:3)]
smp_siz <- smp_siz[-(1:3)]
smp_cnt <- smp_cnt[, -(1:3)]
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04
stp_siz <- (1:itn_num)^(-2 / 3)
apt_rto <- 4e-01

system.time(sel_cof_chn <- cmprunAdaptPMMH(sel_cof, rec_rat, pop_siz, ref_siz, evt_gen, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, stp_siz, apt_rto))

save(sel_cof, rec_rat, pop_siz, ref_siz, evt_gen, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, stp_siz, apt_rto, sel_cof_chn,
     file = "./Output/Output v1.0/REAL v1.3/REAL_COL_RAW_PMMH2n.rda")

load("./Output/Output v1.0/REAL v1.3/REAL_COL_RAW_PMMH2n.rda")

pdf(file = "./Output/Output v1.0/REAL v1.3/REAL_COL_RAW_PMMH2n_Traceplot_SelCoeff.pdf", width = 16, height = 8)
par(mfrow = c(2, 2), mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.5, cex.sub = 1.25, cex.axis = 1.25, cex.lab = 1.25)
plot(1:itn_num, sel_cof_chn[1, 1, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient of black before domestication")

plot(1:itn_num, sel_cof_chn[2, 1, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient of chestnut before domestication")

plot(1:itn_num, sel_cof_chn[1, 2, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient of black after domestication")

plot(1:itn_num, sel_cof_chn[2, 2, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient of chestnut after domestication")
dev.off()

# brn_num <- 1e+04
brn_num <- 0.5 * dim(sel_cof_chn)[3] + 1
sel_cof_chn <- sel_cof_chn[, , brn_num:dim(sel_cof_chn)[3]]

thn_num <- 5e+00
sel_cof_chn <- sel_cof_chn[, , (1:round(dim(sel_cof_chn)[3] / thn_num)) * thn_num]

sel_cof_est <- matrix(NA, nrow = 2, ncol = 2)
sel_cof_est[, 1] <- rowMeans(sel_cof_chn[, 1, ])
sel_cof_est[, 2] <- rowMeans(sel_cof_chn[, 2, ])

sel_cof_hpd <- array(NA, dim = c(2, 2, 2))
sel_cof_hpd[1, , 1] <- HPDinterval(as.mcmc(sel_cof_chn[1, 1, ]), prob = 0.95)
sel_cof_hpd[2, , 1] <- HPDinterval(as.mcmc(sel_cof_chn[2, 1, ]), prob = 0.95)
sel_cof_hpd[1, , 2] <- HPDinterval(as.mcmc(sel_cof_chn[1, 2, ]), prob = 0.95)
sel_cof_hpd[2, , 2] <- HPDinterval(as.mcmc(sel_cof_chn[2, 2, ]), prob = 0.95)

pdf(file = "./Output/Output v1.0/REAL v1.3/REAL_COL_RAW_PMMH2n_Posterior_SelCoeff.pdf", width = 16, height = 16)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
layout(matrix(c(1, 1, 4, 4, 2, 3, 5, 6), nrow = 4, ncol = 2))
grd_num <- 1e+03
sel_cof_pdf <- kde2d(sel_cof_chn[1, 1, ], sel_cof_chn[2, 1, ], n = grd_num)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
      xlab = "Selection coefficient of black", ylab = "Selection coefficient of chestnut",
      main = "Posterior for selection coefficient before domestication")
abline(v = sel_cof_est[1, 1], col = 'black', lty = 2, lwd = 2)
abline(h = sel_cof_est[2, 1], col = 'black', lty = 2, lwd = 2)
# HPDregionplot(as.mcmc(sel_cof_chn), vars = 1:2, n = grd_num, prob = 0.95, col = "blue", lwd = 2, add = TRUE)

hist(sel_cof_chn[1, 1, ], breaks = seq(min(sel_cof_chn[1, 1, ]), max(sel_cof_chn[1, 1, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient of black before domestication")
lines(density(sel_cof_chn[1, 1, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[1, 1], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 1, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 2, 1], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_chn[2, 1, ], breaks = seq(min(sel_cof_chn[2, 1, ]), max(sel_cof_chn[2, 1, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient of chestnut before domestication")
lines(density(sel_cof_chn[2, 1, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[2, 1], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 1, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 2, 1], col = 'blue', lty = 2, lwd = 2)

grd_num <- 1e+03
sel_cof_pdf <- kde2d(sel_cof_chn[1, 2, ], sel_cof_chn[2, 2, ], n = grd_num)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
      xlab = "Selection coefficient of black", ylab = "Selection coefficient of chestnut",
      main = "Posterior for selection coefficient after domestication")
abline(v = sel_cof_est[1, 2], col = 'black', lty = 2, lwd = 2)
abline(h = sel_cof_est[2, 2], col = 'black', lty = 2, lwd = 2)
# HPDregionplot(as.mcmc(sel_cof_chn), vars = 1:2, n = grd_num, prob = 0.95, col = "blue", lwd = 2, add = TRUE)

hist(sel_cof_chn[1, 2, ], breaks = seq(min(sel_cof_chn[1, 2, ]), max(sel_cof_chn[1, 2, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient of black after domestication")
lines(density(sel_cof_chn[1, 2, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[1, 2], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 1, 2], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 2, 2], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_chn[2, 2, ], breaks = seq(min(sel_cof_chn[2, 2, ]), max(sel_cof_chn[2, 2, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient of chestnut after domestication")
lines(density(sel_cof_chn[2, 2, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[2, 2], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 1, 2], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 2, 2], col = 'blue', lty = 2, lwd = 2)
dev.off()

dif_sel_chn <- sel_cof_chn[, 2, ] - sel_cof_chn[, 1, ]

dif_sel_est <- rowMeans(dif_sel_chn)

dif_sel_hpd <- matrix(NA, nrow = 2, ncol = 2)
dif_sel_hpd[1, ] <- HPDinterval(as.mcmc(dif_sel_chn[1, ]), prob = 0.95)
dif_sel_hpd[2, ] <- HPDinterval(as.mcmc(dif_sel_chn[2, ]), prob = 0.95)

pdf(file = "./Output/Output v1.0/REAL v1.3/REAL_COL_RAW_PMMH2n_Posterior_SelChange.pdf", width = 16, height = 16)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
layout(matrix(c(1, 1, 3, 3, 2, 2, 4, 5), nrow = 4, ncol = 2))
grd_num <- 1e+03
sel_cof_pdf <- kde2d(sel_cof_chn[1, 1, ], sel_cof_chn[1, 2, ], n = grd_num)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
      xlab = "Selection coefficient before domestication", ylab = "Selection coefficient after domestication",
      main = "Posterior for selection coefficient of black")
abline(v = sel_cof_est[1, 1], col = 'black', lty = 2, lwd = 2)
abline(h = sel_cof_est[1, 2], col = 'black', lty = 2, lwd = 2)
# HPDregionplot(as.mcmc(sel_cof_chn), vars = 1:2, n = grd_num, prob = 0.95, col = "blue", lwd = 2, add = TRUE)

grd_num <- 1e+03
sel_cof_pdf <- kde2d(sel_cof_chn[2, 1, ], sel_cof_chn[2, 2, ], n = grd_num)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
      xlab = "Selection coefficient before domestication", ylab = "Selection coefficient after domestication",
      main = "Posterior for selection coefficient of chestnut")
abline(v = sel_cof_est[2, 1], col = 'black', lty = 2, lwd = 2)
abline(h = sel_cof_est[2, 2], col = 'black', lty = 2, lwd = 2)
# HPDregionplot(as.mcmc(sel_cof_chn), vars = 1:2, n = grd_num, prob = 0.95, col = "blue", lwd = 2, add = TRUE)

grd_num <- 1e+03
sel_cof_pdf <- kde2d(dif_sel_chn[1, ], dif_sel_chn[2, ], n = grd_num)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
      xlab = "Change in selection coefficient of black", ylab = "Change in selection coefficient of chestnut",
      main = "Posterior for change in selection coefficient before and after domestication")
abline(v = dif_sel_est[1], col = 'black', lty = 2, lwd = 2)
abline(h = dif_sel_est[2], col = 'black', lty = 2, lwd = 2)
# HPDregionplot(as.mcmc(t(dif_sel_chn)), vars = 1:2, n = grd_num, prob = 0.95, col = "blue", lwd = 2, add = TRUE)

hist(dif_sel_chn[1, ], breaks = seq(min(dif_sel_chn[1, ]), max(dif_sel_chn[1, ]), length.out = 50), freq = FALSE,
     xlab = "Change in selection coefficient",
     main = "Posterior for change in selection coefficient of black")
lines(density(dif_sel_chn[1, ]), lwd = 2, col = 'black')
abline(v = dif_sel_est[1], col = 'black', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[1, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[1, 2], col = 'blue', lty = 2, lwd = 2)

hist(dif_sel_chn[2, ], breaks = seq(min(dif_sel_chn[2, ]), max(dif_sel_chn[2, ]), length.out = 50), freq = FALSE,
     xlab = "Change in selection coefficient",
     main = "Posterior for change in selection coefficient of chestnut")
lines(density(dif_sel_chn[2, ]), lwd = 2, col = 'black')
abline(v = dif_sel_est[2], col = 'black', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[2, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[2, 2], col = 'blue', lty = 2, lwd = 2)
dev.off()

########################################

#' Wutke et al. (2016) (Holocene 9700 BC)
load("./Data/REAL_COL_GRP.rda")

set.seed(1)

sel_cof <- matrix(c(0e+00, 0e+00, 0e+00, 0e+00), nrow = 2, ncol = 2)
rec_rat <- 5e-01
pop_siz
ref_siz <- 1.6e+04
evt_gen <- round((-3500 - 2000) / 8) # 3500 BCE
smp_gen
smp_siz
smp_cnt
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04
stp_siz <- (1:itn_num)^(-2 / 3)
apt_rto <- 4e-01

system.time(sel_cof_chn <- cmprunAdaptPMMH(sel_cof, rec_rat, pop_siz, ref_siz, evt_gen, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, stp_siz, apt_rto))

save(sel_cof, rec_rat, pop_siz, ref_siz, evt_gen, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, stp_siz, apt_rto, sel_cof_chn,
     file = "./Output/Output v1.0/REAL v1.3/REAL_COL_GRP_PMMH.rda")

load("./Output/Output v1.0/REAL v1.3/REAL_COL_GRP_PMMH.rda")

pdf(file = "./Output/Output v1.0/REAL v1.3/REAL_COL_GRP_PMMH_Traceplot_SelCoeff.pdf", width = 16, height = 8)
par(mfrow = c(2, 2), mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.5, cex.sub = 1.25, cex.axis = 1.25, cex.lab = 1.25)
plot(1:itn_num, sel_cof_chn[1, 1, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient of black before domestication")

plot(1:itn_num, sel_cof_chn[2, 1, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient of chestnut before domestication")

plot(1:itn_num, sel_cof_chn[1, 2, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient of black after domestication")

plot(1:itn_num, sel_cof_chn[2, 2, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient of chestnut after domestication")
dev.off()

# brn_num <- 1e+04
brn_num <- 0.5 * dim(sel_cof_chn)[3] + 1
sel_cof_chn <- sel_cof_chn[, , brn_num:dim(sel_cof_chn)[3]]

thn_num <- 5e+00
sel_cof_chn <- sel_cof_chn[, , (1:round(dim(sel_cof_chn)[3] / thn_num)) * thn_num]

sel_cof_est <- matrix(NA, nrow = 2, ncol = 2)
sel_cof_est[, 1] <- rowMeans(sel_cof_chn[, 1, ])
sel_cof_est[, 2] <- rowMeans(sel_cof_chn[, 2, ])

sel_cof_hpd <- array(NA, dim = c(2, 2, 2))
sel_cof_hpd[1, , 1] <- HPDinterval(as.mcmc(sel_cof_chn[1, 1, ]), prob = 0.95)
sel_cof_hpd[2, , 1] <- HPDinterval(as.mcmc(sel_cof_chn[2, 1, ]), prob = 0.95)
sel_cof_hpd[1, , 2] <- HPDinterval(as.mcmc(sel_cof_chn[1, 2, ]), prob = 0.95)
sel_cof_hpd[2, , 2] <- HPDinterval(as.mcmc(sel_cof_chn[2, 2, ]), prob = 0.95)

pdf(file = "./Output/Output v1.0/REAL v1.3/REAL_COL_GRP_PMMH_Posterior_SelCoeff.pdf", width = 16, height = 16)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
layout(matrix(c(1, 1, 4, 4, 2, 3, 5, 6), nrow = 4, ncol = 2))
grd_num <- 1e+03
sel_cof_pdf <- kde2d(sel_cof_chn[1, 1, ], sel_cof_chn[2, 1, ], n = grd_num)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
      xlab = "Selection coefficient of black", ylab = "Selection coefficient of chestnut",
      main = "Posterior for selection coefficient before domestication")
abline(v = sel_cof_est[1, 1], col = 'black', lty = 2, lwd = 2)
abline(h = sel_cof_est[2, 1], col = 'black', lty = 2, lwd = 2)
# HPDregionplot(as.mcmc(sel_cof_chn), vars = 1:2, n = grd_num, prob = 0.95, col = "blue", lwd = 2, add = TRUE)

hist(sel_cof_chn[1, 1, ], breaks = seq(min(sel_cof_chn[1, 1, ]), max(sel_cof_chn[1, 1, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient of black before domestication")
lines(density(sel_cof_chn[1, 1, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[1, 1], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 1, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 2, 1], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_chn[2, 1, ], breaks = seq(min(sel_cof_chn[2, 1, ]), max(sel_cof_chn[2, 1, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient of chestnut before domestication")
lines(density(sel_cof_chn[2, 1, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[2, 1], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 1, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 2, 1], col = 'blue', lty = 2, lwd = 2)

grd_num <- 1e+03
sel_cof_pdf <- kde2d(sel_cof_chn[1, 2, ], sel_cof_chn[2, 2, ], n = grd_num)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
      xlab = "Selection coefficient of black", ylab = "Selection coefficient of chestnut",
      main = "Posterior for selection coefficient after domestication")
abline(v = sel_cof_est[1, 2], col = 'black', lty = 2, lwd = 2)
abline(h = sel_cof_est[2, 2], col = 'black', lty = 2, lwd = 2)
# HPDregionplot(as.mcmc(sel_cof_chn), vars = 1:2, n = grd_num, prob = 0.95, col = "blue", lwd = 2, add = TRUE)

hist(sel_cof_chn[1, 2, ], breaks = seq(min(sel_cof_chn[1, 2, ]), max(sel_cof_chn[1, 2, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient of black after domestication")
lines(density(sel_cof_chn[1, 2, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[1, 2], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 1, 2], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 2, 2], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_chn[2, 2, ], breaks = seq(min(sel_cof_chn[2, 2, ]), max(sel_cof_chn[2, 2, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient of chestnut after domestication")
lines(density(sel_cof_chn[2, 2, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[2, 2], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 1, 2], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 2, 2], col = 'blue', lty = 2, lwd = 2)
dev.off()

dif_sel_chn <- sel_cof_chn[, 2, ] - sel_cof_chn[, 1, ]

dif_sel_est <- rowMeans(dif_sel_chn)

dif_sel_hpd <- matrix(NA, nrow = 2, ncol = 2)
dif_sel_hpd[1, ] <- HPDinterval(as.mcmc(dif_sel_chn[1, ]), prob = 0.95)
dif_sel_hpd[2, ] <- HPDinterval(as.mcmc(dif_sel_chn[2, ]), prob = 0.95)

pdf(file = "./Output/Output v1.0/REAL v1.3/REAL_COL_GRP_PMMH_Posterior_SelChange.pdf", width = 16, height = 16)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
layout(matrix(c(1, 1, 3, 3, 2, 2, 4, 5), nrow = 4, ncol = 2))
grd_num <- 1e+03
sel_cof_pdf <- kde2d(sel_cof_chn[1, 1, ], sel_cof_chn[1, 2, ], n = grd_num)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
      xlab = "Selection coefficient before domestication", ylab = "Selection coefficient after domestication",
      main = "Posterior for selection coefficient of black")
abline(v = sel_cof_est[1, 1], col = 'black', lty = 2, lwd = 2)
abline(h = sel_cof_est[1, 2], col = 'black', lty = 2, lwd = 2)
# HPDregionplot(as.mcmc(sel_cof_chn), vars = 1:2, n = grd_num, prob = 0.95, col = "blue", lwd = 2, add = TRUE)

grd_num <- 1e+03
sel_cof_pdf <- kde2d(sel_cof_chn[2, 1, ], sel_cof_chn[2, 2, ], n = grd_num)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
      xlab = "Selection coefficient before domestication", ylab = "Selection coefficient after domestication",
      main = "Posterior for selection coefficient of chestnut")
abline(v = sel_cof_est[2, 1], col = 'black', lty = 2, lwd = 2)
abline(h = sel_cof_est[2, 2], col = 'black', lty = 2, lwd = 2)
# HPDregionplot(as.mcmc(sel_cof_chn), vars = 1:2, n = grd_num, prob = 0.95, col = "blue", lwd = 2, add = TRUE)

grd_num <- 1e+03
sel_cof_pdf <- kde2d(dif_sel_chn[1, ], dif_sel_chn[2, ], n = grd_num)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
      xlab = "Change in selection coefficient of black", ylab = "Change in selection coefficient of chestnut",
      main = "Posterior for change in selection coefficient before and after domestication")
abline(v = dif_sel_est[1], col = 'black', lty = 2, lwd = 2)
abline(h = dif_sel_est[2], col = 'black', lty = 2, lwd = 2)
# HPDregionplot(as.mcmc(t(dif_sel_chn)), vars = 1:2, n = grd_num, prob = 0.95, col = "blue", lwd = 2, add = TRUE)

hist(dif_sel_chn[1, ], breaks = seq(min(dif_sel_chn[1, ]), max(dif_sel_chn[1, ]), length.out = 50), freq = FALSE,
     xlab = "Change in selection coefficient",
     main = "Posterior for change in selection coefficient of black")
lines(density(dif_sel_chn[1, ]), lwd = 2, col = 'black')
abline(v = dif_sel_est[1], col = 'black', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[1, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[1, 2], col = 'blue', lty = 2, lwd = 2)

hist(dif_sel_chn[2, ], breaks = seq(min(dif_sel_chn[2, ]), max(dif_sel_chn[2, ]), length.out = 50), freq = FALSE,
     xlab = "Change in selection coefficient",
     main = "Posterior for change in selection coefficient of chestnut")
lines(density(dif_sel_chn[2, ]), lwd = 2, col = 'black')
abline(v = dif_sel_est[2], col = 'black', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[2, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[2, 2], col = 'blue', lty = 2, lwd = 2)
dev.off()

################################################################################
