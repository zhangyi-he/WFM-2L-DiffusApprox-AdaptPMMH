#' @title Estimating temporally variable selection intensity from ancient DNA data with the flexibility of modelling linkage and epistasis
#' @author Zhangyi He, Xiaoyang Dai, Wenyang Lyu, Mark Beaumont, Feng Yu

#' version 1.0
#' Phenotypes controlled by two gene with epistatic interaction
#' Non-constant natural selection and non-constant demographic histories

#' Fix the linkage disequilibrium to be 0

#' Input: genotype likelihoods
#' Output: posteriors for the selection coefficient and the genotype frequency trajectories of the population

#' Horse base coat colours (ASIP & MC1R)

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
source("./RFUN_COL.R")

################################################################################

#' Raw data of Wutke et al. (2016) from 9322 BC (Holocene 9700 BC)
load("./Data/REAL.rda")

set.seed(11)
ASIP_smp <- ASIP
ASIP_smp <- ASIP_smp[which(rowSums(ASIP_smp[, 4:9]) != 0), ]
int_gen <- -round(max(ASIP_smp$age_mean, ASIP_smp$age_lower, ASIP_smp$age_upper) / 8)
lst_gen <- -round(min(ASIP_smp$age_mean, ASIP_smp$age_lower, ASIP_smp$age_upper) / 8)

MC1R_smp <- MC1R
MC1R_smp <- MC1R_smp[which(rowSums(MC1R_smp[, 4:9]) != 0), ]
int_gen <- -round(max(MC1R_smp$age_mean, MC1R_smp$age_lower, MC1R_smp$age_upper) / 8)
lst_gen <- -round(min(MC1R_smp$age_mean, MC1R_smp$age_lower, MC1R_smp$age_upper) / 8)

min_gen <- 9700 + 2000
ASIP_smp <- ASIP_smp[which(ASIP_smp$age_mean <= min_gen), ]
MC1R_smp <- MC1R_smp[which(MC1R_smp$age_mean <= min_gen), ]

ASIP_smp$age_mean <- -round(ASIP_smp$age_mean / 8)
ASIP_smp$age_lower <- -round(ASIP_smp$age_lower / 8)
ASIP_smp$age_upper <- -round(ASIP_smp$age_upper / 8)
# ASIP_smp[which(ASIP_smp[, 7] == 1), 4:5] <- 1 / 2
# ASIP_smp[which(ASIP_smp[, 8] == 1), 5:6] <- 1 / 2
# ASIP_smp[which(ASIP_smp[, 9] == 1), 4:6] <- 1 / 3
ASIP_smp[which(ASIP_smp[, 7] == 1), 4] <- 1 / 2
ASIP_smp[which(ASIP_smp[, 7] == 1), 5] <- 1 / 2
ASIP_smp[which(ASIP_smp[, 8] == 1), 5] <- 1 / 2
ASIP_smp[which(ASIP_smp[, 8] == 1), 6] <- 1 / 2
ASIP_smp[which(ASIP_smp[, 9] == 1), 4] <- 1 / 4
ASIP_smp[which(ASIP_smp[, 9] == 1), 5] <- 1 / 2
ASIP_smp[which(ASIP_smp[, 9] == 1), 6] <- 1 / 4
ASIP_smp <- ASIP_smp[, -c(2, 3, 7, 8, 9)]

MC1R_smp$age_mean <- -round(MC1R_smp$age_mean / 8)
MC1R_smp$age_lower <- -round(MC1R_smp$age_lower / 8)
MC1R_smp$age_upper <- -round(MC1R_smp$age_upper / 8)
# MC1R_smp[which(MC1R_smp[, 7] == 1), 4:5] <- 1 / 2
# MC1R_smp[which(MC1R_smp[, 8] == 1), 5:6] <- 1 / 2
# MC1R_smp[which(MC1R_smp[, 9] == 1), 4:6] <- 1 / 3
MC1R_smp[which(MC1R_smp[, 7] == 1), 4] <- 1 / 2
MC1R_smp[which(MC1R_smp[, 7] == 1), 5] <- 1 / 2
MC1R_smp[which(MC1R_smp[, 8] == 1), 5] <- 1 / 2
MC1R_smp[which(MC1R_smp[, 8] == 1), 6] <- 1 / 2
MC1R_smp[which(MC1R_smp[, 9] == 1), 4] <- 1 / 4
MC1R_smp[which(MC1R_smp[, 9] == 1), 5] <- 1 / 2
MC1R_smp[which(MC1R_smp[, 9] == 1), 6] <- 1 / 4
MC1R_smp <- MC1R_smp[, -c(2, 3, 7, 8, 9)]

raw_smp <- cbind(ASIP_smp[, 1:4], MC1R_smp[, 2:4])
raw_smp <- raw_smp[order(raw_smp$age_mean), ]
rownames(raw_smp) <- NULL

sel_cof <- matrix(c(0e+00, 0e+00, 0e+00, 0e+00), nrow = 2, ncol = 2)
rec_rat <- 5e-01
pop_siz <- pop_siz[min(raw_smp$age_mean - int_gen + 1):max(raw_smp$age_mean - int_gen + 1)]
ref_siz <- tail(pop_siz, n = 1)
evt_gen <- round((-3500 - 2000) / 8) # 3500 BC (domestication)
raw_smp
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04
stp_siz <- (1:itn_num)^(-2 / 3)
apt_rto <- 4e-01

system.time(PMMH <- cmprunAdaptPMMH(sel_cof, rec_rat, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, stp_siz, apt_rto))

save(sel_cof, rec_rat, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, stp_siz, apt_rto, PMMH,
     file = "./REAL_COL_1.rda")

load("./REAL_COL_1.rda")

sel_cof_chn <- PMMH$sel_cof_chn
frq_pth_chn <- PMMH$frq_pth_chn

pdf(file = "./REAL_COL_1_Traceplot_SelCoeff.pdf", width = 16, height = 12)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:itn_num, sel_cof_chn[1, 1, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient of black before domestication")

plot(1:itn_num, sel_cof_chn[2, 1, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient of chestnut before domestication")

plot(1:itn_num, sel_cof_chn[1, 2, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient of black from domestication")

plot(1:itn_num, sel_cof_chn[2, 2, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient of chestnut from domestication")
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

pdf(file = "./REAL_COL_1_Posterior_SelCoeff.pdf", width = 16, height = 12)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
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

hist(sel_cof_chn[1, 2, ], breaks = seq(min(sel_cof_chn[1, 2, ]), max(sel_cof_chn[1, 2, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient of black from domestication")
lines(density(sel_cof_chn[1, 2, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[1, 2], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 1, 2], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 2, 2], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_chn[2, 2, ], breaks = seq(min(sel_cof_chn[2, 2, ]), max(sel_cof_chn[2, 2, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient of chestnut from domestication")
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

pdf(file = "./REAL_COL_1_Posterior_SelChange.pdf", width = 16, height = 6)
par(mfrow = c(1, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
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

frq_pth_chn <- frq_pth_chn[, , brn_num:dim(frq_pth_chn)[3]]
frq_pth_chn <- frq_pth_chn[, , (1:round(dim(frq_pth_chn)[3] / thn_num)) * thn_num]

frq_pth_est <- matrix(NA, nrow = dim(frq_pth_chn)[1], ncol = dim(frq_pth_chn)[2])
for (i in 1:dim(frq_pth_chn)[1]) {
  frq_pth_est[i, ] <- rowMeans(frq_pth_chn[i, , ])
}

frq_pth_hpd <- array(NA, dim = c(dim(frq_pth_chn)[1], 2, dim(frq_pth_chn)[2]))
for (i in 1:dim(frq_pth_chn)[1]) {
  for (j in 1:dim(frq_pth_chn)[2]) {
    frq_pth_hpd[i, , j] <- HPDinterval(as.mcmc(frq_pth_chn[i, j, ]), prob = 0.95)
  }
}

pdf(file = "./REAL_COL_1_Posterior_HaploTraj.pdf", width = 16, height = 12)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(0, type = 'n', xlim = c(min(raw_smp$age_mean), max(raw_smp$age_mean)), ylim = c(min(frq_pth_chn[1, , ]), max(frq_pth_chn[1, , ])),
     xlab = "Generation", ylab = "Haplotype frequency",
     main = "Posterior for underlying trajectory of AE")
for (i in 1:dim(frq_pth_chn)[3]) {
  lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_chn[1, , i], col = 'grey', lty = 1, lwd = 2)
}
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_est[1, ], col = 'black', lty = 2, lwd = 2)
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_hpd[1, 1, ], col = 'blue', lty = 2, lwd = 2)
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_hpd[1, 2, ], col = 'blue', lty = 2, lwd = 2)

plot(0, type = 'n', xlim = c(min(raw_smp$age_mean), max(raw_smp$age_mean)), ylim = c(min(frq_pth_chn[2, , ]), max(frq_pth_chn[2, , ])),
     xlab = "Generation", ylab = "Haplotype frequency",
     main = "Posterior for underlying trajectory of Ae")
for (i in 1:dim(frq_pth_chn)[3]) {
  lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_chn[2, , i], col = 'grey', lty = 1, lwd = 2)
}
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_est[2, ], col = 'black', lty = 2, lwd = 2)
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_hpd[2, 1, ], col = 'blue', lty = 2, lwd = 2)
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_hpd[2, 2, ], col = 'blue', lty = 2, lwd = 2)

plot(0, type = 'n', xlim = c(min(raw_smp$age_mean), max(raw_smp$age_mean)), ylim = c(min(frq_pth_chn[3, , ]), max(frq_pth_chn[3, , ])),
     xlab = "Generation", ylab = "Haplotype frequency",
     main = "Posterior for underlying trajectory of aE")
for (i in 1:dim(frq_pth_chn)[3]) {
  lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_chn[3, , i], col = 'grey', lty = 1, lwd = 2)
}
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_est[3, ], col = 'black', lty = 2, lwd = 2)
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_hpd[3, 1, ], col = 'blue', lty = 2, lwd = 2)
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_hpd[3, 2, ], col = 'blue', lty = 2, lwd = 2)

plot(0, type = 'n', xlim = c(min(raw_smp$age_mean), max(raw_smp$age_mean)), ylim = c(min(frq_pth_chn[4, , ]), max(frq_pth_chn[4, , ])),
     xlab = "Generation", ylab = "Haplotype frequency",
     main = "Posterior for underlying trajectory of ae")
for (i in 1:dim(frq_pth_chn)[3]) {
  lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_chn[4, , i], col = 'grey', lty = 1, lwd = 2)
}
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_est[4, ], col = 'black', lty = 2, lwd = 2)
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_hpd[4, 1, ], col = 'blue', lty = 2, lwd = 2)
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_hpd[4, 2, ], col = 'blue', lty = 2, lwd = 2)
dev.off()

frq_pth_chn <- cmpconvertHaploFreq(frq_pth_chn, sel_cof_chn, evt_gen, min(raw_smp$age_mean), max(raw_smp$age_mean))$phe_pth

frq_pth_est <- matrix(NA, nrow = dim(frq_pth_chn)[1], ncol = dim(frq_pth_chn)[2])
for (i in 1:dim(frq_pth_chn)[1]) {
  frq_pth_est[i, ] <- rowMeans(frq_pth_chn[i, , ])
}

frq_pth_hpd <- array(NA, dim = c(dim(frq_pth_chn)[1], 2, dim(frq_pth_chn)[2]))
for (i in 1:dim(frq_pth_chn)[1]) {
  for (j in 1:dim(frq_pth_chn)[2]) {
    frq_pth_hpd[i, , j] <- HPDinterval(as.mcmc(frq_pth_chn[i, j, ]), prob = 0.95)
  }
}

pdf(file = "./REAL_COL_1_Posterior_PhenoTraj.pdf", width = 16, height = 12)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(0, type = 'n', xlim = c(min(raw_smp$age_mean), max(raw_smp$age_mean)), ylim = c(min(frq_pth_chn[1, , ]), max(frq_pth_chn[1, , ])),
     xlab = "Generation", ylab = "Phenotype frequency",
     main = "Posterior for underlying trajectory of bay")
for (i in 1:dim(frq_pth_chn)[3]) {
  lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_chn[1, , i], col = 'grey', lty = 1, lwd = 2)
}
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_est[1, ], col = 'black', lty = 2, lwd = 2)
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_hpd[1, 1, ], col = 'blue', lty = 2, lwd = 2)
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_hpd[1, 2, ], col = 'blue', lty = 2, lwd = 2)

plot(0, type = 'n', xlim = c(min(raw_smp$age_mean), max(raw_smp$age_mean)), ylim = c(min(frq_pth_chn[2, , ]), max(frq_pth_chn[2, , ])),
     xlab = "Generation", ylab = "Phenotype frequency",
     main = "Posterior for underlying trajectory of black")
for (i in 1:dim(frq_pth_chn)[3]) {
  lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_chn[2, , i], col = 'grey', lty = 1, lwd = 2)
}
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_est[2, ], col = 'black', lty = 2, lwd = 2)
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_hpd[2, 1, ], col = 'blue', lty = 2, lwd = 2)
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_hpd[2, 2, ], col = 'blue', lty = 2, lwd = 2)

plot(0, type = 'n', xlim = c(min(raw_smp$age_mean), max(raw_smp$age_mean)), ylim = c(min(frq_pth_chn[3, , ]), max(frq_pth_chn[3, , ])),
     xlab = "Generation", ylab = "Phenotype frequency",
     main = "Posterior for underlying trajectory of chestnut")
for (i in 1:dim(frq_pth_chn)[3]) {
  lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_chn[3, , i], col = 'grey', lty = 1, lwd = 2)
}
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_est[3, ], col = 'black', lty = 2, lwd = 2)
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_hpd[3, 1, ], col = 'blue', lty = 2, lwd = 2)
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_hpd[3, 2, ], col = 'blue', lty = 2, lwd = 2)
dev.off()

########################################

#' Raw data of Wutke et al. (2016) from 12500 BC (the sampling time point that the ASIP mutation was first found in the sample)
load("./Data/REAL.rda")

set.seed(21)
ASIP_smp <- ASIP
ASIP_smp <- ASIP_smp[which(rowSums(ASIP_smp[, 4:9]) != 0), ]
int_gen <- -round(max(ASIP_smp$age_mean, ASIP_smp$age_lower, ASIP_smp$age_upper) / 8)
lst_gen <- -round(min(ASIP_smp$age_mean, ASIP_smp$age_lower, ASIP_smp$age_upper) / 8)

MC1R_smp <- MC1R
MC1R_smp <- MC1R_smp[which(rowSums(MC1R_smp[, 4:9]) != 0), ]
int_gen <- -round(max(MC1R_smp$age_mean, MC1R_smp$age_lower, MC1R_smp$age_upper) / 8)
lst_gen <- -round(min(MC1R_smp$age_mean, MC1R_smp$age_lower, MC1R_smp$age_upper) / 8)

min_gen <- max(ASIP_smp$age_mean[which(rowSums(ASIP_smp[, c(5, 6, 8)]) != 0)])
ASIP_smp <- ASIP_smp[which(ASIP_smp$age_mean <= min_gen), ]
MC1R_smp <- MC1R_smp[which(MC1R_smp$age_mean <= min_gen), ]

ASIP_smp$age_mean <- -round(ASIP_smp$age_mean / 8)
ASIP_smp$age_lower <- -round(ASIP_smp$age_lower / 8)
ASIP_smp$age_upper <- -round(ASIP_smp$age_upper / 8)
# ASIP_smp[which(ASIP_smp[, 7] == 1), 4:5] <- 1 / 2
# ASIP_smp[which(ASIP_smp[, 8] == 1), 5:6] <- 1 / 2
# ASIP_smp[which(ASIP_smp[, 9] == 1), 4:6] <- 1 / 3
ASIP_smp[which(ASIP_smp[, 7] == 1), 4] <- 1 / 2
ASIP_smp[which(ASIP_smp[, 7] == 1), 5] <- 1 / 2
ASIP_smp[which(ASIP_smp[, 8] == 1), 5] <- 1 / 2
ASIP_smp[which(ASIP_smp[, 8] == 1), 6] <- 1 / 2
ASIP_smp[which(ASIP_smp[, 9] == 1), 4] <- 1 / 4
ASIP_smp[which(ASIP_smp[, 9] == 1), 5] <- 1 / 2
ASIP_smp[which(ASIP_smp[, 9] == 1), 6] <- 1 / 4
ASIP_smp <- ASIP_smp[, -c(2, 3, 7, 8, 9)]

MC1R_smp$age_mean <- -round(MC1R_smp$age_mean / 8)
MC1R_smp$age_lower <- -round(MC1R_smp$age_lower / 8)
MC1R_smp$age_upper <- -round(MC1R_smp$age_upper / 8)
# MC1R_smp[which(MC1R_smp[, 7] == 1), 4:5] <- 1 / 2
# MC1R_smp[which(MC1R_smp[, 8] == 1), 5:6] <- 1 / 2
# MC1R_smp[which(MC1R_smp[, 9] == 1), 4:6] <- 1 / 3
MC1R_smp[which(MC1R_smp[, 7] == 1), 4] <- 1 / 2
MC1R_smp[which(MC1R_smp[, 7] == 1), 5] <- 1 / 2
MC1R_smp[which(MC1R_smp[, 8] == 1), 5] <- 1 / 2
MC1R_smp[which(MC1R_smp[, 8] == 1), 6] <- 1 / 2
MC1R_smp[which(MC1R_smp[, 9] == 1), 4] <- 1 / 4
MC1R_smp[which(MC1R_smp[, 9] == 1), 5] <- 1 / 2
MC1R_smp[which(MC1R_smp[, 9] == 1), 6] <- 1 / 4
MC1R_smp <- MC1R_smp[, -c(2, 3, 7, 8, 9)]

raw_smp <- cbind(ASIP_smp[, 1:4], MC1R_smp[, 2:4])
raw_smp <- raw_smp[order(raw_smp$age_mean), ]
rownames(raw_smp) <- NULL

sel_cof <- matrix(c(0e+00, 0e+00, 0e+00, 0e+00), nrow = 2, ncol = 2)
rec_rat <- 5e-01
pop_siz <- pop_siz[min(raw_smp$age_mean - int_gen + 1):max(raw_smp$age_mean - int_gen + 1)]
ref_siz <- tail(pop_siz, n = 1)
evt_gen <- round((-3500 - 2000) / 8) # 3500 BC (domestication)
raw_smp
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04
stp_siz <- (1:itn_num)^(-2 / 3)
apt_rto <- 4e-01

system.time(PMMH <- cmprunAdaptPMMH(sel_cof, rec_rat, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, stp_siz, apt_rto))

save(sel_cof, rec_rat, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, stp_siz, apt_rto, PMMH,
     file = "./REAL_COL_2.rda")

load("./REAL_COL_2.rda")

sel_cof_chn <- PMMH$sel_cof_chn
frq_pth_chn <- PMMH$frq_pth_chn

pdf(file = "./REAL_COL_2_Traceplot_SelCoeff.pdf", width = 16, height = 12)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:itn_num, sel_cof_chn[1, 1, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient of black before domestication")

plot(1:itn_num, sel_cof_chn[2, 1, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient of chestnut before domestication")

plot(1:itn_num, sel_cof_chn[1, 2, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient of black from domestication")

plot(1:itn_num, sel_cof_chn[2, 2, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient of chestnut from domestication")
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

pdf(file = "./REAL_COL_2_Posterior_SelCoeff.pdf", width = 16, height = 12)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
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

hist(sel_cof_chn[1, 2, ], breaks = seq(min(sel_cof_chn[1, 2, ]), max(sel_cof_chn[1, 2, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient of black from domestication")
lines(density(sel_cof_chn[1, 2, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[1, 2], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 1, 2], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 2, 2], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_chn[2, 2, ], breaks = seq(min(sel_cof_chn[2, 2, ]), max(sel_cof_chn[2, 2, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient of chestnut from domestication")
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

pdf(file = "./REAL_COL_2_Posterior_SelChange.pdf", width = 16, height = 6)
par(mfrow = c(1, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
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

frq_pth_chn <- frq_pth_chn[, , brn_num:dim(frq_pth_chn)[3]]
frq_pth_chn <- frq_pth_chn[, , (1:round(dim(frq_pth_chn)[3] / thn_num)) * thn_num]

frq_pth_est <- matrix(NA, nrow = dim(frq_pth_chn)[1], ncol = dim(frq_pth_chn)[2])
for (i in 1:dim(frq_pth_chn)[1]) {
  frq_pth_est[i, ] <- rowMeans(frq_pth_chn[i, , ])
}

frq_pth_hpd <- array(NA, dim = c(dim(frq_pth_chn)[1], 2, dim(frq_pth_chn)[2]))
for (i in 1:dim(frq_pth_chn)[1]) {
  for (j in 1:dim(frq_pth_chn)[2]) {
    frq_pth_hpd[i, , j] <- HPDinterval(as.mcmc(frq_pth_chn[i, j, ]), prob = 0.95)
  }
}

pdf(file = "./REAL_COL_2_Posterior_HaploTraj.pdf", width = 16, height = 12)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(0, type = 'n', xlim = c(min(raw_smp$age_mean), max(raw_smp$age_mean)), ylim = c(min(frq_pth_chn[1, , ]), max(frq_pth_chn[1, , ])),
     xlab = "Generation", ylab = "Haplotype frequency",
     main = "Posterior for underlying trajectory of AE")
for (i in 1:dim(frq_pth_chn)[3]) {
  lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_chn[1, , i], col = 'grey', lty = 1, lwd = 2)
}
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_est[1, ], col = 'black', lty = 2, lwd = 2)
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_hpd[1, 1, ], col = 'blue', lty = 2, lwd = 2)
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_hpd[1, 2, ], col = 'blue', lty = 2, lwd = 2)

plot(0, type = 'n', xlim = c(min(raw_smp$age_mean), max(raw_smp$age_mean)), ylim = c(min(frq_pth_chn[2, , ]), max(frq_pth_chn[2, , ])),
     xlab = "Generation", ylab = "Haplotype frequency",
     main = "Posterior for underlying trajectory of Ae")
for (i in 1:dim(frq_pth_chn)[3]) {
  lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_chn[2, , i], col = 'grey', lty = 1, lwd = 2)
}
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_est[2, ], col = 'black', lty = 2, lwd = 2)
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_hpd[2, 1, ], col = 'blue', lty = 2, lwd = 2)
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_hpd[2, 2, ], col = 'blue', lty = 2, lwd = 2)

plot(0, type = 'n', xlim = c(min(raw_smp$age_mean), max(raw_smp$age_mean)), ylim = c(min(frq_pth_chn[3, , ]), max(frq_pth_chn[3, , ])),
     xlab = "Generation", ylab = "Haplotype frequency",
     main = "Posterior for underlying trajectory of aE")
for (i in 1:dim(frq_pth_chn)[3]) {
  lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_chn[3, , i], col = 'grey', lty = 1, lwd = 2)
}
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_est[3, ], col = 'black', lty = 2, lwd = 2)
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_hpd[3, 1, ], col = 'blue', lty = 2, lwd = 2)
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_hpd[3, 2, ], col = 'blue', lty = 2, lwd = 2)

plot(0, type = 'n', xlim = c(min(raw_smp$age_mean), max(raw_smp$age_mean)), ylim = c(min(frq_pth_chn[4, , ]), max(frq_pth_chn[4, , ])),
     xlab = "Generation", ylab = "Haplotype frequency",
     main = "Posterior for underlying trajectory of ae")
for (i in 1:dim(frq_pth_chn)[3]) {
  lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_chn[4, , i], col = 'grey', lty = 1, lwd = 2)
}
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_est[4, ], col = 'black', lty = 2, lwd = 2)
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_hpd[4, 1, ], col = 'blue', lty = 2, lwd = 2)
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_hpd[4, 2, ], col = 'blue', lty = 2, lwd = 2)
dev.off()

frq_pth_chn <- cmpconvertHaploFreq(frq_pth_chn, sel_cof_chn, evt_gen, min(raw_smp$age_mean), max(raw_smp$age_mean))$phe_pth

frq_pth_est <- matrix(NA, nrow = dim(frq_pth_chn)[1], ncol = dim(frq_pth_chn)[2])
for (i in 1:dim(frq_pth_chn)[1]) {
  frq_pth_est[i, ] <- rowMeans(frq_pth_chn[i, , ])
}

frq_pth_hpd <- array(NA, dim = c(dim(frq_pth_chn)[1], 2, dim(frq_pth_chn)[2]))
for (i in 1:dim(frq_pth_chn)[1]) {
  for (j in 1:dim(frq_pth_chn)[2]) {
    frq_pth_hpd[i, , j] <- HPDinterval(as.mcmc(frq_pth_chn[i, j, ]), prob = 0.95)
  }
}

pdf(file = "./REAL_COL_2_Posterior_PhenoTraj.pdf", width = 16, height = 12)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(0, type = 'n', xlim = c(min(raw_smp$age_mean), max(raw_smp$age_mean)), ylim = c(min(frq_pth_chn[1, , ]), max(frq_pth_chn[1, , ])),
     xlab = "Generation", ylab = "Phenotype frequency",
     main = "Posterior for underlying trajectory of bay")
for (i in 1:dim(frq_pth_chn)[3]) {
  lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_chn[1, , i], col = 'grey', lty = 1, lwd = 2)
}
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_est[1, ], col = 'black', lty = 2, lwd = 2)
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_hpd[1, 1, ], col = 'blue', lty = 2, lwd = 2)
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_hpd[1, 2, ], col = 'blue', lty = 2, lwd = 2)

plot(0, type = 'n', xlim = c(min(raw_smp$age_mean), max(raw_smp$age_mean)), ylim = c(min(frq_pth_chn[2, , ]), max(frq_pth_chn[2, , ])),
     xlab = "Generation", ylab = "Phenotype frequency",
     main = "Posterior for underlying trajectory of black")
for (i in 1:dim(frq_pth_chn)[3]) {
  lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_chn[2, , i], col = 'grey', lty = 1, lwd = 2)
}
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_est[2, ], col = 'black', lty = 2, lwd = 2)
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_hpd[2, 1, ], col = 'blue', lty = 2, lwd = 2)
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_hpd[2, 2, ], col = 'blue', lty = 2, lwd = 2)

plot(0, type = 'n', xlim = c(min(raw_smp$age_mean), max(raw_smp$age_mean)), ylim = c(min(frq_pth_chn[3, , ]), max(frq_pth_chn[3, , ])),
     xlab = "Generation", ylab = "Phenotype frequency",
     main = "Posterior for underlying trajectory of chestnut")
for (i in 1:dim(frq_pth_chn)[3]) {
  lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_chn[3, , i], col = 'grey', lty = 1, lwd = 2)
}
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_est[3, ], col = 'black', lty = 2, lwd = 2)
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_hpd[3, 1, ], col = 'blue', lty = 2, lwd = 2)
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_hpd[3, 2, ], col = 'blue', lty = 2, lwd = 2)
dev.off()

########################################

#' Raw data of Wutke et al. (2016) from 4300 BC (the sampling time point that the MC1R mutation was first found in the sample)
load("./Data/REAL.rda")

set.seed(21)
ASIP_smp <- ASIP
ASIP_smp <- ASIP_smp[which(rowSums(ASIP_smp[, 4:9]) != 0), ]
int_gen <- -round(max(ASIP_smp$age_mean, ASIP_smp$age_lower, ASIP_smp$age_upper) / 8)
lst_gen <- -round(min(ASIP_smp$age_mean, ASIP_smp$age_lower, ASIP_smp$age_upper) / 8)

MC1R_smp <- MC1R
MC1R_smp <- MC1R_smp[which(rowSums(MC1R_smp[, 4:9]) != 0), ]
int_gen <- -round(max(MC1R_smp$age_mean, MC1R_smp$age_lower, MC1R_smp$age_upper) / 8)
lst_gen <- -round(min(MC1R_smp$age_mean, MC1R_smp$age_lower, MC1R_smp$age_upper) / 8)

min_gen <- max(MC1R_smp$age_mean[which(rowSums(MC1R_smp[, c(5, 6, 8)]) != 0)])
ASIP_smp <- ASIP_smp[which(ASIP_smp$age_mean <= min_gen), ]
MC1R_smp <- MC1R_smp[which(MC1R_smp$age_mean <= min_gen), ]

ASIP_smp$age_mean <- -round(ASIP_smp$age_mean / 8)
ASIP_smp$age_lower <- -round(ASIP_smp$age_lower / 8)
ASIP_smp$age_upper <- -round(ASIP_smp$age_upper / 8)
# ASIP_smp[which(ASIP_smp[, 7] == 1), 4:5] <- 1 / 2
# ASIP_smp[which(ASIP_smp[, 8] == 1), 5:6] <- 1 / 2
# ASIP_smp[which(ASIP_smp[, 9] == 1), 4:6] <- 1 / 3
ASIP_smp[which(ASIP_smp[, 7] == 1), 4] <- 1 / 2
ASIP_smp[which(ASIP_smp[, 7] == 1), 5] <- 1 / 2
ASIP_smp[which(ASIP_smp[, 8] == 1), 5] <- 1 / 2
ASIP_smp[which(ASIP_smp[, 8] == 1), 6] <- 1 / 2
ASIP_smp[which(ASIP_smp[, 9] == 1), 4] <- 1 / 4
ASIP_smp[which(ASIP_smp[, 9] == 1), 5] <- 1 / 2
ASIP_smp[which(ASIP_smp[, 9] == 1), 6] <- 1 / 4
ASIP_smp <- ASIP_smp[, -c(2, 3, 7, 8, 9)]

MC1R_smp$age_mean <- -round(MC1R_smp$age_mean / 8)
MC1R_smp$age_lower <- -round(MC1R_smp$age_lower / 8)
MC1R_smp$age_upper <- -round(MC1R_smp$age_upper / 8)
# MC1R_smp[which(MC1R_smp[, 7] == 1), 4:5] <- 1 / 2
# MC1R_smp[which(MC1R_smp[, 8] == 1), 5:6] <- 1 / 2
# MC1R_smp[which(MC1R_smp[, 9] == 1), 4:6] <- 1 / 3
MC1R_smp[which(MC1R_smp[, 7] == 1), 4] <- 1 / 2
MC1R_smp[which(MC1R_smp[, 7] == 1), 5] <- 1 / 2
MC1R_smp[which(MC1R_smp[, 8] == 1), 5] <- 1 / 2
MC1R_smp[which(MC1R_smp[, 8] == 1), 6] <- 1 / 2
MC1R_smp[which(MC1R_smp[, 9] == 1), 4] <- 1 / 4
MC1R_smp[which(MC1R_smp[, 9] == 1), 5] <- 1 / 2
MC1R_smp[which(MC1R_smp[, 9] == 1), 6] <- 1 / 4
MC1R_smp <- MC1R_smp[, -c(2, 3, 7, 8, 9)]

raw_smp <- cbind(ASIP_smp[, 1:4], MC1R_smp[, 2:4])
raw_smp <- raw_smp[order(raw_smp$age_mean), ]
rownames(raw_smp) <- NULL

sel_cof <- matrix(c(0e+00, 0e+00, 0e+00, 0e+00), nrow = 2, ncol = 2)
rec_rat <- 5e-01
pop_siz <- pop_siz[min(raw_smp$age_mean - int_gen + 1):max(raw_smp$age_mean - int_gen + 1)]
ref_siz <- tail(pop_siz, n = 1)
evt_gen <- round((-3500 - 2000) / 8) # 3500 BC (domestication)
raw_smp
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04
stp_siz <- (1:itn_num)^(-2 / 3)
apt_rto <- 4e-01

system.time(PMMH <- cmprunAdaptPMMH(sel_cof, rec_rat, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, stp_siz, apt_rto))

save(sel_cof, rec_rat, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, stp_siz, apt_rto, PMMH,
     file = "./REAL_COL_3.rda")

load("./REAL_COL_3.rda")

sel_cof_chn <- PMMH$sel_cof_chn
frq_pth_chn <- PMMH$frq_pth_chn

pdf(file = "./REAL_COL_3_Traceplot_SelCoeff.pdf", width = 16, height = 12)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:itn_num, sel_cof_chn[1, 1, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient of black before domestication")

plot(1:itn_num, sel_cof_chn[2, 1, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient of chestnut before domestication")

plot(1:itn_num, sel_cof_chn[1, 2, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient of black from domestication")

plot(1:itn_num, sel_cof_chn[2, 2, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient of chestnut from domestication")
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

pdf(file = "./REAL_COL_3_Posterior_SelCoeff.pdf", width = 16, height = 12)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
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

hist(sel_cof_chn[1, 2, ], breaks = seq(min(sel_cof_chn[1, 2, ]), max(sel_cof_chn[1, 2, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient of black from domestication")
lines(density(sel_cof_chn[1, 2, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[1, 2], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 1, 2], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 2, 2], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_chn[2, 2, ], breaks = seq(min(sel_cof_chn[2, 2, ]), max(sel_cof_chn[2, 2, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient of chestnut from domestication")
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

pdf(file = "./REAL_COL_3_Posterior_SelChange.pdf", width = 16, height = 6)
par(mfrow = c(1, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
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

frq_pth_chn <- frq_pth_chn[, , brn_num:dim(frq_pth_chn)[3]]
frq_pth_chn <- frq_pth_chn[, , (1:round(dim(frq_pth_chn)[3] / thn_num)) * thn_num]

frq_pth_est <- matrix(NA, nrow = dim(frq_pth_chn)[1], ncol = dim(frq_pth_chn)[2])
for (i in 1:dim(frq_pth_chn)[1]) {
  frq_pth_est[i, ] <- rowMeans(frq_pth_chn[i, , ])
}

frq_pth_hpd <- array(NA, dim = c(dim(frq_pth_chn)[1], 2, dim(frq_pth_chn)[2]))
for (i in 1:dim(frq_pth_chn)[1]) {
  for (j in 1:dim(frq_pth_chn)[2]) {
    frq_pth_hpd[i, , j] <- HPDinterval(as.mcmc(frq_pth_chn[i, j, ]), prob = 0.95)
  }
}

pdf(file = "./REAL_COL_3_Posterior_HaploTraj.pdf", width = 16, height = 12)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(0, type = 'n', xlim = c(min(raw_smp$age_mean), max(raw_smp$age_mean)), ylim = c(min(frq_pth_chn[1, , ]), max(frq_pth_chn[1, , ])),
     xlab = "Generation", ylab = "Haplotype frequency",
     main = "Posterior for underlying trajectory of AE")
for (i in 1:dim(frq_pth_chn)[3]) {
  lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_chn[1, , i], col = 'grey', lty = 1, lwd = 2)
}
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_est[1, ], col = 'black', lty = 2, lwd = 2)
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_hpd[1, 1, ], col = 'blue', lty = 2, lwd = 2)
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_hpd[1, 2, ], col = 'blue', lty = 2, lwd = 2)

plot(0, type = 'n', xlim = c(min(raw_smp$age_mean), max(raw_smp$age_mean)), ylim = c(min(frq_pth_chn[2, , ]), max(frq_pth_chn[2, , ])),
     xlab = "Generation", ylab = "Haplotype frequency",
     main = "Posterior for underlying trajectory of Ae")
for (i in 1:dim(frq_pth_chn)[3]) {
  lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_chn[2, , i], col = 'grey', lty = 1, lwd = 2)
}
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_est[2, ], col = 'black', lty = 2, lwd = 2)
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_hpd[2, 1, ], col = 'blue', lty = 2, lwd = 2)
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_hpd[2, 2, ], col = 'blue', lty = 2, lwd = 2)

plot(0, type = 'n', xlim = c(min(raw_smp$age_mean), max(raw_smp$age_mean)), ylim = c(min(frq_pth_chn[3, , ]), max(frq_pth_chn[3, , ])),
     xlab = "Generation", ylab = "Haplotype frequency",
     main = "Posterior for underlying trajectory of aE")
for (i in 1:dim(frq_pth_chn)[3]) {
  lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_chn[3, , i], col = 'grey', lty = 1, lwd = 2)
}
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_est[3, ], col = 'black', lty = 2, lwd = 2)
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_hpd[3, 1, ], col = 'blue', lty = 2, lwd = 2)
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_hpd[3, 2, ], col = 'blue', lty = 2, lwd = 2)

plot(0, type = 'n', xlim = c(min(raw_smp$age_mean), max(raw_smp$age_mean)), ylim = c(min(frq_pth_chn[4, , ]), max(frq_pth_chn[4, , ])),
     xlab = "Generation", ylab = "Haplotype frequency",
     main = "Posterior for underlying trajectory of ae")
for (i in 1:dim(frq_pth_chn)[3]) {
  lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_chn[4, , i], col = 'grey', lty = 1, lwd = 2)
}
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_est[4, ], col = 'black', lty = 2, lwd = 2)
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_hpd[4, 1, ], col = 'blue', lty = 2, lwd = 2)
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_hpd[4, 2, ], col = 'blue', lty = 2, lwd = 2)
dev.off()

frq_pth_chn <- cmpconvertHaploFreq(frq_pth_chn, sel_cof_chn, evt_gen, min(raw_smp$age_mean), max(raw_smp$age_mean))$phe_pth

frq_pth_est <- matrix(NA, nrow = dim(frq_pth_chn)[1], ncol = dim(frq_pth_chn)[2])
for (i in 1:dim(frq_pth_chn)[1]) {
  frq_pth_est[i, ] <- rowMeans(frq_pth_chn[i, , ])
}

frq_pth_hpd <- array(NA, dim = c(dim(frq_pth_chn)[1], 2, dim(frq_pth_chn)[2]))
for (i in 1:dim(frq_pth_chn)[1]) {
  for (j in 1:dim(frq_pth_chn)[2]) {
    frq_pth_hpd[i, , j] <- HPDinterval(as.mcmc(frq_pth_chn[i, j, ]), prob = 0.95)
  }
}

pdf(file = "./REAL_COL_3_Posterior_PhenoTraj.pdf", width = 16, height = 12)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(0, type = 'n', xlim = c(min(raw_smp$age_mean), max(raw_smp$age_mean)), ylim = c(min(frq_pth_chn[1, , ]), max(frq_pth_chn[1, , ])),
     xlab = "Generation", ylab = "Phenotype frequency",
     main = "Posterior for underlying trajectory of bay")
for (i in 1:dim(frq_pth_chn)[3]) {
  lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_chn[1, , i], col = 'grey', lty = 1, lwd = 2)
}
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_est[1, ], col = 'black', lty = 2, lwd = 2)
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_hpd[1, 1, ], col = 'blue', lty = 2, lwd = 2)
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_hpd[1, 2, ], col = 'blue', lty = 2, lwd = 2)

plot(0, type = 'n', xlim = c(min(raw_smp$age_mean), max(raw_smp$age_mean)), ylim = c(min(frq_pth_chn[2, , ]), max(frq_pth_chn[2, , ])),
     xlab = "Generation", ylab = "Phenotype frequency",
     main = "Posterior for underlying trajectory of black")
for (i in 1:dim(frq_pth_chn)[3]) {
  lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_chn[2, , i], col = 'grey', lty = 1, lwd = 2)
}
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_est[2, ], col = 'black', lty = 2, lwd = 2)
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_hpd[2, 1, ], col = 'blue', lty = 2, lwd = 2)
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_hpd[2, 2, ], col = 'blue', lty = 2, lwd = 2)

plot(0, type = 'n', xlim = c(min(raw_smp$age_mean), max(raw_smp$age_mean)), ylim = c(min(frq_pth_chn[3, , ]), max(frq_pth_chn[3, , ])),
     xlab = "Generation", ylab = "Phenotype frequency",
     main = "Posterior for underlying trajectory of chestnut")
for (i in 1:dim(frq_pth_chn)[3]) {
  lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_chn[3, , i], col = 'grey', lty = 1, lwd = 2)
}
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_est[3, ], col = 'black', lty = 2, lwd = 2)
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_hpd[3, 1, ], col = 'blue', lty = 2, lwd = 2)
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_hpd[3, 2, ], col = 'blue', lty = 2, lwd = 2)
dev.off()

################################################################################
