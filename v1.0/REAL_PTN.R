#' @title Estimating selection coefficients and testing their changes from ancient DNA data with the flexibility of modelling linkage and epistasis
#' @author Xiaoyang Dai, Wenyang Lyu, Mark Beaumont, Feng Yu, Zhangyi He

#' version 1.0
#' Phenotypes controlled by two gene with genetic linkage
#' Non-constant natural selection and non-constant demographic histories

#' Integrate prior knowledge from modern samples (gene polymorphism)

#' Input: genotype likelihoods
#' Output: posteriors for the selection coefficient, the genotype frequency trajectories of the population and the genotypes of the sample

#' Horse pinto coat patterns (KIT13 & KIT116)

# set the directory
setwd("~/Dropbox/Jeffery He/iResearch/Publications/2020/HE2021-WFM-2L-DiffusApprox-PMMHwGibbs-MolEcolResour")

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
source("./Code/Code v1.0/Code v1.0/RFUN_PTN.R")

################################################################################

#' Raw data of Wutke et al. (2016) from 3500 BC
load("./Data/REAL.rda")

set.seed(88)
KIT13_smp <- KIT13
KIT13_smp <- KIT13_smp[which(rowSums(KIT13_smp[, 4:9]) != 0), ]
int_gen <- -round(max(KIT13_smp$age_mean, KIT13_smp$age_lower, KIT13_smp$age_upper) / 8)
lst_gen <- -round(min(KIT13_smp$age_mean, KIT13_smp$age_lower, KIT13_smp$age_upper) / 8)

KIT16_smp <- KIT16
KIT16_smp <- KIT16_smp[which(rowSums(KIT16_smp[, 4:9]) != 0), ]
int_gen <- -round(max(KIT16_smp$age_mean, KIT16_smp$age_lower, KIT16_smp$age_upper) / 8)
lst_gen <- -round(min(KIT16_smp$age_mean, KIT16_smp$age_lower, KIT16_smp$age_upper) / 8)

min_gen <- min(max(KIT13_smp$age_mean[which(rowSums(KIT13_smp[, c(5, 6, 8)]) != 0)]), max(KIT16_smp$age_mean[which(rowSums(KIT16_smp[, c(5, 6, 8)]) != 0)]))
KIT13_smp <- KIT13_smp[which(KIT13_smp$age_mean <= min_gen), ]
KIT16_smp <- KIT16_smp[which(KIT16_smp$age_mean <= min_gen), ]

KIT13_smp$age_mean <- -round(KIT13_smp$age_mean / 8)
KIT13_smp$age_lower <- -round(KIT13_smp$age_lower / 8)
KIT13_smp$age_upper <- -round(KIT13_smp$age_upper / 8)
# KIT13_smp[which(KIT13_smp[, 7] == 1), 4:5] <- 1 / 2
# KIT13_smp[which(KIT13_smp[, 8] == 1), 5:6] <- 1 / 2
# KIT13_smp[which(KIT13_smp[, 9] == 1), 4:6] <- 1 / 3
KIT13_smp[which(KIT13_smp[, 7] == 1), 4] <- 1 / 2
KIT13_smp[which(KIT13_smp[, 7] == 1), 5] <- 1 / 2
KIT13_smp[which(KIT13_smp[, 8] == 1), 5] <- 1 / 2
KIT13_smp[which(KIT13_smp[, 8] == 1), 6] <- 1 / 2
KIT13_smp[which(KIT13_smp[, 9] == 1), 4] <- 1 / 4
KIT13_smp[which(KIT13_smp[, 9] == 1), 5] <- 1 / 2
KIT13_smp[which(KIT13_smp[, 9] == 1), 6] <- 1 / 4
KIT13_smp <- KIT13_smp[, -(7:9)]

KIT16_smp$age_mean <- -round(KIT16_smp$age_mean / 8)
KIT16_smp$age_lower <- -round(KIT16_smp$age_lower / 8)
KIT16_smp$age_upper <- -round(KIT16_smp$age_upper / 8)
# KIT16_smp[which(KIT16_smp[, 7] == 1), 4:5] <- 1 / 2
# KIT16_smp[which(KIT16_smp[, 8] == 1), 5:6] <- 1 / 2
# KIT16_smp[which(KIT16_smp[, 9] == 1), 4:6] <- 1 / 3
KIT16_smp[which(KIT16_smp[, 7] == 1), 4] <- 1 / 2
KIT16_smp[which(KIT16_smp[, 7] == 1), 5] <- 1 / 2
KIT16_smp[which(KIT16_smp[, 8] == 1), 5] <- 1 / 2
KIT16_smp[which(KIT16_smp[, 8] == 1), 6] <- 1 / 2
KIT16_smp[which(KIT16_smp[, 9] == 1), 4] <- 1 / 4
KIT16_smp[which(KIT16_smp[, 9] == 1), 5] <- 1 / 2
KIT16_smp[which(KIT16_smp[, 9] == 1), 6] <- 1 / 4
KIT16_smp <- KIT16_smp[, -(7:9)]

raw_smp <- cbind(KIT13_smp[, 1:3], KIT13_smp[, 4:6], KIT16_smp[, 4:6])
raw_smp <- raw_smp[order(raw_smp$age_mean), ]
rownames(raw_smp) <- NULL

sel_cof <- matrix(c(0e+00, 0e+00, 0e+00, 0e+00, 0e+00, 0e+00), nrow = 3, ncol = 2)
rec_rat <- 1e-08 * 4688
pop_siz <- pop_siz[min(raw_smp$age_mean - int_gen + 1):max(raw_smp$age_mean - int_gen + 1)]
ref_siz <- tail(pop_siz, n = 1)
evt_gen <- round((400 - 2000) / 8) # 400 AD (the Middle Ages)
raw_smp
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04
stp_siz <- (1:itn_num)^(-2 / 3)
apt_rto <- 4e-01

system.time(PMMH <- cmprunAdaptPMMH(sel_cof, rec_rat, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, stp_siz, apt_rto))

save(sel_cof, rec_rat, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, stp_siz, apt_rto, PMMH,
     file = "./Output/Output v1.0/REAL v1.0/REAL_PTN_1.rda")

load("./Output/Output v1.0/REAL v1.0/REAL_PTN_1.rda")

sel_cof_chn <- PMMH$sel_cof_chn[1:2, , ]
frq_pth_chn <- PMMH$frq_pth_chn
imp_smp_chn <- PMMH$imp_smp_chn

pdf(file = "./Output/Output v1.0/REAL v1.0/REAL_PTN_1_Traceplot_SelCoeff.pdf", width = 24, height = 12)
par(mfrow = c(2, 2), mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.5, cex.sub = 1.25, cex.axis = 1.25, cex.lab = 1.25)
plot(1:itn_num, sel_cof_chn[1, 1, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient of tobiano before the Middle Ages")

plot(1:itn_num, sel_cof_chn[2, 1, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient of sabino before the Middle Ages")

plot(1:itn_num, sel_cof_chn[1, 2, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient of tobiano from the Middle Ages")

plot(1:itn_num, sel_cof_chn[2, 2, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient of sabino from the Middle Ages")
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

pdf(file = "./Output/Output v1.0/REAL v1.0/REAL_PTN_1_Posterior_SelCoeff.pdf", width = 24, height = 24)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
layout(matrix(c(1, 1, 4, 4, 2, 3, 5, 6), nrow = 4, ncol = 2))
grd_num <- 1e+03
sel_cof_pdf <- kde2d(sel_cof_chn[1, 1, ], sel_cof_chn[2, 1, ], n = grd_num)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
      xlab = "Selection coefficient of tobiano", ylab = "Selection coefficient of sabino",
      main = "Posterior for selection coefficient before the Middle Ages")
abline(v = sel_cof_est[1, 1], col = 'black', lty = 2, lwd = 2)
abline(h = sel_cof_est[2, 1], col = 'black', lty = 2, lwd = 2)
# HPDregionplot(as.mcmc(sel_cof_chn), vars = 1:2, n = grd_num, prob = 0.95, col = "blue", lwd = 2, add = TRUE)

hist(sel_cof_chn[1, 1, ], breaks = seq(min(sel_cof_chn[1, 1, ]), max(sel_cof_chn[1, 1, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient of tobiano before the Middle Ages")
lines(density(sel_cof_chn[1, 1, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[1, 1], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 1, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 2, 1], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_chn[2, 1, ], breaks = seq(min(sel_cof_chn[2, 1, ]), max(sel_cof_chn[2, 1, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient of sabino before the Middle Ages")
lines(density(sel_cof_chn[2, 1, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[2, 1], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 1, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 2, 1], col = 'blue', lty = 2, lwd = 2)

grd_num <- 1e+03
sel_cof_pdf <- kde2d(sel_cof_chn[1, 2, ], sel_cof_chn[2, 2, ], n = grd_num)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
      xlab = "Selection coefficient of tobiano", ylab = "Selection coefficient of sabino",
      main = "Posterior for selection coefficient from the Middle Ages")
abline(v = sel_cof_est[1, 2], col = 'black', lty = 2, lwd = 2)
abline(h = sel_cof_est[2, 2], col = 'black', lty = 2, lwd = 2)
# HPDregionplot(as.mcmc(sel_cof_chn), vars = 1:2, n = grd_num, prob = 0.95, col = "blue", lwd = 2, add = TRUE)

hist(sel_cof_chn[1, 2, ], breaks = seq(min(sel_cof_chn[1, 2, ]), max(sel_cof_chn[1, 2, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient of tobiano from the Middle Ages")
lines(density(sel_cof_chn[1, 2, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[1, 2], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 1, 2], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 2, 2], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_chn[2, 2, ], breaks = seq(min(sel_cof_chn[2, 2, ]), max(sel_cof_chn[2, 2, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient of sabino from the Middle Ages")
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

pdf(file = "./Output/Output v1.0/REAL v1.0/REAL_PTN_1_Posterior_SelChange.pdf", width = 24, height = 24)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
layout(matrix(c(1, 1, 3, 3, 2, 2, 4, 5), nrow = 4, ncol = 2))
grd_num <- 1e+03
sel_cof_pdf <- kde2d(sel_cof_chn[1, 1, ], sel_cof_chn[1, 2, ], n = grd_num)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
      xlab = "Selection coefficient before the Middle Ages", ylab = "Selection coefficient from the Middle Ages",
      main = "Posterior for selection coefficient of tobiano")
abline(v = sel_cof_est[1, 1], col = 'black', lty = 2, lwd = 2)
abline(h = sel_cof_est[1, 2], col = 'black', lty = 2, lwd = 2)
# HPDregionplot(as.mcmc(sel_cof_chn), vars = 1:2, n = grd_num, prob = 0.95, col = "blue", lwd = 2, add = TRUE)

grd_num <- 1e+03
sel_cof_pdf <- kde2d(sel_cof_chn[2, 1, ], sel_cof_chn[2, 2, ], n = grd_num)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
      xlab = "Selection coefficient before the Middle Ages", ylab = "Selection coefficient from the Middle Ages",
      main = "Posterior for selection coefficient of sabino")
abline(v = sel_cof_est[2, 1], col = 'black', lty = 2, lwd = 2)
abline(h = sel_cof_est[2, 2], col = 'black', lty = 2, lwd = 2)
# HPDregionplot(as.mcmc(sel_cof_chn), vars = 1:2, n = grd_num, prob = 0.95, col = "blue", lwd = 2, add = TRUE)

grd_num <- 1e+03
sel_cof_pdf <- kde2d(dif_sel_chn[1, ], dif_sel_chn[2, ], n = grd_num)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
      xlab = "Change in selection coefficient of tobiano", ylab = "Change in selection coefficient of sabino",
      main = "Posterior for change in selection coefficient before and from the Middle Ages")
abline(v = dif_sel_est[1], col = 'black', lty = 2, lwd = 2)
abline(h = dif_sel_est[2], col = 'black', lty = 2, lwd = 2)
# HPDregionplot(as.mcmc(t(dif_sel_chn)), vars = 1:2, n = grd_num, prob = 0.95, col = "blue", lwd = 2, add = TRUE)

hist(dif_sel_chn[1, ], breaks = seq(min(dif_sel_chn[1, ]), max(dif_sel_chn[1, ]), length.out = 50), freq = FALSE,
     xlab = "Change in selection coefficient",
     main = "Posterior for change in selection coefficient of tobiano")
lines(density(dif_sel_chn[1, ]), lwd = 2, col = 'black')
abline(v = dif_sel_est[1], col = 'black', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[1, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[1, 2], col = 'blue', lty = 2, lwd = 2)

hist(dif_sel_chn[2, ], breaks = seq(min(dif_sel_chn[2, ]), max(dif_sel_chn[2, ]), length.out = 50), freq = FALSE,
     xlab = "Change in selection coefficient",
     main = "Posterior for change in selection coefficient of sabino")
lines(density(dif_sel_chn[2, ]), lwd = 2, col = 'black')
abline(v = dif_sel_est[2], col = 'black', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[2, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[2, 2], col = 'blue', lty = 2, lwd = 2)
dev.off()

frq_pth_chn <- frq_pth_chn[, , brn_num:dim(frq_pth_chn)[3]]
frq_pth_chn <- frq_pth_chn[, , (1:round(dim(frq_pth_chn)[3] / thn_num)) * thn_num]

frq_pth_est <- matrix(NA, nrow = 4, ncol = dim(frq_pth_chn)[2])
frq_pth_est[1, ] <- rowMeans(frq_pth_chn[1, , ])
frq_pth_est[2, ] <- rowMeans(frq_pth_chn[2, , ])
frq_pth_est[3, ] <- rowMeans(frq_pth_chn[3, , ])
frq_pth_est[4, ] <- rowMeans(frq_pth_chn[4, , ])

frq_pth_hpd <- array(NA, dim = c(4, 2, dim(frq_pth_chn)[2]))
for (i in 1:dim(frq_pth_chn)[2]) {
   frq_pth_hpd[1, , i] <- HPDinterval(as.mcmc(frq_pth_chn[1, i, ]), prob = 0.95)
   frq_pth_hpd[2, , i] <- HPDinterval(as.mcmc(frq_pth_chn[2, i, ]), prob = 0.95)
   frq_pth_hpd[3, , i] <- HPDinterval(as.mcmc(frq_pth_chn[3, i, ]), prob = 0.95)
   frq_pth_hpd[4, , i] <- HPDinterval(as.mcmc(frq_pth_chn[4, i, ]), prob = 0.95)
}

pdf(file = "./Output/Output v1.0/REAL v1.0/REAL_PTN_1_Posterior_Traj.pdf", width = 24, height = 12)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
layout(matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2))
plot(0, type = 'n', xlim = c(min(raw_smp$age_mean), max(raw_smp$age_mean)), ylim = c(min(frq_pth_chn[1, , ]), max(frq_pth_chn[1, , ])),
     xlab = "Generation", ylab = "Haplotype frequency",
     main = "Posterior for underlying trajectory of the KM0sb1 haplotype")

for (i in 1:dim(frq_pth_chn)[2]) {
  lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_chn[1, , i], col = 'grey', lty = 1, lwd = 2)
}
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_est[1, ], col = 'black', lty = 2, lwd = 2)
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_hpd[1, 1, ], col = 'blue', lty = 2, lwd = 2)
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_hpd[1, 2, ], col = 'blue', lty = 2, lwd = 2)

plot(0, type = 'n', xlim = c(min(raw_smp$age_mean), max(raw_smp$age_mean)), ylim = c(min(frq_pth_chn[2, , ]), max(frq_pth_chn[2, , ])),
     xlab = "Generation", ylab = "Haplotype frequency",
     main = "Posterior for underlying trajectory of the KM0SB1 haplotype")

for (i in 1:dim(frq_pth_chn)[2]) {
  lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_chn[2, , i], col = 'grey', lty = 1, lwd = 2)
}
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_est[2, ], col = 'black', lty = 2, lwd = 2)
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_hpd[2, 1, ], col = 'blue', lty = 2, lwd = 2)
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_hpd[2, 2, ], col = 'blue', lty = 2, lwd = 2)

plot(0, type = 'n', xlim = c(min(raw_smp$age_mean), max(raw_smp$age_mean)), ylim = c(min(frq_pth_chn[3, , ]), max(frq_pth_chn[3, , ])),
     xlab = "Generation", ylab = "Haplotype frequency",
     main = "Posterior for underlying trajectory of the KM1sb1 haplotype")

for (i in 1:dim(frq_pth_chn)[2]) {
  lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_chn[3, , i], col = 'grey', lty = 1, lwd = 2)
}
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_est[3, ], col = 'black', lty = 2, lwd = 2)
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_hpd[3, 1, ], col = 'blue', lty = 2, lwd = 2)
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_hpd[3, 2, ], col = 'blue', lty = 2, lwd = 2)

plot(0, type = 'n', xlim = c(min(raw_smp$age_mean), max(raw_smp$age_mean)), ylim = c(min(frq_pth_chn[4, , ]), max(frq_pth_chn[4, , ])),
     xlab = "Generation", ylab = "Haplotype frequency",
     main = "Posterior for underlying trajectory of the KM1SB1 haplotype")

for (i in 1:dim(frq_pth_chn)[2]) {
  lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_chn[4, , i], col = 'grey', lty = 1, lwd = 2)
}
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_est[4, ], col = 'black', lty = 2, lwd = 2)
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_hpd[4, 1, ], col = 'blue', lty = 2, lwd = 2)
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_hpd[4, 2, ], col = 'blue', lty = 2, lwd = 2)
dev.off()

imp_smp_chn <- imp_smp_chn[, , brn_num:dim(imp_smp_chn)[3]]
imp_smp_chn <- imp_smp_chn[, , (1:round(dim(imp_smp_chn)[3] / thn_num)) * thn_num]

imp_smp_est <- matrix(NA, nrow = nrow(raw_smp), ncol = 9)
for (i in 1:nrow(imp_smp_est)) {
  imp_smp_est[i, ] <- rowSums(imp_smp_chn[, i, ]) / dim(imp_smp_chn)[3]
}
imp_smp_est <- as.data.frame(cbind(imp_smp_est[, 1] + imp_smp_est[, 2] + imp_smp_est[, 3],
                                   imp_smp_est[, 4] + imp_smp_est[, 5] + imp_smp_est[, 7],
                                   imp_smp_est[, 6] + imp_smp_est[, 8] + imp_smp_est[, 9],
                                   imp_smp_est[, 1] + imp_smp_est[, 4] + imp_smp_est[, 6],
                                   imp_smp_est[, 2] + imp_smp_est[, 5] + imp_smp_est[, 8],
                                   imp_smp_est[, 3] + imp_smp_est[, 7] + imp_smp_est[, 9]))
imp_smp_est

########################################

#' Raw data of Wutke et al. (2016) from 3645 BC
load("./Data/REAL.rda")

set.seed(28)
KIT13_smp <- KIT13
KIT13_smp <- KIT13_smp[which(rowSums(KIT13_smp[, 4:9]) != 0), ]
int_gen <- -round(max(KIT13_smp$age_mean, KIT13_smp$age_lower, KIT13_smp$age_upper) / 8)
lst_gen <- -round(min(KIT13_smp$age_mean, KIT13_smp$age_lower, KIT13_smp$age_upper) / 8)

KIT16_smp <- KIT16
KIT16_smp <- KIT16_smp[which(rowSums(KIT16_smp[, 4:9]) != 0), ]
int_gen <- -round(max(KIT16_smp$age_mean, KIT16_smp$age_lower, KIT16_smp$age_upper) / 8)
lst_gen <- -round(min(KIT16_smp$age_mean, KIT16_smp$age_lower, KIT16_smp$age_upper) / 8)

min_gen <- max(max(KIT13_smp$age_mean[which(rowSums(KIT13_smp[, c(5, 6, 8)]) != 0)]), max(KIT16_smp$age_mean[which(rowSums(KIT16_smp[, c(5, 6, 8)]) != 0)]))
KIT13_smp <- KIT13_smp[which(KIT13_smp$age_mean <= min_gen), ]
KIT16_smp <- KIT16_smp[which(KIT16_smp$age_mean <= min_gen), ]

KIT13_smp$age_mean <- -round(KIT13_smp$age_mean / 8)
KIT13_smp$age_lower <- -round(KIT13_smp$age_lower / 8)
KIT13_smp$age_upper <- -round(KIT13_smp$age_upper / 8)
# KIT13_smp[which(KIT13_smp[, 7] == 1), 4:5] <- 1 / 2
# KIT13_smp[which(KIT13_smp[, 8] == 1), 5:6] <- 1 / 2
# KIT13_smp[which(KIT13_smp[, 9] == 1), 4:6] <- 1 / 3
KIT13_smp[which(KIT13_smp[, 7] == 1), 4] <- 1 / 2
KIT13_smp[which(KIT13_smp[, 7] == 1), 5] <- 1 / 2
KIT13_smp[which(KIT13_smp[, 8] == 1), 5] <- 1 / 2
KIT13_smp[which(KIT13_smp[, 8] == 1), 6] <- 1 / 2
KIT13_smp[which(KIT13_smp[, 9] == 1), 4] <- 1 / 4
KIT13_smp[which(KIT13_smp[, 9] == 1), 5] <- 1 / 2
KIT13_smp[which(KIT13_smp[, 9] == 1), 6] <- 1 / 4
KIT13_smp <- KIT13_smp[, -(7:9)]

KIT16_smp$age_mean <- -round(KIT16_smp$age_mean / 8)
KIT16_smp$age_lower <- -round(KIT16_smp$age_lower / 8)
KIT16_smp$age_upper <- -round(KIT16_smp$age_upper / 8)
# KIT16_smp[which(KIT16_smp[, 7] == 1), 4:5] <- 1 / 2
# KIT16_smp[which(KIT16_smp[, 8] == 1), 5:6] <- 1 / 2
# KIT16_smp[which(KIT16_smp[, 9] == 1), 4:6] <- 1 / 3
KIT16_smp[which(KIT16_smp[, 7] == 1), 4] <- 1 / 2
KIT16_smp[which(KIT16_smp[, 7] == 1), 5] <- 1 / 2
KIT16_smp[which(KIT16_smp[, 8] == 1), 5] <- 1 / 2
KIT16_smp[which(KIT16_smp[, 8] == 1), 6] <- 1 / 2
KIT16_smp[which(KIT16_smp[, 9] == 1), 4] <- 1 / 4
KIT16_smp[which(KIT16_smp[, 9] == 1), 5] <- 1 / 2
KIT16_smp[which(KIT16_smp[, 9] == 1), 6] <- 1 / 4
KIT16_smp <- KIT16_smp[, -(7:9)]

raw_smp <- cbind(KIT13_smp[, 1:3], KIT13_smp[, 4:6], KIT16_smp[, 4:6])
raw_smp <- raw_smp[order(raw_smp$age_mean), ]
rownames(raw_smp) <- NULL

sel_cof <- matrix(c(0e+00, 0e+00, 0e+00, 0e+00, 0e+00, 0e+00), nrow = 3, ncol = 2)
rec_rat <- 1e-08 * 4688
pop_siz <- pop_siz[min(raw_smp$age_mean - int_gen + 1):max(raw_smp$age_mean - int_gen + 1)]
ref_siz <- tail(pop_siz, n = 1)
evt_gen <- round((400 - 2000) / 8) # 400 AD (the Middle Ages)
raw_smp
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04
stp_siz <- (1:itn_num)^(-2 / 3)
apt_rto <- 4e-01

system.time(PMMH <- cmprunAdaptPMMH(sel_cof, rec_rat, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, stp_siz, apt_rto))

save(sel_cof, rec_rat, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, stp_siz, apt_rto, PMMH,
     file = "./Output/Output v1.0/REAL v1.0/REAL_PTN_2.rda")

load("./Output/Output v1.0/REAL v1.0/REAL_PTN_2.rda")

sel_cof_chn <- PMMH$sel_cof_chn[1:2, , ]
frq_pth_chn <- PMMH$frq_pth_chn
imp_smp_chn <- PMMH$imp_smp_chn

pdf(file = "./Output/Output v1.0/REAL v1.0/REAL_PTN_2_Traceplot_SelCoeff.pdf", width = 24, height = 12)
par(mfrow = c(2, 2), mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.5, cex.sub = 1.25, cex.axis = 1.25, cex.lab = 1.25)
plot(1:itn_num, sel_cof_chn[1, 1, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient of tobiano before the Middle Ages")

plot(1:itn_num, sel_cof_chn[2, 1, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient of sabino before the Middle Ages")

plot(1:itn_num, sel_cof_chn[1, 2, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient of tobiano from the Middle Ages")

plot(1:itn_num, sel_cof_chn[2, 2, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient of sabino from the Middle Ages")
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

pdf(file = "./Output/Output v1.0/REAL v1.0/REAL_PTN_2_Posterior_SelCoeff.pdf", width = 24, height = 24)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
layout(matrix(c(1, 1, 4, 4, 2, 3, 5, 6), nrow = 4, ncol = 2))
grd_num <- 1e+03
sel_cof_pdf <- kde2d(sel_cof_chn[1, 1, ], sel_cof_chn[2, 1, ], n = grd_num)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
      xlab = "Selection coefficient of tobiano", ylab = "Selection coefficient of sabino",
      main = "Posterior for selection coefficient before the Middle Ages")
abline(v = sel_cof_est[1, 1], col = 'black', lty = 2, lwd = 2)
abline(h = sel_cof_est[2, 1], col = 'black', lty = 2, lwd = 2)
# HPDregionplot(as.mcmc(sel_cof_chn), vars = 1:2, n = grd_num, prob = 0.95, col = "blue", lwd = 2, add = TRUE)

hist(sel_cof_chn[1, 1, ], breaks = seq(min(sel_cof_chn[1, 1, ]), max(sel_cof_chn[1, 1, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient of tobiano before the Middle Ages")
lines(density(sel_cof_chn[1, 1, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[1, 1], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 1, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 2, 1], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_chn[2, 1, ], breaks = seq(min(sel_cof_chn[2, 1, ]), max(sel_cof_chn[2, 1, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient of sabino before the Middle Ages")
lines(density(sel_cof_chn[2, 1, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[2, 1], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 1, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 2, 1], col = 'blue', lty = 2, lwd = 2)

grd_num <- 1e+03
sel_cof_pdf <- kde2d(sel_cof_chn[1, 2, ], sel_cof_chn[2, 2, ], n = grd_num)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
      xlab = "Selection coefficient of tobiano", ylab = "Selection coefficient of sabino",
      main = "Posterior for selection coefficient from the Middle Ages")
abline(v = sel_cof_est[1, 2], col = 'black', lty = 2, lwd = 2)
abline(h = sel_cof_est[2, 2], col = 'black', lty = 2, lwd = 2)
# HPDregionplot(as.mcmc(sel_cof_chn), vars = 1:2, n = grd_num, prob = 0.95, col = "blue", lwd = 2, add = TRUE)

hist(sel_cof_chn[1, 2, ], breaks = seq(min(sel_cof_chn[1, 2, ]), max(sel_cof_chn[1, 2, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient of tobiano from the Middle Ages")
lines(density(sel_cof_chn[1, 2, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[1, 2], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 1, 2], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 2, 2], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_chn[2, 2, ], breaks = seq(min(sel_cof_chn[2, 2, ]), max(sel_cof_chn[2, 2, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient of sabino from the Middle Ages")
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

pdf(file = "./Output/Output v1.0/REAL v1.0/REAL_PTN_2_Posterior_SelChange.pdf", width = 24, height = 24)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
layout(matrix(c(1, 1, 3, 3, 2, 2, 4, 5), nrow = 4, ncol = 2))
grd_num <- 1e+03
sel_cof_pdf <- kde2d(sel_cof_chn[1, 1, ], sel_cof_chn[1, 2, ], n = grd_num)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
      xlab = "Selection coefficient before the Middle Ages", ylab = "Selection coefficient from the Middle Ages",
      main = "Posterior for selection coefficient of tobiano")
abline(v = sel_cof_est[1, 1], col = 'black', lty = 2, lwd = 2)
abline(h = sel_cof_est[1, 2], col = 'black', lty = 2, lwd = 2)
# HPDregionplot(as.mcmc(sel_cof_chn), vars = 1:2, n = grd_num, prob = 0.95, col = "blue", lwd = 2, add = TRUE)

grd_num <- 1e+03
sel_cof_pdf <- kde2d(sel_cof_chn[2, 1, ], sel_cof_chn[2, 2, ], n = grd_num)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
      xlab = "Selection coefficient before the Middle Ages", ylab = "Selection coefficient from the Middle Ages",
      main = "Posterior for selection coefficient of sabino")
abline(v = sel_cof_est[2, 1], col = 'black', lty = 2, lwd = 2)
abline(h = sel_cof_est[2, 2], col = 'black', lty = 2, lwd = 2)
# HPDregionplot(as.mcmc(sel_cof_chn), vars = 1:2, n = grd_num, prob = 0.95, col = "blue", lwd = 2, add = TRUE)

grd_num <- 1e+03
sel_cof_pdf <- kde2d(dif_sel_chn[1, ], dif_sel_chn[2, ], n = grd_num)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
      xlab = "Change in selection coefficient of tobiano", ylab = "Change in selection coefficient of sabino",
      main = "Posterior for change in selection coefficient before and from the Middle Ages")
abline(v = dif_sel_est[1], col = 'black', lty = 2, lwd = 2)
abline(h = dif_sel_est[2], col = 'black', lty = 2, lwd = 2)
# HPDregionplot(as.mcmc(t(dif_sel_chn)), vars = 1:2, n = grd_num, prob = 0.95, col = "blue", lwd = 2, add = TRUE)

hist(dif_sel_chn[1, ], breaks = seq(min(dif_sel_chn[1, ]), max(dif_sel_chn[1, ]), length.out = 50), freq = FALSE,
     xlab = "Change in selection coefficient",
     main = "Posterior for change in selection coefficient of tobiano")
lines(density(dif_sel_chn[1, ]), lwd = 2, col = 'black')
abline(v = dif_sel_est[1], col = 'black', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[1, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[1, 2], col = 'blue', lty = 2, lwd = 2)

hist(dif_sel_chn[2, ], breaks = seq(min(dif_sel_chn[2, ]), max(dif_sel_chn[2, ]), length.out = 50), freq = FALSE,
     xlab = "Change in selection coefficient",
     main = "Posterior for change in selection coefficient of sabino")
lines(density(dif_sel_chn[2, ]), lwd = 2, col = 'black')
abline(v = dif_sel_est[2], col = 'black', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[2, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[2, 2], col = 'blue', lty = 2, lwd = 2)
dev.off()

frq_pth_chn <- frq_pth_chn[, , brn_num:dim(frq_pth_chn)[3]]
frq_pth_chn <- frq_pth_chn[, , (1:round(dim(frq_pth_chn)[3] / thn_num)) * thn_num]

frq_pth_est <- matrix(NA, nrow = 4, ncol = dim(frq_pth_chn)[2])
frq_pth_est[1, ] <- rowMeans(frq_pth_chn[1, , ])
frq_pth_est[2, ] <- rowMeans(frq_pth_chn[2, , ])
frq_pth_est[3, ] <- rowMeans(frq_pth_chn[3, , ])
frq_pth_est[4, ] <- rowMeans(frq_pth_chn[4, , ])

frq_pth_hpd <- array(NA, dim = c(4, 2, dim(frq_pth_chn)[2]))
for (i in 1:dim(frq_pth_chn)[2]) {
  frq_pth_hpd[1, , i] <- HPDinterval(as.mcmc(frq_pth_chn[1, i, ]), prob = 0.95)
  frq_pth_hpd[2, , i] <- HPDinterval(as.mcmc(frq_pth_chn[2, i, ]), prob = 0.95)
  frq_pth_hpd[3, , i] <- HPDinterval(as.mcmc(frq_pth_chn[3, i, ]), prob = 0.95)
  frq_pth_hpd[4, , i] <- HPDinterval(as.mcmc(frq_pth_chn[4, i, ]), prob = 0.95)
}

pdf(file = "./Output/Output v1.0/REAL v1.0/REAL_PTN_2_Posterior_Traj.pdf", width = 24, height = 12)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
layout(matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2))
plot(0, type = 'n', xlim = c(min(raw_smp$age_mean), max(raw_smp$age_mean)), ylim = c(min(frq_pth_chn[1, , ]), max(frq_pth_chn[1, , ])),
     xlab = "Generation", ylab = "Haplotype frequency",
     main = "Posterior for underlying trajectory of the KM0sb1 haplotype")

for (i in 1:dim(frq_pth_chn)[2]) {
  lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_chn[1, , i], col = 'grey', lty = 1, lwd = 2)
}
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_est[1, ], col = 'black', lty = 2, lwd = 2)
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_hpd[1, 1, ], col = 'blue', lty = 2, lwd = 2)
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_hpd[1, 2, ], col = 'blue', lty = 2, lwd = 2)

plot(0, type = 'n', xlim = c(min(raw_smp$age_mean), max(raw_smp$age_mean)), ylim = c(min(frq_pth_chn[2, , ]), max(frq_pth_chn[2, , ])),
     xlab = "Generation", ylab = "Haplotype frequency",
     main = "Posterior for underlying trajectory of the KM0SB1 haplotype")

for (i in 1:dim(frq_pth_chn)[2]) {
  lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_chn[2, , i], col = 'grey', lty = 1, lwd = 2)
}
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_est[2, ], col = 'black', lty = 2, lwd = 2)
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_hpd[2, 1, ], col = 'blue', lty = 2, lwd = 2)
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_hpd[2, 2, ], col = 'blue', lty = 2, lwd = 2)

plot(0, type = 'n', xlim = c(min(raw_smp$age_mean), max(raw_smp$age_mean)), ylim = c(min(frq_pth_chn[3, , ]), max(frq_pth_chn[3, , ])),
     xlab = "Generation", ylab = "Haplotype frequency",
     main = "Posterior for underlying trajectory of the KM1sb1 haplotype")

for (i in 1:dim(frq_pth_chn)[2]) {
  lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_chn[3, , i], col = 'grey', lty = 1, lwd = 2)
}
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_est[3, ], col = 'black', lty = 2, lwd = 2)
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_hpd[3, 1, ], col = 'blue', lty = 2, lwd = 2)
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_hpd[3, 2, ], col = 'blue', lty = 2, lwd = 2)

plot(0, type = 'n', xlim = c(min(raw_smp$age_mean), max(raw_smp$age_mean)), ylim = c(min(frq_pth_chn[4, , ]), max(frq_pth_chn[4, , ])),
     xlab = "Generation", ylab = "Haplotype frequency",
     main = "Posterior for underlying trajectory of the KM1SB1 haplotype")

for (i in 1:dim(frq_pth_chn)[2]) {
  lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_chn[4, , i], col = 'grey', lty = 1, lwd = 2)
}
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_est[4, ], col = 'black', lty = 2, lwd = 2)
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_hpd[4, 1, ], col = 'blue', lty = 2, lwd = 2)
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_hpd[4, 2, ], col = 'blue', lty = 2, lwd = 2)
dev.off()

imp_smp_chn <- imp_smp_chn[, , brn_num:dim(imp_smp_chn)[3]]
imp_smp_chn <- imp_smp_chn[, , (1:round(dim(imp_smp_chn)[3] / thn_num)) * thn_num]

imp_smp_est <- matrix(NA, nrow = nrow(raw_smp), ncol = 9)
for (i in 1:nrow(imp_smp_est)) {
  imp_smp_est[i, ] <- rowSums(imp_smp_chn[, i, ]) / dim(imp_smp_chn)[3]
}
imp_smp_est <- as.data.frame(cbind(imp_smp_est[, 1] + imp_smp_est[, 2] + imp_smp_est[, 3],
                                   imp_smp_est[, 4] + imp_smp_est[, 5] + imp_smp_est[, 7],
                                   imp_smp_est[, 6] + imp_smp_est[, 8] + imp_smp_est[, 9],
                                   imp_smp_est[, 1] + imp_smp_est[, 4] + imp_smp_est[, 6],
                                   imp_smp_est[, 2] + imp_smp_est[, 5] + imp_smp_est[, 8],
                                   imp_smp_est[, 3] + imp_smp_est[, 7] + imp_smp_est[, 9]))
imp_smp_est

########################################

#' Raw data of Wutke et al. (2016) from 3500 BC (Domestication)
load("./Data/REAL.rda")

set.seed(88)
KIT13_smp <- KIT13
KIT13_smp <- KIT13_smp[which(rowSums(KIT13_smp[, 4:9]) != 0), ]
int_gen <- -round(max(KIT13_smp$age_mean, KIT13_smp$age_lower, KIT13_smp$age_upper) / 8)
lst_gen <- -round(min(KIT13_smp$age_mean, KIT13_smp$age_lower, KIT13_smp$age_upper) / 8)

KIT16_smp <- KIT16
KIT16_smp <- KIT16_smp[which(rowSums(KIT16_smp[, 4:9]) != 0), ]
int_gen <- -round(max(KIT16_smp$age_mean, KIT16_smp$age_lower, KIT16_smp$age_upper) / 8)
lst_gen <- -round(min(KIT16_smp$age_mean, KIT16_smp$age_lower, KIT16_smp$age_upper) / 8)

min_gen <- 3500 + 2000
KIT13_smp <- KIT13_smp[which(KIT13_smp$age_mean <= min_gen), ]
KIT16_smp <- KIT16_smp[which(KIT16_smp$age_mean <= min_gen), ]

KIT13_smp$age_mean <- -round(KIT13_smp$age_mean / 8)
KIT13_smp$age_lower <- -round(KIT13_smp$age_lower / 8)
KIT13_smp$age_upper <- -round(KIT13_smp$age_upper / 8)
# KIT13_smp[which(KIT13_smp[, 7] == 1), 4:5] <- 1 / 2
# KIT13_smp[which(KIT13_smp[, 8] == 1), 5:6] <- 1 / 2
# KIT13_smp[which(KIT13_smp[, 9] == 1), 4:6] <- 1 / 3
KIT13_smp[which(KIT13_smp[, 7] == 1), 4] <- 1 / 2
KIT13_smp[which(KIT13_smp[, 7] == 1), 5] <- 1 / 2
KIT13_smp[which(KIT13_smp[, 8] == 1), 5] <- 1 / 2
KIT13_smp[which(KIT13_smp[, 8] == 1), 6] <- 1 / 2
KIT13_smp[which(KIT13_smp[, 9] == 1), 4] <- 1 / 4
KIT13_smp[which(KIT13_smp[, 9] == 1), 5] <- 1 / 2
KIT13_smp[which(KIT13_smp[, 9] == 1), 6] <- 1 / 4
KIT13_smp <- KIT13_smp[, -(7:9)]

KIT16_smp$age_mean <- -round(KIT16_smp$age_mean / 8)
KIT16_smp$age_lower <- -round(KIT16_smp$age_lower / 8)
KIT16_smp$age_upper <- -round(KIT16_smp$age_upper / 8)
# KIT16_smp[which(KIT16_smp[, 7] == 1), 4:5] <- 1 / 2
# KIT16_smp[which(KIT16_smp[, 8] == 1), 5:6] <- 1 / 2
# KIT16_smp[which(KIT16_smp[, 9] == 1), 4:6] <- 1 / 3
KIT16_smp[which(KIT16_smp[, 7] == 1), 4] <- 1 / 2
KIT16_smp[which(KIT16_smp[, 7] == 1), 5] <- 1 / 2
KIT16_smp[which(KIT16_smp[, 8] == 1), 5] <- 1 / 2
KIT16_smp[which(KIT16_smp[, 8] == 1), 6] <- 1 / 2
KIT16_smp[which(KIT16_smp[, 9] == 1), 4] <- 1 / 4
KIT16_smp[which(KIT16_smp[, 9] == 1), 5] <- 1 / 2
KIT16_smp[which(KIT16_smp[, 9] == 1), 6] <- 1 / 4
KIT16_smp <- KIT16_smp[, -(7:9)]

raw_smp <- cbind(KIT13_smp[, 1:3], KIT13_smp[, 4:6], KIT16_smp[, 4:6])
raw_smp <- raw_smp[order(raw_smp$age_mean), ]
rownames(raw_smp) <- NULL

sel_cof <- matrix(c(0e+00, 0e+00, 0e+00, 0e+00, 0e+00, 0e+00), nrow = 3, ncol = 2)
rec_rat <- 1e-08 * 4688
pop_siz <- pop_siz[min(raw_smp$age_mean - int_gen + 1):max(raw_smp$age_mean - int_gen + 1)]
ref_siz <- tail(pop_siz, n = 1)
evt_gen <- round((400 - 2000) / 8) # 400 AD (the Middle Ages)
raw_smp
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04
stp_siz <- (1:itn_num)^(-2 / 3)
apt_rto <- 4e-01

system.time(PMMH <- cmprunAdaptPMMH(sel_cof, rec_rat, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, stp_siz, apt_rto))

save(sel_cof, rec_rat, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, stp_siz, apt_rto, PMMH,
     file = "./Output/Output v1.0/REAL v1.0/REAL_PTN_3.rda")

load("./Output/Output v1.0/REAL v1.0/REAL_PTN_3.rda")

sel_cof_chn <- PMMH$sel_cof_chn[1:2, , ]
frq_pth_chn <- PMMH$frq_pth_chn
imp_smp_chn <- PMMH$imp_smp_chn

pdf(file = "./Output/Output v1.0/REAL v1.0/REAL_PTN_3_Traceplot_SelCoeff.pdf", width = 24, height = 12)
par(mfrow = c(2, 2), mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.5, cex.sub = 1.25, cex.axis = 1.25, cex.lab = 1.25)
plot(1:itn_num, sel_cof_chn[1, 1, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient of tobiano before the Middle Ages")

plot(1:itn_num, sel_cof_chn[2, 1, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient of sabino before the Middle Ages")

plot(1:itn_num, sel_cof_chn[1, 2, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient of tobiano from the Middle Ages")

plot(1:itn_num, sel_cof_chn[2, 2, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient of sabino from the Middle Ages")
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

pdf(file = "./Output/Output v1.0/REAL v1.0/REAL_PTN_3_Posterior_SelCoeff.pdf", width = 24, height = 24)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
layout(matrix(c(1, 1, 4, 4, 2, 3, 5, 6), nrow = 4, ncol = 2))
grd_num <- 1e+03
sel_cof_pdf <- kde2d(sel_cof_chn[1, 1, ], sel_cof_chn[2, 1, ], n = grd_num)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
      xlab = "Selection coefficient of tobiano", ylab = "Selection coefficient of sabino",
      main = "Posterior for selection coefficient before the Middle Ages")
abline(v = sel_cof_est[1, 1], col = 'black', lty = 2, lwd = 2)
abline(h = sel_cof_est[2, 1], col = 'black', lty = 2, lwd = 2)
# HPDregionplot(as.mcmc(sel_cof_chn), vars = 1:2, n = grd_num, prob = 0.95, col = "blue", lwd = 2, add = TRUE)

hist(sel_cof_chn[1, 1, ], breaks = seq(min(sel_cof_chn[1, 1, ]), max(sel_cof_chn[1, 1, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient of tobiano before the Middle Ages")
lines(density(sel_cof_chn[1, 1, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[1, 1], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 1, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 2, 1], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_chn[2, 1, ], breaks = seq(min(sel_cof_chn[2, 1, ]), max(sel_cof_chn[2, 1, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient of sabino before the Middle Ages")
lines(density(sel_cof_chn[2, 1, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[2, 1], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 1, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 2, 1], col = 'blue', lty = 2, lwd = 2)

grd_num <- 1e+03
sel_cof_pdf <- kde2d(sel_cof_chn[1, 2, ], sel_cof_chn[2, 2, ], n = grd_num)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
      xlab = "Selection coefficient of tobiano", ylab = "Selection coefficient of sabino",
      main = "Posterior for selection coefficient from the Middle Ages")
abline(v = sel_cof_est[1, 2], col = 'black', lty = 2, lwd = 2)
abline(h = sel_cof_est[2, 2], col = 'black', lty = 2, lwd = 2)
# HPDregionplot(as.mcmc(sel_cof_chn), vars = 1:2, n = grd_num, prob = 0.95, col = "blue", lwd = 2, add = TRUE)

hist(sel_cof_chn[1, 2, ], breaks = seq(min(sel_cof_chn[1, 2, ]), max(sel_cof_chn[1, 2, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient of tobiano from the Middle Ages")
lines(density(sel_cof_chn[1, 2, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[1, 2], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 1, 2], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 2, 2], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_chn[2, 2, ], breaks = seq(min(sel_cof_chn[2, 2, ]), max(sel_cof_chn[2, 2, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient of sabino from the Middle Ages")
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

pdf(file = "./Output/Output v1.0/REAL v1.0/REAL_PTN_3_Posterior_SelChange.pdf", width = 24, height = 24)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
layout(matrix(c(1, 1, 3, 3, 2, 2, 4, 5), nrow = 4, ncol = 2))
grd_num <- 1e+03
sel_cof_pdf <- kde2d(sel_cof_chn[1, 1, ], sel_cof_chn[1, 2, ], n = grd_num)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
      xlab = "Selection coefficient before the Middle Ages", ylab = "Selection coefficient from the Middle Ages",
      main = "Posterior for selection coefficient of tobiano")
abline(v = sel_cof_est[1, 1], col = 'black', lty = 2, lwd = 2)
abline(h = sel_cof_est[1, 2], col = 'black', lty = 2, lwd = 2)
# HPDregionplot(as.mcmc(sel_cof_chn), vars = 1:2, n = grd_num, prob = 0.95, col = "blue", lwd = 2, add = TRUE)

grd_num <- 1e+03
sel_cof_pdf <- kde2d(sel_cof_chn[2, 1, ], sel_cof_chn[2, 2, ], n = grd_num)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
      xlab = "Selection coefficient before the Middle Ages", ylab = "Selection coefficient from the Middle Ages",
      main = "Posterior for selection coefficient of sabino")
abline(v = sel_cof_est[2, 1], col = 'black', lty = 2, lwd = 2)
abline(h = sel_cof_est[2, 2], col = 'black', lty = 2, lwd = 2)
# HPDregionplot(as.mcmc(sel_cof_chn), vars = 1:2, n = grd_num, prob = 0.95, col = "blue", lwd = 2, add = TRUE)

grd_num <- 1e+03
sel_cof_pdf <- kde2d(dif_sel_chn[1, ], dif_sel_chn[2, ], n = grd_num)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
      xlab = "Change in selection coefficient of tobiano", ylab = "Change in selection coefficient of sabino",
      main = "Posterior for change in selection coefficient before and from the Middle Ages")
abline(v = dif_sel_est[1], col = 'black', lty = 2, lwd = 2)
abline(h = dif_sel_est[2], col = 'black', lty = 2, lwd = 2)
# HPDregionplot(as.mcmc(t(dif_sel_chn)), vars = 1:2, n = grd_num, prob = 0.95, col = "blue", lwd = 2, add = TRUE)

hist(dif_sel_chn[1, ], breaks = seq(min(dif_sel_chn[1, ]), max(dif_sel_chn[1, ]), length.out = 50), freq = FALSE,
     xlab = "Change in selection coefficient",
     main = "Posterior for change in selection coefficient of tobiano")
lines(density(dif_sel_chn[1, ]), lwd = 2, col = 'black')
abline(v = dif_sel_est[1], col = 'black', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[1, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[1, 2], col = 'blue', lty = 2, lwd = 2)

hist(dif_sel_chn[2, ], breaks = seq(min(dif_sel_chn[2, ]), max(dif_sel_chn[2, ]), length.out = 50), freq = FALSE,
     xlab = "Change in selection coefficient",
     main = "Posterior for change in selection coefficient of sabino")
lines(density(dif_sel_chn[2, ]), lwd = 2, col = 'black')
abline(v = dif_sel_est[2], col = 'black', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[2, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[2, 2], col = 'blue', lty = 2, lwd = 2)
dev.off()

frq_pth_chn <- frq_pth_chn[, , brn_num:dim(frq_pth_chn)[3]]
frq_pth_chn <- frq_pth_chn[, , (1:round(dim(frq_pth_chn)[3] / thn_num)) * thn_num]

frq_pth_est <- matrix(NA, nrow = 4, ncol = dim(frq_pth_chn)[2])
frq_pth_est[1, ] <- rowMeans(frq_pth_chn[1, , ])
frq_pth_est[2, ] <- rowMeans(frq_pth_chn[2, , ])
frq_pth_est[3, ] <- rowMeans(frq_pth_chn[3, , ])
frq_pth_est[4, ] <- rowMeans(frq_pth_chn[4, , ])

frq_pth_hpd <- array(NA, dim = c(4, 2, dim(frq_pth_chn)[2]))
for (i in 1:dim(frq_pth_chn)[2]) {
   frq_pth_hpd[1, , i] <- HPDinterval(as.mcmc(frq_pth_chn[1, i, ]), prob = 0.95)
   frq_pth_hpd[2, , i] <- HPDinterval(as.mcmc(frq_pth_chn[2, i, ]), prob = 0.95)
   frq_pth_hpd[3, , i] <- HPDinterval(as.mcmc(frq_pth_chn[3, i, ]), prob = 0.95)
   frq_pth_hpd[4, , i] <- HPDinterval(as.mcmc(frq_pth_chn[4, i, ]), prob = 0.95)
}

pdf(file = "./Output/Output v1.0/REAL v1.0/REAL_PTN_3_Posterior_Traj.pdf", width = 24, height = 12)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
layout(matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2))
plot(0, type = 'n', xlim = c(min(raw_smp$age_mean), max(raw_smp$age_mean)), ylim = c(min(frq_pth_chn[1, , ]), max(frq_pth_chn[1, , ])),
     xlab = "Generation", ylab = "Haplotype frequency",
     main = "Posterior for underlying trajectory of the KM0sb1 haplotype")

for (i in 1:dim(frq_pth_chn)[2]) {
  lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_chn[1, , i], col = 'grey', lty = 1, lwd = 2)
}
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_est[1, ], col = 'black', lty = 2, lwd = 2)
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_hpd[1, 1, ], col = 'blue', lty = 2, lwd = 2)
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_hpd[1, 2, ], col = 'blue', lty = 2, lwd = 2)

plot(0, type = 'n', xlim = c(min(raw_smp$age_mean), max(raw_smp$age_mean)), ylim = c(min(frq_pth_chn[2, , ]), max(frq_pth_chn[2, , ])),
     xlab = "Generation", ylab = "Haplotype frequency",
     main = "Posterior for underlying trajectory of the KM0SB1 haplotype")

for (i in 1:dim(frq_pth_chn)[2]) {
  lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_chn[2, , i], col = 'grey', lty = 1, lwd = 2)
}
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_est[2, ], col = 'black', lty = 2, lwd = 2)
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_hpd[2, 1, ], col = 'blue', lty = 2, lwd = 2)
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_hpd[2, 2, ], col = 'blue', lty = 2, lwd = 2)

plot(0, type = 'n', xlim = c(min(raw_smp$age_mean), max(raw_smp$age_mean)), ylim = c(min(frq_pth_chn[3, , ]), max(frq_pth_chn[3, , ])),
     xlab = "Generation", ylab = "Haplotype frequency",
     main = "Posterior for underlying trajectory of the KM1sb1 haplotype")

for (i in 1:dim(frq_pth_chn)[2]) {
  lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_chn[3, , i], col = 'grey', lty = 1, lwd = 2)
}
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_est[3, ], col = 'black', lty = 2, lwd = 2)
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_hpd[3, 1, ], col = 'blue', lty = 2, lwd = 2)
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_hpd[3, 2, ], col = 'blue', lty = 2, lwd = 2)

plot(0, type = 'n', xlim = c(min(raw_smp$age_mean), max(raw_smp$age_mean)), ylim = c(min(frq_pth_chn[4, , ]), max(frq_pth_chn[4, , ])),
     xlab = "Generation", ylab = "Haplotype frequency",
     main = "Posterior for underlying trajectory of the KM1SB1 haplotype")

for (i in 1:dim(frq_pth_chn)[2]) {
  lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_chn[4, , i], col = 'grey', lty = 1, lwd = 2)
}
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_est[4, ], col = 'black', lty = 2, lwd = 2)
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_hpd[4, 1, ], col = 'blue', lty = 2, lwd = 2)
lines(min(raw_smp$age_mean):max(raw_smp$age_mean), frq_pth_hpd[4, 2, ], col = 'blue', lty = 2, lwd = 2)
dev.off()

imp_smp_chn <- imp_smp_chn[, , brn_num:dim(imp_smp_chn)[3]]
imp_smp_chn <- imp_smp_chn[, , (1:round(dim(imp_smp_chn)[3] / thn_num)) * thn_num]

imp_smp_est <- matrix(NA, nrow = nrow(raw_smp), ncol = 9)
for (i in 1:nrow(imp_smp_est)) {
  imp_smp_est[i, ] <- rowSums(imp_smp_chn[, i, ]) / dim(imp_smp_chn)[3]
}
imp_smp_est <- as.data.frame(cbind(imp_smp_est[, 1] + imp_smp_est[, 2] + imp_smp_est[, 3],
                                   imp_smp_est[, 4] + imp_smp_est[, 5] + imp_smp_est[, 7],
                                   imp_smp_est[, 6] + imp_smp_est[, 8] + imp_smp_est[, 9],
                                   imp_smp_est[, 1] + imp_smp_est[, 4] + imp_smp_est[, 6],
                                   imp_smp_est[, 2] + imp_smp_est[, 5] + imp_smp_est[, 8],
                                   imp_smp_est[, 3] + imp_smp_est[, 7] + imp_smp_est[, 9]))
imp_smp_est

################################################################################
