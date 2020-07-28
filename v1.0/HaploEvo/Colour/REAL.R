#' @title Inferring fluctuating selection acting on horse coat colours and patterns from ancient DNA data
#' @author Zhangyi He, Xiaoyang Dai, Mark Beaumont and Feng Yu

#' Haplotype evolution (Wright-Fisher diffusion)

#' Horse base coat colour: ASIP and MC1R (fluctuating selection caused by domestication)

# set the directory
setwd("~/Dropbox/Jeffery He/iResearch/Publications/2017/HE2019-2L-WFD-PMMH-Horse-MBE")

source("./Code/Code v1.0/HaploEvo/Horse coat colour/HE2017_rfun.R")

#install.packages("RColorBrewer")
library("RColorBrewer")

#install.packages("ggplot2")
library("ggplot2")

#install.packages("plot3D")
library("plot3D")

################################################################################

load("./Output/Output v1.0/HaploEvo/Horse coat colour/REAL_data.rda")

int_sel_cof <- c(0e+00, 0e+00, 0e+00, 0e+00)
rec_rat <- 5e-01
pop_siz <- 1.6e+04
dom_gen <- -5500 / 8
smp_gen <- -smp_ybp / 8
smp_gen <- round(smp_gen[-(1:6)])
smp_siz <- colSums(smp_gen_cnt) + colSums(mis_gen_cnt)
smp_siz <- smp_siz[-(1:6)]
smp_cnt <- smp_gen_cnt[, -(1:6)]
mis_cnt <- mis_gen_cnt[, -(1:6)]
ptn_num <- 5e+00
pcl_num <- 2e+03
# gap_num <- 5e+02
# 
# system.time(OptNum <- calculateOptimalParticleNum(int_sel_cof, rec_rat, pop_siz, dom_gen, smp_gen, smp_siz, smp_cnt, mis_cnt, ptn_num, pcl_num, gap_num))
# save(int_sel_cof, rec_rat, pop_siz, dom_gen, smp_gen, smp_siz, smp_cnt, mis_cnt, ptn_num, pcl_num, gap_num, OptNum,
#      file = "./Output/Output v1.0/HaploEvo/Horse coat colour/REAL_OptNum.rda")
# 
# load("./Output/Output v1.0/HaploEvo/Horse coat colour/REAL_OptNum.rda")
# 
# opt_pcl_num <- OptNum$opt_pcl_num
# log_lik_sdv <- OptNum$log_lik_sdv
# 
# pdf(file = "./Output/Output v1.0/HaploEvo/Horse coat colour/REAL_OptNum.pdf", width = 10, height = 10)
# par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
# plot(opt_pcl_num, log_lik_sdv, type = 'b', lwd = 2, 
#      xlab = "Particle number", ylab = "Log-likelihood standard deviation", main = "Optimal particle number in the PMMH")
# abline(h = 1.7, col = 'red', lty = 2, lwd = 2)
# abline(h = 1.0, col = 'red', lty = 2, lwd = 2)
# dev.off()
# 
# 
# pcl_num <- head(opt_pcl_num[which(log_lik_sdv > 1.0 & log_lik_sdv < 1.7)], n = 1)
itn_num <- 5e+04
ER <- TRUE
PA <- TRUE
nap_num <- itn_num * 0.1

system.time(PMMH <- cmprunPMMH(int_sel_cof, rec_rat, pop_siz, dom_gen, smp_gen, smp_siz, smp_cnt, mis_cnt, ptn_num, pcl_num, itn_num, ER, PA, nap_num))

save(rec_rat, pop_siz, dom_gen, smp_gen, smp_siz, smp_cnt, mis_cnt, ptn_num, pcl_num, itn_num, ER, PA, nap_num, PMMH,
     file = "./Output/Output v1.0/HaploEvo/Horse coat colour/REAL_PMMH.rda")

load("./Output/Output v1.0/HaploEvo/Horse coat colour/REAL_PMMH.rda")

sel_cof_b_chn <- PMMH$sel_cof_b_chn
sel_cof_c_chn <- PMMH$sel_cof_c_chn

pdf(file = "./Output/Output v1.0/HaploEvo/Horse coat colour/REAL_PMMH_traceplot.pdf", width = 20, height = 10)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
plot(1:itn_num, sel_cof_b_chn[1, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient", main = "Trace plot of s of black coat before domestication")
plot(1:itn_num, sel_cof_b_chn[2, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient", main = "Trace plot of s of black coat after domestication")
plot(1:itn_num, sel_cof_c_chn[1, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient", main = "Trace plot of s of chestnut coat before domestication")
plot(1:itn_num, sel_cof_c_chn[2, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient", main = "Trace plot of s of chestnut coat after domestication")
dev.off()

brn_num <- itn_num * 0.2
sel_cof_b_chn <- sel_cof_b_chn[, brn_num:ncol(sel_cof_b_chn)]
sel_cof_c_chn <- sel_cof_c_chn[, brn_num:ncol(sel_cof_c_chn)]

thn_num <- 4e+00
sel_cof_b_chn <- sel_cof_b_chn[, (1:round(ncol(sel_cof_b_chn) / thn_num)) * thn_num]
sel_cof_c_chn <- sel_cof_c_chn[, (1:round(ncol(sel_cof_c_chn) / thn_num)) * thn_num]

sel_cof_b_mmse <- c(mean(sel_cof_b_chn[1, ]), mean(sel_cof_b_chn[2, ]))
sel_cof_c_mmse <- c(mean(sel_cof_c_chn[1, ]), mean(sel_cof_c_chn[2, ]))
    
sel_cof_b_hpd <- rbind(HPDinterval(as.mcmc(sel_cof_b_chn[1, ]), prob = 0.95), HPDinterval(as.mcmc(sel_cof_b_chn[2, ]), prob = 0.95))
sel_cof_c_hpd <- rbind(HPDinterval(as.mcmc(sel_cof_c_chn[1, ]), prob = 0.95), HPDinterval(as.mcmc(sel_cof_c_chn[2, ]), prob = 0.95))

pdf(file = "./Output/Output v1.0/HaploEvo/Horse coat colour/REAL_PMMH_posterior.pdf", width = 30, height = 20)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
layout(matrix(c(1, 1, 2, 2, 3, 3, 4, 4, 5, 6, 7, 8), nrow = 4, ncol = 3))
grd_num <- 1e+03
sel_cof_pdf <- kde2d(sel_cof_b_chn[1, ], sel_cof_c_chn[1, ], n = grd_num)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
      xlab = "Selection coefficient of black coat", ylab = "Selection coefficient of chestnut coat",
      main = "Posterior for s of black and chestnut coat before domestication")
abline(v = sel_cof_b_mmse[1], col = 'black', lty = 2, lwd = 2)
abline(h = sel_cof_c_mmse[1], col = 'black', lty = 2, lwd = 2)

grd_num <- 1e+03
sel_cof_pdf <- kde2d(sel_cof_b_chn[2, ], sel_cof_c_chn[2, ], n = grd_num)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
      xlab = "Selection coefficient of black coat", ylab = "Selection coefficient of chestnut coat",
      main = "Posterior for s of black and chestnut coat after domestication")
abline(v = sel_cof_b_mmse[2], col = 'black', lty = 2, lwd = 2)
abline(h = sel_cof_c_mmse[2], col = 'black', lty = 2, lwd = 2)

grd_num <- 1e+03
sel_cof_pdf <- kde2d(sel_cof_b_chn[1, ], sel_cof_b_chn[2, ], n = grd_num)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
      xlab = "Selection coefficient before domestication", ylab = "Selection coefficient after domestication",
      main = "Posterior for s of black coat before and after domestication")
abline(v = sel_cof_b_mmse[1], col = 'black', lty = 2, lwd = 2)
abline(h = sel_cof_b_mmse[2], col = 'black', lty = 2, lwd = 2)

grd_num <- 1e+03
sel_cof_pdf <- kde2d(sel_cof_c_chn[1, ], sel_cof_c_chn[2, ], n = grd_num)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
      xlab = "Selection coefficient before domestication", ylab = "Selection coefficient after domestication",
      main = "Posterior for s of chestnut coat before and after domestication")
abline(v = sel_cof_c_mmse[1], col = 'black', lty = 2, lwd = 2)
abline(h = sel_cof_c_mmse[2], col = 'black', lty = 2, lwd = 2)

hist(sel_cof_b_chn[1, ], breaks = seq(min(sel_cof_b_chn[1, ]), max(sel_cof_b_chn[1, ]), length.out = 30), freq = FALSE,
     xlab = "Selection coefficient", main = "Marginal posterior for s of black coat before domestication")
lines(density(sel_cof_b_chn[1, ]), lwd = 2, col = 'black')
abline(v = sel_cof_b_mmse[1], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_b_hpd[1, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_b_hpd[1, 2], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_b_chn[2, ], breaks = seq(min(sel_cof_b_chn[2, ]), max(sel_cof_b_chn[2, ]), length.out = 30), freq = FALSE,
     xlab = "Selection coefficient", main = "Marginal posterior for s of black coat after domestication")
lines(density(sel_cof_b_chn[2, ]), lwd = 2, col = 'black')
abline(v = sel_cof_b_mmse[2], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_b_hpd[2, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_b_hpd[2, 2], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_c_chn[1, ], breaks = seq(min(sel_cof_c_chn[1, ]), max(sel_cof_c_chn[1, ]), length.out = 30), freq = FALSE,
     xlab = "Selection coefficient", main = "Marginal posterior for s of chestnut coat before domestication")
lines(density(sel_cof_c_chn[1, ]), lwd = 2, col = 'black')
abline(v = sel_cof_c_mmse[1], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_c_hpd[1, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_c_hpd[1, 2], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_c_chn[2, ], breaks = seq(min(sel_cof_c_chn[2, ]), max(sel_cof_c_chn[2, ]), length.out = 30), freq = FALSE,
     xlab = "Selection coefficient", main = "Marginal posterior for s of chestnut coat after domestication")
lines(density(sel_cof_c_chn[2, ]), lwd = 2, col = 'black')
abline(v = sel_cof_c_mmse[2], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_c_hpd[2, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_c_hpd[2, 2], col = 'blue', lty = 2, lwd = 2)
dev.off()
