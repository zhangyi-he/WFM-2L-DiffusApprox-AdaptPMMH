#' @title Inferring fluctuating selection acting on horse coat colours and patterns from ancient DNA data
#' @author Zhangyi He, Xiaoyang Dai, Mark Beaumont and Feng Yu

#' Haplotype evolution (Wright-Fisher diffusion)

#' Horse white spotting pattern: KIT13 and KIT16 (constant selection)

# set the directory
setwd("~/Dropbox/Jeffery He/iResearch/Publications/2017/HE2019-2L-WFD-PMMH-Horse-MBE")

source("./Code/Code v1.0/HaploEvo/Horse coat pattern/HE2017_rfun.R")

#install.packages("RColorBrewer")
library("RColorBrewer")

#install.packages("ggplot2")
library("ggplot2")

#install.packages("plot3D")
library("plot3D")

################################################################################

load("./Output/Output v1.0/HaploEvo/Horse coat pattern/REAL_data.rda")

int_sel_cof <- c(0e-00, 0e-00, 0e-00)
rec_rat <- 5e-05
pop_siz <- 1.6e+04
smp_gen <- -smp_ybp / 8
smp_gen <- round(smp_gen[-(1:21)])
smp_siz <- colSums(smp_gen_cnt) + colSums(mis_gen_cnt)
smp_siz <- smp_siz[-(1:21)]
smp_cnt <- smp_gen_cnt[, -(1:21)]
mis_cnt <- mis_gen_cnt[, -(1:21)]
ptn_num <- 5e+00
pcl_num <- 2e+03
# gap_num <- 5e+02
# 
# system.time(OptNum <- calculateOptimalParticleNum(int_sel_cof, rec_rat, pop_siz, smp_gen, smp_siz, smp_cnt, mis_cnt, ptn_num, pcl_num, gap_num))
# save(int_sel_cof, rec_rat, pop_siz, smp_gen, smp_siz, smp_cnt, mis_cnt, ptn_num, pcl_num, gap_num, OptNum,
#      file = "./Output/Output v1.0/HaploEvo/Horse coat pattern/REAL_OptNum.rda")
# 
# load("./Output/Output v1.0/HaploEvo/Horse coat pattern/REAL_OptNum.rda")
# 
# opt_pcl_num <- OptNum$opt_pcl_num
# log_lik_sdv <- OptNum$log_lik_sdv
# 
# pdf(file = "./Output/Output v1.0/HaploEvo/Horse coat pattern/REAL_OptNum.pdf", width = 10, height = 10)
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

system.time(PMMH <- cmprunPMMH(int_sel_cof, rec_rat, pop_siz, smp_gen, smp_siz, smp_cnt, mis_cnt, ptn_num, pcl_num, itn_num, ER, PA, nap_num))

save(rec_rat, pop_siz, smp_gen, smp_siz, smp_cnt, mis_cnt, ptn_num, pcl_num, itn_num, ER, PA, nap_num, PMMH,
     file = "./Output/Output v1.0/HaploEvo/Horse coat pattern/REAL_PMMH.rda")

load("./Output/Output v1.0/HaploEvo/Horse coat pattern/REAL_PMMH.rda")

sel_cof_to_chn <- PMMH$sel_cof_to_chn
sel_cof_sb_chn <- PMMH$sel_cof_sb_chn
sel_cof_ts_chn <- PMMH$sel_cof_ts_chn

pdf(file = "./Output/Output v1.0/HaploEvo/Horse coat pattern/REAL_PMMH_traceplot.pdf", width = 10, height = 15)
par(mfrow = c(3, 1), mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
plot(1:itn_num, sel_cof_to_chn[1:itn_num], type = 'l', 
     xlab = "Iteration", ylab = "Selection coefficient", main = "Trace plot of selection coefficient of tobiano coat")
plot(1:itn_num, sel_cof_sb_chn[1:itn_num], type = 'l', 
     xlab = "Iteration", ylab = "Selection coefficient", main = "Trace plot of selection coefficient of sabino coat")
plot(1:itn_num, sel_cof_ts_chn[1:itn_num], type = 'l', 
     xlab = "Iteration", ylab = "Selection coefficient", main = "Trace plot of selection coefficient of mixed coat")
dev.off()

brn_num <- itn_num * 0.2
sel_cof_to_chn <- sel_cof_to_chn[brn_num:length(sel_cof_to_chn)]
sel_cof_sb_chn <- sel_cof_sb_chn[brn_num:length(sel_cof_sb_chn)]
sel_cof_ts_chn <- sel_cof_ts_chn[brn_num:length(sel_cof_ts_chn)]

thn_num <- 4e+00
sel_cof_to_chn <- sel_cof_to_chn[(1:round(length(sel_cof_to_chn) / thn_num)) * thn_num]
sel_cof_sb_chn <- sel_cof_sb_chn[(1:round(length(sel_cof_sb_chn) / thn_num)) * thn_num]
sel_cof_ts_chn <- sel_cof_ts_chn[(1:round(length(sel_cof_ts_chn) / thn_num)) * thn_num]

sel_cof_to_mmse <- mean(sel_cof_to_chn)
sel_cof_sb_mmse <- mean(sel_cof_sb_chn)
sel_cof_ts_mmse <- mean(sel_cof_ts_chn)
    
sel_cof_to_hpd <- HPDinterval(as.mcmc(sel_cof_to_chn), prob = 0.95)
sel_cof_sb_hpd <- HPDinterval(as.mcmc(sel_cof_sb_chn), prob = 0.95)
sel_cof_ts_hpd <- HPDinterval(as.mcmc(sel_cof_ts_chn), prob = 0.95)

pdf(file = "./Output/Output v1.0/HaploEvo/Horse coat pattern/REAL_PMMH_posterior.pdf", width = 20, height = 20)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
layout(matrix(c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 5, 6), nrow = 6, ncol = 2))
grd_num <- 1e+03
sel_cof_pdf <- kde2d(sel_cof_to_chn, sel_cof_sb_chn, n = grd_num)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
      xlab = "Selection coefficient of tobiano coat", ylab = "Selection coefficient of sabino coat",
      main = "Posterior for selection coefficients of tobiano and sabino coat")
abline(v = sel_cof_to_mmse, col = 'black', lty = 2, lwd = 2)
abline(h = sel_cof_sb_mmse, col = 'black', lty = 2, lwd = 2)

grd_num <- 1e+03
sel_cof_pdf <- kde2d(sel_cof_to_chn, sel_cof_ts_chn, n = grd_num)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
      xlab = "Selection coefficient of tobiano coat", ylab = "Selection coefficient of mixed coat",
      main = "Posterior for selection coefficients of tobiano and mixed coat")
abline(v = sel_cof_to_mmse, col = 'black', lty = 2, lwd = 2)
abline(h = sel_cof_ts_mmse, col = 'black', lty = 2, lwd = 2)

grd_num <- 1e+03
sel_cof_pdf <- kde2d(sel_cof_ts_chn, sel_cof_sb_chn, n = grd_num)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
      xlab = "Selection coefficient of sabino coat", ylab = "Selection coefficient of sabino coat",
      main = "Posterior for selection coefficients of sabino and mixed coat")
abline(v = sel_cof_ts_mmse, col = 'black', lty = 2, lwd = 2)
abline(h = sel_cof_sb_mmse, col = 'black', lty = 2, lwd = 2)

hist(sel_cof_to_chn, breaks = seq(min(sel_cof_to_chn), max(sel_cof_to_chn), length.out = 30), freq = FALSE, 
     xlab = "Selection coefficient", main = "Marginal posterior for selection coefficient of tobiano coat")
lines(density(sel_cof_to_chn), lwd = 2, col = 'black')
abline(v = sel_cof_to_mmse, col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_to_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_to_hpd[2], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_sb_chn, breaks = seq(min(sel_cof_sb_chn), max(sel_cof_sb_chn), length.out = 30), freq = FALSE, 
     xlab = "Selection coefficient", main = "Marginal posterior for selection coefficient of sabino coat")
lines(density(sel_cof_sb_chn), lwd = 2, col = 'black')
abline(v = sel_cof_sb_mmse, col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_sb_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_sb_hpd[2], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_ts_chn, breaks = seq(min(sel_cof_ts_chn), max(sel_cof_ts_chn), length.out = 30), freq = FALSE, 
     xlab = "Selection coefficient", main = "Marginal posterior for selection coefficient of mixed coat")
lines(density(sel_cof_ts_chn), lwd = 2, col = 'black')
abline(v = sel_cof_ts_mmse, col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_ts_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_ts_hpd[2], col = 'blue', lty = 2, lwd = 2)
dev.off()
