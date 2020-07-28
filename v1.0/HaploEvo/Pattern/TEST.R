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

test_seed <- 1
set.seed(test_seed)

model <- "WFM"
missing <- FALSE
sel_cof <- c(5e-03, 1e-03, -5e-03)
rec_rat <- 5e-05
pop_siz <- 5e+03
int_frq <- rmultinom(1, size = pop_siz, prob = rep(1, 4) / sum(rep(1, 4))) / pop_siz
smp_gen <- (0:10) * 50
smp_siz <- rep(50, 11)

sim_HMM_WFM <- cmpsimulateHMM(model, missing, sel_cof, rec_rat, pop_siz, int_frq, smp_gen, smp_siz)
smp_gen <- sim_HMM_WFM$smp_gen
smp_siz <- sim_HMM_WFM$smp_siz
smp_cnt <- sim_HMM_WFM$smp_cnt
mis_cnt <- sim_HMM_WFM$mis_cnt
smp_frq <- sim_HMM_WFM$smp_cnt / sim_HMM_WFM$smp_siz
pop_frq <- sim_HMM_WFM$pop_frq

save(sel_cof, rec_rat, pop_siz, smp_gen, smp_siz, smp_cnt, mis_cnt, smp_frq, pop_frq, 
     file = "./Output/Output v1.0/HaploEvo/Horse coat pattern/TEST_data.rda")

load("./Output/Output v1.0/HaploEvo/Horse coat pattern/TEST_data.rda")

pdf(file = "./Output/Output v1.0/HaploEvo/Horse coat pattern/TEST_data.pdf", width = 30, height = 15)
par(mfrow = c(3, 3), oma = c(0, 0, 3, 0), mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
k <- min(smp_gen):max(smp_gen)
plot(k, pop_frq[1, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_frq[1, ], pop_frq[1, ]), max(smp_frq[1, ], pop_frq[1, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "Genotype KM0KM0/sb1sb1")
points(smp_gen, smp_frq[1, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_frq[2, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_frq[2, ], pop_frq[2, ]), max(smp_frq[2, ], pop_frq[2, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "Genotype KM0KM0/sb1SB1")
points(smp_gen, smp_frq[2, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_frq[3, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_frq[3, ], pop_frq[3, ]), max(smp_frq[3, ], pop_frq[3, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "Genotype KM0KM0/SB1SB1")
points(smp_gen, smp_frq[3, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_frq[4, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_frq[4, ], pop_frq[4, ]), max(smp_frq[4, ], pop_frq[4, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "Genotype KM0KM1/sb1sb1")
points(smp_gen, smp_frq[4, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_frq[5, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_frq[5, ], pop_frq[5, ]), max(smp_frq[5, ], pop_frq[5, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "Genotype KM0KM1/sb1SB1")
points(smp_gen, smp_frq[5, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_frq[6, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_frq[6, ], pop_frq[6, ]), max(smp_frq[6, ], pop_frq[6, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "Genotype KM0KM1/SB1SB1")
points(smp_gen, smp_frq[6, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_frq[7, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_frq[7, ], pop_frq[7, ]), max(smp_frq[7, ], pop_frq[7, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "Genotype KM1KM1/sb1sb1")
points(smp_gen, smp_frq[7, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_frq[8, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_frq[8, ], pop_frq[8, ]), max(smp_frq[8, ], pop_frq[8, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "Genotype KM1KM1/sb1SB1")
points(smp_gen, smp_frq[8, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_frq[9, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_frq[9, ], pop_frq[9, ]), max(smp_frq[9, ], pop_frq[9, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "Genotype KM1KM1/SB1SB1")
points(smp_gen, smp_frq[9, ], col = 'red', pch = 17, cex = 1)
title("A simulated dataset generated with the Wright-Fisher model", outer = TRUE)
dev.off()

########################################

load("./Output/Output v1.0/HaploEvo/Horse coat pattern/TEST_data.rda")

int_sel_cof <- c(0e-00, 0e-00, 0e-00)
rec_rat
pop_siz
smp_gen
smp_siz
smp_cnt
mis_cnt
ptn_num <- 5e+00
pcl_num <- 1e+03
gap_num <- 1e+02

system.time(OptNum <- calculateOptimalParticleNum(int_sel_cof, rec_rat, pop_siz, smp_gen, smp_siz, smp_cnt, mis_cnt, ptn_num, pcl_num, gap_num))

save(sel_cof, rec_rat, pop_siz, smp_gen, smp_siz, smp_cnt, mis_cnt, ptn_num, pcl_num, gap_num, OptNum,
     file = "./Output/Output v1.0/HaploEvo/Horse coat pattern/TEST_OptNum.rda")

load("./Output/Output v1.0/HaploEvo/Horse coat pattern/TEST_OptNum.rda")

opt_pcl_num <- OptNum$opt_pcl_num
log_lik_sdv <- OptNum$log_lik_sdv

pdf(file = "./Output/Output v1.0/HaploEvo/Horse coat pattern/TEST_OptNum.pdf", width = 10, height = 10)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
plot(opt_pcl_num, log_lik_sdv, type = 'b', lwd = 2, 
     xlab = "Particle number", ylab = "Log-likelihood standard deviation", main = "Optimal particle number in the PMMH")
abline(h = 1.7, col = 'red', lty = 2, lwd = 2)
abline(h = 1.0, col = 'red', lty = 2, lwd = 2)
dev.off()

########################################

load("./Output/Output v1.0/HaploEvo/Horse coat pattern/TEST_data.rda")

int_sel_cof <- c(0e-00, 0e-00, 0e-00)
rec_rat
pop_siz
smp_gen
smp_siz
smp_cnt
mis_cnt
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04
ER <- TRUE
PA <- FALSE
nap_num <- itn_num * 0.1

system.time(PMMH <- cmprunPMMH(int_sel_cof, rec_rat, pop_siz, smp_gen, smp_siz, smp_cnt, mis_cnt, ptn_num, pcl_num, itn_num, ER, PA, nap_num))

save(sel_cof, rec_rat, pop_siz, smp_gen, smp_siz, smp_cnt, mis_cnt, smp_frq, pop_frq, ptn_num, pcl_num, itn_num, ER, PA, nap_num, PMMH,
     file = "./Output/Output v1.0/HaploEvo/Horse coat pattern/TEST_PMMH.rda")

load("./Output/Output v1.0/HaploEvo/Horse coat pattern/TEST_PMMH.rda")

sel_cof_to_chn <- PMMH$sel_cof_to_chn
sel_cof_sb_chn <- PMMH$sel_cof_sb_chn
sel_cof_ts_chn <- PMMH$sel_cof_ts_chn

pdf(file = "./Output/Output v1.0/HaploEvo/Horse coat pattern/TEST_PMMH_traceplot.pdf", width = 10, height = 15)
par(mfrow = c(3, 1), mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
plot(1:itn_num, sel_cof_to_chn[1:itn_num], type = 'l', 
     xlab = "Iteration", ylab = "Selection coefficient", main = "Trace plot of selection coefficient of tobiano coat")
abline(h = sel_cof[1], col = 'red', lty = 2, lwd = 2)
plot(1:itn_num, sel_cof_sb_chn[1:itn_num], type = 'l', 
     xlab = "Iteration", ylab = "Selection coefficient", main = "Trace plot of selection coefficient of sabino coat")
abline(h = sel_cof[2], col = 'red', lty = 2, lwd = 2)
plot(1:itn_num, sel_cof_ts_chn[1:itn_num], type = 'l', 
     xlab = "Iteration", ylab = "Selection coefficient", main = "Trace plot of selection coefficient of mixed coat")
abline(h = sel_cof[3], col = 'red', lty = 2, lwd = 2)
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

pdf(file = "./Output/Output v1.0/HaploEvo/Horse coat pattern/TEST_PMMH_posterior.pdf", width = 20, height = 20)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
layout(matrix(c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 5, 6), nrow = 6, ncol = 2))
grd_num <- 1e+03
sel_cof_pdf <- kde2d(sel_cof_to_chn, sel_cof_sb_chn, n = grd_num)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
      xlab = "Selection coefficient of tobiano coat", ylab = "Selection coefficient of sabino coat",
      main = "Posterior for selection coefficients of tobiano and sabino coat")
abline(v = sel_cof[1], col = 'red', lty = 2, lwd = 2)
abline(h = sel_cof[2], col = 'red', lty = 2, lwd = 2)
abline(v = sel_cof_to_mmse, col = 'black', lty = 2, lwd = 2)
abline(h = sel_cof_sb_mmse, col = 'black', lty = 2, lwd = 2)

grd_num <- 1e+03
sel_cof_pdf <- kde2d(sel_cof_to_chn, sel_cof_ts_chn, n = grd_num)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
      xlab = "Selection coefficient of tobiano coat", ylab = "Selection coefficient of mixed coat",
      main = "Posterior for selection coefficients of tobiano and mixed coat")
abline(v = sel_cof[1], col = 'red', lty = 2, lwd = 2)
abline(h = sel_cof[3], col = 'red', lty = 2, lwd = 2)
abline(v = sel_cof_to_mmse, col = 'black', lty = 2, lwd = 2)
abline(h = sel_cof_ts_mmse, col = 'black', lty = 2, lwd = 2)

grd_num <- 1e+03
sel_cof_pdf <- kde2d(sel_cof_ts_chn, sel_cof_sb_chn, n = grd_num)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
      xlab = "Selection coefficient of sabino coat", ylab = "Selection coefficient of sabino coat",
      main = "Posterior for selection coefficients of sabino and mixed coat")
abline(v = sel_cof[3], col = 'red', lty = 2, lwd = 2)
abline(h = sel_cof[2], col = 'red', lty = 2, lwd = 2)
abline(v = sel_cof_ts_mmse, col = 'black', lty = 2, lwd = 2)
abline(h = sel_cof_sb_mmse, col = 'black', lty = 2, lwd = 2)

hist(sel_cof_to_chn, breaks = seq(min(sel_cof_to_chn), max(sel_cof_to_chn), length.out = 30), freq = FALSE, 
     xlab = "Selection coefficient", main = "Marginal posterior for selection coefficient of tobiano coat")
lines(density(sel_cof_to_chn), lwd = 2, col = 'black')
abline(v = sel_cof[1], col = 'red', lty = 2, lwd = 2)
abline(v = sel_cof_to_mmse, col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_to_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_to_hpd[2], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_sb_chn, breaks = seq(min(sel_cof_sb_chn), max(sel_cof_sb_chn), length.out = 30), freq = FALSE, 
     xlab = "Selection coefficient", main = "Marginal posterior for selection coefficient of sabino coat")
lines(density(sel_cof_sb_chn), lwd = 2, col = 'black')
abline(v = sel_cof[2], col = 'red', lty = 2, lwd = 2)
abline(v = sel_cof_sb_mmse, col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_sb_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_sb_hpd[2], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_ts_chn, breaks = seq(min(sel_cof_ts_chn), max(sel_cof_ts_chn), length.out = 30), freq = FALSE, 
     xlab = "Selection coefficient", main = "Marginal posterior for selection coefficient of mixed coat")
lines(density(sel_cof_ts_chn), lwd = 2, col = 'black')
abline(v = sel_cof[3], col = 'red', lty = 2, lwd = 2)
abline(v = sel_cof_ts_mmse, col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_ts_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_ts_hpd[2], col = 'blue', lty = 2, lwd = 2)
dev.off()
