#' @title A Bayesian approach for estimating selection coefficients and testing their changes from ancient DNA data
#' @author Xiaoyang Dai, Mark Beaumont, Feng Yu, Ludovic Orlando, Zhangyi He

#' version 1.2
#' Horse coat patterns (KIT13 & KIT116) under non-constant natural selection and non-constant demographic histories (N/A is not allowed)

# set the directory
setwd("~/Dropbox/Jeffery He/iResearch/Publications/2018/HE2021-WFM-2L-DiffusApprox-PMMH-MolEcol")

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
source("./Code/Code v1.0/Code v1.2/RFUN_PTN.R")

################################################################################

#' Simulate the haplotype frequency trajectories according to the two-locus Wright-Fisher model with selection
#' Parameter setting
#' @param sel_cof the selection coefficients of the tobiano, sabino and mixed phenotypes
#' @param rec_rat the recombination rate between the KIT13 and KIT16 loci
#' @param pop_siz the size of the horse population (non-constant)
#' @param int_frq the initial haplotype frequencies of the population
#' @param evt_gen the generation that the event of interest occurred
#' @param int_gen the generation that the simulated haplotype frequency trajectories started
#' @param lst_gen the generation that the simulated haplotype frequency trajectories ended

# sel_cof <- matrix(c(1e-02, 5e-03, 1e-03, 5e-03, 1e-02, 1e-03), nrow = 3, ncol = 2)
# rec_rat <- 5e-05
# pop_siz <- c(rep(1e+04, length.out = 201), rep(5e+03, length.out = 200), rep(1e+04, length.out = 100))
# int_frq <- c(5e-01, 2e-01, 2e-01, 1e-01)
# evt_gen <- 240
# int_gen <- 0
# lst_gen <- 500
#
# frq_pth <- cmpsimulateWFM(sel_cof, rec_rat, pop_siz, int_frq, evt_gen, int_gen, lst_gen)$hap_frq_pth
#
# k <- int_gen:lst_gen
# plot(k, frq_pth[1, ], type = "l", lwd = 1.5,
#      xlab = "Generation", ylab = "Haplotype frequency",
#      main = "WFM: the KM0sb1 haplotype frequency trajectory")
# plot(k, frq_pth[2, ], type = "l", lwd = 1.5,
#      xlab = "Generation", ylab = "Haplotype frequency",
#      main = "WFM: the KM0SB1 haplotype frequency trajectory")
# plot(k, frq_pth[3, ], type = "l", lwd = 1.5,
#      xlab = "Generation", ylab = "Haplotype frequency",
#      main = "WFM: the KM1sb1 haplotype frequency trajectory")
# plot(k, frq_pth[4, ], type = "l", lwd = 1.5,
#      xlab = "Generation", ylab = "Haplotype frequency",
#      main = "WFM: the KM1SB1 haplotype frequency trajectory")

########################################

#' Simulate the haplotype frequency trajectories according to the two-locus Wright-Fisher diffusion with selection using the Euler-Maruyama method
#' Parameter setting
#' @param sel_cof the selection coefficients of the tobiano, sabino and mixed phenotypes
#' @param rec_rat the recombination rate between the KIT13 and KIT16 loci
#' @param pop_siz the size of the horse population (non-constant)
#' @param ref_siz the reference size of the horse population
#' @param int_frq the initial haplotype frequencies of the population
#' @param evt_gen the generation that the event of interest occurred
#' @param int_gen the generation that the simulated haplotype frequency trajectories started
#' @param lst_gen the generation that the simulated haplotype frequency trajectories ended
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param dat_aug = TRUE/FALSE (return the simulated sample trajectory with data augmentation or not)

# sel_cof <- matrix(c(1e-02, 5e-03, 1e-03, 5e-03, 1e-02, 1e-03), nrow = 3, ncol = 2)
# rec_rat <- 5e-05
# pop_siz <- c(rep(1e+04, length.out = 201), rep(5e+03, length.out = 200), rep(1e+04, length.out = 100))
# ref_siz <- 1e+04
# int_frq <- c(5e-01, 2e-01, 2e-01, 1e-01)
# evt_gen <- 240
# int_gen <- 0
# lst_gen <- 500
# ptn_num <- 5e+00
#
# frq_pth <- cmpsimulateWFD(sel_cof, rec_rat, pop_siz, ref_siz, int_frq, evt_gen, int_gen, lst_gen, ptn_num, dat_aug = TRUE)
#
# t <- (int_gen:(int_gen + (lst_gen - int_gen) * ptn_num)) / 2 / ref_siz
# plot(t, frq_pth[1, ], type = "l", lwd = 1.5,
#      xlab = "Time", ylab = "Haplotype frequency",
#      main = "WFD: the KM0sb1 haplotype frequency trajectory")
# plot(t, frq_pth[2, ], type = "l", lwd = 1.5,
#      xlab = "Time", ylab = "Haplotype frequency",
#      main = "WFD: the KM0SB1 haplotype frequency trajectory")
# plot(t, frq_pth[3, ], type = "l", lwd = 1.5,
#      xlab = "Time", ylab = "Haplotype frequency",
#      main = "WFD: the KM1sb1 haplotype frequency trajectory")
# plot(t, frq_pth[4, ], type = "l", lwd = 1.5,
#      xlab = "Time", ylab = "Haplotype frequency",
#      main = "WFD: the KM1SB1 haplotype frequency trajectory")

########################################

#' Compare the simulation generated with the Wright-Fisher model and the Wright-Fisher diffusion
# sel_cof <- matrix(c(1e-02, 5e-03, 1e-03, 5e-03, 1e-02, 1e-03), nrow = 3, ncol = 2)
# rec_rat <- 5e-05
# pop_siz <- c(rep(1e+04, length.out = 201), rep(5e+03, length.out = 200), rep(1e+04, length.out = 100))
# ref_siz <- 1e+04
# int_frq <- c(5e-01, 2e-01, 2e-01, 1e-01)
# evt_gen <- 240
# int_gen <- 0
# lst_gen <- 500
# ptn_num <- 5e+00
# sim_num <- 1e+06
#
# smp_WFM <- matrix(NA, nrow = 4, ncol = sim_num)
# smp_WFD <- matrix(NA, nrow = 4, ncol = sim_num)
# for (i in 1:sim_num) {
#   print(i)
#   smp_WFM[, i] <- cmpsimulateWFM(sel_cof, rec_rat, pop_siz, int_frq, evt_gen, int_gen, lst_gen)$hap_frq_pth[, (lst_gen - int_gen) + 1]
#   smp_WFD[, i] <- cmpsimulateWFD(sel_cof, rec_rat, pop_siz, ref_siz, int_frq, evt_gen, int_gen, lst_gen, ptn_num, dat_aug = FALSE)[, (lst_gen - int_gen) + 1]
# }
#
# hist(smp_WFM[1, ], breaks = seq(min(smp_WFM[1, ], smp_WFD[1, ]), max(smp_WFM[1, ], smp_WFD[1, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
#      xlim = c(min(smp_WFM[1, ], smp_WFD[1, ]), max(smp_WFM[1, ], smp_WFD[1, ])),
#      xlab = "Haplotype frequency", main = paste("Histogram of the KM0sb1 haplotype at generation", lst_gen))
# hist(smp_WFD[1, ], breaks = seq(min(smp_WFM[1, ], smp_WFD[1, ]), max(smp_WFM[1, ], smp_WFD[1, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
# hist(smp_WFM[2, ], breaks = seq(min(smp_WFM[2, ], smp_WFD[2, ]), max(smp_WFM[2, ], smp_WFD[2, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
#      xlim = c(min(smp_WFM[2, ], smp_WFD[2, ]), max(smp_WFM[2, ], smp_WFD[2, ])),
#      xlab = "Haplotype frequency", main = paste("Histogram of the KM0SB1 haplotype at generation", lst_gen))
# hist(smp_WFD[2, ], breaks = seq(min(smp_WFM[2, ], smp_WFD[2, ]), max(smp_WFM[2, ], smp_WFD[2, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
# hist(smp_WFM[3, ], breaks = seq(min(smp_WFM[3, ], smp_WFD[3, ]), max(smp_WFM[3, ], smp_WFD[3, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
#      xlim = c(min(smp_WFM[3, ], smp_WFD[3, ]), max(smp_WFM[3, ], smp_WFD[3, ])),
#      xlab = "Haplotype frequency", main = paste("Histogram of the KM1sb1 haplotype at generation", lst_gen))
# hist(smp_WFD[3, ], breaks = seq(min(smp_WFM[3, ], smp_WFD[3, ]), max(smp_WFM[3, ], smp_WFD[3, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
# hist(smp_WFM[4, ], breaks = seq(min(smp_WFM[4, ], smp_WFD[4, ]), max(smp_WFM[4, ], smp_WFD[4, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
#      xlim = c(min(smp_WFM[4, ], smp_WFD[4, ]), max(smp_WFM[4, ], smp_WFD[4, ])),
#      xlab = "Haplotype frequency", main = paste("Histogram of the KM1SB1 haplotype at generation", lst_gen))
# hist(smp_WFD[4, ], breaks = seq(min(smp_WFM[4, ], smp_WFD[4, ]), max(smp_WFM[4, ], smp_WFD[4, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)

################################################################################

#' Simulate the hidden Markov model
#' Parameter setting
#' @param model = "WFM"/"WFD" (return the observations from the underlying population evolving according to the WFM or the WFD)
#' @param sel_cof the selection coefficients of the tobiano, sabino and mixed phenotypes
#' @param rec_rat the recombination rate between the KIT13 and KIT16 loci
#' @param pop_siz the size of the horse population (non-constant)
#' @param int_con the initial haplotype frequencies of the population / the initial mutant allele frequencies and the linkage disequilibrium of the population
#' @param evt_gen the generation that the event of interest occurred
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_siz the count of the horses drawn from the population at all sampling time points
#' @param obs_hap = TRUE/FALSE (return the simulated sample genotypes with haplotype information or not)
#' @param ref_siz the reference size of the horse population
#' @param ptn_num the number of the subintervals divided per generation in the Euler-Maruyama method for the WFD

#' Simulate the dataset under the Wright-Fisher model
# model <- "WFM"
# sel_cof <- matrix(c(1e-02, 5e-03, 1e-03, 5e-03, 1e-02, 1e-03), nrow = 3, ncol = 2)
# rec_rat <- 5e-05
# pop_siz <- c(rep(1e+04, length.out = 201), rep(5e+03, length.out = 200), rep(1e+04, length.out = 100))
# int_con <- c(5e-01, 2e-01, 2e-01, 1e-01)
# evt_gen <- 240
# smp_gen <- (0:10) * 50
# smp_siz <- rep(100, 11)
# obs_hap <- FALSE
#
# sim_HMM_WFM <- cmpsimulateHMM(model, sel_cof, rec_rat, pop_siz, int_con, evt_gen, smp_gen, smp_siz, obs_hap)
# smp_gen <- sim_HMM_WFM$smp_gen
# smp_siz <- sim_HMM_WFM$smp_siz
# smp_gen_cnt <- sim_HMM_WFM$smp_gen_cnt
# smp_gen_frq <- sim_HMM_WFM$smp_gen_frq
# pop_gen_frq <- sim_HMM_WFM$pop_gen_frq
#
# k <- min(smp_gen):max(smp_gen)
# plot(k, pop_gen_frq[1, ], type = 'l', lwd = 1.5,
#      xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_gen_frq[1, ], pop_gen_frq[1, ]), max(smp_gen_frq[1, ], pop_gen_frq[1, ])),
#      xlab = "Generation", ylab = "Genotype frequency",
#      main = "WFM-HMM: the KM0KM0/sb1sb1 genotype")
# points(smp_gen, smp_gen_frq[1, ], col = 'red', pch = 17, cex = 1)
# plot(k, pop_gen_frq[2, ], type = 'l', lwd = 1.5,
#      xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_gen_frq[2, ], pop_gen_frq[2, ]), max(smp_gen_frq[2, ], pop_gen_frq[2, ])),
#      xlab = "Generation", ylab = "Genotype frequency",
#      main = "WFM-HMM: the KM0KM0/sb1SB1 genotype")
# points(smp_gen, smp_gen_frq[2, ], col = 'red', pch = 17, cex = 1)
# plot(k, pop_gen_frq[3, ], type = 'l', lwd = 1.5,
#      xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_gen_frq[3, ], pop_gen_frq[3, ]), max(smp_gen_frq[3, ], pop_gen_frq[3, ])),
#      xlab = "Generation", ylab = "Genotype frequency",
#      main = "WFM-HMM: the KM0KM0/SB1SB1 genotype")
# points(smp_gen, smp_gen_frq[3, ], col = 'red', pch = 17, cex = 1)
# plot(k, pop_gen_frq[4, ], type = 'l', lwd = 1.5,
#      xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_gen_frq[4, ], pop_gen_frq[4, ]), max(smp_gen_frq[4, ], pop_gen_frq[4, ])),
#      xlab = "Generation", ylab = "Genotype frequency",
#      main = "WFM-HMM: the KM0KM1/sb1sb1 genotype")
# points(smp_gen, smp_gen_frq[4, ], col = 'red', pch = 17, cex = 1)
# plot(k, pop_gen_frq[5, ], type = 'l', lwd = 1.5,
#      xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_gen_frq[5, ], pop_gen_frq[5, ]), max(smp_gen_frq[5, ], pop_gen_frq[5, ])),
#      xlab = "Generation", ylab = "Genotype frequency",
#      main = "WFM-HMM: the KM0KM1/sb1SB1 genotype")
# points(smp_gen, smp_gen_frq[5, ], col = 'red', pch = 17, cex = 1)
# plot(k, pop_gen_frq[7, ], type = 'l', lwd = 1.5,
#      xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_gen_frq[7, ], pop_gen_frq[7, ]), max(smp_gen_frq[7, ], pop_gen_frq[7, ])),
#      xlab = "Generation", ylab = "Genotype frequency",
#      main = "WFM-HMM: the KM0KM1/SB1SB1 genotype")
# points(smp_gen, smp_gen_frq[7, ], col = 'red', pch = 17, cex = 1)
# plot(k, pop_gen_frq[6, ], type = 'l', lwd = 1.5,
#      xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_gen_frq[6, ], pop_gen_frq[6, ]), max(smp_gen_frq[6, ], pop_gen_frq[6, ])),
#      xlab = "Generation", ylab = "Genotype frequency",
#      main = "WFM-HMM: the KM1KM1/sb1sb1 genotype")
# points(smp_gen, smp_gen_frq[6, ], col = 'red', pch = 17, cex = 1)
# plot(k, pop_gen_frq[8, ], type = 'l', lwd = 1.5,
#      xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_gen_frq[8, ], pop_gen_frq[8, ]), max(smp_gen_frq[8, ], pop_gen_frq[8, ])),
#      xlab = "Generation", ylab = "Genotype frequency",
#      main = "WFM-HMM: the KM1KM1/SB1sb1 genotype")
# points(smp_gen, smp_gen_frq[8, ], col = 'red', pch = 17, cex = 1)
# plot(k, pop_gen_frq[9, ], type = 'l', lwd = 1.5,
#      xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_gen_frq[9, ], pop_gen_frq[9, ]), max(smp_gen_frq[9, ], pop_gen_frq[9, ])),
#      xlab = "Generation", ylab = "Genotype frequency",
#      main = "WFM-HMM: the KM1KM1/SB1SB1 genotype")
# points(smp_gen, smp_gen_frq[9, ], col = 'red', pch = 17, cex = 1)

####################

#' Simulate the dataset under the Wright-Fisher diffusion
# model <- "WFD"
# sel_cof <- matrix(c(1e-02, 5e-03, 1e-03, 5e-03, 1e-02, 1e-03), nrow = 3, ncol = 2)
# rec_rat <- 5e-05
# pop_siz <- c(rep(1e+04, length.out = 201), rep(5e+03, length.out = 200), rep(1e+04, length.out = 100))
# int_con <- c(5e-01, 2e-01, 2e-01, 1e-01)
# evt_gen <- 240
# smp_gen <- (0:10) * 50
# smp_siz <- rep(100, 11)
# ref_siz <- 1e+04
# obs_hap <- FALSE
# ptn_num <- 5e+00
#
# sim_HMM_WFD <- cmpsimulateHMM(model, sel_cof, rec_rat, pop_siz, int_con, evt_gen, smp_gen, smp_siz, ref_siz, obs_hap, ptn_num)
# smp_gen <- sim_HMM_WFD$smp_gen
# smp_siz <- sim_HMM_WFD$smp_siz
# smp_gen_cnt <- sim_HMM_WFD$smp_gen_cnt
# smp_gen_frq <- sim_HMM_WFD$smp_gen_frq
# pop_gen_frq <- sim_HMM_WFD$pop_gen_frq
#
# k <- min(smp_gen):max(smp_gen)
# plot(k, pop_gen_frq[1, ], type = 'l', lwd = 1.5,
#      xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_gen_frq[1, ], pop_gen_frq[1, ]), max(smp_gen_frq[1, ], pop_gen_frq[1, ])),
#      xlab = "Generation", ylab = "Genotype frequency",
#      main = "WFD-HMM: the KM0KM0/sb1sb1 genotype")
# points(smp_gen, smp_gen_frq[1, ], col = 'red', pch = 17, cex = 1)
# plot(k, pop_gen_frq[2, ], type = 'l', lwd = 1.5,
#      xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_gen_frq[2, ], pop_gen_frq[2, ]), max(smp_gen_frq[2, ], pop_gen_frq[2, ])),
#      xlab = "Generation", ylab = "Genotype frequency",
#      main = "WFD-HMM: the KM0KM0/sb1SB1 genotype")
# points(smp_gen, smp_gen_frq[2, ], col = 'red', pch = 17, cex = 1)
# plot(k, pop_gen_frq[3, ], type = 'l', lwd = 1.5,
#      xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_gen_frq[3, ], pop_gen_frq[3, ]), max(smp_gen_frq[3, ], pop_gen_frq[3, ])),
#      xlab = "Generation", ylab = "Genotype frequency",
#      main = "WFD-HMM: the KM0KM0/SB1SB1 genotype")
# points(smp_gen, smp_gen_frq[3, ], col = 'red', pch = 17, cex = 1)
# plot(k, pop_gen_frq[4, ], type = 'l', lwd = 1.5,
#      xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_gen_frq[4, ], pop_gen_frq[4, ]), max(smp_gen_frq[4, ], pop_gen_frq[4, ])),
#      xlab = "Generation", ylab = "Genotype frequency",
#      main = "WFD-HMM: the KM0KM1/sb1sb1 genotype")
# points(smp_gen, smp_gen_frq[4, ], col = 'red', pch = 17, cex = 1)
# plot(k, pop_gen_frq[5, ], type = 'l', lwd = 1.5,
#      xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_gen_frq[5, ], pop_gen_frq[5, ]), max(smp_gen_frq[5, ], pop_gen_frq[5, ])),
#      xlab = "Generation", ylab = "Genotype frequency",
#      main = "WFD-HMM: the KM0KM1/sb1SB1 genotype")
# points(smp_gen, smp_gen_frq[5, ], col = 'red', pch = 17, cex = 1)
# plot(k, pop_gen_frq[7, ], type = 'l', lwd = 1.5,
#      xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_gen_frq[7, ], pop_gen_frq[7, ]), max(smp_gen_frq[7, ], pop_gen_frq[7, ])),
#      xlab = "Generation", ylab = "Genotype frequency",
#      main = "WFD-HMM: the KM0KM1/SB1SB1 genotype")
# points(smp_gen, smp_gen_frq[7, ], col = 'red', pch = 17, cex = 1)
# plot(k, pop_gen_frq[6, ], type = 'l', lwd = 1.5,
#      xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_gen_frq[6, ], pop_gen_frq[6, ]), max(smp_gen_frq[6, ], pop_gen_frq[6, ])),
#      xlab = "Generation", ylab = "Genotype frequency",
#      main = "WFD-HMM: the KM1KM1/sb1sb1 genotype")
# points(smp_gen, smp_gen_frq[6, ], col = 'red', pch = 17, cex = 1)
# plot(k, pop_gen_frq[8, ], type = 'l', lwd = 1.5,
#      xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_gen_frq[8, ], pop_gen_frq[8, ]), max(smp_gen_frq[8, ], pop_gen_frq[8, ])),
#      xlab = "Generation", ylab = "Genotype frequency",
#      main = "WFD-HMM: the KM1KM1/SB1sb1 genotype")
# points(smp_gen, smp_gen_frq[8, ], col = 'red', pch = 17, cex = 1)
# plot(k, pop_gen_frq[9, ], type = 'l', lwd = 1.5,
#      xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_gen_frq[9, ], pop_gen_frq[9, ]), max(smp_gen_frq[9, ], pop_gen_frq[9, ])),
#      xlab = "Generation", ylab = "Genotype frequency",
#      main = "WFD-HMM: the KM1KM1/SB1SB1 genotype")
# points(smp_gen, smp_gen_frq[9, ], col = 'red', pch = 17, cex = 1)

################################################################################

#' Generate a simulated dataset under the Wright-Fisher model
test_seed <- 17
set.seed(test_seed)

model <- "WFM"
sel_cof <- matrix(c(1e-02, 5e-03, 1e-03, 5e-03, 1e-02, 1e-03), nrow = 3, ncol = 2)
rec_rat <- 5e-05
pop_siz <- c(rep(1e+04, length.out = 201), rep(5e+03, length.out = 200), rep(1e+04, length.out = 100))
int_con <- c(5e-01, 2e-01, 2e-01, 1e-01)
evt_gen <- 240
smp_gen <- (0:10) * 50
smp_siz <- rep(100, 11)
obs_hap <- FALSE

sim_HMM_WFM <- cmpsimulateHMM(model, sel_cof, rec_rat, pop_siz, int_con, evt_gen, smp_gen, smp_siz, obs_hap)
smp_gen <- sim_HMM_WFM$smp_gen
smp_siz <- sim_HMM_WFM$smp_siz
smp_cnt <- sim_HMM_WFM$smp_gen_cnt
smp_frq <- sim_HMM_WFM$smp_gen_frq
pop_frq <- sim_HMM_WFM$pop_gen_frq
# pop_frq[5, ] <- pop_frq[5, ] + pop_frq[7, ]
# pop_frq <- pop_frq[-7, ]

save(model, sel_cof, rec_rat, pop_siz, int_con, evt_gen, smp_gen, smp_siz, obs_hap, smp_cnt, smp_frq, pop_frq,
     file = "./Output/Output v1.0/Test v1.2/TEST_PTN_SimData.rda")

load("./Output/Output v1.0/Test v1.2/TEST_PTN_SimData.rda")

pdf(file = "./Output/Output v1.0/Test v1.2/TEST_PTN_SimData.pdf", width = 24, height = 18)
par(mfrow = c(3, 3), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
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

plot(k, pop_frq[7, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_frq[7, ], pop_frq[7, ]), max(smp_frq[7, ], pop_frq[7, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "Genotype KM0KM1/SB1SB1")
points(smp_gen, smp_frq[7, ], col = 'red', pch = 17, cex = 1)

plot(k, pop_frq[6, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_frq[6, ], pop_frq[6, ]), max(smp_frq[6, ], pop_frq[6, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "Genotype KM1KM1/sb1sb1")
points(smp_gen, smp_frq[6, ], col = 'red', pch = 17, cex = 1)

plot(k, pop_frq[8, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_frq[8, ], pop_frq[8, ]), max(smp_frq[8, ], pop_frq[8, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "Genotype KM1KM1/SB1sb1")
points(smp_gen, smp_frq[8, ], col = 'red', pch = 17, cex = 1)

plot(k, pop_frq[9, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_frq[9, ], pop_frq[9, ]), max(smp_frq[9, ], pop_frq[9, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "Genotype KM1KM1/SB1SB1")
points(smp_gen, smp_frq[9, ], col = 'red', pch = 17, cex = 1)
dev.off()

########################################

#' Run the bootstrap particle filter (BPF) with the two-locus Wright-Fisher diffusion with selection
#' Parameter setting
#' @param sel_cof the selection coefficients of the tobiano, sabino and mixed phenotypes
#' @param rec_rat the recombination rate between the KIT13 and KIT16 loci
#' @param pop_siz the size of the horse population (non-constant)
#' @param ref_siz the reference size of the horse population
#' @param evt_gen the generation that the event of interest occurred
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_siz the count of the horses drawn from the population at all sampling time points
#' @param smp_cnt the count of the genotypes observed in the sample at all sampling time points
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter

load("./Output/Output v1.0/Test v1.2/TEST_PTN_SimData.rda")

set.seed(test_seed)

sel_cof
rec_rat
pop_siz
ref_siz <- 1e+04
evt_gen
smp_gen
smp_siz
smp_cnt
ptn_num <- 5e+00
pcl_num <- 1e+05

system.time(BPF <- cmprunBPF(sel_cof, rec_rat, pop_siz, ref_siz, evt_gen, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num))

save(sel_cof, rec_rat, pop_siz, ref_siz, evt_gen, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, BPF,
     file = "./Output/Output v1.0/Test v1.2/TEST_PTN_BPF.rda")

load("./Output/Output v1.0/Test v1.2/TEST_PTN_BPF.rda")

lik <- rep(1, pcl_num)
wght <- BPF$wght
for (k in 1:length(smp_gen)) {
  lik <- lik * (cumsum(wght[, k]) / (1:pcl_num))
}

pdf(file = "./Output/Output v1.0/Test v1.2/TEST_PTN_BPF_Likelihood.pdf", width = 8, height = 6)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:pcl_num, log(lik), type = 'l',
     xlab = "Number of particles", ylab = "Log likelihood",
     main = "Log likelihood through the bootstrap particle filter")
dev.off()

pop_frq_pre_resmp <- BPF$gen_frq_pre_resmp
pop_frq_pst_resmp <- BPF$gen_frq_pst_resmp

pdf(file = "./Output/Output v1.0/Test v1.2/TEST_PTN_BPF_Particle.pdf", width = 72, height = 66)
par(mfrow = c(11, 9), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
for (k in 1:length(smp_gen)) {
  hist_pst_resmp <- hist(pop_frq_pst_resmp[1, , k], breaks = seq(min(pop_frq_pst_resmp[1, , k], pop_frq_pre_resmp[1, , k]), max(pop_frq_pst_resmp[1, , k], pop_frq_pre_resmp[1, , k]), length.out = 50), plot = FALSE)
  hist_pre_resmp <- hist(pop_frq_pre_resmp[1, , k], breaks = seq(min(pop_frq_pst_resmp[1, , k], pop_frq_pre_resmp[1, , k]), max(pop_frq_pst_resmp[1, , k], pop_frq_pre_resmp[1, , k]), length.out = 50), plot = FALSE)
  hist(pop_frq_pst_resmp[1, , k], breaks = seq(min(pop_frq_pst_resmp[1, , k], pop_frq_pre_resmp[1, , k]), max(pop_frq_pst_resmp[1, , k], pop_frq_pre_resmp[1, , k]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
       xlim = c(min(pop_frq_pst_resmp[1, , k], pop_frq_pre_resmp[1, , k], smp_frq[1, k]), max(pop_frq_pst_resmp[1, , k], pop_frq_pre_resmp[1, , k], smp_frq[1, k])), ylim = c(0, max(hist_pst_resmp$density, hist_pre_resmp$density)),
       xlab = "Genotype frequency",
       main = paste("Genotype KM0KM0/sb1sb1 at generation", smp_gen[k]))
  hist(pop_frq_pre_resmp[1, , k], breaks = seq(min(pop_frq_pst_resmp[1, , k], pop_frq_pre_resmp[1, , k]), max(pop_frq_pst_resmp[1, , k], pop_frq_pre_resmp[1, , k]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
  abline(v = smp_frq[1, k], col = 'red', lty = 2, lwd = 2)

  hist_pst_resmp <- hist(pop_frq_pst_resmp[2, , k], breaks = seq(min(pop_frq_pst_resmp[2, , k], pop_frq_pre_resmp[2, , k]), max(pop_frq_pst_resmp[2, , k], pop_frq_pre_resmp[2, , k]), length.out = 50), plot = FALSE)
  hist_pre_resmp <- hist(pop_frq_pre_resmp[2, , k], breaks = seq(min(pop_frq_pst_resmp[2, , k], pop_frq_pre_resmp[2, , k]), max(pop_frq_pst_resmp[2, , k], pop_frq_pre_resmp[2, , k]), length.out = 50), plot = FALSE)
  hist(pop_frq_pst_resmp[2, , k], breaks = seq(min(pop_frq_pst_resmp[2, , k], pop_frq_pre_resmp[2, , k]), max(pop_frq_pst_resmp[2, , k], pop_frq_pre_resmp[2, , k]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
       xlim = c(min(pop_frq_pst_resmp[2, , k], pop_frq_pre_resmp[2, , k], smp_frq[2, k]), max(pop_frq_pst_resmp[2, , k], pop_frq_pre_resmp[2, , k], smp_frq[2, k])), ylim = c(0, max(hist_pst_resmp$density, hist_pre_resmp$density)),
       xlab = "Genotype frequency",
       main = paste("Genotype KM0KM0/sb1SB1 at generation", smp_gen[k]))
  hist(pop_frq_pre_resmp[2, , k], breaks = seq(min(pop_frq_pst_resmp[2, , k], pop_frq_pre_resmp[2, , k]), max(pop_frq_pst_resmp[2, , k], pop_frq_pre_resmp[2, , k]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
  abline(v = smp_frq[2, k], col = 'red', lty = 2, lwd = 2)

  hist_pst_resmp <- hist(pop_frq_pst_resmp[3, , k], breaks = seq(min(pop_frq_pst_resmp[3, , k], pop_frq_pre_resmp[3, , k]), max(pop_frq_pst_resmp[3, , k], pop_frq_pre_resmp[3, , k]), length.out = 50), plot = FALSE)
  hist_pre_resmp <- hist(pop_frq_pre_resmp[3, , k], breaks = seq(min(pop_frq_pst_resmp[3, , k], pop_frq_pre_resmp[3, , k]), max(pop_frq_pst_resmp[3, , k], pop_frq_pre_resmp[3, , k]), length.out = 50), plot = FALSE)
  hist(pop_frq_pst_resmp[3, , k], breaks = seq(min(pop_frq_pst_resmp[3, , k], pop_frq_pre_resmp[3, , k]), max(pop_frq_pst_resmp[3, , k], pop_frq_pre_resmp[3, , k]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
       xlim = c(min(pop_frq_pst_resmp[3, , k], pop_frq_pre_resmp[3, , k], smp_frq[3, k]), max(pop_frq_pst_resmp[3, , k], pop_frq_pre_resmp[3, , k], smp_frq[3, k])), ylim = c(0, max(hist_pst_resmp$density, hist_pre_resmp$density)),
       xlab = "Genotype frequency",
       main = paste("Genotype KM0KM0/SB1SB1 at generation", smp_gen[k]))
  hist(pop_frq_pre_resmp[3, , k], breaks = seq(min(pop_frq_pst_resmp[3, , k], pop_frq_pre_resmp[3, , k]), max(pop_frq_pst_resmp[3, , k], pop_frq_pre_resmp[3, , k]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
  abline(v = smp_frq[3, k], col = 'red', lty = 2, lwd = 2)

  hist_pst_resmp <- hist(pop_frq_pst_resmp[4, , k], breaks = seq(min(pop_frq_pst_resmp[4, , k], pop_frq_pre_resmp[4, , k]), max(pop_frq_pst_resmp[4, , k], pop_frq_pre_resmp[4, , k]), length.out = 50), plot = FALSE)
  hist_pre_resmp <- hist(pop_frq_pre_resmp[4, , k], breaks = seq(min(pop_frq_pst_resmp[4, , k], pop_frq_pre_resmp[4, , k]), max(pop_frq_pst_resmp[4, , k], pop_frq_pre_resmp[4, , k]), length.out = 50), plot = FALSE)
  hist(pop_frq_pst_resmp[4, , k], breaks = seq(min(pop_frq_pst_resmp[4, , k], pop_frq_pre_resmp[4, , k]), max(pop_frq_pst_resmp[4, , k], pop_frq_pre_resmp[4, , k]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
       xlim = c(min(pop_frq_pst_resmp[4, , k], pop_frq_pre_resmp[4, , k], smp_frq[4, k]), max(pop_frq_pst_resmp[4, , k], pop_frq_pre_resmp[4, , k], smp_frq[4, k])), ylim = c(0, max(hist_pst_resmp$density, hist_pre_resmp$density)),
       xlab = "Genotype frequency",
       main = paste("Genotype KM0KM1/sb1sb1 at generation", smp_gen[k]))
  hist(pop_frq_pre_resmp[4, , k], breaks = seq(min(pop_frq_pst_resmp[4, , k], pop_frq_pre_resmp[4, , k]), max(pop_frq_pst_resmp[4, , k], pop_frq_pre_resmp[4, , k]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
  abline(v = smp_frq[4, k], col = 'red', lty = 2, lwd = 2)

  hist_pst_resmp <- hist(pop_frq_pst_resmp[5, , k], breaks = seq(min(pop_frq_pst_resmp[5, , k], pop_frq_pre_resmp[5, , k]), max(pop_frq_pst_resmp[5, , k], pop_frq_pre_resmp[5, , k]), length.out = 50), plot = FALSE)
  hist_pre_resmp <- hist(pop_frq_pre_resmp[5, , k], breaks = seq(min(pop_frq_pst_resmp[5, , k], pop_frq_pre_resmp[5, , k]), max(pop_frq_pst_resmp[5, , k], pop_frq_pre_resmp[5, , k]), length.out = 50), plot = FALSE)
  hist(pop_frq_pst_resmp[5, , k], breaks = seq(min(pop_frq_pst_resmp[5, , k], pop_frq_pre_resmp[5, , k]), max(pop_frq_pst_resmp[5, , k], pop_frq_pre_resmp[5, , k]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
       xlim = c(min(pop_frq_pst_resmp[5, , k], pop_frq_pre_resmp[5, , k], smp_frq[5, k]), max(pop_frq_pst_resmp[5, , k], pop_frq_pre_resmp[5, , k], smp_frq[5, k])), ylim = c(0, max(hist_pst_resmp$density, hist_pre_resmp$density)),
       xlab = "Genotype frequency",
       main = paste("Genotype KM0KM1/sb1SB1 at generation", smp_gen[k]))
  hist(pop_frq_pre_resmp[5, , k], breaks = seq(min(pop_frq_pst_resmp[5, , k], pop_frq_pre_resmp[5, , k]), max(pop_frq_pst_resmp[5, , k], pop_frq_pre_resmp[5, , k]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
  abline(v = smp_frq[5, k], col = 'red', lty = 2, lwd = 2)

  hist_pst_resmp <- hist(pop_frq_pst_resmp[7, , k], breaks = seq(min(pop_frq_pst_resmp[7, , k], pop_frq_pre_resmp[7, , k]), max(pop_frq_pst_resmp[7, , k], pop_frq_pre_resmp[7, , k]), length.out = 50), plot = FALSE)
  hist_pre_resmp <- hist(pop_frq_pre_resmp[7, , k], breaks = seq(min(pop_frq_pst_resmp[7, , k], pop_frq_pre_resmp[7, , k]), max(pop_frq_pst_resmp[7, , k], pop_frq_pre_resmp[7, , k]), length.out = 50), plot = FALSE)
  hist(pop_frq_pst_resmp[7, , k], breaks = seq(min(pop_frq_pst_resmp[7, , k], pop_frq_pre_resmp[7, , k]), max(pop_frq_pst_resmp[7, , k], pop_frq_pre_resmp[7, , k]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
       xlim = c(min(pop_frq_pst_resmp[7, , k], pop_frq_pre_resmp[7, , k], smp_frq[7, k]), max(pop_frq_pst_resmp[7, , k], pop_frq_pre_resmp[7, , k], smp_frq[7, k])), ylim = c(0, max(hist_pst_resmp$density, hist_pre_resmp$density)),
       xlab = "Genotype frequency",
       main = paste("Genotype KM0KM1/SB1SB1 at generation", smp_gen[k]))
  hist(pop_frq_pre_resmp[7, , k], breaks = seq(min(pop_frq_pst_resmp[7, , k], pop_frq_pre_resmp[7, , k]), max(pop_frq_pst_resmp[7, , k], pop_frq_pre_resmp[7, , k]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
  abline(v = smp_frq[7, k], col = 'red', lty = 2, lwd = 2)

  hist_pst_resmp <- hist(pop_frq_pst_resmp[6, , k], breaks = seq(min(pop_frq_pst_resmp[6, , k], pop_frq_pre_resmp[6, , k]), max(pop_frq_pst_resmp[6, , k], pop_frq_pre_resmp[6, , k]), length.out = 50), plot = FALSE)
  hist_pre_resmp <- hist(pop_frq_pre_resmp[6, , k], breaks = seq(min(pop_frq_pst_resmp[6, , k], pop_frq_pre_resmp[6, , k]), max(pop_frq_pst_resmp[6, , k], pop_frq_pre_resmp[6, , k]), length.out = 50), plot = FALSE)
  hist(pop_frq_pst_resmp[6, , k], breaks = seq(min(pop_frq_pst_resmp[6, , k], pop_frq_pre_resmp[6, , k]), max(pop_frq_pst_resmp[6, , k], pop_frq_pre_resmp[6, , k]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
       xlim = c(min(pop_frq_pst_resmp[6, , k], pop_frq_pre_resmp[6, , k], smp_frq[6, k]), max(pop_frq_pst_resmp[6, , k], pop_frq_pre_resmp[6, , k], smp_frq[6, k])), ylim = c(0, max(hist_pst_resmp$density, hist_pre_resmp$density)),
       xlab = "Genotype frequency",
       main = paste("Genotype KM1KM1/sb1sb1 at generation", smp_gen[k]))
  hist(pop_frq_pre_resmp[6, , k], breaks = seq(min(pop_frq_pst_resmp[6, , k], pop_frq_pre_resmp[6, , k]), max(pop_frq_pst_resmp[6, , k], pop_frq_pre_resmp[6, , k]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
  abline(v = smp_frq[6, k], col = 'red', lty = 2, lwd = 2)

  hist_pst_resmp <- hist(pop_frq_pst_resmp[8, , k], breaks = seq(min(pop_frq_pst_resmp[8, , k], pop_frq_pre_resmp[8, , k]), max(pop_frq_pst_resmp[8, , k], pop_frq_pre_resmp[8, , k]), length.out = 50), plot = FALSE)
  hist_pre_resmp <- hist(pop_frq_pre_resmp[8, , k], breaks = seq(min(pop_frq_pst_resmp[8, , k], pop_frq_pre_resmp[8, , k]), max(pop_frq_pst_resmp[8, , k], pop_frq_pre_resmp[8, , k]), length.out = 50), plot = FALSE)
  hist(pop_frq_pst_resmp[8, , k], breaks = seq(min(pop_frq_pst_resmp[8, , k], pop_frq_pre_resmp[8, , k]), max(pop_frq_pst_resmp[8, , k], pop_frq_pre_resmp[8, , k]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
       xlim = c(min(pop_frq_pst_resmp[8, , k], pop_frq_pre_resmp[8, , k], smp_frq[8, k]), max(pop_frq_pst_resmp[8, , k], pop_frq_pre_resmp[8, , k], smp_frq[8, k])), ylim = c(0, max(hist_pst_resmp$density, hist_pre_resmp$density)),
       xlab = "Genotype frequency",
       main = paste("Genotype KM1KM1/SB1sb1 at generation", smp_gen[k]))
  hist(pop_frq_pre_resmp[8, , k], breaks = seq(min(pop_frq_pst_resmp[8, , k], pop_frq_pre_resmp[8, , k]), max(pop_frq_pst_resmp[8, , k], pop_frq_pre_resmp[8, , k]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
  abline(v = smp_frq[8, k], col = 'red', lty = 2, lwd = 2)

  hist_pst_resmp <- hist(pop_frq_pst_resmp[9, , k], breaks = seq(min(pop_frq_pst_resmp[9, , k], pop_frq_pre_resmp[9, , k]), max(pop_frq_pst_resmp[9, , k], pop_frq_pre_resmp[9, , k]), length.out = 50), plot = FALSE)
  hist_pre_resmp <- hist(pop_frq_pre_resmp[9, , k], breaks = seq(min(pop_frq_pst_resmp[9, , k], pop_frq_pre_resmp[9, , k]), max(pop_frq_pst_resmp[9, , k], pop_frq_pre_resmp[9, , k]), length.out = 50), plot = FALSE)
  hist(pop_frq_pst_resmp[9, , k], breaks = seq(min(pop_frq_pst_resmp[9, , k], pop_frq_pre_resmp[9, , k]), max(pop_frq_pst_resmp[9, , k], pop_frq_pre_resmp[9, , k]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
       xlim = c(min(pop_frq_pst_resmp[9, , k], pop_frq_pre_resmp[9, , k], smp_frq[9, k]), max(pop_frq_pst_resmp[9, , k], pop_frq_pre_resmp[9, , k], smp_frq[9, k])), ylim = c(0, max(hist_pst_resmp$density, hist_pre_resmp$density)),
       xlab = "Genotype frequency",
       main = paste("Genotype KM1KM1/SB1SB1 at generation", smp_gen[k]))
  hist(pop_frq_pre_resmp[9, , k], breaks = seq(min(pop_frq_pst_resmp[9, , k], pop_frq_pre_resmp[9, , k]), max(pop_frq_pst_resmp[9, , k], pop_frq_pre_resmp[9, , k]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
  abline(v = smp_frq[9, k], col = 'red', lty = 2, lwd = 2)
}
dev.off()

########################################

#' Calculate the optimal particle number in the particle marginal Metropolis-Hastings (PMMH)
#' Parameter settings
#' @param sel_cof the selection coefficients of the tobiano, sabino and mixed phenotypes
#' @param rec_rat the recombination rate between the KIT13 and KIT16 loci
#' @param pop_siz the size of the horse population (non-constant)
#' @param ref_siz the reference size of the horse population
#' @param evt_gen the generation that the event of interest occurred
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_siz the count of the horses drawn from the population at all sampling time points
#' @param smp_cnt the count of the genotypes observed in the sample at all sampling time points
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter
#' @param gap_num the number of particles increased or decreased in the optimal particle number search

load("./Output/Output v1.0/Test v1.2/TEST_PTN_SimData.rda")

set.seed(test_seed)

sel_cof
rec_rat
pop_siz
ref_siz <- 1e+04
evt_gen
smp_gen
smp_siz
smp_cnt
ptn_num <- 5e+00
pcl_num <- 1e+03
gap_num <- 1e+02

system.time(OptNum <- calculateOptimalParticleNum(sel_cof, rec_rat, pop_siz, ref_siz, evt_gen, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, gap_num))

save(sel_cof, rec_rat, pop_siz, ref_siz, evt_gen, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, gap_num, OptNum,
     file = "./Output/Output v1.0/Test v1.2/TEST_PTN_OptNum.rda")

load("./Output/Output v1.0/Test v1.2/TEST_PTN_OptNum.rda")

opt_pcl_num <- OptNum$opt_pcl_num
log_lik_sdv <- OptNum$log_lik_sdv

pdf(file = "./Output/Output v1.0/Test v1.2/TEST_PTN_OptNum.pdf", width = 8, height = 6)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(opt_pcl_num, log_lik_sdv, type = 'b', lwd = 2,
     xlab = "Particle number", ylab = "Log-likelihood standard deviation",
     main = "Optimal particle number in the PMMH")
abline(h = 1.7, col = 'red', lty = 2, lwd = 2)
abline(h = 1.0, col = 'red', lty = 2, lwd = 2)
dev.off()

########################################

#' Run the particle marginal Metropolis-Hastings (PMMH)
#' Parameter settings
#' @param sel_cof the selection coefficients of the tobiano, sabino and mixed phenotypes
#' @param rec_rat the recombination rate between the KIT13 and KIT16 loci
#' @param pop_siz the size of the horse population (non-constant)
#' @param ref_siz the reference size of the horse population
#' @param evt_gen the generation that the event of interest occurred
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_siz the count of the horses drawn from the population at all sampling time points
#' @param smp_cnt the count of the genotypes observed in the sample at all sampling time points
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter
#' @param itn_num the number of the iterations carried out in the particle marginal Metropolis-Hastings

load("./Output/Output v1.0/Test v1.2/TEST_PTN_SimData.rda")

set.seed(test_seed)

sel_cof <- matrix(c(0e+00, 0e+00, 0e+00, 0e+00, 0e+00, 0e+00), nrow = 3, ncol = 2)
rec_rat
pop_siz
ref_siz <- 1e+04
evt_gen
smp_gen
smp_siz
smp_cnt
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 5e+04

system.time(sel_cof_chn <- cmprunPMMH(sel_cof, rec_rat, pop_siz, ref_siz, evt_gen, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num))

load("./Output/Output v1.0/Test v1.2/TEST_PTN_SimData.rda")

save(sel_cof, rec_rat, pop_siz, ref_siz, evt_gen, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, sel_cof_chn,
     file = "./Output/Output v1.0/Test v1.2/TEST_PTN_PMMH.rda")

load("./Output/Output v1.0/Test v1.2/TEST_PTN_PMMH.rda")

pdf(file = "./Output/Output v1.0/Test v1.2/TEST_PTN_PMMH_Traceplot.pdf", width = 24, height = 12)
par(mfrow = c(2, 3), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:itn_num, sel_cof_chn[1, 1, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for sel coeff of tobiano pre-event")
abline(h = sel_cof[1, 1], col = 'red', lty = 2, lwd = 2)

plot(1:itn_num, sel_cof_chn[2, 1, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for sel coeff of sabino pre-event")
abline(h = sel_cof[2, 1], col = 'red', lty = 2, lwd = 2)

plot(1:itn_num, sel_cof_chn[3, 1, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for sel coeff of mixed pre-event")
abline(h = sel_cof[3 ,1], col = 'red', lty = 2, lwd = 2)

plot(1:itn_num, sel_cof_chn[1, 2, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for sel coeff of tobiano post-event")
abline(h = sel_cof[1, 2], col = 'red', lty = 2, lwd = 2)

plot(1:itn_num, sel_cof_chn[2, 2, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for sel coeff of sabino post-event")
abline(h = sel_cof[2, 2], col = 'red', lty = 2, lwd = 2)

plot(1:itn_num, sel_cof_chn[3, 2, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for sel coeff of mixed post-event")
abline(h = sel_cof[3, 2], col = 'red', lty = 2, lwd = 2)
dev.off()

pdf(file = "./Output/Output v1.0/Test v1.2/TEST_PTN_PMMH_Autocorrplot.pdf", width = 24, height = 24)
par(mfrow = c(4, 3), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
# effectiveSize(as.mcmc(t(sel_cof_chn)))

brn_num <- 1e+04
sel_cof_chn <- sel_cof_chn[, , brn_num:dim(sel_cof_chn)[3]]
# effectiveSize(as.mcmc(t(sel_cof_chn)))

autocorr.plot(sel_cof_chn[1, 1, ], lag.max = 50, auto.layout = FALSE, main = "Autocorrelation plot for sel coeff of tobiano pre-event before thinning")
autocorr.plot(sel_cof_chn[2, 1, ], lag.max = 50, auto.layout = FALSE, main = "Autocorrelation plot for sel coeff of sabino pre-event before thinning")
autocorr.plot(sel_cof_chn[3, 1, ], lag.max = 50, auto.layout = FALSE, main = "Autocorrelation plot for sel coeff of mixed pre-event before thinning")
autocorr.plot(sel_cof_chn[1, 2, ], lag.max = 50, auto.layout = FALSE, main = "Autocorrelation plot for sel coeff of tobiano post-event before thinning")
autocorr.plot(sel_cof_chn[2, 2, ], lag.max = 50, auto.layout = FALSE, main = "Autocorrelation plot for sel coeff of sabino post-event before thinning")
autocorr.plot(sel_cof_chn[3, 2, ], lag.max = 50, auto.layout = FALSE, main = "Autocorrelation plot for sel coeff of mixed post-event before thinning")

thn_num <- 8e+00
sel_cof_chn <- sel_cof_chn[, , (1:round(dim(sel_cof_chn)[3] / thn_num)) * thn_num]

autocorr.plot(sel_cof_chn[1, 1, ], lag.max = 50, auto.layout = FALSE, main = "Autocorrelation plot for sel coeff of tobiano pre-event after thinning")
autocorr.plot(sel_cof_chn[2, 1, ], lag.max = 50, auto.layout = FALSE, main = "Autocorrelation plot for sel coeff of sabino pre-event after thinning")
autocorr.plot(sel_cof_chn[3, 1, ], lag.max = 50, auto.layout = FALSE, main = "Autocorrelation plot for sel coeff of mixed pre-event after thinning")
autocorr.plot(sel_cof_chn[1, 2, ], lag.max = 50, auto.layout = FALSE, main = "Autocorrelation plot for sel coeff of tobiano post-event after thinning")
autocorr.plot(sel_cof_chn[2, 2, ], lag.max = 50, auto.layout = FALSE, main = "Autocorrelation plot for sel coeff of sabino post-event after thinning")
autocorr.plot(sel_cof_chn[3, 2, ], lag.max = 50, auto.layout = FALSE, main = "Autocorrelation plot for sel coeff of mixed post-event after thinning")
dev.off()

sel_cof_est <- matrix(NA, nrow = 3, ncol = 2)
sel_cof_est[, 1] <- rowMeans(sel_cof_chn[, 1, ])
sel_cof_est[, 2] <- rowMeans(sel_cof_chn[, 2, ])

sel_cof_hpd <- array(NA, dim = c(3, 2, 2))
sel_cof_hpd[1, , 1] <- HPDinterval(as.mcmc(sel_cof_chn[1, 1, ]), prob = 0.95)
sel_cof_hpd[2, , 1] <- HPDinterval(as.mcmc(sel_cof_chn[2, 1, ]), prob = 0.95)
sel_cof_hpd[3, , 1] <- HPDinterval(as.mcmc(sel_cof_chn[3, 1, ]), prob = 0.95)
sel_cof_hpd[1, , 2] <- HPDinterval(as.mcmc(sel_cof_chn[1, 2, ]), prob = 0.95)
sel_cof_hpd[2, , 2] <- HPDinterval(as.mcmc(sel_cof_chn[2, 2, ]), prob = 0.95)
sel_cof_hpd[3, , 2] <- HPDinterval(as.mcmc(sel_cof_chn[3, 2, ]), prob = 0.95)

pdf(file = "./Output/Output v1.0/Test v1.2/TEST_PTN_PMMH_Posterior.pdf", width = 24, height = 12)
par(mfrow = c(2, 3), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(sel_cof_chn[1, 1, ], breaks = seq(min(sel_cof_chn[1, 1, ]), max(sel_cof_chn[1, 1, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for sel coeff of tobiano pre-event")
lines(density(sel_cof_chn[1, 1, ]), lwd = 2, col = 'black')
abline(v = sel_cof[1, 1], col = 'red', lty = 2, lwd = 2)
abline(v = sel_cof_est[1, 1], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 1, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 2, 1], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_chn[2, 1, ], breaks = seq(min(sel_cof_chn[2, 1, ]), max(sel_cof_chn[2, 1, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for sel coeff of sabino pre-event")
lines(density(sel_cof_chn[2, 1, ]), lwd = 2, col = 'black')
abline(v = sel_cof[2, 1], col = 'red', lty = 2, lwd = 2)
abline(v = sel_cof_est[2, 1], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 1, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 2, 1], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_chn[3, 1, ], breaks = seq(min(sel_cof_chn[3, 1, ]), max(sel_cof_chn[3, 1, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for sel coeff of mixed pre-event")
lines(density(sel_cof_chn[3, 1, ]), lwd = 2, col = 'black')
abline(v = sel_cof[3, 1], col = 'red', lty = 2, lwd = 2)
abline(v = sel_cof_est[3, 1], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[3, 1, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[3, 2, 1], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_chn[1, 2, ], breaks = seq(min(sel_cof_chn[1, 2, ]), max(sel_cof_chn[1, 2, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for sel coeff of tobiano post-event")
lines(density(sel_cof_chn[1, 2, ]), lwd = 2, col = 'black')
abline(v = sel_cof[1, 2], col = 'red', lty = 2, lwd = 2)
abline(v = sel_cof_est[1, 2], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 1, 2], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 2, 2], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_chn[2, 2, ], breaks = seq(min(sel_cof_chn[2, 2, ]), max(sel_cof_chn[2, 2, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for sel coeff of sabino post-event")
lines(density(sel_cof_chn[2, 2, ]), lwd = 2, col = 'black')
abline(v = sel_cof[2, 2], col = 'red', lty = 2, lwd = 2)
abline(v = sel_cof_est[2, 2], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 1, 2], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 2, 2], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_chn[3, 2, ], breaks = seq(min(sel_cof_chn[3, 2, ]), max(sel_cof_chn[3, 2, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for sel coeff of mixed post-event")
lines(density(sel_cof_chn[3, 2, ]), lwd = 2, col = 'black')
abline(v = sel_cof[3, 2], col = 'red', lty = 2, lwd = 2)
abline(v = sel_cof_est[3, 2], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[3, 1, 2], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[3, 2, 2], col = 'blue', lty = 2, lwd = 2)
dev.off()

dif_sel_chn <- sel_cof_chn[, 2, ] - sel_cof_chn[, 1, ]

dif_sel_est <- rowMeans(dif_sel_chn)

dif_sel_hpd <- matrix(NA, nrow = 3, ncol = 2)
dif_sel_hpd[1, ] <- HPDinterval(as.mcmc(dif_sel_chn[1, ]), prob = 0.95)
dif_sel_hpd[2, ] <- HPDinterval(as.mcmc(dif_sel_chn[2, ]), prob = 0.95)
dif_sel_hpd[3, ] <- HPDinterval(as.mcmc(dif_sel_chn[3, ]), prob = 0.95)

pdf(file = "./Output/Output v1.0/Test v1.2/TEST_PTN_PMMH_Posterior_SelChange.pdf", width = 8, height = 18)
par(mfrow = c(3, 1), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(dif_sel_chn[1, ], breaks = seq(min(dif_sel_chn[1, ]), max(dif_sel_chn[1, ]), length.out = 50), freq = FALSE,
     xlab = "Change in the selection coefficient",
     main = "Posterior for change in sel coeff of tobiano")
lines(density(dif_sel_chn[1, ]), lwd = 2, col = 'black')
abline(v = 0, col = 'red', lty = 2, lwd = 2)
abline(v = dif_sel_est[1], col = 'black', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[1, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[1, 2], col = 'blue', lty = 2, lwd = 2)

hist(dif_sel_chn[2, ], breaks = seq(min(dif_sel_chn[2, ]), max(dif_sel_chn[2, ]), length.out = 50), freq = FALSE,
     xlab = "Change in the selection coefficient",
     main = "Posterior for change in sel coeff of sabino")
lines(density(dif_sel_chn[2, ]), lwd = 2, col = 'black')
abline(v = 0, col = 'red', lty = 2, lwd = 2)
abline(v = dif_sel_est[2], col = 'black', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[2, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[2, 2], col = 'blue', lty = 2, lwd = 2)

hist(dif_sel_chn[3, ], breaks = seq(min(dif_sel_chn[3, ]), max(dif_sel_chn[3, ]), length.out = 50), freq = FALSE,
     xlab = "Change in the selection coefficient",
     main = "Posterior for change in sel coeff of mixed")
lines(density(dif_sel_chn[3, ]), lwd = 2, col = 'black')
abline(v = 0, col = 'red', lty = 2, lwd = 2)
abline(v = dif_sel_est[3], col = 'black', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[3, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[3, 2], col = 'blue', lty = 2, lwd = 2)
dev.off()

########################################

#' Run the adaptive particle marginal Metropolis-Hastings (PMMH)
#' Parameter settings
#' @param sel_cof the selection coefficients of the tobiano, sabino and mixed phenotypes
#' @param rec_rat the recombination rate between the KIT13 and KIT16 loci
#' @param pop_siz the size of the horse population (non-constant)
#' @param ref_siz the reference size of the horse population
#' @param evt_gen the generation that the event of interest occurred
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_siz the count of the horses drawn from the population at all sampling time points
#' @param smp_cnt the count of the genotypes observed in the sample at all sampling time points
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter
#' @param itn_num the number of the iterations carried out in the PMMH
#' @param stp_siz the step size sequence in the adaptive setting (decaying to zero)
#' @param apt_rto the target mean acceptance probability of the adaptive setting

load("./Output/Output v1.0/Test v1.2/TEST_PTN_SimData.rda")

set.seed(test_seed)

sel_cof <- matrix(c(0e+00, 0e+00, 0e+00, 0e+00, 0e+00, 0e+00), nrow = 3, ncol = 2)
rec_rat
pop_siz
ref_siz <- 1e+04
evt_gen
smp_gen
smp_siz
smp_cnt
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 5e+04
stp_siz <- (1:itn_num)^(-2 / 3)
apt_rto <- 4e-01

system.time(sel_cof_chn <- cmprunAdaptPMMH(sel_cof, rec_rat, pop_siz, ref_siz, evt_gen, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, stp_siz, apt_rto))

load("./Output/Output v1.0/Test v1.2/TEST_PTN_SimData.rda")

save(sel_cof, rec_rat, pop_siz, ref_siz, evt_gen, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, stp_siz, apt_rto, sel_cof_chn,
     file = "./Output/Output v1.0/Test v1.2/TEST_PTN_AdaptPMMH.rda")

load("./Output/Output v1.0/Test v1.2/TEST_PTN_AdaptPMMH.rda")

pdf(file = "./Output/Output v1.0/Test v1.2/TEST_PTN_AdaptPMMH_Traceplot.pdf", width = 24, height = 12)
par(mfrow = c(2, 3), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:itn_num, sel_cof_chn[1, 1, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for sel coeff of tobiano pre-event")
abline(h = sel_cof[1, 1], col = 'red', lty = 2, lwd = 2)

plot(1:itn_num, sel_cof_chn[2, 1, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for sel coeff of sabino pre-event")
abline(h = sel_cof[2, 1], col = 'red', lty = 2, lwd = 2)

plot(1:itn_num, sel_cof_chn[3, 1, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for sel coeff of mixed pre-event")
abline(h = sel_cof[3 ,1], col = 'red', lty = 2, lwd = 2)

plot(1:itn_num, sel_cof_chn[1, 2, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for sel coeff of tobiano post-event")
abline(h = sel_cof[1, 2], col = 'red', lty = 2, lwd = 2)

plot(1:itn_num, sel_cof_chn[2, 2, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for sel coeff of sabino post-event")
abline(h = sel_cof[2, 2], col = 'red', lty = 2, lwd = 2)

plot(1:itn_num, sel_cof_chn[3, 2, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for sel coeff of mixed post-event")
abline(h = sel_cof[3, 2], col = 'red', lty = 2, lwd = 2)
dev.off()

pdf(file = "./Output/Output v1.0/Test v1.2/TEST_PTN_AdaptPMMH_Autocorrplot.pdf", width = 24, height = 24)
par(mfrow = c(4, 3), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
# effectiveSize(as.mcmc(t(sel_cof_chn)))

brn_num <- 1e+04
sel_cof_chn <- sel_cof_chn[, , brn_num:dim(sel_cof_chn)[3]]
# effectiveSize(as.mcmc(t(sel_cof_chn)))

autocorr.plot(sel_cof_chn[1, 1, ], lag.max = 50, auto.layout = FALSE, main = "Autocorrelation plot for sel coeff of tobiano pre-event before thinning")
autocorr.plot(sel_cof_chn[2, 1, ], lag.max = 50, auto.layout = FALSE, main = "Autocorrelation plot for sel coeff of sabino pre-event before thinning")
autocorr.plot(sel_cof_chn[3, 1, ], lag.max = 50, auto.layout = FALSE, main = "Autocorrelation plot for sel coeff of mixed pre-event before thinning")
autocorr.plot(sel_cof_chn[1, 2, ], lag.max = 50, auto.layout = FALSE, main = "Autocorrelation plot for sel coeff of tobiano post-event before thinning")
autocorr.plot(sel_cof_chn[2, 2, ], lag.max = 50, auto.layout = FALSE, main = "Autocorrelation plot for sel coeff of sabino post-event before thinning")
autocorr.plot(sel_cof_chn[3, 2, ], lag.max = 50, auto.layout = FALSE, main = "Autocorrelation plot for sel coeff of mixed post-event before thinning")

thn_num <- 8e+00
sel_cof_chn <- sel_cof_chn[, , (1:round(dim(sel_cof_chn)[3] / thn_num)) * thn_num]

autocorr.plot(sel_cof_chn[1, 1, ], lag.max = 50, auto.layout = FALSE, main = "Autocorrelation plot for sel coeff of tobiano pre-event after thinning")
autocorr.plot(sel_cof_chn[2, 1, ], lag.max = 50, auto.layout = FALSE, main = "Autocorrelation plot for sel coeff of sabino pre-event after thinning")
autocorr.plot(sel_cof_chn[3, 1, ], lag.max = 50, auto.layout = FALSE, main = "Autocorrelation plot for sel coeff of mixed pre-event after thinning")
autocorr.plot(sel_cof_chn[1, 2, ], lag.max = 50, auto.layout = FALSE, main = "Autocorrelation plot for sel coeff of tobiano post-event after thinning")
autocorr.plot(sel_cof_chn[2, 2, ], lag.max = 50, auto.layout = FALSE, main = "Autocorrelation plot for sel coeff of sabino post-event after thinning")
autocorr.plot(sel_cof_chn[3, 2, ], lag.max = 50, auto.layout = FALSE, main = "Autocorrelation plot for sel coeff of mixed post-event after thinning")
dev.off()

sel_cof_est <- matrix(NA, nrow = 3, ncol = 2)
sel_cof_est[, 1] <- rowMeans(sel_cof_chn[, 1, ])
sel_cof_est[, 2] <- rowMeans(sel_cof_chn[, 2, ])

sel_cof_hpd <- array(NA, dim = c(3, 2, 2))
sel_cof_hpd[1, , 1] <- HPDinterval(as.mcmc(sel_cof_chn[1, 1, ]), prob = 0.95)
sel_cof_hpd[2, , 1] <- HPDinterval(as.mcmc(sel_cof_chn[2, 1, ]), prob = 0.95)
sel_cof_hpd[3, , 1] <- HPDinterval(as.mcmc(sel_cof_chn[3, 1, ]), prob = 0.95)
sel_cof_hpd[1, , 2] <- HPDinterval(as.mcmc(sel_cof_chn[1, 2, ]), prob = 0.95)
sel_cof_hpd[2, , 2] <- HPDinterval(as.mcmc(sel_cof_chn[2, 2, ]), prob = 0.95)
sel_cof_hpd[3, , 2] <- HPDinterval(as.mcmc(sel_cof_chn[3, 2, ]), prob = 0.95)

pdf(file = "./Output/Output v1.0/Test v1.2/TEST_PTN_AdaptPMMH_Posterior.pdf", width = 24, height = 12)
par(mfrow = c(2, 3), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(sel_cof_chn[1, 1, ], breaks = seq(min(sel_cof_chn[1, 1, ]), max(sel_cof_chn[1, 1, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for sel coeff of tobiano pre-event")
lines(density(sel_cof_chn[1, 1, ]), lwd = 2, col = 'black')
abline(v = sel_cof[1, 1], col = 'red', lty = 2, lwd = 2)
abline(v = sel_cof_est[1, 1], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 1, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 2, 1], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_chn[2, 1, ], breaks = seq(min(sel_cof_chn[2, 1, ]), max(sel_cof_chn[2, 1, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for sel coeff of sabino pre-event")
lines(density(sel_cof_chn[2, 1, ]), lwd = 2, col = 'black')
abline(v = sel_cof[2, 1], col = 'red', lty = 2, lwd = 2)
abline(v = sel_cof_est[2, 1], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 1, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 2, 1], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_chn[3, 1, ], breaks = seq(min(sel_cof_chn[3, 1, ]), max(sel_cof_chn[3, 1, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for sel coeff of mixed pre-event")
lines(density(sel_cof_chn[3, 1, ]), lwd = 2, col = 'black')
abline(v = sel_cof[3, 1], col = 'red', lty = 2, lwd = 2)
abline(v = sel_cof_est[3, 1], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[3, 1, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[3, 2, 1], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_chn[1, 2, ], breaks = seq(min(sel_cof_chn[1, 2, ]), max(sel_cof_chn[1, 2, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for sel coeff of tobiano post-event")
lines(density(sel_cof_chn[1, 2, ]), lwd = 2, col = 'black')
abline(v = sel_cof[1, 2], col = 'red', lty = 2, lwd = 2)
abline(v = sel_cof_est[1, 2], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 1, 2], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 2, 2], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_chn[2, 2, ], breaks = seq(min(sel_cof_chn[2, 2, ]), max(sel_cof_chn[2, 2, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for sel coeff of sabino post-event")
lines(density(sel_cof_chn[2, 2, ]), lwd = 2, col = 'black')
abline(v = sel_cof[2, 2], col = 'red', lty = 2, lwd = 2)
abline(v = sel_cof_est[2, 2], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 1, 2], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 2, 2], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_chn[3, 2, ], breaks = seq(min(sel_cof_chn[3, 2, ]), max(sel_cof_chn[3, 2, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for sel coeff of mixed post-event")
lines(density(sel_cof_chn[3, 2, ]), lwd = 2, col = 'black')
abline(v = sel_cof[3, 2], col = 'red', lty = 2, lwd = 2)
abline(v = sel_cof_est[3, 2], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[3, 1, 2], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[3, 2, 2], col = 'blue', lty = 2, lwd = 2)
dev.off()

dif_sel_chn <- sel_cof_chn[, 2, ] - sel_cof_chn[, 1, ]

dif_sel_est <- rowMeans(dif_sel_chn)

dif_sel_hpd <- matrix(NA, nrow = 3, ncol = 2)
dif_sel_hpd[1, ] <- HPDinterval(as.mcmc(dif_sel_chn[1, ]), prob = 0.95)
dif_sel_hpd[2, ] <- HPDinterval(as.mcmc(dif_sel_chn[2, ]), prob = 0.95)
dif_sel_hpd[3, ] <- HPDinterval(as.mcmc(dif_sel_chn[3, ]), prob = 0.95)

pdf(file = "./Output/Output v1.0/Test v1.2/TEST_PTN_AdaptPMMH_Posterior_SelChange.pdf", width = 8, height = 18)
par(mfrow = c(3, 1), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(dif_sel_chn[1, ], breaks = seq(min(dif_sel_chn[1, ]), max(dif_sel_chn[1, ]), length.out = 50), freq = FALSE,
     xlab = "Change in the selection coefficient",
     main = "Posterior for change in sel coeff of tobiano")
lines(density(dif_sel_chn[1, ]), lwd = 2, col = 'black')
abline(v = 0, col = 'red', lty = 2, lwd = 2)
abline(v = dif_sel_est[1], col = 'black', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[1, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[1, 2], col = 'blue', lty = 2, lwd = 2)

hist(dif_sel_chn[2, ], breaks = seq(min(dif_sel_chn[2, ]), max(dif_sel_chn[2, ]), length.out = 50), freq = FALSE,
     xlab = "Change in the selection coefficient",
     main = "Posterior for change in sel coeff of sabino")
lines(density(dif_sel_chn[2, ]), lwd = 2, col = 'black')
abline(v = 0, col = 'red', lty = 2, lwd = 2)
abline(v = dif_sel_est[2], col = 'black', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[2, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[2, 2], col = 'blue', lty = 2, lwd = 2)

hist(dif_sel_chn[3, ], breaks = seq(min(dif_sel_chn[3, ]), max(dif_sel_chn[3, ]), length.out = 50), freq = FALSE,
     xlab = "Change in the selection coefficient",
     main = "Posterior for change in sel coeff of mixed")
lines(density(dif_sel_chn[3, ]), lwd = 2, col = 'black')
abline(v = 0, col = 'red', lty = 2, lwd = 2)
abline(v = dif_sel_est[3], col = 'black', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[3, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[3, 2], col = 'blue', lty = 2, lwd = 2)
dev.off()

########################################

#' Run the Bayesian procedure for the inference of natural selection
#' Parameter settings
#' @param sel_cof the selection coefficients of the tobiano, sabino and mixed phenotypes
#' @param rec_rat the recombination rate between the KIT13 and KIT16 loci
#' @param pop_siz the size of the horse population (non-constant)
#' @param ref_siz the reference size of the horse population
#' @param evt_gen the generation that the event of interest occurred
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_siz the count of the chromosomes drawn from the population at all sampling time points
#' @param smp_cnt the count of the mutant alleles observed in the sample at all sampling time points
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter
#' @param itn_num the number of the iterations carried out in the particle marginal Metropolis-Hastings
#' @param brn_num the number of the iterations for burn-in
#' @param thn_num the number of the iterations for thinning
#' @param adp_set = TRUE/FALSE (return the result with the adaptive setting or not)
#' @param stp_siz the step size sequence in the adaptive setting (decaying to zero)
#' @param apt_rto the target mean acceptance probability of the adaptive setting

load("./Output/Output v1.0/Test v1.2/TEST_PTN_SimData.rda")

set.seed(test_seed)

sel_cof <- matrix(c(0e+00, 0e+00, 0e+00, 0e+00, 0e+00, 0e+00), nrow = 3, ncol = 2)
rec_rat
pop_siz
ref_siz <- 1e+04
evt_gen
smp_gen
smp_siz
smp_cnt
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04
brn_num <- 5e+03
thn_num <- 3e+00
adp_set <- TRUE
stp_siz <- (1:itn_num)^(-2 / 3)
apt_rto <- 4e-01

system.time(BayesianProcedure <- cmprunBayesianProcedure(sel_cof, rec_rat, pop_siz, ref_siz, evt_gen, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, brn_num, thn_num, adp_set, stp_siz, apt_rto))

load("./Output/Output v1.0/Test v1.2/TEST_PTN_SimData.rda")

save(sel_cof, rec_rat, pop_siz, ref_siz, evt_gen, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, brn_num, thn_num, adp_set, stp_siz, apt_rto, BayesianProcedure,
     file = "./Output/Output v1.0/Test v1.2/TEST_PTN_BayesProc.rda")

load("./Output/Output v1.0/Test v1.2/TEST_PTN_BayesProc.rda")

sel_cof_chn <- BayesianProcedure$sel_cof_chn

sel_cof_est <- BayesianProcedure$sel_cof_est

sel_cof_hpd <- BayesianProcedure$sel_cof_hpd

pdf(file = "./Output/Output v1.0/Test v1.2/TEST_PTN_BayesProc_Posterior_SelCoeff.pdf", width = 24, height = 12)
par(mfrow = c(2, 3), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(sel_cof_chn[1, 1, ], breaks = seq(min(sel_cof_chn[1, 1, ]), max(sel_cof_chn[1, 1, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for sel coeff of tobiano pre-event")
lines(density(sel_cof_chn[1, 1, ]), lwd = 2, col = 'black')
abline(v = sel_cof[1, 1], col = 'red', lty = 2, lwd = 2)
abline(v = sel_cof_est[1, 1], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 1, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 2, 1], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_chn[2, 1, ], breaks = seq(min(sel_cof_chn[2, 1, ]), max(sel_cof_chn[2, 1, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for sel coeff of sabino pre-event")
lines(density(sel_cof_chn[2, 1, ]), lwd = 2, col = 'black')
abline(v = sel_cof[2, 1], col = 'red', lty = 2, lwd = 2)
abline(v = sel_cof_est[2, 1], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 1, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 2, 1], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_chn[3, 1, ], breaks = seq(min(sel_cof_chn[3, 1, ]), max(sel_cof_chn[3, 1, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for sel coeff of mixed pre-event")
lines(density(sel_cof_chn[3, 1, ]), lwd = 2, col = 'black')
abline(v = sel_cof[3, 1], col = 'red', lty = 2, lwd = 2)
abline(v = sel_cof_est[3, 1], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[3, 1, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[3, 2, 1], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_chn[1, 2, ], breaks = seq(min(sel_cof_chn[1, 2, ]), max(sel_cof_chn[1, 2, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for sel coeff of tobiano post-event")
lines(density(sel_cof_chn[1, 2, ]), lwd = 2, col = 'black')
abline(v = sel_cof[1, 2], col = 'red', lty = 2, lwd = 2)
abline(v = sel_cof_est[1, 2], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 1, 2], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 2, 2], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_chn[2, 2, ], breaks = seq(min(sel_cof_chn[2, 2, ]), max(sel_cof_chn[2, 2, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for sel coeff of sabino post-event")
lines(density(sel_cof_chn[2, 2, ]), lwd = 2, col = 'black')
abline(v = sel_cof[2, 2], col = 'red', lty = 2, lwd = 2)
abline(v = sel_cof_est[2, 2], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 1, 2], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 2, 2], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_chn[3, 2, ], breaks = seq(min(sel_cof_chn[3, 2, ]), max(sel_cof_chn[3, 2, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for sel coeff of mixed post-event")
lines(density(sel_cof_chn[3, 2, ]), lwd = 2, col = 'black')
abline(v = sel_cof[3, 2], col = 'red', lty = 2, lwd = 2)
abline(v = sel_cof_est[3, 2], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[3, 1, 2], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[3, 2, 2], col = 'blue', lty = 2, lwd = 2)
dev.off()

dif_sel_chn <- BayesianProcedure$dif_sel_chn

dif_sel_est <- BayesianProcedure$dif_sel_est

dif_sel_hpd <- BayesianProcedure$dif_sel_hpd

pdf(file = "./Output/Output v1.0/Test v1.2/TEST_PTN_BayesProc_Posterior_SelChange.pdf", width = 8, height = 18)
par(mfrow = c(3, 1), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(dif_sel_chn[1, ], breaks = seq(min(dif_sel_chn[1, ]), max(dif_sel_chn[1, ]), length.out = 50), freq = FALSE,
     xlab = "Change in the selection coefficient",
     main = "Posterior for change in sel coeff of tobiano")
lines(density(dif_sel_chn[1, ]), lwd = 2, col = 'black')
abline(v = 0, col = 'red', lty = 2, lwd = 2)
abline(v = dif_sel_est[1], col = 'black', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[1, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[1, 2], col = 'blue', lty = 2, lwd = 2)

hist(dif_sel_chn[2, ], breaks = seq(min(dif_sel_chn[2, ]), max(dif_sel_chn[2, ]), length.out = 50), freq = FALSE,
     xlab = "Change in the selection coefficient",
     main = "Posterior for change in sel coeff of sabino")
lines(density(dif_sel_chn[2, ]), lwd = 2, col = 'black')
abline(v = 0, col = 'red', lty = 2, lwd = 2)
abline(v = dif_sel_est[2], col = 'black', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[2, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[2, 2], col = 'blue', lty = 2, lwd = 2)

hist(dif_sel_chn[3, ], breaks = seq(min(dif_sel_chn[3, ]), max(dif_sel_chn[3, ]), length.out = 50), freq = FALSE,
     xlab = "Change in the selection coefficient",
     main = "Posterior for change in sel coeff of mixed")
lines(density(dif_sel_chn[3, ]), lwd = 2, col = 'black')
abline(v = 0, col = 'red', lty = 2, lwd = 2)
abline(v = dif_sel_est[3], col = 'black', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[3, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[3, 2], col = 'blue', lty = 2, lwd = 2)
dev.off()

################################################################################
