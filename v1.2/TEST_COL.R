#' @title Inferring natural selection acting on horse coat colours and patterns during the process of domestication from ancient DNA data
#' @author Xiaoyang Dai, Sile Hu, Mark Beaumont, Feng Yu, Zhangyi He

#' version 1.2

#' Horse coat colours (ASIP & MC1R) under non-constant natural selection and non-constant demographic histories (N/A is not allowed)

#install.packages("RColorBrewer")
library("RColorBrewer")

#install.packages("viridis")
library("viridis")

#install.packages("ggplot2")
library("ggplot2")

#install.packages("plot3D")
library("plot3D")

# set the directory
setwd("~/Dropbox/Jeffery He/iResearch/Publications/2018/HE2021-WFM-2L-DiffusApprox-PMMH-Horse-MolEcol")

# call R functions
source("./Code/Code v1.2/RFUN_COL.R")

################################################################################

#' Simulate the haplotype frequency trajectories according to the two-locus Wright-Fisher model with selection
#' Parameter setting
#' @param sel_cof the selection coefficients of the black and chestnut phenotypes
#' @param rec_rat the recombination rate between the ASIP and MC1R loci
#' @param pop_siz the size of the horse population (non-constant)
#' @param int_frq the initial haplotype frequencies of the population
#' @param evt_gen the generation that the event of interest occurred
#' @param int_gen the generation that the simulated haplotype frequency trajectories started
#' @param lst_gen the generation that the simulated haplotype frequency trajectories ended

sel_cof <- matrix(c(0e+00, 5e-03, 1e-02, 5e-03), nrow = 2, ncol = 2)
rec_rat <- 5e-01
pop_siz <- c(rep(1e+04, length.out = 201), rep(5e+03, length.out = 200), rep(1e+04, length.out = 100))
int_frq <- c(6e-01, 2e-01, 1e-01, 1e-01)
evt_gen <- 210
int_gen <- 0
lst_gen <- 500

frq_pth <- cmpsimulateWFM(sel_cof, rec_rat, pop_siz, int_frq, evt_gen, int_gen, lst_gen)$hap_frq_pth

k <- int_gen:lst_gen
plot(k, frq_pth[1, ], type = "l", lwd = 1.5,
     xlab = "Generation", ylab = "Haplotype frequency",
     main = "WFM: the AE haplotype frequency trajectory")
plot(k, frq_pth[2, ], type = "l", lwd = 1.5,
     xlab = "Generation", ylab = "Haplotype frequency",
     main = "WFM: the Ae haplotype frequency trajectory")
plot(k, frq_pth[3, ], type = "l", lwd = 1.5,
     xlab = "Generation", ylab = "Haplotype frequency",
     main = "WFM: the aE haplotype frequency trajectory")
plot(k, frq_pth[4, ], type = "l", lwd = 1.5,
     xlab = "Generation", ylab = "Haplotype frequency",
     main = "WFM: the ae haplotype frequency trajectory")

########################################

#' Simulate the haplotype frequency trajectories according to the two-locus Wright-Fisher diffusion with selection using the Euler-Maruyama method
#' Parameter setting
#' @param sel_cof the selection coefficients of the black and chestnut phenotypes
#' @param rec_rat the recombination rate between the ASIP and MC1R loci
#' @param pop_siz the size of the horse population (non-constant)
#' @param ref_siz the reference size of the horse population
#' @param int_frq the initial haplotype frequencies of the population
#' @param evt_gen the generation that the event of interest occurred
#' @param int_gen the generation that the simulated haplotype frequency trajectories started
#' @param lst_gen the generation that the simulated haplotype frequency trajectories ended
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param dat_aug = TRUE/FALSE (return the simulated sample trajectory with data augmentation or not)

sel_cof <- matrix(c(0e+00, 5e-03, 1e-02, 5e-03), nrow = 2, ncol = 2)
rec_rat <- 5e-01
pop_siz <- c(rep(1e+04, length.out = 201), rep(5e+03, length.out = 200), rep(1e+04, length.out = 100))
ref_siz <- 1e+04
int_frq <- c(6e-01, 2e-01, 1e-01, 1e-01)
evt_gen <- 210
int_gen <- 0
lst_gen <- 500
ptn_num <- 5e+00

frq_pth <- cmpsimulateWFD(sel_cof, rec_rat, pop_siz, ref_siz, int_frq, evt_gen, int_gen, lst_gen, ptn_num, dat_aug = TRUE)

t <- (int_gen:(int_gen + (lst_gen - int_gen) * ptn_num)) / 2 / ref_siz
plot(t, frq_pth[1, ], type = "l", lwd = 1.5,
     xlab = "Time", ylab = "Haplotype frequency",
     main = "WFD: the AE haplotype frequency trajectory")
plot(t, frq_pth[2, ], type = "l", lwd = 1.5,
     xlab = "Time", ylab = "Haplotype frequency",
     main = "WFD: the Ae haplotype frequency trajectory")
plot(t, frq_pth[3, ], type = "l", lwd = 1.5,
     xlab = "Time", ylab = "Haplotype frequency",
     main = "WFD: the aE haplotype frequency trajectory")
plot(t, frq_pth[4, ], type = "l", lwd = 1.5,
     xlab = "Time", ylab = "Haplotype frequency",
     main = "WFD: the ae haplotype frequency trajectory")

########################################

#' Compare the simulation generated with the Wright-Fisher model and the Wright-Fisher diffusion
sel_cof <- matrix(c(0e+00, 5e-03, 1e-02, 5e-03), nrow = 2, ncol = 2)
rec_rat <- 5e-01
pop_siz <- c(rep(1e+04, length.out = 201), rep(5e+03, length.out = 200), rep(1e+04, length.out = 100))
ref_siz <- 1e+04
int_frq <- c(6e-01, 2e-01, 1e-01, 1e-01)
evt_gen <- 210
int_gen <- 0
lst_gen <- 500
ptn_num <- 5e+00
sim_num <- 1e+06

smp_WFM <- matrix(NA, nrow = 4, ncol = sim_num)
smp_WFD <- matrix(NA, nrow = 4, ncol = sim_num)
for (i in 1:sim_num) {
  print(i)
  smp_WFM[, i] <- cmpsimulateWFM(sel_cof, rec_rat, pop_siz, int_frq, evt_gen, int_gen, lst_gen)$hap_frq_pth[, (lst_gen - int_gen) + 1]
  smp_WFD[, i] <- cmpsimulateWFD(sel_cof, rec_rat, pop_siz, ref_siz, int_frq, evt_gen, int_gen, lst_gen, ptn_num, dat_aug = FALSE)[, (lst_gen - int_gen) + 1]
}

hist(smp_WFM[1, ], breaks = seq(min(smp_WFM[1, ], smp_WFD[1, ]), max(smp_WFM[1, ], smp_WFD[1, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(smp_WFM[1, ], smp_WFD[1, ]), max(smp_WFM[1, ], smp_WFD[1, ])),
     xlab = "Haplotype frequency", main = paste("Histogram of the AE haplotype at generation", lst_gen))
hist(smp_WFD[1, ], breaks = seq(min(smp_WFM[1, ], smp_WFD[1, ]), max(smp_WFM[1, ], smp_WFD[1, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)

hist(smp_WFM[2, ], breaks = seq(min(smp_WFM[2, ], smp_WFD[2, ]), max(smp_WFM[2, ], smp_WFD[2, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(smp_WFM[2, ], smp_WFD[2, ]), max(smp_WFM[2, ], smp_WFD[2, ])),
     xlab = "Haplotype frequency", main = paste("Histogram of the Ae haplotype at generation", lst_gen))
hist(smp_WFD[2, ], breaks = seq(min(smp_WFM[2, ], smp_WFD[2, ]), max(smp_WFM[2, ], smp_WFD[2, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)

hist(smp_WFM[3, ], breaks = seq(min(smp_WFM[3, ], smp_WFD[3, ]), max(smp_WFM[3, ], smp_WFD[3, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(smp_WFM[3, ], smp_WFD[3, ]), max(smp_WFM[3, ], smp_WFD[3, ])),
     xlab = "Haplotype frequency", main = paste("Histogram of the aE haplotype at generation", lst_gen))
hist(smp_WFD[3, ], breaks = seq(min(smp_WFM[3, ], smp_WFD[3, ]), max(smp_WFM[3, ], smp_WFD[3, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)

hist(smp_WFM[4, ], breaks = seq(min(smp_WFM[4, ], smp_WFD[4, ]), max(smp_WFM[4, ], smp_WFD[4, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(smp_WFM[4, ], smp_WFD[4, ]), max(smp_WFM[4, ], smp_WFD[4, ])),
     xlab = "Haplotype frequency", main = paste("Histogram of the ae haplotype at generation", lst_gen))
hist(smp_WFD[4, ], breaks = seq(min(smp_WFM[4, ], smp_WFD[4, ]), max(smp_WFM[4, ], smp_WFD[4, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)

################################################################################

#' Simulate the hidden Markov model
#' Parameter setting
#' @param model = "WFM"/"WFD" (return the observations from the underlying population evolving according to the WFM or the WFD)
#' @param sel_cof the selection coefficients of the black and chestnut phenotypes
#' @param rec_rat the recombination rate between the ASIP and MC1R loci
#' @param pop_siz the size of the horse population (non-constant)
#' @param int_frq the initial haplotype frequencies of the population
#' @param evt_gen the generation that the event of interest occurred
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_siz the count of the horses drawn from the population at all sampling time points
#' @param ref_siz the reference size of the horse population
#' @param ptn_num the number of the subintervals divided per generation in the Euler-Maruyama method for the WFD

#' Simulate the dataset under the Wright-Fisher model
model <- "WFM"
sel_cof <- matrix(c(0e+00, 5e-03, 1e-02, 5e-03), nrow = 2, ncol = 2)
rec_rat <- 5e-01
pop_siz <- c(rep(1e+04, length.out = 201), rep(5e+03, length.out = 200), rep(1e+04, length.out = 100))
int_frq <- c(6e-01, 2e-01, 1e-01, 1e-01)
evt_gen <- 210
smp_gen <- (0:10) * 50
smp_siz <- rep(100, 11)

sim_HMM_WFM <- cmpsimulateHMM(model, sel_cof, rec_rat, pop_siz, int_frq, evt_gen, smp_gen, smp_siz)
smp_gen <- sim_HMM_WFM$smp_gen
smp_siz <- sim_HMM_WFM$smp_siz
smp_gen_cnt <- sim_HMM_WFM$smp_gen_cnt
smp_gen_frq <- sim_HMM_WFM$smp_gen_frq
pop_gen_frq <- sim_HMM_WFM$pop_gen_frq
pop_hap_frq <- sim_HMM_WFM$pop_hap_frq

pop_gen_frq[5, ] <- pop_gen_frq[5, ] + pop_gen_frq[7, ]
pop_gen_frq <- pop_gen_frq[-7, ]

k <- min(smp_gen):max(smp_gen)
plot(k, pop_gen_frq[1, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_gen_frq[1, ], pop_gen_frq[1, ]), max(smp_gen_frq[1, ], pop_gen_frq[1, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "WFM-HMM: the AA/EE genotype")
points(smp_gen, smp_gen_frq[1, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_gen_frq[2, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_gen_frq[2, ], pop_gen_frq[2, ]), max(smp_gen_frq[2, ], pop_gen_frq[2, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "WFM-HMM: the AA/Ee genotype")
points(smp_gen, smp_gen_frq[2, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_gen_frq[3, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_gen_frq[3, ], pop_gen_frq[3, ]), max(smp_gen_frq[3, ], pop_gen_frq[3, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "WFM-HMM: the AA/ee genotype")
points(smp_gen, smp_gen_frq[3, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_gen_frq[4, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_gen_frq[4, ], pop_gen_frq[4, ]), max(smp_gen_frq[4, ], pop_gen_frq[4, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "WFM-HMM: the Aa/EE genotype")
points(smp_gen, smp_gen_frq[4, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_gen_frq[5, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_gen_frq[5, ], pop_gen_frq[5, ]), max(smp_gen_frq[5, ], pop_gen_frq[5, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "WFM-HMM: the Aa/Ee genotype")
points(smp_gen, smp_gen_frq[5, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_gen_frq[6, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_gen_frq[6, ], pop_gen_frq[6, ]), max(smp_gen_frq[6, ], pop_gen_frq[6, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "WFM-HMM: the aa/EE genotype")
points(smp_gen, smp_gen_frq[6, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_gen_frq[7, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_gen_frq[7, ], pop_gen_frq[7, ]), max(smp_gen_frq[7, ], pop_gen_frq[7, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "WFM-HMM: the Aa/ee genotype")
points(smp_gen, smp_gen_frq[7, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_gen_frq[8, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_gen_frq[8, ], pop_gen_frq[8, ]), max(smp_gen_frq[8, ], pop_gen_frq[8, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "WFM-HMM: the aa/Ee genotype")
points(smp_gen, smp_gen_frq[8, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_gen_frq[9, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_gen_frq[9, ], pop_gen_frq[9, ]), max(smp_gen_frq[9, ], pop_gen_frq[9, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "WFM-HMM: the aa/ee genotype")
points(smp_gen, smp_gen_frq[9, ], col = 'red', pch = 17, cex = 1)

####################

#' Simulate the dataset under the Wright-Fisher diffusion
model <- "WFD"
sel_cof <- matrix(c(0e+00, 5e-03, 1e-02, 5e-03), nrow = 2, ncol = 2)
rec_rat <- 5e-01
pop_siz <- c(rep(1e+04, length.out = 201), rep(5e+03, length.out = 200), rep(1e+04, length.out = 100))
int_frq <- c(6e-01, 2e-01, 1e-01, 1e-01)
evt_gen <- 210
smp_gen <- (0:10) * 50
smp_siz <- rep(100, 11)
ref_siz <- 1e+04
ptn_num <- 5e+00

sim_HMM_WFD <- cmpsimulateHMM(model, sel_cof, rec_rat, pop_siz, int_frq, evt_gen, smp_gen, smp_siz, ref_siz, ptn_num)
smp_gen <- sim_HMM_WFD$smp_gen
smp_siz <- sim_HMM_WFD$smp_siz
smp_gen_cnt <- sim_HMM_WFD$smp_gen_cnt
smp_gen_frq <- sim_HMM_WFD$smp_gen_frq
pop_gen_frq <- sim_HMM_WFD$pop_gen_frq
pop_hap_frq <- sim_HMM_WFD$pop_hap_frq

pop_gen_frq[5, ] <- pop_gen_frq[5, ] + pop_gen_frq[7, ]
pop_gen_frq <- pop_gen_frq[-7, ]

k <- min(smp_gen):max(smp_gen)
plot(k, pop_gen_frq[1, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_gen_frq[1, ], pop_gen_frq[1, ]), max(smp_gen_frq[1, ], pop_gen_frq[1, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "WFD-HMM: the AA/EE genotype")
points(smp_gen, smp_gen_frq[1, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_gen_frq[2, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_gen_frq[2, ], pop_gen_frq[2, ]), max(smp_gen_frq[2, ], pop_gen_frq[2, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "WFD-HMM: the AA/Ee genotype")
points(smp_gen, smp_gen_frq[2, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_gen_frq[3, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_gen_frq[3, ], pop_gen_frq[3, ]), max(smp_gen_frq[3, ], pop_gen_frq[3, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "WFD-HMM: the AA/ee genotype")
points(smp_gen, smp_gen_frq[3, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_gen_frq[4, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_gen_frq[4, ], pop_gen_frq[4, ]), max(smp_gen_frq[4, ], pop_gen_frq[4, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "WFD-HMM: the Aa/EE genotype")
points(smp_gen, smp_gen_frq[4, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_gen_frq[5, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_gen_frq[5, ], pop_gen_frq[5, ]), max(smp_gen_frq[5, ], pop_gen_frq[5, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "WFD-HMM: the Aa/Ee genotype")
points(smp_gen, smp_gen_frq[5, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_gen_frq[6, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_gen_frq[6, ], pop_gen_frq[6, ]), max(smp_gen_frq[6, ], pop_gen_frq[6, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "WFD-HMM: the aa/EE genotype")
points(smp_gen, smp_gen_frq[6, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_gen_frq[7, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_gen_frq[7, ], pop_gen_frq[7, ]), max(smp_gen_frq[7, ], pop_gen_frq[7, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "WFD-HMM: the Aa/ee genotype")
points(smp_gen, smp_gen_frq[7, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_gen_frq[8, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_gen_frq[8, ], pop_gen_frq[8, ]), max(smp_gen_frq[8, ], pop_gen_frq[8, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "WFD-HMM: the aa/Ee genotype")
points(smp_gen, smp_gen_frq[8, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_gen_frq[9, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_gen_frq[9, ], pop_gen_frq[9, ]), max(smp_gen_frq[9, ], pop_gen_frq[9, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "WFD-HMM: the aa/ee genotype")
points(smp_gen, smp_gen_frq[9, ], col = 'red', pch = 17, cex = 1)

################################################################################

#' Generate a simulated dataset under the Wright-Fisher model
test_seed <- 27
set.seed(test_seed)

model <- "WFM"
sel_cof <- matrix(c(0e+00, 5e-03, 1e-02, 5e-03), nrow = 2, ncol = 2)
rec_rat <- 5e-01
pop_siz <- c(rep(1e+04, length.out = 201), rep(5e+03, length.out = 200), rep(1e+04, length.out = 100))
int_frq <- c(6e-01, 2e-01, 1e-01, 1e-01)
evt_gen <- 210
smp_gen <- (0:10) * 50
smp_siz <- rep(100, 11)

sim_HMM_WFM <- cmpsimulateHMM(model, sel_cof, rec_rat, pop_siz, int_frq, evt_gen, smp_gen, smp_siz)
smp_gen <- sim_HMM_WFM$smp_gen
smp_siz <- sim_HMM_WFM$smp_siz
smp_cnt <- sim_HMM_WFM$smp_gen_cnt
smp_gen_frq <- sim_HMM_WFM$smp_gen_frq
pop_gen_frq <- sim_HMM_WFM$pop_gen_frq
pop_hap_frq <- sim_HMM_WFM$pop_hap_frq

save(sel_cof, rec_rat, pop_siz, int_frq, evt_gen, smp_gen, smp_siz, smp_cnt, smp_gen_frq, pop_hap_frq, pop_gen_frq, pop_hap_frq,
     file = "./Output/Output v1.2/TEST_SimData.rda")

load("./Output/Output v1.2/TEST_SimData.rda")

pop_gen_frq[5, ] <- pop_gen_frq[5, ] + pop_gen_frq[7, ]
pop_gen_frq <- pop_gen_frq[-7, ]

pdf(file = "./Output/Output v1.2/TEST_SimData.pdf", width = 16, height = 12)
par(mfrow = c(3, 3), oma = c(0, 0, 3, 0), mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
k <- min(smp_gen):max(smp_gen)
plot(k, pop_gen_frq[1, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_gen_frq[1, ], pop_gen_frq[1, ]), max(smp_gen_frq[1, ], pop_gen_frq[1, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "Genotype AA/EE")
points(smp_gen, smp_gen_frq[1, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_gen_frq[2, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_gen_frq[2, ], pop_gen_frq[2, ]), max(smp_gen_frq[2, ], pop_gen_frq[2, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "Genotype AA/Ee")
points(smp_gen, smp_gen_frq[2, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_gen_frq[3, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_gen_frq[3, ], pop_gen_frq[3, ]), max(smp_gen_frq[3, ], pop_gen_frq[3, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "Genotype AA/ee")
points(smp_gen, smp_gen_frq[3, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_gen_frq[4, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_gen_frq[4, ], pop_gen_frq[4, ]), max(smp_gen_frq[4, ], pop_gen_frq[4, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "Genotype Aa/EE")
points(smp_gen, smp_gen_frq[4, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_gen_frq[5, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_gen_frq[5, ], pop_gen_frq[5, ]), max(smp_gen_frq[5, ], pop_gen_frq[5, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "Genotype Aa/Ee")
points(smp_gen, smp_gen_frq[5, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_gen_frq[6, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_gen_frq[6, ], pop_gen_frq[6, ]), max(smp_gen_frq[6, ], pop_gen_frq[6, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "Genotype aa/EE")
points(smp_gen, smp_gen_frq[6, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_gen_frq[7, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_gen_frq[7, ], pop_gen_frq[7, ]), max(smp_gen_frq[7, ], pop_gen_frq[7, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "Genotype Aa/ee")
points(smp_gen, smp_gen_frq[7, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_gen_frq[8, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_gen_frq[8, ], pop_gen_frq[8, ]), max(smp_gen_frq[8, ], pop_gen_frq[8, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "Genotype aa/Ee")
points(smp_gen, smp_gen_frq[8, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_gen_frq[9, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_gen_frq[9, ], pop_gen_frq[9, ]), max(smp_gen_frq[9, ], pop_gen_frq[9, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "Genotype aa/ee")
points(smp_gen, smp_gen_frq[9, ], col = 'red', pch = 17, cex = 1)
dev.off()

########################################

#' Run the bootstrap particle filter (BPF) with the two-locus Wright-Fisher diffusion with selection
#' Parameter setting
#' @param sel_cof the selection coefficients of the black and chestnut phenotypes
#' @param rec_rat the recombination rate between the ASIP and MC1R loci
#' @param pop_siz the size of the horse population (non-constant)
#' @param ref_siz the reference size of the horse population
#' @param evt_gen the generation that the event of interest occurred
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_siz the count of the horses drawn from the population at all sampling time points
#' @param smp_cnt the count of the genotypes observed in the sample at all sampling time points
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter

load("./Output/Output v1.2/TEST_SimData.rda")

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
     file = "./Output/Output v1.2/TEST_BPF.rda")

load("./Output/Output v1.2/TEST_BPF.rda")

lik <- rep(1, pcl_num)
wght <- BPF$wght
for (k in 1:length(smp_gen)) {
  lik <- lik * (cumsum(wght[, k]) / (1:pcl_num))
}

pdf(file = "./Output/Output v1.2/TEST_BPF_Likelihood.pdf", width = 12, height = 9)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:pcl_num, log(lik), type = 'l',
     xlab = "Number of particles", ylab = "Log likelihood",
     main = "Log likelihood estimated with the bootstrap particle filter")
dev.off()

pop_gen_frq_pre_resmp <- BPF$gen_frq_pre_resmp
pop_gen_frq_pst_resmp <- BPF$gen_frq_pst_resmp

pdf(file = "./Output/Output v1.2/TEST_BPF_Particle.pdf", width = 72, height = 66)
par(mfrow = c(11, 9), oma = c(0, 0, 3, 0), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
for (k in 1:length(smp_gen)) {
  hist_pst_resmp <- hist(pop_gen_frq_pst_resmp[1, , k], breaks = seq(min(pop_gen_frq_pst_resmp[1, , k], pop_gen_frq_pre_resmp[1, , k]), max(pop_gen_frq_pst_resmp[1, , k], pop_gen_frq_pre_resmp[1, , k]), length.out = 50), plot = FALSE)
  hist_pre_resmp <- hist(pop_gen_frq_pre_resmp[1, , k], breaks = seq(min(pop_gen_frq_pst_resmp[1, , k], pop_gen_frq_pre_resmp[1, , k]), max(pop_gen_frq_pst_resmp[1, , k], pop_gen_frq_pre_resmp[1, , k]), length.out = 50), plot = FALSE)
  hist(pop_gen_frq_pst_resmp[1, , k], breaks = seq(min(pop_gen_frq_pst_resmp[1, , k], pop_gen_frq_pre_resmp[1, , k]), max(pop_gen_frq_pst_resmp[1, , k], pop_gen_frq_pre_resmp[1, , k]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
       xlim = c(min(pop_gen_frq_pst_resmp[1, , k], pop_gen_frq_pre_resmp[1, , k], smp_gen_frq[1, k]), max(pop_gen_frq_pst_resmp[1, , k], pop_gen_frq_pre_resmp[1, , k], smp_gen_frq[1, k])), ylim = c(0, max(hist_pst_resmp$density, hist_pre_resmp$density)),
       xlab = "Genotype frequency",
       main = paste("Genotype AA/EE at generation", smp_gen[k]))
  hist(pop_gen_frq_pre_resmp[1, , k], breaks = seq(min(pop_gen_frq_pst_resmp[1, , k], pop_gen_frq_pre_resmp[1, , k]), max(pop_gen_frq_pst_resmp[1, , k], pop_gen_frq_pre_resmp[1, , k]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
  abline(v = smp_gen_frq[1, k], col = 'red', lty = 2, lwd = 2)

  hist_pst_resmp <- hist(pop_gen_frq_pst_resmp[2, , k], breaks = seq(min(pop_gen_frq_pst_resmp[2, , k], pop_gen_frq_pre_resmp[2, , k]), max(pop_gen_frq_pst_resmp[2, , k], pop_gen_frq_pre_resmp[2, , k]), length.out = 50), plot = FALSE)
  hist_pre_resmp <- hist(pop_gen_frq_pre_resmp[2, , k], breaks = seq(min(pop_gen_frq_pst_resmp[2, , k], pop_gen_frq_pre_resmp[2, , k]), max(pop_gen_frq_pst_resmp[2, , k], pop_gen_frq_pre_resmp[2, , k]), length.out = 50), plot = FALSE)
  hist(pop_gen_frq_pst_resmp[2, , k], breaks = seq(min(pop_gen_frq_pst_resmp[2, , k], pop_gen_frq_pre_resmp[2, , k]), max(pop_gen_frq_pst_resmp[2, , k], pop_gen_frq_pre_resmp[2, , k]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
       xlim = c(min(pop_gen_frq_pst_resmp[2, , k], pop_gen_frq_pre_resmp[2, , k], smp_gen_frq[2, k]), max(pop_gen_frq_pst_resmp[2, , k], pop_gen_frq_pre_resmp[2, , k], smp_gen_frq[2, k])), ylim = c(0, max(hist_pst_resmp$density, hist_pre_resmp$density)),
       xlab = "Genotype frequency",
       main = paste("Genotype AA/Ee at generation", smp_gen[k]))
  hist(pop_gen_frq_pre_resmp[2, , k], breaks = seq(min(pop_gen_frq_pst_resmp[2, , k], pop_gen_frq_pre_resmp[2, , k]), max(pop_gen_frq_pst_resmp[2, , k], pop_gen_frq_pre_resmp[2, , k]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
  abline(v = smp_gen_frq[2, k], col = 'red', lty = 2, lwd = 2)

  hist_pst_resmp <- hist(pop_gen_frq_pst_resmp[3, , k], breaks = seq(min(pop_gen_frq_pst_resmp[3, , k], pop_gen_frq_pre_resmp[3, , k]), max(pop_gen_frq_pst_resmp[3, , k], pop_gen_frq_pre_resmp[3, , k]), length.out = 50), plot = FALSE)
  hist_pre_resmp <- hist(pop_gen_frq_pre_resmp[3, , k], breaks = seq(min(pop_gen_frq_pst_resmp[3, , k], pop_gen_frq_pre_resmp[3, , k]), max(pop_gen_frq_pst_resmp[3, , k], pop_gen_frq_pre_resmp[3, , k]), length.out = 50), plot = FALSE)
  hist(pop_gen_frq_pst_resmp[3, , k], breaks = seq(min(pop_gen_frq_pst_resmp[3, , k], pop_gen_frq_pre_resmp[3, , k]), max(pop_gen_frq_pst_resmp[3, , k], pop_gen_frq_pre_resmp[3, , k]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
       xlim = c(min(pop_gen_frq_pst_resmp[3, , k], pop_gen_frq_pre_resmp[3, , k], smp_gen_frq[3, k]), max(pop_gen_frq_pst_resmp[3, , k], pop_gen_frq_pre_resmp[3, , k], smp_gen_frq[3, k])), ylim = c(0, max(hist_pst_resmp$density, hist_pre_resmp$density)),
       xlab = "Genotype frequency",
       main = paste("Genotype AA/ee at generation", smp_gen[k]))
  hist(pop_gen_frq_pre_resmp[3, , k], breaks = seq(min(pop_gen_frq_pst_resmp[3, , k], pop_gen_frq_pre_resmp[3, , k]), max(pop_gen_frq_pst_resmp[3, , k], pop_gen_frq_pre_resmp[3, , k]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
  abline(v = smp_gen_frq[3, k], col = 'red', lty = 2, lwd = 2)

  hist_pst_resmp <- hist(pop_gen_frq_pst_resmp[4, , k], breaks = seq(min(pop_gen_frq_pst_resmp[4, , k], pop_gen_frq_pre_resmp[4, , k]), max(pop_gen_frq_pst_resmp[4, , k], pop_gen_frq_pre_resmp[4, , k]), length.out = 50), plot = FALSE)
  hist_pre_resmp <- hist(pop_gen_frq_pre_resmp[4, , k], breaks = seq(min(pop_gen_frq_pst_resmp[4, , k], pop_gen_frq_pre_resmp[4, , k]), max(pop_gen_frq_pst_resmp[4, , k], pop_gen_frq_pre_resmp[4, , k]), length.out = 50), plot = FALSE)
  hist(pop_gen_frq_pst_resmp[4, , k], breaks = seq(min(pop_gen_frq_pst_resmp[4, , k], pop_gen_frq_pre_resmp[4, , k]), max(pop_gen_frq_pst_resmp[4, , k], pop_gen_frq_pre_resmp[4, , k]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
       xlim = c(min(pop_gen_frq_pst_resmp[4, , k], pop_gen_frq_pre_resmp[4, , k], smp_gen_frq[4, k]), max(pop_gen_frq_pst_resmp[4, , k], pop_gen_frq_pre_resmp[4, , k], smp_gen_frq[4, k])), ylim = c(0, max(hist_pst_resmp$density, hist_pre_resmp$density)),
       xlab = "Genotype frequency",
       main = paste("Genotype Aa/EE at generation", smp_gen[k]))
  hist(pop_gen_frq_pre_resmp[4, , k], breaks = seq(min(pop_gen_frq_pst_resmp[4, , k], pop_gen_frq_pre_resmp[4, , k]), max(pop_gen_frq_pst_resmp[4, , k], pop_gen_frq_pre_resmp[4, , k]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
  abline(v = smp_gen_frq[4, k], col = 'red', lty = 2, lwd = 2)
  hist_pst_resmp <- hist(pop_gen_frq_pst_resmp[5, , k], breaks = seq(min(pop_gen_frq_pst_resmp[5, , k], pop_gen_frq_pre_resmp[5, , k]), max(pop_gen_frq_pst_resmp[5, , k], pop_gen_frq_pre_resmp[5, , k]), length.out = 50), plot = FALSE)
  hist_pre_resmp <- hist(pop_gen_frq_pre_resmp[5, , k], breaks = seq(min(pop_gen_frq_pst_resmp[5, , k], pop_gen_frq_pre_resmp[5, , k]), max(pop_gen_frq_pst_resmp[5, , k], pop_gen_frq_pre_resmp[5, , k]), length.out = 50), plot = FALSE)
  hist(pop_gen_frq_pst_resmp[5, , k], breaks = seq(min(pop_gen_frq_pst_resmp[5, , k], pop_gen_frq_pre_resmp[5, , k]), max(pop_gen_frq_pst_resmp[5, , k], pop_gen_frq_pre_resmp[5, , k]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
       xlim = c(min(pop_gen_frq_pst_resmp[5, , k], pop_gen_frq_pre_resmp[5, , k], smp_gen_frq[5, k]), max(pop_gen_frq_pst_resmp[5, , k], pop_gen_frq_pre_resmp[5, , k], smp_gen_frq[5, k])), ylim = c(0, max(hist_pst_resmp$density, hist_pre_resmp$density)),
       xlab = "Genotype frequency",
       main = paste("Genotype Aa/Ee at generation", smp_gen[k]))
  hist(pop_gen_frq_pre_resmp[5, , k], breaks = seq(min(pop_gen_frq_pst_resmp[5, , k], pop_gen_frq_pre_resmp[5, , k]), max(pop_gen_frq_pst_resmp[5, , k], pop_gen_frq_pre_resmp[5, , k]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
  abline(v = smp_gen_frq[5, k], col = 'red', lty = 2, lwd = 2)

  hist_pst_resmp <- hist(pop_gen_frq_pst_resmp[6, , k], breaks = seq(min(pop_gen_frq_pst_resmp[6, , k], pop_gen_frq_pre_resmp[6, , k]), max(pop_gen_frq_pst_resmp[6, , k], pop_gen_frq_pre_resmp[6, , k]), length.out = 50), plot = FALSE)
  hist_pre_resmp <- hist(pop_gen_frq_pre_resmp[6, , k], breaks = seq(min(pop_gen_frq_pst_resmp[6, , k], pop_gen_frq_pre_resmp[6, , k]), max(pop_gen_frq_pst_resmp[6, , k], pop_gen_frq_pre_resmp[6, , k]), length.out = 50), plot = FALSE)
  hist(pop_gen_frq_pst_resmp[6, , k], breaks = seq(min(pop_gen_frq_pst_resmp[6, , k], pop_gen_frq_pre_resmp[6, , k]), max(pop_gen_frq_pst_resmp[6, , k], pop_gen_frq_pre_resmp[6, , k]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
       xlim = c(min(pop_gen_frq_pst_resmp[6, , k], pop_gen_frq_pre_resmp[6, , k], smp_gen_frq[6, k]), max(pop_gen_frq_pst_resmp[6, , k], pop_gen_frq_pre_resmp[6, , k], smp_gen_frq[6, k])), ylim = c(0, max(hist_pst_resmp$density, hist_pre_resmp$density)),
       xlab = "Genotype frequency",
       main = paste("Genotype aa/EE at generation", smp_gen[k]))
  hist(pop_gen_frq_pre_resmp[6, , k], breaks = seq(min(pop_gen_frq_pst_resmp[6, , k], pop_gen_frq_pre_resmp[6, , k]), max(pop_gen_frq_pst_resmp[6, , k], pop_gen_frq_pre_resmp[6, , k]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
  abline(v = smp_gen_frq[6, k], col = 'red', lty = 2, lwd = 2)

  hist_pst_resmp <- hist(pop_gen_frq_pst_resmp[7, , k], breaks = seq(min(pop_gen_frq_pst_resmp[7, , k], pop_gen_frq_pre_resmp[7, , k]), max(pop_gen_frq_pst_resmp[7, , k], pop_gen_frq_pre_resmp[7, , k]), length.out = 50), plot = FALSE)
  hist_pre_resmp <- hist(pop_gen_frq_pre_resmp[7, , k], breaks = seq(min(pop_gen_frq_pst_resmp[7, , k], pop_gen_frq_pre_resmp[7, , k]), max(pop_gen_frq_pst_resmp[7, , k], pop_gen_frq_pre_resmp[7, , k]), length.out = 50), plot = FALSE)
  hist(pop_gen_frq_pst_resmp[7, , k], breaks = seq(min(pop_gen_frq_pst_resmp[7, , k], pop_gen_frq_pre_resmp[7, , k]), max(pop_gen_frq_pst_resmp[7, , k], pop_gen_frq_pre_resmp[7, , k]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
       xlim = c(min(pop_gen_frq_pst_resmp[7, , k], pop_gen_frq_pre_resmp[7, , k], smp_gen_frq[7, k]), max(pop_gen_frq_pst_resmp[7, , k], pop_gen_frq_pre_resmp[7, , k], smp_gen_frq[7, k])), ylim = c(0, max(hist_pst_resmp$density, hist_pre_resmp$density)),
       xlab = "Genotype frequency",
       main = paste("Genotype Aa/ee at generation", smp_gen[k]))
  hist(pop_gen_frq_pre_resmp[7, , k], breaks = seq(min(pop_gen_frq_pst_resmp[7, , k], pop_gen_frq_pre_resmp[7, , k]), max(pop_gen_frq_pst_resmp[7, , k], pop_gen_frq_pre_resmp[7, , k]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
  abline(v = smp_gen_frq[7, k], col = 'red', lty = 2, lwd = 2)

  hist_pst_resmp <- hist(pop_gen_frq_pst_resmp[8, , k], breaks = seq(min(pop_gen_frq_pst_resmp[8, , k], pop_gen_frq_pre_resmp[8, , k]), max(pop_gen_frq_pst_resmp[8, , k], pop_gen_frq_pre_resmp[8, , k]), length.out = 50), plot = FALSE)
  hist_pre_resmp <- hist(pop_gen_frq_pre_resmp[8, , k], breaks = seq(min(pop_gen_frq_pst_resmp[8, , k], pop_gen_frq_pre_resmp[8, , k]), max(pop_gen_frq_pst_resmp[8, , k], pop_gen_frq_pre_resmp[8, , k]), length.out = 50), plot = FALSE)
  hist(pop_gen_frq_pst_resmp[8, , k], breaks = seq(min(pop_gen_frq_pst_resmp[8, , k], pop_gen_frq_pre_resmp[8, , k]), max(pop_gen_frq_pst_resmp[8, , k], pop_gen_frq_pre_resmp[8, , k]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
       xlim = c(min(pop_gen_frq_pst_resmp[8, , k], pop_gen_frq_pre_resmp[8, , k], smp_gen_frq[8, k]), max(pop_gen_frq_pst_resmp[8, , k], pop_gen_frq_pre_resmp[8, , k], smp_gen_frq[8, k])), ylim = c(0, max(hist_pst_resmp$density, hist_pre_resmp$density)),
       xlab = "Genotype frequency",
       main = paste("Genotype aa/Ee at generation", smp_gen[k]))
  hist(pop_gen_frq_pre_resmp[8, , k], breaks = seq(min(pop_gen_frq_pst_resmp[8, , k], pop_gen_frq_pre_resmp[8, , k]), max(pop_gen_frq_pst_resmp[8, , k], pop_gen_frq_pre_resmp[8, , k]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
  abline(v = smp_gen_frq[8, k], col = 'red', lty = 2, lwd = 2)

  hist_pst_resmp <- hist(pop_gen_frq_pst_resmp[9, , k], breaks = seq(min(pop_gen_frq_pst_resmp[9, , k], pop_gen_frq_pre_resmp[9, , k]), max(pop_gen_frq_pst_resmp[9, , k], pop_gen_frq_pre_resmp[9, , k]), length.out = 50), plot = FALSE)
  hist_pre_resmp <- hist(pop_gen_frq_pre_resmp[9, , k], breaks = seq(min(pop_gen_frq_pst_resmp[9, , k], pop_gen_frq_pre_resmp[9, , k]), max(pop_gen_frq_pst_resmp[9, , k], pop_gen_frq_pre_resmp[9, , k]), length.out = 50), plot = FALSE)
  hist(pop_gen_frq_pst_resmp[9, , k], breaks = seq(min(pop_gen_frq_pst_resmp[9, , k], pop_gen_frq_pre_resmp[9, , k]), max(pop_gen_frq_pst_resmp[9, , k], pop_gen_frq_pre_resmp[9, , k]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
       xlim = c(min(pop_gen_frq_pst_resmp[9, , k], pop_gen_frq_pre_resmp[9, , k], smp_gen_frq[9, k]), max(pop_gen_frq_pst_resmp[9, , k], pop_gen_frq_pre_resmp[9, , k], smp_gen_frq[9, k])), ylim = c(0, max(hist_pst_resmp$density, hist_pre_resmp$density)),
       xlab = "Genotype frequency",
       main = paste("Genotype aa/ee at generation", smp_gen[k]))
  hist(pop_gen_frq_pre_resmp[9, , k], breaks = seq(min(pop_gen_frq_pst_resmp[9, , k], pop_gen_frq_pre_resmp[9, , k]), max(pop_gen_frq_pst_resmp[9, , k], pop_gen_frq_pre_resmp[9, , k]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
  abline(v = smp_gen_frq[9, k], col = 'red', lty = 2, lwd = 2)
}
dev.off()

########################################

#' Calculate the optimal particle number in the particle marginal Metropolis-Hastings (PMMH)
#' Parameter settings
#' @param sel_cof the selection coefficients of the black and chestnut phenotypes
#' @param rec_rat the recombination rate between the ASIP and MC1R loci
#' @param pop_siz the size of the horse population (non-constant)
#' @param ref_siz the reference size of the horse population
#' @param evt_gen the generation that the event of interest occurred
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_siz the count of the horses drawn from the population at all sampling time points
#' @param smp_cnt the count of the genotypes observed in the sample at all sampling time points
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter
#' @param gap_num the number of particles increased or decreased in the optimal particle number search

load("./Output/Output v1.2/TEST_SimData.rda")

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
     file = "./Output/Output v1.2/TEST_OptNum.rda")

load("./Output/Output v1.2/TEST_OptNum.rda")

opt_pcl_num <- OptNum$opt_pcl_num
log_lik_sdv <- OptNum$log_lik_sdv

pdf(file = "./Output/Output v1.2/TEST_OptNum.rda.pdf", width = 12, height = 9)
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
#' @param sel_cof the selection coefficients of the black and chestnut phenotypes
#' @param rec_rat the recombination rate between the ASIP and MC1R loci
#' @param pop_siz the size of the horse population (non-constant)
#' @param ref_siz the reference size of the horse population
#' @param evt_gen the generation that the event of interest occurred
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_siz the count of the horses drawn from the population at all sampling time points
#' @param smp_cnt the count of the genotypes observed in the sample at all sampling time points
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter
#' @param itn_num the number of the iterations carried out in the particle marginal Metropolis-Hastings

load("./Output/Output v1.2/TEST_SimData.rda")

set.seed(test_seed)

sel_cof <- matrix(c(0e+00, 0e+00, 0e+00, 0e+00), nrow = 2, ncol = 2)
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

load("./Output/Output v1.2/TEST_SimData.rda")

save(sel_cof, rec_rat, pop_siz, ref_siz, evt_gen, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, sel_cof_chn,
     file = "./Output/Output v1.2/TEST_PMMH.rda")

load("./Output/Output v1.2/TEST_PMMH.rda")

pdf(file = "./Output/Output v1.2/TEST_PMMH_Traceplot.pdf", width = 12, height = 9)
par(mfrow = c(4, 1), mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
plot(1:itn_num, sel_cof_chn[1, 1, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot of the black before the event")
abline(h = sel_cof[1, 1], col = 'red', lty = 2, lwd = 2)
plot(1:itn_num, sel_cof_chn[2, 1, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot of the chestnut before the event")
abline(h = sel_cof[2, 1], col = 'red', lty = 2, lwd = 2)
plot(1:itn_num, sel_cof_chn[1, 2, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot of the black after the event")
abline(h = sel_cof[1, 2], col = 'red', lty = 2, lwd = 2)
plot(1:itn_num, sel_cof_chn[2, 2, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot of the chestnut after the event")
abline(h = sel_cof[2, 2], col = 'red', lty = 2, lwd = 2)
dev.off()

brn_num <- 1e+04
sel_cof_chn <- sel_cof_chn[, , brn_num:dim(sel_cof_chn)[3]]

thn_num <- 8e+00
sel_cof_chn <- sel_cof_chn[, , (1:round(dim(sel_cof_chn)[3] / thn_num)) * thn_num]

sel_cof_est <- matrix(NA, nrow = 2, ncol = 2)
sel_cof_est[, 1] <- rowMeans(sel_cof_chn[, 1, ])
sel_cof_est[, 2] <- rowMeans(sel_cof_chn[, 2, ])

sel_cof_hpd <- array(NA, dim = c(2, 2, 2))
sel_cof_hpd[1, , 1] <- HPDinterval(as.mcmc(sel_cof_chn[1, 1, ]), prob = 0.95)
sel_cof_hpd[2, , 1] <- HPDinterval(as.mcmc(sel_cof_chn[2, 1, ]), prob = 0.95)
sel_cof_hpd[1, , 2] <- HPDinterval(as.mcmc(sel_cof_chn[1, 2, ]), prob = 0.95)
sel_cof_hpd[2, , 2] <- HPDinterval(as.mcmc(sel_cof_chn[2, 2, ]), prob = 0.95)

pdf(file = "./Output/Output v1.2/TEST_PMMH_Posterior.pdf", width = 12, height = 9)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
hist(sel_cof_chn[1, 1, ], breaks = seq(min(sel_cof_chn[1, 1, ]), max(sel_cof_chn[1, 1, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the black before the event")
lines(density(sel_cof_chn[1, 1, ]), lwd = 2, col = 'black')
abline(v = sel_cof[1, 1], col = 'red', lty = 2, lwd = 2)
abline(v = sel_cof_est[1, 1], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 1, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 2, 1], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_chn[2, 1, ], breaks = seq(min(sel_cof_chn[2, 1, ]), max(sel_cof_chn[2, 1, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the chestnut before the event")
lines(density(sel_cof_chn[2, 1, ]), lwd = 2, col = 'black')
abline(v = sel_cof[2, 1], col = 'red', lty = 2, lwd = 2)
abline(v = sel_cof_est[2, 1], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 1, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 2, 1], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_chn[1, 2, ], breaks = seq(min(sel_cof_chn[1, 2, ]), max(sel_cof_chn[1, 2, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the black after the event")
lines(density(sel_cof_chn[1, 2, ]), lwd = 2, col = 'black')
abline(v = sel_cof[1, 2], col = 'red', lty = 2, lwd = 2)
abline(v = sel_cof_est[1, 2], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 1, 2], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 2, 2], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_chn[2, 2, ], breaks = seq(min(sel_cof_chn[2, 2, ]), max(sel_cof_chn[2, 2, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the chestnut after the event")
lines(density(sel_cof_chn[2, 2, ]), lwd = 2, col = 'black')
abline(v = sel_cof[2, 2], col = 'red', lty = 2, lwd = 2)
abline(v = sel_cof_est[2, 2], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 1, 2], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 2, 2], col = 'blue', lty = 2, lwd = 2)
dev.off()

########################################

#' Run the Bayesian procedure for the inference of natural selection
#' Parameter settings
#' @param sel_cof the selection coefficients of the black and chestnut phenotypes
#' @param rec_rat the recombination rate between the ASIP and MC1R loci
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

load("./Output/Output v1.2/TEST_SimData.rda")

set.seed(test_seed)

sel_cof <- matrix(c(0e+00, 0e+00, 0e+00, 0e+00), nrow = 2, ncol = 2)
rec_rat
pop_siz
ref_siz <- 1e+04
evt_gen
smp_gen
smp_siz
smp_cnt
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 1e+04
brn_num <- 4e+03
thn_num <- 3e+00

system.time(BayesianProcedure <- cmprunBayesianProcedure(sel_cof, rec_rat, pop_siz, ref_siz, evt_gen, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, brn_num, thn_num))

load("./Output/Output v1.2/TEST_SimData.rda")

save(sel_cof, rec_rat, pop_siz, ref_siz, evt_gen, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, brn_num, thn_num, BayesianProcedure,
     file = "./Output/Output v1.2/TEST_BayesianProcedure.rda")

load("./Output/Output v1.2/TEST_BayesianProcedure.rda")

sel_cof_chn <- BayesianProcedure$sel_cof_chn

sel_cof_est <- BayesianProcedure$sel_cof_est

sel_cof_hpd <- BayesianProcedure$sel_cof_hpd

pdf(file = "./Output/Output v1.2/TEST_BayesianProcedure_Posterior.pdf", width = 12, height = 9)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
hist(sel_cof_chn[1, 1, ], breaks = seq(min(sel_cof_chn[1, 1, ]), max(sel_cof_chn[1, 1, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the black before the event")
lines(density(sel_cof_chn[1, 1, ]), lwd = 2, col = 'black')
abline(v = sel_cof[1, 1], col = 'red', lty = 2, lwd = 2)
abline(v = sel_cof_est[1, 1], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 1, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 2, 1], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_chn[2, 1, ], breaks = seq(min(sel_cof_chn[2, 1, ]), max(sel_cof_chn[2, 1, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the chestnut before the event")
lines(density(sel_cof_chn[2, 1, ]), lwd = 2, col = 'black')
abline(v = sel_cof[2, 1], col = 'red', lty = 2, lwd = 2)
abline(v = sel_cof_est[2, 1], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 1, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 2, 1], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_chn[1, 2, ], breaks = seq(min(sel_cof_chn[1, 2, ]), max(sel_cof_chn[1, 2, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the black after the event")
lines(density(sel_cof_chn[1, 2, ]), lwd = 2, col = 'black')
abline(v = sel_cof[1, 2], col = 'red', lty = 2, lwd = 2)
abline(v = sel_cof_est[1, 2], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 1, 2], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 2, 2], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_chn[2, 2, ], breaks = seq(min(sel_cof_chn[2, 2, ]), max(sel_cof_chn[2, 2, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the chestnut after the event")
lines(density(sel_cof_chn[2, 2, ]), lwd = 2, col = 'black')
abline(v = sel_cof[2, 2], col = 'red', lty = 2, lwd = 2)
abline(v = sel_cof_est[2, 2], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 1, 2], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 2, 2], col = 'blue', lty = 2, lwd = 2)
dev.off()

################################################################################
