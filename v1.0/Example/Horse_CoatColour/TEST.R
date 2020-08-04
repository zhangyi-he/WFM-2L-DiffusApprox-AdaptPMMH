#' @title Inferring natural selection acting on horse coat colours and patterns during the process of domestication
#' @author Zhangyi He, Xiaoyang Dai, Mark Beaumont and Feng Yu

#' version 1.0

# set the directory
setwd("~/Documents/GitHub/WFM-2L-DiffusApprox-PMMH-Horse/v1.0/Example/Horse_CoatColour")

# call R functions
source("./RFUN.R")

#install.packages("RColorBrewer")
library("RColorBrewer")

#install.packages("viridis")
library("viridis")

#install.packages("ggplot2")
library("ggplot2")

#install.packages("plot3D")
library("plot3D")

setwd("~/Dropbox/Jeffery He/iResearch/Publications/2018/HE2020-WFM-2L-DiffusApprox-PMMH-Horse")

################################################################################

#' Simulate the haplotype frequency trajectories according to the Wright-Fisher model
#' parameter settings
#' @param s_b the selection coefficient of the black phenotype
#' @param s_c the selection coefficient of the chestnut phenotype
#' @param r_cb the recombination rate between ASIP and MC1R (r_bu = 0.5 is predetermined since ASIP and MC1R are located on a separate chromosome)
#' @param N the number of individuals in the population
#' @param int_frq the initial haplotype frequencies of the population
#' @param int_gen the first generation of the simulated haplotype frequency trajectories
#' @param lst_gen the last generation of the simulated haplotype frequency trajectories

s_b <- 5e-03
s_c <- 1e-03
r <- 5e-01
N <- 1e+04
int_frq <- rep(1, length.out = 4) / 4
int_gen <- 0
lst_gen <- 500

frq_pth <- cmpsimulateTLWFMS(s_b, s_c, r, N, int_frq, int_gen, lst_gen)

k <- int_gen:lst_gen
plot(k, frq_pth[1, ], type = "l", lwd = 1.5,
     xlab = "Generation", ylab = "Haplotype frequency",
     main = "A frequency trajectory of the AE haplotype")
plot(k, frq_pth[2, ], type = "l", lwd = 1.5,
     xlab = "Generation", ylab = "Haplotype frequency",
     main = "A frequency trajectory of the Ae haplotype")
plot(k, frq_pth[3, ], type = "l", lwd = 1.5,
     xlab = "Generation", ylab = "Haplotype frequency",
     main = "A frequency trajectory of the aE haplotype")
plot(k, frq_pth[4, ], type = "l", lwd = 1.5,
     xlab = "Generation", ylab = "Haplotype frequency",
     main = "A frequency trajectory of the ae haplotype")

########################################

#' Simulate the haplotype frequency trajectories according to the Wright-Fisher diffusion using the Euler-Maruyama method
#' Parameter setting
#' @param s_b the selection coefficient of the black phenotype
#' @param s_c the selection coefficient of the chestnut phenotype
#' @param r the recombination rate between ASIP and MC1R (r_bu = 0.5 is predetermined since ASIP and MC1R are located on a separate chromosome)
#' @param N the number of individuals in the population
#' @param int_frq the initial haplotype frequencies of the population
#' @param int_gen the first generation of the simulated haplotype frequency trajectories
#' @param lst_gen the last generation of the simulated haplotype frequency trajectories
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param data_augmentation = TRUE/FALSE (return the simulated sample trajectory with data augmentation or not)

s_b <- 5e-03
s_c <- 1e-03
r <- 5e-01
N <- 1e+04
int_frq <- rep(1, length.out = 4) / 4
int_gen <- 0
lst_gen <- 500
ptn_num <- 5e+00

frq_pth <- cmpsimulateTLWFDS(s_b, s_c, r, N, int_frq, int_gen, lst_gen, ptn_num, data_augmentation = TRUE)

t <- (int_gen:(int_gen + (lst_gen - int_gen) * ptn_num)) / 2 / N
plot(t, frq_pth[1, ], type = "l", lwd = 1.5,
     xlab = "Generation", ylab = "Haplotype frequency",
     main = "A frequency trajectory of the AE haplotype")
plot(t, frq_pth[2, ], type = "l", lwd = 1.5,
     xlab = "Generation", ylab = "Haplotype frequency",
     main = "A frequency trajectory of the Ae haplotype")
plot(t, frq_pth[3, ], type = "l", lwd = 1.5,
     xlab = "Generation", ylab = "Haplotype frequency",
     main = "A frequency trajectory of the aE haplotype")
plot(t, frq_pth[4, ], type = "l", lwd = 1.5,
     xlab = "Generation", ylab = "Haplotype frequency",
     main = "A frequency trajectory of the ae haplotype")

########################################

#' Compare the simulation generated with the Wright-Fisher model and the Wright-Fisher diffusion
s_b <- 5e-03
s_c <- 1e-03
r <- 5e-01
N <- 1e+04
int_frq <- rep(1, length.out = 4) / 4
int_gen <- 0
lst_gen <- 100
ptn_num <- 5e+00
sim_num <- 1e+06

sim_frq_WFM <- matrix(NA, nrow = 4, ncol = sim_num)
sim_frq_WFD <- matrix(NA, nrow = 4, ncol = sim_num)
for (i in 1:sim_num) {
  print(i)
  sim_frq_WFM[, i] <- cmpsimulateTLWFMS(s_b, s_c, r, N, int_frq, int_gen, lst_gen)[, (lst_gen - int_gen) + 1]
  sim_frq_WFD[, i] <- cmpsimulateTLWFDS(s_b, s_c, r, N, int_frq, int_gen, lst_gen, ptn_num, data_augmentation = FALSE)[, (lst_gen - int_gen) + 1]
}

save(s_b, s_c, r, N, int_frq, int_gen, lst_gen, ptn_num, sim_num, sim_frq_WFM, sim_frq_WFD,
     file = "./Output/Output v1.0/Test/WFM_vs_WFD.rda")

load("./Output/Output v1.0/Test/WFM_vs_WFD.rda")

pdf(file = "./Output/Output v1.0/Test/WFM_vs_WFD.pdf", width = 16, height = 12)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), oma = c(0, 0, 3, 0), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(sim_frq_WFM[1, ], breaks = seq(min(sim_frq_WFM[1, ], sim_frq_WFD[1, ]), max(sim_frq_WFM[1, ], sim_frq_WFD[1, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(sim_frq_WFM[1, ], sim_frq_WFD[1, ]), max(sim_frq_WFM[1, ], sim_frq_WFD[1, ])),
     xlab = "Haplotype frequency", main = "Haplotype AE")
hist(sim_frq_WFD[1, ], breaks = seq(min(sim_frq_WFM[1, ], sim_frq_WFD[1, ]), max(sim_frq_WFM[1, ], sim_frq_WFD[1, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)

hist(sim_frq_WFM[2, ], breaks = seq(min(sim_frq_WFM[2, ], sim_frq_WFD[2, ]), max(sim_frq_WFM[2, ], sim_frq_WFD[2, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(sim_frq_WFM[2, ], sim_frq_WFD[2, ]), max(sim_frq_WFM[2, ], sim_frq_WFD[2, ])),
     xlab = "Haplotype frequency", main = "Haplotype Ae")
hist(sim_frq_WFD[2, ], breaks = seq(min(sim_frq_WFM[2, ], sim_frq_WFD[2, ]), max(sim_frq_WFM[2, ], sim_frq_WFD[2, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)

hist(sim_frq_WFM[3, ], breaks = seq(min(sim_frq_WFM[3, ], sim_frq_WFD[3, ]), max(sim_frq_WFM[3, ], sim_frq_WFD[3, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(sim_frq_WFM[3, ], sim_frq_WFD[3, ]), max(sim_frq_WFM[3, ], sim_frq_WFD[3, ])),
     xlab = "Haplotype frequency", main = "Haplotype aE")
hist(sim_frq_WFD[3, ], breaks = seq(min(sim_frq_WFM[3, ], sim_frq_WFD[3, ]), max(sim_frq_WFM[3, ], sim_frq_WFD[3, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)

hist(sim_frq_WFM[4, ], breaks = seq(min(sim_frq_WFM[4, ], sim_frq_WFD[4, ]), max(sim_frq_WFM[4, ], sim_frq_WFD[4, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(sim_frq_WFM[4, ], sim_frq_WFD[4, ]), max(sim_frq_WFM[4, ], sim_frq_WFD[4, ])),
     xlab = "Haplotype frequency", main = "Haplotype ae")
hist(sim_frq_WFD[4, ], breaks = seq(min(sim_frq_WFM[4, ], sim_frq_WFD[4, ]), max(sim_frq_WFM[4, ], sim_frq_WFD[4, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
title(paste("Histograms of the haplotype frequencies in generation", lst_gen, "under the Wright-Fisher model and the Wright-Fisher diffusion"), outer = TRUE)
dev.off()

################################################################################

#' Simulate the hidden Markov model
#' parameter settings
#' @param model = "WFM"/"WFD" (return the observations from the underlying population under the WFM or WFD)
#' @param s_b the selection coefficient of the black phenotype
#' @param s_c the selection coefficient of the chestnut phenotype
#' @param r the recombination rate between ASIP and MC1R (r_bu = 0.5 is predetermined since ASIP and MC1R are located on a separate chromosome)
#' @param N the number of individuals in the population
#' @param int_frq the initial haplotype frequencies of the population
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_siz the count of the chromosomes drawn from the population at all sampling time points
#' @param ptn_num the number of the subintervals divided per generation in the Euler-Maruyama method for the WFD

#' Simulate the dataset under the Wright-Fisher model
model <- "WFM"
sel_cof <- c(1e-02, 5e-03)
dom_par <- c(5e-01, 5e-01)
rec_rat <- 1e-03
pop_siz <- 5e+03
int_frq <- rep(1, length.out = 4) / 4
smp_gen <- (0:10) * 50
smp_siz <- rep(50, 11)

sim_HMM_WFM <- cmpsimulateHMM(model, sel_cof, dom_par, rec_rat, pop_siz, int_frq, smp_gen, smp_siz)
smp_gen <- sim_HMM_WFM$smp_gen
smp_siz <- sim_HMM_WFM$smp_siz
smp_cnt <- sim_HMM_WFM$smp_cnt
smp_frq <- sim_HMM_WFM$smp_frq
pop_frq <- sim_HMM_WFM$pop_frq

k <- min(smp_gen):max(smp_gen)
smp_ale_frq <- smp_ale_cnt %*% diag(1 / smp_siz)
plot(k, pop_ale_frq[1, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_ale_frq[1, ], pop_ale_frq[1, ]), max(smp_ale_frq[1, ], pop_ale_frq[1, ])),
     xlab = "Generation", ylab = "Allele frequency",
     main = "A simulated dataset of the mutant allele at locus A generated with the Wright-Fisher model")
points(smp_gen, smp_ale_frq[1, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_ale_frq[2, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_ale_frq[2, ], pop_ale_frq[2, ]), max(smp_ale_frq[2, ], pop_ale_frq[2, ])),
     xlab = "Generation", ylab = "Allele frequency",
     main = "A simulated dataset of the mutant allele at locus B generated with the Wright-Fisher model")
points(smp_gen, smp_ale_frq[2, ], col = 'red', pch = 17, cex = 1)

####################

#' Simulate the dataset under the Wright-Fisher diffusion
model <- "WFD"
sel_cof <- c(1e-02, 5e-03)
dom_par <- c(5e-01, 5e-01)
rec_rat <- 1e-03
pop_siz <- 5e+03
int_frq <- c(1e-01, 2e-01, 3e-01, 4e-01)
smp_gen <- (0:10) * 50
smp_siz <- rep(50, 11)
ptn_num <- 5e+00

sim_HMM_WFD <- cmpsimulateHMM(model, sel_cof, dom_par, rec_rat, pop_siz, int_frq, smp_gen, smp_siz, ptn_num)
smp_gen <- sim_HMM_WFM$smp_gen
smp_siz <- sim_HMM_WFM$smp_siz
smp_cnt <- sim_HMM_WFM$smp_cnt
smp_frq <- sim_HMM_WFM$smp_frq
pop_frq <- sim_HMM_WFM$pop_frq

k <- min(smp_gen):max(smp_gen)
smp_ale_frq <- smp_ale_cnt %*% diag(1 / smp_siz)
plot(k, pop_ale_frq[1, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_ale_frq[1, ], pop_ale_frq[1, ]), max(smp_ale_frq[1, ], pop_ale_frq[1, ])),
     xlab = "Generation", ylab = "Allele frequency",
     main = "A simulated dataset of the mutant allele at locus A generated with the Wright-Fisher model")
points(smp_gen, smp_ale_frq[1, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_ale_frq[2, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_ale_frq[2, ], pop_ale_frq[2, ]), max(smp_ale_frq[2, ], pop_ale_frq[2, ])),
     xlab = "Generation", ylab = "Allele frequency",
     main = "A simulated dataset of the mutant allele at locus B generated with the Wright-Fisher model")
points(smp_gen, smp_ale_frq[2, ], col = 'red', pch = 17, cex = 1)

################################################################################

#' Generate a simulated dataset under the Wright-Fisher model
test_seed <- 7
set.seed(test_seed)

model <- "WFM"
sel_cof <- c(1e-02, 5e-03)
dom_par <- c(5e-01, 5e-01)
rec_rat <- 1e-03
pop_siz <- 5e+03
int_frq <- c(1e-01, 2e-01, 3e-01, 4e-01)
smp_gen <- (0:10) * 50
smp_siz <- rep(50, 11)

SimData <- cmpsimulateHMM(model, sel_cof, dom_par, rec_rat, pop_siz, int_frq, smp_gen, smp_siz)
smp_gen <- SimData$smp_gen
smp_siz <- SimData$smp_siz
smp_hap_cnt <- SimData$smp_hap_cnt
pop_hap_frq <- SimData$pop_hap_frq
smp_ale_cnt <- SimData$smp_ale_cnt
pop_ale_frq <- SimData$pop_ale_frq

save(sel_cof, dom_par, rec_rat, pop_siz, smp_gen, smp_siz, smp_hap_cnt, pop_hap_frq, smp_ale_cnt, pop_ale_frq,
     file = "./Output/Output v2.1/Test v2.1/TEST_2L_SimData.rda")

load("./Output/Output v2.1/Test v2.1/TEST_2L_SimData.rda")

k <- min(smp_gen):max(smp_gen)
smp_ale_frq <- smp_ale_cnt %*% diag(1 / smp_siz)

pdf(file = "./Output/Output v2.1/Test v2.1/TEST_2L_SimData.pdf", width = 20, height = 10)
par(mfrow = c(1, 2), oma = c(0, 0, 3, 0), mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
plot(k, pop_ale_frq[1, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(pop_ale_frq[1, ], smp_ale_frq[1, ]), max(pop_ale_frq[1, ], smp_ale_frq[1, ])),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Mutant allele at locus A")
points(smp_gen, smp_ale_frq[1, ], col = 'red', pch = 17, cex = 1)

plot(k, pop_ale_frq[2, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(pop_ale_frq[2, ], smp_ale_frq[2, ]), max(pop_ale_frq[2, ], smp_ale_frq[2, ])),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Mutant allele at locus B")
points(smp_gen, smp_ale_frq[2, ], col = 'red', pch = 17, cex = 1)
title("A simulated dataset without missing values generated with the Wright-Fisher model", outer = TRUE)
dev.off()

########################################

#' Run the bootstrap particle filter (BPF) with the Wright-Fisher diffusion
#' Parameter setting
#' @param s_b the selection coefficient of the black phenotype
#' @param s_c the selection coefficient of the chestnut phenotype
#' @param r the recombination rate between ASIP and MC1R (r_bu = 0.5 is predetermined since ASIP and MC1R are located on a separate chromosome)
#' @param N the number of individuals in the population
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_siz the count of the chromosomes drawn from the population at all sampling time points
#' @param smp_cnt 
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter

load("./Output/Output v2.1/Test v2.1/TEST_2L_SimData.rda")

set.seed(test_seed)

s_b
s_c
r
N
smp_gen
smp_siz
smp_cnt
ptn_num <- 5e+00
pcl_num <- 1e+05

system.time(BPF <- cmprunBPF(s_b, s_c, r, N, smp_gen, smp_siz, smp_ale_cnt, ptn_num, pcl_num))

save(s_b, s_c, r, N, smp_gen, smp_siz, smp_ale_cnt, ptn_num, pcl_num, BPF,
     file = "./Output/Output v1.0/Test/BPF.rda")

load("./Output/Output v1.0/Test/BPF.rda")

lik <- rep(1, pcl_num)
wght <- BPF$wght
for (k in 1:length(smp_gen)) {
  lik <- lik * (cumsum(wght[, k]) / (1:pcl_num))
}

pdf(file = "./Output/Output v1.0/Test/BPF_Likelihood.pdf", width = 12, height = 9)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:pcl_num, log(lik), type = 'l',
     xlab = "Number of particles", ylab = "Log likelihood",
     main = "Log likelihood estimated with the bootstrap particle filter")
dev.off()

pop_frq_pre_resmp <- BPF$pop_frq_pre_resmp
pop_frq_pst_resmp <- BPF$pop_frq_pst_resmp

pdf(file = "./Output/Output v2.1/Test v2.1/TEST_2L_BPF_Particle.pdf", width = 20, height = 55)
par(mfrow = c(11, 4), oma = c(0, 0, 3, 0), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
for (k in 1:length(smp_gen)) {
  hist_pst_resmp <- hist(pop_hap_frq_pst_resmp[1, , k], breaks = seq(min(pop_hap_frq_pst_resmp[1, , k], pop_hap_frq_pre_resmp[1, , k]), max(pop_hap_frq_pst_resmp[1, , k], pop_hap_frq_pre_resmp[1, , k]), length.out = 50), plot = FALSE)
  hist_pre_resmp <- hist(pop_hap_frq_pre_resmp[1, , k], breaks = seq(min(pop_hap_frq_pst_resmp[1, , k], pop_hap_frq_pre_resmp[1, , k]), max(pop_hap_frq_pst_resmp[1, , k], pop_hap_frq_pre_resmp[1, , k]), length.out = 50), plot = FALSE)
  hist(pop_hap_frq_pst_resmp[1, , k], breaks = seq(min(pop_hap_frq_pst_resmp[1, , k], pop_hap_frq_pre_resmp[1, , k]), max(pop_hap_frq_pst_resmp[1, , k], pop_hap_frq_pre_resmp[1, , k]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
       xlim = c(min(pop_hap_frq_pst_resmp[1, , k], pop_hap_frq_pre_resmp[1, , k], smp_hap_frq[1, k]), max(pop_hap_frq_pst_resmp[1, , k], pop_hap_frq_pre_resmp[1, , k], smp_hap_frq[1, k])), ylim = c(0, max(hist_pst_resmp$density, hist_pre_resmp$density)),
       xlab = "Haplotype frequency",
       main = paste("Haplotype A1B1 in generation", smp_gen[k]))
  hist(pop_hap_frq_pre_resmp[1, , k], breaks = seq(min(pop_hap_frq_pst_resmp[1, , k], pop_hap_frq_pre_resmp[1, , k]), max(pop_hap_frq_pst_resmp[1, , k], pop_hap_frq_pre_resmp[1, , k]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
  abline(v = smp_hap_frq[1, k], col = 'red', lty = 2, lwd = 2)

  hist_pst_resmp <- hist(pop_hap_frq_pst_resmp[2, , k], breaks = seq(min(pop_hap_frq_pst_resmp[2, , k], pop_hap_frq_pre_resmp[2, , k]), max(pop_hap_frq_pst_resmp[2, , k], pop_hap_frq_pre_resmp[2, , k]), length.out = 50), plot = FALSE)
  hist_pre_resmp <- hist(pop_hap_frq_pre_resmp[2, , k], breaks = seq(min(pop_hap_frq_pst_resmp[2, , k], pop_hap_frq_pre_resmp[2, , k]), max(pop_hap_frq_pst_resmp[2, , k], pop_hap_frq_pre_resmp[2, , k]), length.out = 50), plot = FALSE)
  hist(pop_hap_frq_pst_resmp[2, , k], breaks = seq(min(pop_hap_frq_pst_resmp[2, , k], pop_hap_frq_pre_resmp[2, , k]), max(pop_hap_frq_pst_resmp[2, , k], pop_hap_frq_pre_resmp[2, , k]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
       xlim = c(min(pop_hap_frq_pst_resmp[2, , k], pop_hap_frq_pre_resmp[2, , k], smp_hap_frq[2, k]), max(pop_hap_frq_pst_resmp[2, , k], pop_hap_frq_pre_resmp[2, , k], smp_hap_frq[2, k])), ylim = c(0, max(hist_pst_resmp$density, hist_pre_resmp$density)),
       xlab = "Haplotype frequency",
       main = paste("Haplotype A1B2 in generation", smp_gen[k]))
  hist(pop_hap_frq_pre_resmp[2, , k], breaks = seq(min(pop_hap_frq_pst_resmp[2, , k], pop_hap_frq_pre_resmp[2, , k]), max(pop_hap_frq_pst_resmp[2, , k], pop_hap_frq_pre_resmp[2, , k]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
  abline(v = smp_hap_frq[2, k], col = 'red', lty = 2, lwd = 2)

  hist_pst_resmp <- hist(pop_hap_frq_pst_resmp[3, , k], breaks = seq(min(pop_hap_frq_pst_resmp[3, , k], pop_hap_frq_pre_resmp[3, , k]), max(pop_hap_frq_pst_resmp[3, , k], pop_hap_frq_pre_resmp[3, , k]), length.out = 50), plot = FALSE)
  hist_pre_resmp <- hist(pop_hap_frq_pre_resmp[3, , k], breaks = seq(min(pop_hap_frq_pst_resmp[3, , k], pop_hap_frq_pre_resmp[3, , k]), max(pop_hap_frq_pst_resmp[3, , k], pop_hap_frq_pre_resmp[3, , k]), length.out = 50), plot = FALSE)
  hist(pop_hap_frq_pst_resmp[3, , k], breaks = seq(min(pop_hap_frq_pst_resmp[3, , k], pop_hap_frq_pre_resmp[3, , k]), max(pop_hap_frq_pst_resmp[3, , k], pop_hap_frq_pre_resmp[3, , k]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
       xlim = c(min(pop_hap_frq_pst_resmp[3, , k], pop_hap_frq_pre_resmp[3, , k], smp_hap_frq[3, k]), max(pop_hap_frq_pst_resmp[3, , k], pop_hap_frq_pre_resmp[3, , k], smp_hap_frq[3, k])), ylim = c(0, max(hist_pst_resmp$density, hist_pre_resmp$density)),
       xlab = "Haplotype frequency",
       main = paste("Haplotype A2B1 in generation", smp_gen[k]))
  hist(pop_hap_frq_pre_resmp[3, , k], breaks = seq(min(pop_hap_frq_pst_resmp[3, , k], pop_hap_frq_pre_resmp[3, , k]), max(pop_hap_frq_pst_resmp[3, , k], pop_hap_frq_pre_resmp[3, , k]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
  abline(v = smp_hap_frq[3, k], col = 'red', lty = 2, lwd = 2)

  hist_pst_resmp <- hist(pop_hap_frq_pst_resmp[4, , k], breaks = seq(min(pop_hap_frq_pst_resmp[4, , k], pop_hap_frq_pre_resmp[4, , k]), max(pop_hap_frq_pst_resmp[4, , k], pop_hap_frq_pre_resmp[4, , k]), length.out = 50), plot = FALSE)
  hist_pre_resmp <- hist(pop_hap_frq_pre_resmp[4, , k], breaks = seq(min(pop_hap_frq_pst_resmp[4, , k], pop_hap_frq_pre_resmp[4, , k]), max(pop_hap_frq_pst_resmp[4, , k], pop_hap_frq_pre_resmp[4, , k]), length.out = 50), plot = FALSE)
  hist(pop_hap_frq_pst_resmp[4, , k], breaks = seq(min(pop_hap_frq_pst_resmp[4, , k], pop_hap_frq_pre_resmp[4, , k]), max(pop_hap_frq_pst_resmp[4, , k], pop_hap_frq_pre_resmp[4, , k]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
       xlim = c(min(pop_hap_frq_pst_resmp[4, , k], pop_hap_frq_pre_resmp[4, , k], smp_hap_frq[4, k]), max(pop_hap_frq_pst_resmp[4, , k], pop_hap_frq_pre_resmp[4, , k], smp_hap_frq[4, k])), ylim = c(0, max(hist_pst_resmp$density, hist_pre_resmp$density)),
       xlab = "Haplotype frequency",
       main = paste("Haplotype A2B2 in generation", smp_gen[k]))
  hist(pop_hap_frq_pre_resmp[4, , k], breaks = seq(min(pop_hap_frq_pst_resmp[4, , k], pop_hap_frq_pre_resmp[4, , k]), max(pop_hap_frq_pst_resmp[4, , k], pop_hap_frq_pre_resmp[4, , k]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
  abline(v = smp_hap_frq[4, k], col = 'red', lty = 2, lwd = 2)
}
title("Histograms of the pre- and post-resampling particles", outer = TRUE)
dev.off()

########################################

#' Calculate the optimal particle number in the particle marginal Metropolis-Hastings (PMMH)
#' Parameter settings
#' @param s_b the selection coefficient of the black phenotype
#' @param s_c the selection coefficient of the chestnut phenotype
#' @param r the recombination rate between ASIP and MC1R (r_bu = 0.5 is predetermined since ASIP and MC1R are located on a separate chromosome)
#' @param N the number of individuals in the population
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_siz the count of the chromosomes drawn from the population at all sampling time points
#' @param smp_cnt 
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter
#' @param gap_num the number of particles increased or decreased in the optimal particle number search

load("./Output/Output v1.0/Test/SimData.rda")

set.seed(test_seed)

s_b
s_c
r
N
smp_gen
smp_siz
smp_cnt
ptn_num <- 5e+00
pcl_num <- 1e+03
gap_num <- 1e+02

system.time(OptNum <- calculateOptimalParticleNum(s_b, s_c, r, N, smp_gen, smp_siz, smp_ale_cnt, ptn_num, pcl_num, gap_num))

save(s_b, s_c, r, N, smp_gen, smp_siz, smp_ale_cnt, ptn_num, pcl_num, gap_num, OptNum,
     file = "./Output/Output v1.0/Test/OptNum.rda")

load("./Output/Output v1.0/Test/OptNum.rda")

opt_pcl_num <- OptNum$opt_pcl_num
log_lik_sdv <- OptNum$log_lik_sdv

pdf(file = "./Output/Output v1.0/Test/OptNum.pdf", width = 12, height = 9)
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
#' @param s_b the selection coefficient of the black phenotype
#' @param s_c the selection coefficient of the chestnut phenotype
#' @param r the recombination rate between ASIP and MC1R (r_bu = 0.5 is predetermined since ASIP and MC1R are located on a separate chromosome)
#' @param N the number of individuals in the population
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_siz the count of the chromosomes drawn from the population at all sampling time points
#' @param smp_cnt 
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter
#' @param itn_num the number of the iterations carried out in the particle marginal Metropolis-Hastings

load("./Output/Output v1.0/Test/SimData.rda")

set.seed(test_seed)

s_b <- 0e+00
s_c <- 0e+00
r <- 5e-01
N <- 1e+04
smp_gen
smp_siz
smp_cnt
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04

system.time(PMMH <- cmprunPMMH(s_b, s_c, r, N, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num))

save(s_b, s_c, r, N, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, PMMH,
     file = "./Output/Output v1.0/Test/PMMH.rda")

load("./Output/Output v1.0/Test/PMMH.rda")

load("./Output/Output v1.0/Test/SimData.rda")

s_b_chn <- PMMH$s_b_chn
s_c_chn <- PMMH$s_c_chn

pdf(file = "./Output/Output v1.0/Test/PMMH_Traceplot.pdf", width = 12, height = 9)
par(mfrow = c(2, 1), mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
plot(1:itn_num, s_b_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot of the selection coefficient of the black phenotype")
abline(h = s_b, col = 'red', lty = 2, lwd = 2)
plot(1:itn_num, s_c_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot of the selection coefficient of the chestnut phenotype")
abline(h = s_c, col = 'red', lty = 2, lwd = 2)
dev.off()

brn_num <- 4e+03
s_b_chn <- s_b_chn[brn_num:length(s_b_chn)]
s_c_chn <- s_c_chn[brn_num:length(s_c_chn)]

thn_num <- 8e+00
s_b_chn <- s_b_chn[(1:round(length(s_b_chn) / thn_num)) * thn_num]
s_c_chn <- s_c_chn[(1:round(length(s_c_chn) / thn_num)) * thn_num]

s_b_mmse <- mean(s_b_chn)
s_c_mmse <- mean(s_c_chn)

pdf(file = "./Output/Output v1.0/Test/PMMH_MarginalPosterior.pdf", width = 16, height = 6)
par(mfrow = c(3, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
hist(s_b_chn, breaks = seq(min(s_b_chn), max(s_b_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the selection coefficient of the black phenotype")
lines(density(s_b_chn), lwd = 2, col = 'black')
abline(v = s_b, col = 'red', lty = 2, lwd = 2)
abline(v = s_b_mmse, col = 'black', lty = 2, lwd = 2)

hist(s_c_chn, breaks = seq(min(s_c_chn), max(s_c_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the selection coefficient of the chestnut phenotype")
lines(density(s_c_chn), lwd = 2, col = 'black')
abline(v = s_c, col = 'red', lty = 2, lwd = 2)
abline(v = s_c_mmse, col = 'black', lty = 2, lwd = 2)
dev.off()

grd_num <- 1e+03
sel_cof_pdf <- kde2d(s_b_chn, s_c_chn, n = grd_num)
pdf(file = "./Output/Output v1.0/Test/PMMH_JointPosterior.pdf", width = 8, height = 6)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
      xlab = "Selection coefficient of the black phenotype", ylab = "Selection coefficient of the chestnut phenotype",
      main = "Posterior for the selection coefficients of the black and chestnut phenotypes")
abline(v = s_b, col = 'red', lty = 2, lwd = 2)
abline(h = s_c, col = 'red', lty = 2, lwd = 2)
dev.off()

########################################

#' Run the Bayesian procedure for the inference of natural selection
#' Parameter settings
#' @param s_b the selection coefficient of the black phenotype
#' @param s_c the selection coefficient of the chestnut phenotype
#' @param r the recombination rate between ASIP and MC1R (r_bu = 0.5 is predetermined since ASIP and MC1R are located on a separate chromosome)
#' @param N the number of individuals in the population
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_siz the count of the chromosomes drawn from the population at all sampling time points
#' @param smp_cnt 
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter
#' @param itn_num the number of the iterations carried out in the particle marginal Metropolis-Hastings
#' @param brn_num the number of the iterations for burn-in
#' @param thn_num the number of the iterations for thinning
#' @param grd_num the number of the grids in the kernel density estimation

load("./Output/Output v1.0/Test/SimData.rda")

set.seed(test_seed)

s_b <- 0e+00
s_c <- 0e+00
r <- 5e-01
N <- 1e+04
smp_gen
smp_siz
smp_cnt
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04
brn_num <- 4e+03
thn_num <- 8e+00

system.time(BayesianProcedure <- cmprunBayesianProcedure(ss_b, s_c, r, N, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, brn_num, thn_num))

save(s_b, s_c, r, N, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, brn_num, thn_num, BayesianProcedure,
     file = "./Output/Output v1.0/Test/BayesianProcedure.rda")

load("./Output/Output v1.0/Test/BayesianProcedure.rda")

load("./Output/Output v1.0/Test/SimData.rda")

s_b_chn <- BayesianProcedure$s_b_chn
s_c_chn <- BayesianProcedure$s_c_chn

s_b_mmse <- BayesianProcedure$s_b_mmse
s_c_mmse <- BayesianProcedure$s_c_mmse

s_c_hpd <- BayesianProcedure$s_c_hpd
s_b_hpd <- BayesianProcedure$s_b_hpd

pdf(file = "./Output/Output v1.0/Test/BayesianProcedure_Posterior.pdf", width = 16, height = 6)
par(mfrow = c(1, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
hist(s_b_chn, breaks = seq(min(s_b_chn), max(s_b_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the selection coefficients of the black phenotype")
lines(density(s_b_chn), lwd = 2, col = 'black')
abline(v = s_b, col = 'red', lty = 2, lwd = 2)
abline(v = s_b_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = s_b_hpd[2], col = 'blue', lty = 2, lwd = 2)

hist(s_c_chn, breaks = seq(min(s_c_chn), max(s_c_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the selection coefficients of the chestnut phenotype")
lines(density(s_c_chn), lwd = 2, col = 'black')
abline(v = s_c, col = 'red', lty = 2, lwd = 2)
abline(v = s_c_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = s_c_hpd[2], col = 'blue', lty = 2, lwd = 2)
dev.off()

################################################################################
