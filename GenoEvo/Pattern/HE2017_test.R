#' @title Inferring fluctuating selection acting on horse coat colours and patterns from ancient DNA data
#' @author Zhangyi He, Xiaoyang Dai, Mark Beaumont and Feng Yu

#' Genotype evolution (Wright-Fisher diffusion)

#' Horse white spotting pattern: KIT13 and KIT16 (constant selection)

# set the directory
setwd("~/Dropbox/Jeffery He/iResearch/Publications/2017/HE2019-2L-WFD-PMMH-Horse-MBE")

source("./Code/Code v1.0/GenoEvo/Horse coat pattern/HE2017_rfun.R")

#install.packages("RColorBrewer")
library("RColorBrewer")

#install.packages("ggplot2")
library("ggplot2")

#install.packages("plot3D")
library("plot3D")

################################################################################

#' Simulate the allele/haplotype/genotype frequency trajectories according to the two-locus Wright-Fisher model with selection
#' Parameter setting
#' @param type = "allele"/"haplotype"/"genotype" (return the simulated allele/haplotype/genotype frequency trajectories)
#' @param phased = TRUE/FALSE (return the simulated phased/unphased gentoype frequency trajectories)
#' @param sel_cof the selection coefficients of the tobiano, sabino and mixed phenotypes
#' @param rec_rat the recombination rate between the KIT13 and KIT16 loci
#' @param pop_siz the number of the horses in the population
#' @param int_frq the initial haplotype frequencies of the population
#' @param int_gen the first generation of the simulated haplotype frequency trajectories
#' @param lst_gen the last generation of the simulated haplotype frequency trajectories

sel_cof <- c(5e-03, 1e-03, -5e-03)
rec_rat <- 5e-05
pop_siz <- 5e+03
int_frq <- rmultinom(1, size = pop_siz, prob = rep(1, 10) / sum(rep(1, 10))) / pop_siz
int_gen <- 0
lst_gen <- 500

#' simulate the allele frequency trajectories
ale_frq_pth <- cmpsimulateTLWFMS(type = "allele", phased = FALSE, sel_cof, rec_rat, pop_siz, int_frq, int_gen, lst_gen)

k <- int_gen:lst_gen
plot(k, ale_frq_pth[1, ], type = "l", lwd = 1.5,
     xlab = "Generation", ylab = "Allele frequency",
     main = "A frequency trajectory of the KM1 allele generated with the Wright-Fisher model")
plot(k, ale_frq_pth[2, ], type = "l", lwd = 1.5,
     xlab = "Generation", ylab = "Allele frequency",
     main = "A frequency trajectory of the SB1 allele generated with the Wright-Fisher model")

#' simulate the haplotype frequency trajectories
hap_frq_pth <- cmpsimulateTLWFMS(type = "haplotype", phased = FALSE, sel_cof, rec_rat, pop_siz, int_frq, int_gen, lst_gen)

k <- int_gen:lst_gen
plot(k, hap_frq_pth[1, ], type = "l", lwd = 1.5, 
     xlab = "Generation", ylab = "Haplotype frequency", 
     main = "A frequency trajectory of the KM0sb1 haplotype generated with the Wright-Fisher model")
plot(k, hap_frq_pth[2, ], type = "l", lwd = 1.5, 
     xlab = "Generation", ylab = "Haplotype frequency", 
     main = "A frequency trajectory of the KM0SB1 haplotype generated with the Wright-Fisher model")
plot(k, hap_frq_pth[3, ], type = "l", lwd = 1.5, 
     xlab = "Generation", ylab = "Haplotype frequency", 
     main = "A frequency trajectory of the KM1sb1 haplotype generated with the Wright-Fisher model")
plot(k, hap_frq_pth[4, ], type = "l", lwd = 1.5, 
     xlab = "Generation", ylab = "Haplotype frequency", 
     main = "A frequency trajectory of the KM1SB1 haplotype generated with the Wright-Fisher model")

#' simulate the unphased genotype frequency trajectories
gen_frq_pth <- cmpsimulateTLWFMS(type = "genotype", phased = FALSE, sel_cof, rec_rat, pop_siz, int_frq, int_gen, lst_gen)

k <- int_gen:lst_gen
plot(k, gen_frq_pth[1, ], type = "l", lwd = 1.5,
     xlab = "Generation", ylab = "Genotype frequency",
     main = "A frequency trajectory of the KM0KM0/sb1sb1 genotype generated with the Wright-Fisher model")
plot(k, gen_frq_pth[2, ], type = "l", lwd = 1.5,
     xlab = "Generation", ylab = "Genotype frequency",
     main = "A frequency trajectory of the KM0KM0/sb1SB1 genotype generated with the Wright-Fisher model")
plot(k, gen_frq_pth[3, ], type = "l", lwd = 1.5,
     xlab = "Generation", ylab = "Genotype frequency",
     main = "A frequency trajectory of the KM0KM0/SB1SB1 genotype generated with the Wright-Fisher model")
plot(k, gen_frq_pth[4, ], type = "l", lwd = 1.5,
     xlab = "Generation", ylab = "Genotype frequency",
     main = "A frequency trajectory of the KM0KM1/sb1sb1 genotype generated with the Wright-Fisher model")
plot(k, gen_frq_pth[5, ], type = "l", lwd = 1.5,
     xlab = "Generation", ylab = "Genotype frequency",
     main = "A frequency trajectory of the KM0KM1/sb1SB1 genotype generated with the Wright-Fisher model")
plot(k, gen_frq_pth[6, ], type = "l", lwd = 1.5,
     xlab = "Generation", ylab = "Genotype frequency",
     main = "A frequency trajectory of the KM0KM1/SB1SB1 genotype generated with the Wright-Fisher model")
plot(k, gen_frq_pth[7, ], type = "l", lwd = 1.5,
     xlab = "Generation", ylab = "Genotype frequency",
     main = "A frequency trajectory of the KM1KM1/sb1sb1 genotype generated with the Wright-Fisher model")
plot(k, gen_frq_pth[8, ], type = "l", lwd = 1.5,
     xlab = "Generation", ylab = "Genotype frequency",
     main = "A frequency trajectory of the KM1KM1/sb1SB1 genotype generated with the Wright-Fisher model")
plot(k, gen_frq_pth[9, ], type = "l", lwd = 1.5,
     xlab = "Generation", ylab = "Genotype frequency",
     main = "A frequency trajectory of the KM1KM1/SB1SB1 genotype generated with the Wright-Fisher model")

#' simulate the phased genotype frequency trajectories
gen_frq_pth <- cmpsimulateTLWFMS(type = "genotype", phased = TRUE, sel_cof, rec_rat, pop_siz, int_frq, int_gen, lst_gen)

k <- int_gen:lst_gen
plot(k, gen_frq_pth[1, ], type = "l", lwd = 1.5,
     xlab = "Generation", ylab = "Genotype frequency",
     main = "A frequency trajectory of the KM0sb1/KM0sb1 genotype generated with the Wright-Fisher model")
plot(k, gen_frq_pth[2, ], type = "l", lwd = 1.5,
     xlab = "Generation", ylab = "Genotype frequency",
     main = "A frequency trajectory of the KM0sb1/KM0SB1 genotype generated with the Wright-Fisher model")
plot(k, gen_frq_pth[3, ], type = "l", lwd = 1.5,
     xlab = "Generation", ylab = "Genotype frequency",
     main = "A frequency trajectory of the KM0SB1/KM0SB1 genotype generated with the Wright-Fisher model")
plot(k, gen_frq_pth[4, ], type = "l", lwd = 1.5,
     xlab = "Generation", ylab = "Genotype frequency",
     main = "A frequency trajectory of the KM0sb1/KM1sb1 genotype generated with the Wright-Fisher model")
plot(k, gen_frq_pth[5, ], type = "l", lwd = 1.5,
     xlab = "Generation", ylab = "Genotype frequency",
     main = "A frequency trajectory of the KM0SB1/KM1sb1 genotype generated with the Wright-Fisher model")
plot(k, gen_frq_pth[6, ], type = "l", lwd = 1.5,
     xlab = "Generation", ylab = "Genotype frequency",
     main = "A frequency trajectory of the KM0SB1/KM1SB1 genotype generated with the Wright-Fisher model")
plot(k, gen_frq_pth[7, ], type = "l", lwd = 1.5,
     xlab = "Generation", ylab = "Genotype frequency",
     main = "A frequency trajectory of the KM1sb1/KM1sb1 genotype generated with the Wright-Fisher model")
plot(k, gen_frq_pth[8, ], type = "l", lwd = 1.5,
     xlab = "Generation", ylab = "Genotype frequency",
     main = "A frequency trajectory of the KM1sb1/KM1SB1 genotype generated with the Wright-Fisher model")
plot(k, gen_frq_pth[9, ], type = "l", lwd = 1.5,
     xlab = "Generation", ylab = "Genotype frequency",
     main = "A frequency trajectory of the KM1SB1/KM1SB1 genotype generated with the Wright-Fisher model")
plot(k, gen_frq_pth[10, ], type = "l", lwd = 1.5,
     xlab = "Generation", ylab = "Genotype frequency",
     main = "A frequency trajectory of the KM0sb1/KM1SB1 genotype generated with the Wright-Fisher model")

########################################

#' Simulate the allele/haplotype/genotype frequency trajectories according to the two-locus Wright-Fisher diffusion with selection using the Euler-Maruyama method
#' Parameter setting
#' @param type = "allele"/"haplotype"/"genotype" (return the simulated allele/haplotype/genotype frequency trajectories)
#' @param phased = TRUE/FALSE (return the simulated phased/unphased gentoype frequency trajectories)
#' @param sel_cof the selection coefficients of the tobiano, sabino and mixed phenotypes
#' @param rec_rat the recombination rate between the KIT13 and KIT16 loci
#' @param pop_siz the number of the horses in the population
#' @param int_frq the initial haplotype frequencies of the population
#' @param int_gen the first generation of the simulated haplotype frequency trajectories
#' @param lst_gen the last generation of the simulated haplotype frequency trajectories
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param data_augmentation = TRUE/FALSE (return the simulated sample trajectory with data augmentation or not)

sel_cof <- c(5e-03, 1e-03, -5e-03)
rec_rat <- 5e-05
pop_siz <- 5e+03
int_frq <- rmultinom(1, size = pop_siz, prob = rep(1, 10) / sum(rep(1, 10))) / pop_siz
int_gen <- 0
lst_gen <- 500
ptn_num <- 5e+00

#' simulate the allele frequency trajectories
ale_frq_pth <- cmpsimulateTLWFDS(type = "allele", phased = FALSE, sel_cof, rec_rat, pop_siz, int_frq, int_gen, lst_gen, ptn_num)

k <- int_gen:lst_gen
plot(k, ale_frq_pth[1, ], type = "l", lwd = 1.5,
     xlab = "Generation", ylab = "Allele frequency",
     main = "A frequency trajectory of the KM1 allele generated with the Wright-Fisher model")
plot(k, ale_frq_pth[2, ], type = "l", lwd = 1.5,
     xlab = "Generation", ylab = "Allele frequency",
     main = "A frequency trajectory of the SB1 allele generated with the Wright-Fisher model")

#' simulate the haplotype frequency trajectories
hap_frq_pth <- cmpsimulateTLWFDS(type = "haplotype", phased = FALSE, sel_cof, rec_rat, pop_siz, int_frq, int_gen, lst_gen, ptn_num)

k <- int_gen:lst_gen
plot(k, hap_frq_pth[1, ], type = "l", lwd = 1.5, 
     xlab = "Generation", ylab = "Haplotype frequency", 
     main = "A frequency trajectory of the KM0sb1 haplotype generated with the Wright-Fisher model")
plot(k, hap_frq_pth[2, ], type = "l", lwd = 1.5, 
     xlab = "Generation", ylab = "Haplotype frequency", 
     main = "A frequency trajectory of the KM0SB1 haplotype generated with the Wright-Fisher model")
plot(k, hap_frq_pth[3, ], type = "l", lwd = 1.5, 
     xlab = "Generation", ylab = "Haplotype frequency", 
     main = "A frequency trajectory of the KM1sb1 haplotype generated with the Wright-Fisher model")
plot(k, hap_frq_pth[4, ], type = "l", lwd = 1.5, 
     xlab = "Generation", ylab = "Haplotype frequency", 
     main = "A frequency trajectory of the KM1SB1 haplotype generated with the Wright-Fisher model")

#' simulate the unphased genotype frequency trajectories
gen_frq_pth <- cmpsimulateTLWFDS(type = "genotype", phased = FALSE, sel_cof, rec_rat, pop_siz, int_frq, int_gen, lst_gen, ptn_num)

k <- int_gen:lst_gen
plot(k, gen_frq_pth[1, ], type = "l", lwd = 1.5,
     xlab = "Generation", ylab = "Genotype frequency",
     main = "A frequency trajectory of the KM0KM0/sb1sb1 genotype generated with the Wright-Fisher model")
plot(k, gen_frq_pth[2, ], type = "l", lwd = 1.5,
     xlab = "Generation", ylab = "Genotype frequency",
     main = "A frequency trajectory of the KM0KM0/sb1SB1 genotype generated with the Wright-Fisher model")
plot(k, gen_frq_pth[3, ], type = "l", lwd = 1.5,
     xlab = "Generation", ylab = "Genotype frequency",
     main = "A frequency trajectory of the KM0KM0/SB1SB1 genotype generated with the Wright-Fisher model")
plot(k, gen_frq_pth[4, ], type = "l", lwd = 1.5,
     xlab = "Generation", ylab = "Genotype frequency",
     main = "A frequency trajectory of the KM0KM1/sb1sb1 genotype generated with the Wright-Fisher model")
plot(k, gen_frq_pth[5, ], type = "l", lwd = 1.5,
     xlab = "Generation", ylab = "Genotype frequency",
     main = "A frequency trajectory of the KM0KM1/sb1SB1 genotype generated with the Wright-Fisher model")
plot(k, gen_frq_pth[6, ], type = "l", lwd = 1.5,
     xlab = "Generation", ylab = "Genotype frequency",
     main = "A frequency trajectory of the KM0KM1/SB1SB1 genotype generated with the Wright-Fisher model")
plot(k, gen_frq_pth[7, ], type = "l", lwd = 1.5,
     xlab = "Generation", ylab = "Genotype frequency",
     main = "A frequency trajectory of the KM1KM1/sb1sb1 genotype generated with the Wright-Fisher model")
plot(k, gen_frq_pth[8, ], type = "l", lwd = 1.5,
     xlab = "Generation", ylab = "Genotype frequency",
     main = "A frequency trajectory of the KM1KM1/sb1SB1 genotype generated with the Wright-Fisher model")
plot(k, gen_frq_pth[9, ], type = "l", lwd = 1.5,
     xlab = "Generation", ylab = "Genotype frequency",
     main = "A frequency trajectory of the KM1KM1/SB1SB1 genotype generated with the Wright-Fisher model")

#' simulate the phased genotype frequency trajectories
gen_frq_pth <- cmpsimulateTLWFDS(type = "genotype", phased = TRUE, sel_cof, rec_rat, pop_siz, int_frq, int_gen, lst_gen, ptn_num)

k <- int_gen:lst_gen
plot(k, gen_frq_pth[1, ], type = "l", lwd = 1.5,
     xlab = "Generation", ylab = "Genotype frequency",
     main = "A frequency trajectory of the KM0sb1/KM0sb1 genotype generated with the Wright-Fisher model")
plot(k, gen_frq_pth[2, ], type = "l", lwd = 1.5,
     xlab = "Generation", ylab = "Genotype frequency",
     main = "A frequency trajectory of the KM0sb1/KM0SB1 genotype generated with the Wright-Fisher model")
plot(k, gen_frq_pth[3, ], type = "l", lwd = 1.5,
     xlab = "Generation", ylab = "Genotype frequency",
     main = "A frequency trajectory of the KM0SB1/KM0SB1 genotype generated with the Wright-Fisher model")
plot(k, gen_frq_pth[4, ], type = "l", lwd = 1.5,
     xlab = "Generation", ylab = "Genotype frequency",
     main = "A frequency trajectory of the KM0sb1/KM1sb1 genotype generated with the Wright-Fisher model")
plot(k, gen_frq_pth[5, ], type = "l", lwd = 1.5,
     xlab = "Generation", ylab = "Genotype frequency",
     main = "A frequency trajectory of the KM0SB1/KM1sb1 genotype generated with the Wright-Fisher model")
plot(k, gen_frq_pth[6, ], type = "l", lwd = 1.5,
     xlab = "Generation", ylab = "Genotype frequency",
     main = "A frequency trajectory of the KM0SB1/KM1SB1 genotype generated with the Wright-Fisher model")
plot(k, gen_frq_pth[7, ], type = "l", lwd = 1.5,
     xlab = "Generation", ylab = "Genotype frequency",
     main = "A frequency trajectory of the KM1sb1/KM1sb1 genotype generated with the Wright-Fisher model")
plot(k, gen_frq_pth[8, ], type = "l", lwd = 1.5,
     xlab = "Generation", ylab = "Genotype frequency",
     main = "A frequency trajectory of the KM1sb1/KM1SB1 genotype generated with the Wright-Fisher model")
plot(k, gen_frq_pth[9, ], type = "l", lwd = 1.5,
     xlab = "Generation", ylab = "Genotype frequency",
     main = "A frequency trajectory of the KM1SB1/KM1SB1 genotype generated with the Wright-Fisher model")
plot(k, gen_frq_pth[10, ], type = "l", lwd = 1.5,
     xlab = "Generation", ylab = "Genotype frequency",
     main = "A frequency trajectory of the KM0sb1/KM1SB1 genotype generated with the Wright-Fisher model")

########################################

#' Compare the simulation generated with the Wright-Fisher model and the Wright-Fisher diffusion
sel_cof <- c(5e-03, 1e-03, -5e-03)
rec_rat <- 5e-05
pop_siz <- 5e+03
int_frq <- rmultinom(1, size = pop_siz, prob = rep(1, 10) / sum(rep(1, 10))) / pop_siz
int_gen <- 0
lst_gen <- 500
ptn_num <- 5e+00
sim_num <- 1e+06

#' compare the probability distributions for haplotype frequencies
sim_WFM_gen_frq <- matrix(NA, nrow = 10, ncol = sim_num)
sim_WFD_gen_frq <- matrix(NA, nrow = 10, ncol = sim_num)
for (i in 1:sim_num) {
  print(i)
  sim_WFM_gen_frq[, i] <- cmpsimulateTLWFMS(type = "genotype", phased = TRUE, sel_cof, rec_rat, pop_siz, int_frq, int_gen, lst_gen)[, (lst_gen - int_gen) + 1]
  sim_WFD_gen_frq[, i] <- cmpsimulateTLWFDS(type = "genotype", phased = TRUE, sel_cof, rec_rat, pop_siz, int_frq, int_gen, lst_gen, ptn_num)[, (lst_gen - int_gen) + 1]
}

save(sel_cof, rec_rat, pop_siz, int_frq, int_gen, lst_gen, ptn_num, sim_num, sim_WFM_gen_frq, sim_WFD_gen_frq,
     file = "./Output/Output v1.0/GenoEvo/Horse coat pattern/TEST_WFM_vs_WFD.rda")

load("./Output/Output v1.0/GenoEvo/Horse coat pattern/TEST_WFM_vs_WFD.rda")

pdf(file = "./Output/Output v1.0/GenoEvo/Horse coat pattern/TEST_WFM_vs_WFD.pdf", width = 20, height = 15)
par(mfrow = c(3, 4), mar = c(5.5, 5, 5.5, 2.5), oma = c(0, 0, 3, 0), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
hist(sim_WFM_gen_frq[1, ], breaks = seq(min(sim_WFM_gen_frq[1, ], sim_WFD_gen_frq[1, ]), max(sim_WFM_gen_frq[1, ], sim_WFD_gen_frq[1, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(sim_WFM_gen_frq[1, ], sim_WFD_gen_frq[1, ]), max(sim_WFM_gen_frq[1, ], sim_WFD_gen_frq[1, ])),
     xlab = "Genotype frequency", main = "Genotype KM0sb1/KM0sb1")
hist(sim_WFD_gen_frq[1, ], breaks = seq(min(sim_WFM_gen_frq[1, ], sim_WFD_gen_frq[1, ]), max(sim_WFM_gen_frq[1, ], sim_WFD_gen_frq[1, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
hist(sim_WFM_gen_frq[2, ], breaks = seq(min(sim_WFM_gen_frq[2, ], sim_WFD_gen_frq[2, ]), max(sim_WFM_gen_frq[2, ], sim_WFD_gen_frq[2, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(sim_WFM_gen_frq[2, ], sim_WFD_gen_frq[2, ]), max(sim_WFM_gen_frq[2, ], sim_WFD_gen_frq[2, ])),
     xlab = "Genotype frequency", main = "Genotype KM0sb1/KM0SB1")
hist(sim_WFD_gen_frq[2, ], breaks = seq(min(sim_WFM_gen_frq[2, ], sim_WFD_gen_frq[2, ]), max(sim_WFM_gen_frq[2, ], sim_WFD_gen_frq[2, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
hist(sim_WFM_gen_frq[3, ], breaks = seq(min(sim_WFM_gen_frq[3, ], sim_WFD_gen_frq[3, ]), max(sim_WFM_gen_frq[3, ], sim_WFD_gen_frq[3, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(sim_WFM_gen_frq[3, ], sim_WFD_gen_frq[3, ]), max(sim_WFM_gen_frq[3, ], sim_WFD_gen_frq[3, ])),
     xlab = "Genotype frequency", main = "Genotype KM0SB1/KM0SB1")
hist(sim_WFD_gen_frq[3, ], breaks = seq(min(sim_WFM_gen_frq[3, ], sim_WFD_gen_frq[3, ]), max(sim_WFM_gen_frq[3, ], sim_WFD_gen_frq[3, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
hist(sim_WFM_gen_frq[4, ], breaks = seq(min(sim_WFM_gen_frq[4, ], sim_WFD_gen_frq[4, ]), max(sim_WFM_gen_frq[4, ], sim_WFD_gen_frq[4, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(sim_WFM_gen_frq[4, ], sim_WFD_gen_frq[4, ]), max(sim_WFM_gen_frq[4, ], sim_WFD_gen_frq[4, ])),
     xlab = "Genotype frequency", main = "Genotype KM0sb1/KM1sb1")
hist(sim_WFD_gen_frq[4, ], breaks = seq(min(sim_WFM_gen_frq[4, ], sim_WFD_gen_frq[4, ]), max(sim_WFM_gen_frq[4, ], sim_WFD_gen_frq[4, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
hist(sim_WFM_gen_frq[5, ], breaks = seq(min(sim_WFM_gen_frq[5, ], sim_WFD_gen_frq[5, ]), max(sim_WFM_gen_frq[5, ], sim_WFD_gen_frq[5, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(sim_WFM_gen_frq[5, ], sim_WFD_gen_frq[5, ]), max(sim_WFM_gen_frq[5, ], sim_WFD_gen_frq[5, ])),
     xlab = "Genotype frequency", main = "Genotype KM0SB1/KM1sb1")
hist(sim_WFD_gen_frq[5, ], breaks = seq(min(sim_WFM_gen_frq[5, ], sim_WFD_gen_frq[5, ]), max(sim_WFM_gen_frq[5, ], sim_WFD_gen_frq[5, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
hist(sim_WFM_gen_frq[6, ], breaks = seq(min(sim_WFM_gen_frq[6, ], sim_WFD_gen_frq[6, ]), max(sim_WFM_gen_frq[6, ], sim_WFD_gen_frq[6, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(sim_WFM_gen_frq[6, ], sim_WFD_gen_frq[6, ]), max(sim_WFM_gen_frq[6, ], sim_WFD_gen_frq[6, ])),
     xlab = "Genotype frequency", main = "Genotype KM0SB1/KM1SB1")
hist(sim_WFD_gen_frq[6, ], breaks = seq(min(sim_WFM_gen_frq[6, ], sim_WFD_gen_frq[6, ]), max(sim_WFM_gen_frq[6, ], sim_WFD_gen_frq[6, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
hist(sim_WFM_gen_frq[7, ], breaks = seq(min(sim_WFM_gen_frq[7, ], sim_WFD_gen_frq[7, ]), max(sim_WFM_gen_frq[7, ], sim_WFD_gen_frq[7, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(sim_WFM_gen_frq[7, ], sim_WFD_gen_frq[7, ]), max(sim_WFM_gen_frq[7, ], sim_WFD_gen_frq[7, ])),
     xlab = "Genotype frequency", main = "Genotype KM1sb1/KM1sb1")
hist(sim_WFD_gen_frq[7, ], breaks = seq(min(sim_WFM_gen_frq[7, ], sim_WFD_gen_frq[7, ]), max(sim_WFM_gen_frq[7, ], sim_WFD_gen_frq[7, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
hist(sim_WFM_gen_frq[8, ], breaks = seq(min(sim_WFM_gen_frq[8, ], sim_WFD_gen_frq[8, ]), max(sim_WFM_gen_frq[8, ], sim_WFD_gen_frq[8, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(sim_WFM_gen_frq[8, ], sim_WFD_gen_frq[8, ]), max(sim_WFM_gen_frq[8, ], sim_WFD_gen_frq[8, ])),
     xlab = "Genotype frequency", main = "Genotype KM1sb1/KM1SB1")
hist(sim_WFD_gen_frq[8, ], breaks = seq(min(sim_WFM_gen_frq[8, ], sim_WFD_gen_frq[8, ]), max(sim_WFM_gen_frq[8, ], sim_WFD_gen_frq[8, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
hist(sim_WFM_gen_frq[9, ], breaks = seq(min(sim_WFM_gen_frq[9, ], sim_WFD_gen_frq[9, ]), max(sim_WFM_gen_frq[9, ], sim_WFD_gen_frq[9, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(sim_WFM_gen_frq[9, ], sim_WFD_gen_frq[9, ]), max(sim_WFM_gen_frq[9, ], sim_WFD_gen_frq[9, ])),
     xlab = "Genotype frequency", main = "Genotype KM1SB1/KM1SB1")
hist(sim_WFD_gen_frq[9, ], breaks = seq(min(sim_WFM_gen_frq[9, ], sim_WFD_gen_frq[9, ]), max(sim_WFM_gen_frq[9, ], sim_WFD_gen_frq[9, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
hist(sim_WFM_gen_frq[10, ], breaks = seq(min(sim_WFM_gen_frq[10, ], sim_WFD_gen_frq[10, ]), max(sim_WFM_gen_frq[10, ], sim_WFD_gen_frq[10, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(sim_WFM_gen_frq[10, ], sim_WFD_gen_frq[10, ]), max(sim_WFM_gen_frq[10, ], sim_WFD_gen_frq[10, ])),
     xlab = "Genotype frequency", main = "Genotype KM0sb1/KM1SB1")
hist(sim_WFD_gen_frq[10, ], breaks = seq(min(sim_WFM_gen_frq[10, ], sim_WFD_gen_frq[10, ]), max(sim_WFM_gen_frq[10, ], sim_WFD_gen_frq[10, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
title(paste("Histograms of the genotype frequencies in generation", lst_gen, "under the Wright-Fisher model and the Wright-Fisher diffusion"), outer = TRUE)
dev.off()

################################################################################

#' Simulate the hidden Markov model
#' Parameter setting
#' @param model = "WFM"/"WFD" (return the observations from the underlying population evolving according to the WFM/WFD)
#' @param missing = TRUE/FALSE (return the observations with missing values or not)
#' @param sel_cof the selection coefficients of the tobiano, sabino and mixed phenotypes
#' @param rec_rat the recombination rate between the KIT13 and KIT16 loci
#' @param pop_siz the number of the horses in the population
#' @param int_frq the initial frequencies of the four haplotypes in the population
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_siz the count of the diploid individuals drawn from the population at all sampling time points
#' @param ptn_num the number of the subintervals divided per generation in the Euler-Maruyama method for the WFD

#' Simulate the dataset under the Wright-Fisher model 
model <- "WFM"
sel_cof <- c(5e-03, 1e-03, -5e-03)
rec_rat <- 5e-05
pop_siz <- 5e+03
int_frq <- rmultinom(1, size = pop_siz, prob = rep(1, 10) / sum(rep(1, 10))) / pop_siz
smp_gen <- (0:10) * 50
smp_siz <- rep(50, 11)

#' without uncertain genotypes
missing <- FALSE
sim_HMM_WFM <- cmpsimulateHMM(model, missing, sel_cof, rec_rat, pop_siz, int_frq, smp_gen, smp_siz)
smp_gen <- sim_HMM_WFM$smp_gen
smp_siz <- sim_HMM_WFM$smp_siz
smp_cnt <- sim_HMM_WFM$smp_cnt
mis_cnt <- sim_HMM_WFM$mis_cnt
smp_frq <- sim_HMM_WFM$smp_cnt / sim_HMM_WFM$smp_siz
pop_frq <- sim_HMM_WFM$pop_frq

k <- min(smp_gen):max(smp_gen)
plot(k, pop_frq[1, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_frq[1, ], pop_frq[1, ]), max(smp_frq[1, ], pop_frq[1, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "A simulated dataset of the KM0KM0/sb1sb1 genotype generated with the Wright-Fisher model")
points(smp_gen, smp_frq[1, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_frq[2, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_frq[2, ], pop_frq[2, ]), max(smp_frq[2, ], pop_frq[2, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "A simulated dataset of the KM0KM0/sb1SB1 genotype generated with the Wright-Fisher model")
points(smp_gen, smp_frq[2, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_frq[3, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_frq[3, ], pop_frq[3, ]), max(smp_frq[3, ], pop_frq[3, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "A simulated dataset of the KM0KM0/SB1SB1 genotype generated with the Wright-Fisher model")
points(smp_gen, smp_frq[3, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_frq[4, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_frq[4, ], pop_frq[4, ]), max(smp_frq[4, ], pop_frq[4, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "A simulated dataset of the KM0KM1/sb1sb1 genotype generated with the Wright-Fisher model")
points(smp_gen, smp_frq[4, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_frq[5, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_frq[5, ], pop_frq[5, ]), max(smp_frq[5, ], pop_frq[5, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "A simulated dataset of the KM0KM1/sb1SB1 genotype generated with the Wright-Fisher model")
points(smp_gen, smp_frq[5, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_frq[6, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_frq[6, ], pop_frq[6, ]), max(smp_frq[6, ], pop_frq[6, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "A simulated dataset of the KM0KM1/SB1SB1 genotype generated with the Wright-Fisher model")
points(smp_gen, smp_frq[6, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_frq[7, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_frq[7, ], pop_frq[7, ]), max(smp_frq[7, ], pop_frq[7, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "A simulated dataset of the KM1KM1/sb1sb1 genotype generated with the Wright-Fisher model")
points(smp_gen, smp_frq[7, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_frq[8, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_frq[8, ], pop_frq[8, ]), max(smp_frq[8, ], pop_frq[8, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "A simulated dataset of the KM1KM1/sb1SB1 genotype generated with the Wright-Fisher model")
points(smp_gen, smp_frq[8, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_frq[9, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_frq[9, ], pop_frq[9, ]), max(smp_frq[9, ], pop_frq[9, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "A simulated dataset of the KM1KM1/SB1SB1 genotype generated with the Wright-Fisher model")
points(smp_gen, smp_frq[9, ], col = 'red', pch = 17, cex = 1)

#' with uncertain genotypes
missing <- TRUE
sim_HMM_WFM <- cmpsimulateHMM(model, missing, sel_cof, rec_rat, pop_siz, int_frq, smp_gen, smp_siz)
smp_gen <- sim_HMM_WFM$smp_gen
smp_siz <- sim_HMM_WFM$smp_siz
smp_cnt <- sim_HMM_WFM$smp_cnt
mis_cnt <- sim_HMM_WFM$mis_cnt
smp_frq <- sim_HMM_WFM$smp_cnt / sim_HMM_WFM$smp_siz
pop_frq <- sim_HMM_WFM$pop_frq

k <- min(smp_gen):max(smp_gen)
plot(k, pop_frq[1, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_frq[1, ], pop_frq[1, ]), max(smp_frq[1, ], pop_frq[1, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "A simulated dataset of the KM0KM0/sb1sb1 genotype generated with the Wright-Fisher model")
points(smp_gen, smp_frq[1, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_frq[2, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_frq[2, ], pop_frq[2, ]), max(smp_frq[2, ], pop_frq[2, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "A simulated dataset of the KM0KM0/sb1SB1 genotype generated with the Wright-Fisher model")
points(smp_gen, smp_frq[2, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_frq[3, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_frq[3, ], pop_frq[3, ]), max(smp_frq[3, ], pop_frq[3, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "A simulated dataset of the KM0KM0/SB1SB1 genotype generated with the Wright-Fisher model")
points(smp_gen, smp_frq[3, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_frq[4, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_frq[4, ], pop_frq[4, ]), max(smp_frq[4, ], pop_frq[4, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "A simulated dataset of the KM0KM1/sb1sb1 genotype generated with the Wright-Fisher model")
points(smp_gen, smp_frq[4, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_frq[5, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_frq[5, ], pop_frq[5, ]), max(smp_frq[5, ], pop_frq[5, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "A simulated dataset of the KM0KM1/sb1SB1 genotype generated with the Wright-Fisher model")
points(smp_gen, smp_frq[5, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_frq[6, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_frq[6, ], pop_frq[6, ]), max(smp_frq[6, ], pop_frq[6, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "A simulated dataset of the KM0KM1/SB1SB1 genotype generated with the Wright-Fisher model")
points(smp_gen, smp_frq[6, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_frq[7, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_frq[7, ], pop_frq[7, ]), max(smp_frq[7, ], pop_frq[7, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "A simulated dataset of the KM1KM1/sb1sb1 genotype generated with the Wright-Fisher model")
points(smp_gen, smp_frq[7, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_frq[8, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_frq[8, ], pop_frq[8, ]), max(smp_frq[8, ], pop_frq[8, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "A simulated dataset of the KM1KM1/sb1SB1 genotype generated with the Wright-Fisher model")
points(smp_gen, smp_frq[8, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_frq[9, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_frq[9, ], pop_frq[9, ]), max(smp_frq[9, ], pop_frq[9, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "A simulated dataset of the KM1KM1/SB1SB1 genotype generated with the Wright-Fisher model")
points(smp_gen, smp_frq[9, ], col = 'red', pch = 17, cex = 1)

####################

#' Simulate the dataset under the Wright-Fisher diffusion
model <- "WFD"
sel_cof <- c(5e-03, 1e-03, -5e-03)
rec_rat <- 5e-05
pop_siz <- 5e+03
int_frq <- rmultinom(1, size = pop_siz, prob = rep(1, 10) / sum(rep(1, 10))) / pop_siz
smp_gen <- (0:10) * 50
smp_siz <- rep(50, 11)
ptn_num <- 5e+00

#' without uncertain genotypes
missing <- FALSE
sim_HMM_WFD <- cmpsimulateHMM(model, missing, sel_cof, rec_rat, pop_siz, int_frq, smp_gen, smp_siz, ptn_num)
smp_gen <- sim_HMM_WFD$smp_gen
smp_siz <- sim_HMM_WFD$smp_siz
smp_cnt <- sim_HMM_WFD$smp_cnt
mis_cnt <- sim_HMM_WFD$mis_cnt
smp_frq <- sim_HMM_WFD$smp_cnt / sim_HMM_WFD$smp_siz
pop_frq <- sim_HMM_WFD$pop_frq

k <- min(smp_gen):max(smp_gen)
plot(k, pop_frq[1, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_frq[1, ], pop_frq[1, ]), max(smp_frq[1, ], pop_frq[1, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "A simulated dataset of the KM0KM0/sb1sb1 genotype generated with the Wright-Fisher diffusion")
points(smp_gen, smp_frq[1, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_frq[2, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_frq[2, ], pop_frq[2, ]), max(smp_frq[2, ], pop_frq[2, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "A simulated dataset of the KM0KM0/sb1SB1 genotype generated with the Wright-Fisher diffusion")
points(smp_gen, smp_frq[2, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_frq[3, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_frq[3, ], pop_frq[3, ]), max(smp_frq[3, ], pop_frq[3, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "A simulated dataset of the KM0KM0/SB1SB1 genotype generated with the Wright-Fisher diffusion")
points(smp_gen, smp_frq[3, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_frq[4, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_frq[4, ], pop_frq[4, ]), max(smp_frq[4, ], pop_frq[4, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "A simulated dataset of the KM0KM1/sb1sb1 genotype generated with the Wright-Fisher diffusion")
points(smp_gen, smp_frq[4, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_frq[5, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_frq[5, ], pop_frq[5, ]), max(smp_frq[5, ], pop_frq[5, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "A simulated dataset of the KM0KM1/sb1SB1 genotype generated with the Wright-Fisher diffusion")
points(smp_gen, smp_frq[5, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_frq[6, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_frq[6, ], pop_frq[6, ]), max(smp_frq[6, ], pop_frq[6, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "A simulated dataset of the KM0KM1/SB1SB1 genotype generated with the Wright-Fisher diffusion")
points(smp_gen, smp_frq[6, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_frq[7, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_frq[7, ], pop_frq[7, ]), max(smp_frq[7, ], pop_frq[7, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "A simulated dataset of the KM1KM1/sb1sb1 genotype generated with the Wright-Fisher diffusion")
points(smp_gen, smp_frq[7, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_frq[8, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_frq[8, ], pop_frq[8, ]), max(smp_frq[8, ], pop_frq[8, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "A simulated dataset of the KM1KM1/sb1SB1 genotype generated with the Wright-Fisher diffusion")
points(smp_gen, smp_frq[8, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_frq[9, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_frq[9, ], pop_frq[9, ]), max(smp_frq[9, ], pop_frq[9, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "A simulated dataset of the KM1KM1/SB1SB1 genotype generated with the Wright-Fisher diffusion")
points(smp_gen, smp_frq[9, ], col = 'red', pch = 17, cex = 1)

#' with uncertain genotypes
missing <- TRUE
sim_HMM_WFD <- cmpsimulateHMM(model, missing, sel_cof, rec_rat, pop_siz, int_frq, smp_gen, smp_siz, ptn_num)
smp_gen <- sim_HMM_WFD$smp_gen
smp_siz <- sim_HMM_WFD$smp_siz
smp_cnt <- sim_HMM_WFD$smp_cnt
mis_cnt <- sim_HMM_WFD$mis_cnt
smp_frq <- sim_HMM_WFD$smp_cnt / sim_HMM_WFD$smp_siz
pop_frq <- sim_HMM_WFD$pop_frq

k <- min(smp_gen):max(smp_gen)
plot(k, pop_frq[1, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_frq[1, ], pop_frq[1, ]), max(smp_frq[1, ], pop_frq[1, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "A simulated dataset of the KM0KM0/sb1sb1 genotype generated with the Wright-Fisher diffusion")
points(smp_gen, smp_frq[1, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_frq[2, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_frq[2, ], pop_frq[2, ]), max(smp_frq[2, ], pop_frq[2, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "A simulated dataset of the KM0KM0/sb1SB1 genotype generated with the Wright-Fisher diffusion")
points(smp_gen, smp_frq[2, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_frq[3, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_frq[3, ], pop_frq[3, ]), max(smp_frq[3, ], pop_frq[3, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "A simulated dataset of the KM0KM0/SB1SB1 genotype generated with the Wright-Fisher diffusion")
points(smp_gen, smp_frq[3, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_frq[4, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_frq[4, ], pop_frq[4, ]), max(smp_frq[4, ], pop_frq[4, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "A simulated dataset of the KM0KM1/sb1sb1 genotype generated with the Wright-Fisher diffusion")
points(smp_gen, smp_frq[4, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_frq[5, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_frq[5, ], pop_frq[5, ]), max(smp_frq[5, ], pop_frq[5, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "A simulated dataset of the KM0KM1/sb1SB1 genotype generated with the Wright-Fisher diffusion")
points(smp_gen, smp_frq[5, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_frq[6, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_frq[6, ], pop_frq[6, ]), max(smp_frq[6, ], pop_frq[6, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "A simulated dataset of the KM0KM1/SB1SB1 genotype generated with the Wright-Fisher diffusion")
points(smp_gen, smp_frq[6, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_frq[7, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_frq[7, ], pop_frq[7, ]), max(smp_frq[7, ], pop_frq[7, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "A simulated dataset of the KM1KM1/sb1sb1 genotype generated with the Wright-Fisher diffusion")
points(smp_gen, smp_frq[7, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_frq[8, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_frq[8, ], pop_frq[8, ]), max(smp_frq[8, ], pop_frq[8, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "A simulated dataset of the KM1KM1/sb1SB1 genotype generated with the Wright-Fisher diffusion")
points(smp_gen, smp_frq[8, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_frq[9, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_frq[9, ], pop_frq[9, ]), max(smp_frq[9, ], pop_frq[9, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "A simulated dataset of the KM1KM1/SB1SB1 genotype generated with the Wright-Fisher diffusion")
points(smp_gen, smp_frq[9, ], col = 'red', pch = 17, cex = 1)

################################################################################

#' Generate a simulated dataset without uncertain phenotypes under the Wright-Fisher model 
test_seed <- 48
set.seed(test_seed)

model <- "WFM"
missing <- FALSE
sel_cof <- c(5e-03, 1e-03, -5e-03)
rec_rat <- 5e-05
pop_siz <- 5e+03
int_frq <- rmultinom(1, size = pop_siz, prob = rep(1, 10) / sum(rep(1, 10))) / pop_siz
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
     file = "./Output/Output v1.0/GenoEvo/Horse coat pattern/TEST_withoutNA_simulated_dataset.rda")

load("./Output/Output v1.0/GenoEvo/Horse coat pattern/TEST_withoutNA_simulated_dataset.rda")

pdf(file = "./Output/Output v1.0/GenoEvo/Horse coat pattern/TEST_withoutNA_simulated_dataset.pdf", width = 30, height = 15)
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
title("A simulated dataset without uncertain genotypes generated with the Wright-Fisher model", outer = TRUE)
dev.off()

####################

#' Generate a simulated dataset with uncertain phenotypes under the Wright-Fisher model 
test_seed <- 48
set.seed(test_seed)

model <- "WFM"
missing <- TRUE
sel_cof <- c(5e-03, 1e-03, -5e-03)
rec_rat <- 5e-05
pop_siz <- 5e+03
int_frq <- rmultinom(1, size = pop_siz, prob = rep(1, 10) / sum(rep(1, 10))) / pop_siz
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
     file = "./Output/Output v1.0/GenoEvo/Horse coat pattern/TEST_withNA_simulated_dataset.rda")

load("./Output/Output v1.0/GenoEvo/Horse coat pattern/TEST_withNA_simulated_dataset.rda")

pdf(file = "./Output/Output v1.0/GenoEvo/Horse coat pattern/TEST_withNA_simulated_dataset.pdf", width = 30, height = 15)
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
title("A simulated dataset without uncertain genotypes generated with the Wright-Fisher model", outer = TRUE)
dev.off()

########################################

#' Run the bootstrap particle filter (BPF) with the two-locus Wright-Fisher diffusion with selection
#' Parameter setting
#' @param sel_cof the selection coefficients of the tobiano, sabino and mixed phenotypes
#' @param rec_rat the recombination rate between the KIT13 and KIT16 loci
#' @param pop_siz the number of the horses in the population
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_siz the count of the horses drawn from the population at all sampling time points
#' @param smp_cnt the count of the certain genotypes observed in the sample at all sampling time points
#' @param mis_cnt the count of the uncertain genotypes observed in the sample at all sampling time points
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter

#' without uncertain genotypes
load("./Output/Output v1.0/GenoEvo/Horse coat pattern/TEST_withoutNA_simulated_dataset.rda")

set.seed(test_seed)

sel_cof
rec_rat
pop_siz
smp_gen
smp_siz
smp_cnt
mis_cnt
ptn_num <- 5e+00
pcl_num <- 1e+06

system.time(BPF <- cmprunBPF(sel_cof, rec_rat, pop_siz, smp_gen, smp_siz, smp_cnt, mis_cnt, ptn_num, pcl_num))

save(sel_cof, rec_rat, pop_siz, smp_gen, smp_siz, smp_cnt, mis_cnt, smp_frq, pop_frq, ptn_num, pcl_num, BPF,
     file = "./Output/Output v1.0/GenoEvo/Horse coat pattern/TEST_withoutNA_BPF.rda")

load("./Output/Output v1.0/GenoEvo/Horse coat pattern/TEST_withoutNA_BPF.rda")

lik <- rep(1, pcl_num)
wght <- BPF$wght
for (k in 1:length(smp_gen)) {
  lik <- lik * (cumsum(wght[, k]) / (1:pcl_num))
}

pdf(file = "./Output/Output v1.0/GenoEvo/Horse coat pattern/TEST_withoutNA_BPF_likelihood.pdf", width = 20, height = 10)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
plot(1:pcl_num, log(lik), type = 'l', 
     xlab = "Number of particles", ylab = "Log likelihood", main = "Log likelihood estimated with the bootstrap particle filter")
dev.off()

pop_frq_pre_resmp <- BPF$pop_frq_pre_resmp
pop_frq_pst_resmp <- BPF$pop_frq_pst_resmp

pdf(file = "./Output/Output v1.0/GenoEvo/Horse coat pattern/TEST_withoutNA_BPF_particle.pdf", width = 50, height = 50)
par(mfrow = c(11, 9), oma = c(0, 0, 3, 0), mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
for (k in 1:length(smp_gen)) {
  hist(pop_frq_pst_resmp[1, , k], breaks = seq(min(pop_frq_pst_resmp[1, , k], pop_frq_pre_resmp[1, , k]), max(pop_frq_pst_resmp[1, , k], pop_frq_pre_resmp[1, , k]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
       xlim = c(min(pop_frq_pst_resmp[1, , k], pop_frq_pre_resmp[1, , k], smp_frq[1, k]), max(pop_frq_pst_resmp[1, , k], pop_frq_pre_resmp[1, , k], smp_frq[1, k])),
       xlab = "Genotype frequency", main = paste("Genotype KM0KM0/sb1sb1 in generation", smp_gen[k]))
  hist(pop_frq_pre_resmp[1, , k], breaks = seq(min(pop_frq_pst_resmp[1, , k], pop_frq_pre_resmp[1, , k]), max(pop_frq_pst_resmp[1, , k], pop_frq_pre_resmp[1, , k]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
  abline(v = smp_frq[1, k], col = 'red', lty = 2, lwd = 2)
  hist(pop_frq_pst_resmp[2, , k], breaks = seq(min(pop_frq_pst_resmp[2, , k], pop_frq_pre_resmp[2, , k]), max(pop_frq_pst_resmp[2, , k], pop_frq_pre_resmp[2, , k]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
       xlim = c(min(pop_frq_pst_resmp[2, , k], pop_frq_pre_resmp[2, , k], smp_frq[2, k]), max(pop_frq_pst_resmp[2, , k], pop_frq_pre_resmp[2, , k], smp_frq[2, k])),
       xlab = "Genotype frequency", main = paste("Genotype KM0KM0/sb1SB1 in generation", smp_gen[k]))
  hist(pop_frq_pre_resmp[2, , k], breaks = seq(min(pop_frq_pst_resmp[2, , k], pop_frq_pre_resmp[2, , k]), max(pop_frq_pst_resmp[2, , k], pop_frq_pre_resmp[2, , k]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
  abline(v = smp_frq[2, k], col = 'red', lty = 2, lwd = 2)
  hist(pop_frq_pst_resmp[3, , k], breaks = seq(min(pop_frq_pst_resmp[3, , k], pop_frq_pre_resmp[3, , k]), max(pop_frq_pst_resmp[3, , k], pop_frq_pre_resmp[3, , k]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
       xlim = c(min(pop_frq_pst_resmp[3, , k], pop_frq_pre_resmp[3, , k], smp_frq[3, k]), max(pop_frq_pst_resmp[3, , k], pop_frq_pre_resmp[3, , k], smp_frq[3, k])),
       xlab = "Genotype frequency", main = paste("Genotype KM0KM0/SB1SB1 in generation", smp_gen[k]))
  hist(pop_frq_pre_resmp[3, , k], breaks = seq(min(pop_frq_pst_resmp[3, , k], pop_frq_pre_resmp[3, , k]), max(pop_frq_pst_resmp[3, , k], pop_frq_pre_resmp[3, , k]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
  abline(v = smp_frq[3, k], col = 'red', lty = 2, lwd = 2)
  hist(pop_frq_pst_resmp[4, , k], breaks = seq(min(pop_frq_pst_resmp[4, , k], pop_frq_pre_resmp[4, , k]), max(pop_frq_pst_resmp[4, , k], pop_frq_pre_resmp[4, , k]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
       xlim = c(min(pop_frq_pst_resmp[4, , k], pop_frq_pre_resmp[4, , k], smp_frq[4, k]), max(pop_frq_pst_resmp[4, , k], pop_frq_pre_resmp[4, , k], smp_frq[4, k])),
       xlab = "Genotype frequency", main = paste("Genotype KM0KM1/sb1sb1 in generation", smp_gen[k]))
  hist(pop_frq_pre_resmp[4, , k], breaks = seq(min(pop_frq_pst_resmp[4, , k], pop_frq_pre_resmp[4, , k]), max(pop_frq_pst_resmp[4, , k], pop_frq_pre_resmp[4, , k]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
  abline(v = smp_frq[4, k], col = 'red', lty = 2, lwd = 2)
  hist(pop_frq_pst_resmp[5, , k], breaks = seq(min(pop_frq_pst_resmp[5, , k], pop_frq_pre_resmp[5, , k]), max(pop_frq_pst_resmp[5, , k], pop_frq_pre_resmp[5, , k]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
       xlim = c(min(pop_frq_pst_resmp[5, , k], pop_frq_pre_resmp[5, , k], smp_frq[5, k]), max(pop_frq_pst_resmp[5, , k], pop_frq_pre_resmp[5, , k], smp_frq[5, k])),
       xlab = "Genotype frequency", main = paste("Genotype KM0KM1/sb1SB1 in generation", smp_gen[k]))
  hist(pop_frq_pre_resmp[5, , k], breaks = seq(min(pop_frq_pst_resmp[5, , k], pop_frq_pre_resmp[5, , k]), max(pop_frq_pst_resmp[5, , k], pop_frq_pre_resmp[5, , k]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
  abline(v = smp_frq[5, k], col = 'red', lty = 2, lwd = 2)
  hist(pop_frq_pst_resmp[6, , k], breaks = seq(min(pop_frq_pst_resmp[6, , k], pop_frq_pre_resmp[6, , k]), max(pop_frq_pst_resmp[6, , k], pop_frq_pre_resmp[6, , k]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
       xlim = c(min(pop_frq_pst_resmp[6, , k], pop_frq_pre_resmp[6, , k], smp_frq[6, k]), max(pop_frq_pst_resmp[6, , k], pop_frq_pre_resmp[6, , k], smp_frq[6, k])),
       xlab = "Genotype frequency", main = paste("Genotype KM0KM1/SB1SB1 in generation", smp_gen[k]))
  hist(pop_frq_pre_resmp[6, , k], breaks = seq(min(pop_frq_pst_resmp[6, , k], pop_frq_pre_resmp[6, , k]), max(pop_frq_pst_resmp[6, , k], pop_frq_pre_resmp[6, , k]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
  abline(v = smp_frq[6, k], col = 'red', lty = 2, lwd = 2)
  hist(pop_frq_pst_resmp[7, , k], breaks = seq(min(pop_frq_pst_resmp[7, , k], pop_frq_pre_resmp[7, , k]), max(pop_frq_pst_resmp[7, , k], pop_frq_pre_resmp[7, , k]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
       xlim = c(min(pop_frq_pst_resmp[7, , k], pop_frq_pre_resmp[7, , k], smp_frq[7, k]), max(pop_frq_pst_resmp[7, , k], pop_frq_pre_resmp[7, , k], smp_frq[7, k])),
       xlab = "Genotype frequency", main = paste("Genotype KM1KM1/sb1sb1 in generation", smp_gen[k]))
  hist(pop_frq_pre_resmp[7, , k], breaks = seq(min(pop_frq_pst_resmp[7, , k], pop_frq_pre_resmp[7, , k]), max(pop_frq_pst_resmp[7, , k], pop_frq_pre_resmp[7, , k]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
  abline(v = smp_frq[7, k], col = 'red', lty = 2, lwd = 2)
  hist(pop_frq_pst_resmp[8, , k], breaks = seq(min(pop_frq_pst_resmp[8, , k], pop_frq_pre_resmp[8, , k]), max(pop_frq_pst_resmp[8, , k], pop_frq_pre_resmp[8, , k]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
       xlim = c(min(pop_frq_pst_resmp[8, , k], pop_frq_pre_resmp[8, , k], smp_frq[8, k]), max(pop_frq_pst_resmp[8, , k], pop_frq_pre_resmp[8, , k], smp_frq[8, k])),
       xlab = "Genotype frequency", main = paste("Genotype KM1KM1/sb1SB1 in generation", smp_gen[k]))
  hist(pop_frq_pre_resmp[8, , k], breaks = seq(min(pop_frq_pst_resmp[8, , k], pop_frq_pre_resmp[8, , k]), max(pop_frq_pst_resmp[8, , k], pop_frq_pre_resmp[8, , k]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
  abline(v = smp_frq[8, k], col = 'red', lty = 2, lwd = 2)
  hist(pop_frq_pst_resmp[9, , k], breaks = seq(min(pop_frq_pst_resmp[9, , k], pop_frq_pre_resmp[9, , k]), max(pop_frq_pst_resmp[9, , k], pop_frq_pre_resmp[9, , k]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
       xlim = c(min(pop_frq_pst_resmp[9, , k], pop_frq_pre_resmp[9, , k], smp_frq[9, k]), max(pop_frq_pst_resmp[9, , k], pop_frq_pre_resmp[9, , k], smp_frq[9, k])),
       xlab = "Genotype frequency", main = paste("Genotype KM1KM1/SB1SB1 in generation", smp_gen[k]))
  hist(pop_frq_pre_resmp[9, , k], breaks = seq(min(pop_frq_pst_resmp[9, , k], pop_frq_pre_resmp[9, , k]), max(pop_frq_pst_resmp[9, , k], pop_frq_pre_resmp[9, , k]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
  abline(v = smp_frq[9, k], col = 'red', lty = 2, lwd = 2)
}
title("Histograms of the particles before and after resampling", outer = TRUE)
dev.off()

####################

#' with uncertain genotypes
load("./Output/Output v1.0/GenoEvo/Horse coat pattern/TEST_withNA_simulated_dataset.rda")

set.seed(test_seed)

sel_cof
rec_rat
pop_siz
smp_gen
smp_siz
smp_cnt
mis_cnt
ptn_num <- 5e+00
pcl_num <- 1e+06

system.time(BPF <- cmprunBPF(sel_cof, rec_rat, pop_siz, smp_gen, smp_siz, smp_cnt, mis_cnt, ptn_num, pcl_num))

save(sel_cof, rec_rat, pop_siz, smp_gen, smp_siz, smp_cnt, mis_cnt, smp_frq, pop_frq, ptn_num, pcl_num, BPF,
     file = "./Output/Output v1.0/GenoEvo/Horse coat pattern/TEST_withNA_BPF.rda")

load("./Output/Output v1.0/GenoEvo/Horse coat pattern/TEST_withNA_BPF.rda")

lik <- rep(1, pcl_num)
wght <- BPF$wght
for (k in 1:length(smp_gen)) {
  lik <- lik * (cumsum(wght[, k]) / (1:pcl_num))
}

pdf(file = "./Output/Output v1.0/GenoEvo/Horse coat pattern/TEST_withNA_BPF_likelihood.pdf", width = 20, height = 10)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
plot(1:pcl_num, log(lik), type = 'l', 
     xlab = "Number of particles", ylab = "Log likelihood", main = "Log likelihood estimated with the bootstrap particle filter")
dev.off()

pop_frq_pre_resmp <- BPF$pop_frq_pre_resmp
pop_frq_pst_resmp <- BPF$pop_frq_pst_resmp

pdf(file = "./Output/Output v1.0/GenoEvo/Horse coat pattern/TEST_withNA_BPF_particle.pdf", width = 50, height = 50)
par(mfrow = c(11, 9), oma = c(0, 0, 3, 0), mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
for (k in 1:length(smp_gen)) {
  hist(pop_frq_pst_resmp[1, , k], breaks = seq(min(pop_frq_pst_resmp[1, , k], pop_frq_pre_resmp[1, , k]), max(pop_frq_pst_resmp[1, , k], pop_frq_pre_resmp[1, , k]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
       xlim = c(min(pop_frq_pst_resmp[1, , k], pop_frq_pre_resmp[1, , k], smp_frq[1, k]), max(pop_frq_pst_resmp[1, , k], pop_frq_pre_resmp[1, , k], smp_frq[1, k])),
       xlab = "Genotype frequency", main = paste("Genotype KM0KM0/sb1sb1 in generation", smp_gen[k]))
  hist(pop_frq_pre_resmp[1, , k], breaks = seq(min(pop_frq_pst_resmp[1, , k], pop_frq_pre_resmp[1, , k]), max(pop_frq_pst_resmp[1, , k], pop_frq_pre_resmp[1, , k]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
  abline(v = smp_frq[1, k], col = 'red', lty = 2, lwd = 2)
  hist(pop_frq_pst_resmp[2, , k], breaks = seq(min(pop_frq_pst_resmp[2, , k], pop_frq_pre_resmp[2, , k]), max(pop_frq_pst_resmp[2, , k], pop_frq_pre_resmp[2, , k]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
       xlim = c(min(pop_frq_pst_resmp[2, , k], pop_frq_pre_resmp[2, , k], smp_frq[2, k]), max(pop_frq_pst_resmp[2, , k], pop_frq_pre_resmp[2, , k], smp_frq[2, k])),
       xlab = "Genotype frequency", main = paste("Genotype KM0KM0/sb1SB1 in generation", smp_gen[k]))
  hist(pop_frq_pre_resmp[2, , k], breaks = seq(min(pop_frq_pst_resmp[2, , k], pop_frq_pre_resmp[2, , k]), max(pop_frq_pst_resmp[2, , k], pop_frq_pre_resmp[2, , k]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
  abline(v = smp_frq[2, k], col = 'red', lty = 2, lwd = 2)
  hist(pop_frq_pst_resmp[3, , k], breaks = seq(min(pop_frq_pst_resmp[3, , k], pop_frq_pre_resmp[3, , k]), max(pop_frq_pst_resmp[3, , k], pop_frq_pre_resmp[3, , k]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
       xlim = c(min(pop_frq_pst_resmp[3, , k], pop_frq_pre_resmp[3, , k], smp_frq[3, k]), max(pop_frq_pst_resmp[3, , k], pop_frq_pre_resmp[3, , k], smp_frq[3, k])),
       xlab = "Genotype frequency", main = paste("Genotype KM0KM0/SB1SB1 in generation", smp_gen[k]))
  hist(pop_frq_pre_resmp[3, , k], breaks = seq(min(pop_frq_pst_resmp[3, , k], pop_frq_pre_resmp[3, , k]), max(pop_frq_pst_resmp[3, , k], pop_frq_pre_resmp[3, , k]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
  abline(v = smp_frq[3, k], col = 'red', lty = 2, lwd = 2)
  hist(pop_frq_pst_resmp[4, , k], breaks = seq(min(pop_frq_pst_resmp[4, , k], pop_frq_pre_resmp[4, , k]), max(pop_frq_pst_resmp[4, , k], pop_frq_pre_resmp[4, , k]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
       xlim = c(min(pop_frq_pst_resmp[4, , k], pop_frq_pre_resmp[4, , k], smp_frq[4, k]), max(pop_frq_pst_resmp[4, , k], pop_frq_pre_resmp[4, , k], smp_frq[4, k])),
       xlab = "Genotype frequency", main = paste("Genotype KM0KM1/sb1sb1 in generation", smp_gen[k]))
  hist(pop_frq_pre_resmp[4, , k], breaks = seq(min(pop_frq_pst_resmp[4, , k], pop_frq_pre_resmp[4, , k]), max(pop_frq_pst_resmp[4, , k], pop_frq_pre_resmp[4, , k]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
  abline(v = smp_frq[4, k], col = 'red', lty = 2, lwd = 2)
  hist(pop_frq_pst_resmp[5, , k], breaks = seq(min(pop_frq_pst_resmp[5, , k], pop_frq_pre_resmp[5, , k]), max(pop_frq_pst_resmp[5, , k], pop_frq_pre_resmp[5, , k]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
       xlim = c(min(pop_frq_pst_resmp[5, , k], pop_frq_pre_resmp[5, , k], smp_frq[5, k]), max(pop_frq_pst_resmp[5, , k], pop_frq_pre_resmp[5, , k], smp_frq[5, k])),
       xlab = "Genotype frequency", main = paste("Genotype KM0KM1/sb1SB1 in generation", smp_gen[k]))
  hist(pop_frq_pre_resmp[5, , k], breaks = seq(min(pop_frq_pst_resmp[5, , k], pop_frq_pre_resmp[5, , k]), max(pop_frq_pst_resmp[5, , k], pop_frq_pre_resmp[5, , k]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
  abline(v = smp_frq[5, k], col = 'red', lty = 2, lwd = 2)
  hist(pop_frq_pst_resmp[6, , k], breaks = seq(min(pop_frq_pst_resmp[6, , k], pop_frq_pre_resmp[6, , k]), max(pop_frq_pst_resmp[6, , k], pop_frq_pre_resmp[6, , k]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
       xlim = c(min(pop_frq_pst_resmp[6, , k], pop_frq_pre_resmp[6, , k], smp_frq[6, k]), max(pop_frq_pst_resmp[6, , k], pop_frq_pre_resmp[6, , k], smp_frq[6, k])),
       xlab = "Genotype frequency", main = paste("Genotype KM0KM1/SB1SB1 in generation", smp_gen[k]))
  hist(pop_frq_pre_resmp[6, , k], breaks = seq(min(pop_frq_pst_resmp[6, , k], pop_frq_pre_resmp[6, , k]), max(pop_frq_pst_resmp[6, , k], pop_frq_pre_resmp[6, , k]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
  abline(v = smp_frq[6, k], col = 'red', lty = 2, lwd = 2)
  hist(pop_frq_pst_resmp[7, , k], breaks = seq(min(pop_frq_pst_resmp[7, , k], pop_frq_pre_resmp[7, , k]), max(pop_frq_pst_resmp[7, , k], pop_frq_pre_resmp[7, , k]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
       xlim = c(min(pop_frq_pst_resmp[7, , k], pop_frq_pre_resmp[7, , k], smp_frq[7, k]), max(pop_frq_pst_resmp[7, , k], pop_frq_pre_resmp[7, , k], smp_frq[7, k])),
       xlab = "Genotype frequency", main = paste("Genotype KM1KM1/sb1sb1 in generation", smp_gen[k]))
  hist(pop_frq_pre_resmp[7, , k], breaks = seq(min(pop_frq_pst_resmp[7, , k], pop_frq_pre_resmp[7, , k]), max(pop_frq_pst_resmp[7, , k], pop_frq_pre_resmp[7, , k]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
  abline(v = smp_frq[7, k], col = 'red', lty = 2, lwd = 2)
  hist(pop_frq_pst_resmp[8, , k], breaks = seq(min(pop_frq_pst_resmp[8, , k], pop_frq_pre_resmp[8, , k]), max(pop_frq_pst_resmp[8, , k], pop_frq_pre_resmp[8, , k]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
       xlim = c(min(pop_frq_pst_resmp[8, , k], pop_frq_pre_resmp[8, , k], smp_frq[8, k]), max(pop_frq_pst_resmp[8, , k], pop_frq_pre_resmp[8, , k], smp_frq[8, k])),
       xlab = "Genotype frequency", main = paste("Genotype KM1KM1/sb1SB1 in generation", smp_gen[k]))
  hist(pop_frq_pre_resmp[8, , k], breaks = seq(min(pop_frq_pst_resmp[8, , k], pop_frq_pre_resmp[8, , k]), max(pop_frq_pst_resmp[8, , k], pop_frq_pre_resmp[8, , k]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
  abline(v = smp_frq[8, k], col = 'red', lty = 2, lwd = 2)
  hist(pop_frq_pst_resmp[9, , k], breaks = seq(min(pop_frq_pst_resmp[9, , k], pop_frq_pre_resmp[9, , k]), max(pop_frq_pst_resmp[9, , k], pop_frq_pre_resmp[9, , k]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
       xlim = c(min(pop_frq_pst_resmp[9, , k], pop_frq_pre_resmp[9, , k], smp_frq[9, k]), max(pop_frq_pst_resmp[9, , k], pop_frq_pre_resmp[9, , k], smp_frq[9, k])),
       xlab = "Genotype frequency", main = paste("Genotype KM1KM1/SB1SB1 in generation", smp_gen[k]))
  hist(pop_frq_pre_resmp[9, , k], breaks = seq(min(pop_frq_pst_resmp[9, , k], pop_frq_pre_resmp[9, , k]), max(pop_frq_pst_resmp[9, , k], pop_frq_pre_resmp[9, , k]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
  abline(v = smp_frq[9, k], col = 'red', lty = 2, lwd = 2)
}
title("Histograms of the particles before and after resampling", outer = TRUE)
dev.off()

########################################

#' Calculate the optimal particle number in the particle marginal Metropolis-Hastings (PMMH)
#' Parameter settings
#' @param sel_cof the selection coefficients of the tobiano, sabino and mixed phenotypes
#' @param rec_rat the recombination rate between the KIT13 and KIT16 loci
#' @param pop_siz the number of the horses in the population
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_siz the count of the horses drawn from the population at all sampling time points
#' @param smp_cnt the count of the certain genotypes observed in the sample at all sampling time points
#' @param mis_cnt the count of the uncertain genotypes observed in the sample at all sampling time points
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the particle marginal Metropolis-Hastings
#' @param gap_num the number of particles increased or decreased in the optimal particle number search

#' without uncertain genotypes
load("./Output/Output v1.0/GenoEvo/Horse coat pattern/TEST_withoutNA_simulated_dataset.rda")

set.seed(test_seed)

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
     file = "./Output/Output v1.0/GenoEvo/Horse coat pattern/TEST_withoutNA_OptNum.rda")

load("./Output/Output v1.0/GenoEvo/Horse coat pattern/TEST_withoutNA_OptNum.rda")

opt_pcl_num <- OptNum$opt_pcl_num
log_lik_sdv <- OptNum$log_lik_sdv

pdf(file = "./Output/Output v1.0/GenoEvo/Horse coat pattern/TEST_withoutNA_OptNum.pdf", width = 10, height = 10)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
plot(opt_pcl_num, log_lik_sdv, type = 'b', lwd = 2, 
     xlab = "Particle number", ylab = "Log-likelihood standard deviation", main = "Optimal particle number in the PMMH")
abline(h = 1.7, col = 'red', lty = 2, lwd = 2)
abline(h = 1.0, col = 'red', lty = 2, lwd = 2)
dev.off()

####################

#' with uncertain genotypes
load("./Output/Output v1.0/GenoEvo/Horse coat pattern/TEST_withNA_simulated_dataset.rda")

set.seed(test_seed)

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
     file = "./Output/Output v1.0/GenoEvo/Horse coat pattern/TEST_withNA_OptNum.rda")

load("./Output/Output v1.0/GenoEvo/Horse coat pattern/TEST_withNA_OptNum.rda")

opt_pcl_num <- OptNum$opt_pcl_num
log_lik_sdv <- OptNum$log_lik_sdv

pdf(file = "./Output/Output v1.0/GenoEvo/Horse coat pattern/TEST_withNA_OptNum.pdf", width = 10, height = 10)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
plot(opt_pcl_num, log_lik_sdv, type = 'b', lwd = 2, 
     xlab = "Particle number", ylab = "Log-likelihood standard deviation", main = "Optimal particle number in the PMMH")
abline(h = 1.7, col = 'red', lty = 2, lwd = 2)
abline(h = 1.0, col = 'red', lty = 2, lwd = 2)
dev.off()

########################################

#' Run the particle marginal Metropolis-Hastings (PMMH)
#' Parameter settings
#' @param sel_cof the initial values of the selection coefficients of the tobiano, sabino and mixed phenotypes
#' @param rec_rat the recombination rate between the KIT13 and KIT16 loci
#' @param pop_siz the number of the horses in the population
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_siz the count of the horses drawn from the population at all sampling time points
#' @param smp_cnt the count of the certain genotypes observed in the sample at all sampling time points
#' @param mis_cnt the count of the uncertain genotypes observed in the sample at all sampling time points
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter
#' @param itn_num the number of the iterations carried out in the particle marginal Metropolis-Hastings
#' @param ER = TRUE/FALSE (use the early rejection technique or not in the particle marginal Metropolis-Hastings)
#' @param PA = TRUE/FALSE (use the parallel adaptive technique or not in the particle marginal Metropolis-Hastings)
#' @param nap_num the number of the iterations carried out as a non-adaptation period in the parallel adaptive particle marginal Metropolis-Hastings

#' without uncertain genotypes
load("./Output/Output v1.0/GenoEvo/Horse coat pattern/TEST_withoutNA_simulated_dataset.rda")

set.seed(test_seed)

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
PA <- TRUE
nap_num <- itn_num * 0.1

system.time(PMMH <- cmprunPMMH(int_sel_cof, rec_rat, pop_siz, smp_gen, smp_siz, smp_cnt, mis_cnt, ptn_num, pcl_num, itn_num, ER, PA, nap_num))

save(sel_cof, rec_rat, pop_siz, smp_gen, smp_siz, smp_cnt, mis_cnt, smp_frq, pop_frq, ptn_num, pcl_num, itn_num, ER, PA, nap_num, PMMH,
     file = "./Output/Output v1.0/GenoEvo/Horse coat pattern/TEST_withoutNA_PMMH.rda")

load("./Output/Output v1.0/GenoEvo/Horse coat pattern/TEST_withoutNA_PMMH.rda")

sel_cof_to_chn <- PMMH$sel_cof_to_chn
sel_cof_sb_chn <- PMMH$sel_cof_sb_chn
sel_cof_ts_chn <- PMMH$sel_cof_ts_chn

pdf(file = "./Output/Output v1.0/GenoEvo/Horse coat pattern/TEST_withoutNA_PMMH_traceplot.pdf", width = 10, height = 15)
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

pdf(file = "./Output/Output v1.0/GenoEvo/Horse coat pattern/TEST_withoutNA_PMMH_posterior.pdf", width = 20, height = 20)
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

####################

#' with uncertain genotypes
load("./Output/Output v1.0/GenoEvo/Horse coat pattern/TEST_withNA_simulated_dataset.rda")

set.seed(test_seed)

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
PA <- TRUE
nap_num <- itn_num * 0.1

system.time(PMMH <- cmprunPMMH(int_sel_cof, rec_rat, pop_siz, smp_gen, smp_siz, smp_cnt, mis_cnt, ptn_num, pcl_num, itn_num, ER, PA, nap_num))

save(sel_cof, rec_rat, pop_siz, smp_gen, smp_siz, smp_cnt, mis_cnt, smp_frq, pop_frq, ptn_num, pcl_num, itn_num, ER, PA, nap_num, PMMH,
     file = "./Output/Output v1.0/GenoEvo/Horse coat pattern/TEST_withNA_PMMH.rda")

load("./Output/Output v1.0/GenoEvo/Horse coat pattern/TEST_withNA_PMMH.rda")

sel_cof_to_chn <- PMMH$sel_cof_to_chn
sel_cof_sb_chn <- PMMH$sel_cof_sb_chn
sel_cof_ts_chn <- PMMH$sel_cof_ts_chn

pdf(file = "./Output/Output v1.0/GenoEvo/Horse coat pattern/TEST_withNA_PMMH_traceplot.pdf", width = 10, height = 15)
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

pdf(file = "./Output/Output v1.0/GenoEvo/Horse coat pattern/TEST_withNA_PMMH_posterior.pdf", width = 20, height = 20)
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

########################################

#' Run the Bayesian procedure for the inference of natural selection
#' Parameter settings
#' @param sel_cof the initial values of the selection coefficients of the tobiano, sabino and mixed phenotypes
#' @param rec_rat the recombination rate between the KIT13 and KIT16 loci
#' @param pop_siz the number of the horses in the population
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_siz the count of the horses drawn from the population at all sampling time points
#' @param smp_cnt the count of the certain genotypes observed in the sample at all sampling time points
#' @param mis_cnt the count of the uncertain genotypes observed in the sample at all sampling time points
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of the particles generated in the bootstrap particle filter
#' @param itn_num the number of the iterations carried out in the particle marginal Metropolis-Hastings
#' @param ER = TRUE/FALSE (use the early rejection technique or not in the particle marginal Metropolis-Hastings)
#' @param PA = TRUE/FALSE (use the parallel adaptive technique or not in the particle marginal Metropolis-Hastings)
#' @param brn_num the number of the iterations for burn-in
#' @param thn_num the number of the iterations for thinning
#' @param est = the Bayesian estimator (i.e., MAP estimator or MMSE estimator)
#' @param nap_num the number of the iterations carried out as a non-adaptation period in the parallel adaptive particle marginal Metropolis-Hastings
#' @param grd_num the number of the grids in the kernel density estimation

#' without uncertain genotypes
load("./Output/Output v1.0/GenoEvo/Horse coat pattern/TEST_withoutNA_simulated_dataset.rda")

set.seed(test_seed)

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
PA <- TRUE
brn_num <- itn_num * 0.2
thn_num <- 4e+00
est <- "MMSE"
nap_num <- itn_num * 0.1
grd_num <- 1e+03

system.time(BayesianProcedure <- cmprunBayesianProcedure(int_sel_cof, rec_rat, pop_siz, smp_gen, smp_siz, smp_cnt, mis_cnt, ptn_num, pcl_num, itn_num, ER, PA, brn_num, thn_num, est, nap_num, grd_num))

save(sel_cof, rec_rat, pop_siz, smp_gen, smp_siz, smp_cnt, mis_cnt, ptn_num, pcl_num, itn_num, ER, PA, brn_num, thn_num, est, nap_num, grd_num, BayesianProcedure,
     file = "./Output/Output v1.0/GenoEvo/Horse coat pattern/TEST_withoutNA_BayesianProcedure.rda")

load("./Output/Output v1.0/GenoEvo/Horse coat pattern/TEST_withoutNA_BayesianProcedure.rda")

sel_cof_to_chn <- BayesianProcedure$sel_cof_to_chn
sel_cof_sb_chn <- BayesianProcedure$sel_cof_sb_chn
sel_cof_ts_chn <- BayesianProcedure$sel_cof_ts_chn

sel_cof_to_mmse <- BayesianProcedure$sel_cof_to_mmse
sel_cof_sb_mmse <- BayesianProcedure$sel_cof_sb_mmse
sel_cof_ts_mmse <- BayesianProcedure$sel_cof_ts_mmse

sel_cof_to_hpd <- BayesianProcedure$sel_cof_to_hpd
sel_cof_sb_hpd <- BayesianProcedure$sel_cof_sb_hpd
sel_cof_ts_hpd <- BayesianProcedure$sel_cof_ts_hpd

pdf(file = "./Output/Output v1.0/GenoEvo/Horse coat pattern/TEST_withoutNA_BayesianProcedure_posterior.pdf", width = 10, height = 15)
par(mfrow = c(3, 1), mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
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

####################

#' with uncertain genotypes
load("./Output/Output v1.0/GenoEvo/Horse coat pattern/TEST_withNA_simulated_dataset.rda")

set.seed(test_seed)

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
PA <- TRUE
brn_num <- itn_num * 0.2
thn_num <- 4e+00
est <- "MMSE"
nap_num <- itn_num * 0.1
grd_num <- 1e+03

system.time(BayesianProcedure <- cmprunBayesianProcedure(int_sel_cof, rec_rat, pop_siz, smp_gen, smp_siz, smp_cnt, mis_cnt, ptn_num, pcl_num, itn_num, ER, PA, brn_num, thn_num, est, nap_num, grd_num))

save(sel_cof, rec_rat, pop_siz, smp_gen, smp_siz, smp_cnt, mis_cnt, ptn_num, pcl_num, itn_num, ER, PA, brn_num, thn_num, est, nap_num, grd_num, BayesianProcedure,
     file = "./Output/Output v1.0/GenoEvo/Horse coat pattern/TEST_withNA_BayesianProcedure.rda")

load("./Output/Output v1.0/GenoEvo/Horse coat pattern/TEST_withNA_BayesianProcedure.rda")

sel_cof_to_chn <- BayesianProcedure$sel_cof_to_chn
sel_cof_sb_chn <- BayesianProcedure$sel_cof_sb_chn
sel_cof_ts_chn <- BayesianProcedure$sel_cof_ts_chn

sel_cof_to_mmse <- BayesianProcedure$sel_cof_to_mmse
sel_cof_sb_mmse <- BayesianProcedure$sel_cof_sb_mmse
sel_cof_ts_mmse <- BayesianProcedure$sel_cof_ts_mmse

sel_cof_to_hpd <- BayesianProcedure$sel_cof_to_hpd
sel_cof_sb_hpd <- BayesianProcedure$sel_cof_sb_hpd
sel_cof_ts_hpd <- BayesianProcedure$sel_cof_ts_hpd

pdf(file = "./Output/Output v1.0/GenoEvo/Horse coat pattern/TEST_withNA_BayesianProcedure_posterior.pdf", width = 10, height = 15)
par(mfrow = c(3, 1), mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
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

################################################################################
