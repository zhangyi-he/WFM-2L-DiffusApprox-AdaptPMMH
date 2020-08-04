#' @title Inferring natural selection acting on horse coat colours and patterns during the process of domestication
#' @author Zhangyi He, Xiaoyang Dai, Mark Beaumont and Feng Yu

#' version 1.0

#' R functions

#install.packages("MASS")
library("MASS")

#install.packages("coda")
library("coda")

#install.packages("inline")
library("inline")
#install.packages("Rcpp")
library("Rcpp")
#install.packages("RcppArmadillo")
library("RcppArmadillo")

#install.packages("compiler")
library("compiler")
#enableJIT(1)

# call C++ functions
sourceCpp("./CFUN.cpp")

################################################################################

#' Simulate the haplotype frequency trajectories according to the Wright-Fisher model
#' parameter settings
#' @param s_b the selection coefficient of the black phenotype
#' @param s_c the selection coefficient of the chestnut phenotype
#' @param r the recombination rate between ASIP and MC1R (r_bu = 0.5 is predetermined since ASIP and MC1R are located on a separate chromosome)
#' @param N the number of individuals in the population
#' @param int_frq the initial haplotype frequencies of the population
#' @param int_gen the first generation of the simulated haplotype frequency trajectories
#' @param lst_gen the last generation of the simulated haplotype frequency trajectories

#' Standard version
simulateTLWFMS <- function(s_b, s_c, r, N, int_frq, int_gen, lst_gen) {
  fts_mat <- calculateFitnessMat_arma(s_b, s_c)
  # print(fts_mat)

  frq_pth <- simulateTLWFMS_arma(fts_mat, r, N, int_frq, int_gen, lst_gen)

  return(frq_pth)
}
#' Compiled version
cmpsimulateTLWFMS <- cmpfun(simulateTLWFMS)

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

#' Standard version
simulateTLWFDS <- function(s_b, s_c, r, N, int_frq, int_gen, lst_gen, ptn_num, data_augmentation = TRUE) {
  alpha_b <- 2 * N * s_b
  alpha_c <- 2 * N * s_c
  rho <- 4 * N * r

  frq_pth <- simulateTLWFDS_arma(alpha_b, alpha_c, rho, N, int_frq, int_gen, lst_gen, ptn_num)

  if (data_augmentation == FALSE) {
    return(frq_pth[, (0:(lst_gen - int_gen)) * ptn_num + 1])
  } else {
    return(frq_pth)
  }
}
#' Compiled version
cmpsimulateTLWFDS <- cmpfun(simulateTLWFDS)

########################################

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

#' Standard version
simulateHMM <- function(model, s_b, s_c, r, N, int_frq, smp_gen, smp_siz, ...) {
  int_gen <- min(smp_gen)
  lst_gen <- max(smp_gen)

  # generate the population haplotype frequency trajectories and the population allele frequency trajectories
  if (model == "WFM") {
    pop_hap_frq <- cmpsimulateTLWFMS(s_b, s_c, r, N, int_frq, int_gen, lst_gen)
  }
  if (model == "WFD") {
    pop_hap_frq <- cmpsimulateTLWFDS(s_b, s_c, r, N, int_frq, int_gen, lst_gen, ptn_num, data_augmentation = FALSE)
  }
  pop_ale_frq <- matrix(NA, nrow = 2, ncol = (lst_gen - int_gen) + 1)
  pop_ale_frq[1, ] <- pop_hap_frq[1, ] + pop_hap_frq[2, ]
  pop_ale_frq[2, ] <- pop_hap_frq[1, ] + pop_hap_frq[3, ]

  # generate the sample haplotype counts and the sample allele counts at all sampling time points
  smp_hap_cnt <- matrix(NA, nrow = 4, ncol = length(smp_gen))
  smp_ale_cnt <- matrix(NA, nrow = 2, ncol = length(smp_gen))
  for (k in 1:length(smp_gen)) {
    smp_hap_cnt[, k] <- rmultinom(1, size = smp_siz[k], prob = pop_hap_frq[, smp_gen[k] - int_gen + 1])
    smp_ale_cnt[1, k] <- smp_hap_cnt[1, k] + smp_hap_cnt[2, k]
    smp_ale_cnt[2, k] <- smp_hap_cnt[1, k] + smp_hap_cnt[3, k]
  }

  return(list(smp_gen = smp_gen,
              smp_siz = smp_siz,
              smp_hap_cnt = smp_hap_cnt,
              smp_ale_cnt = smp_ale_cnt,
              pop_hap_frq = pop_hap_frq,
              pop_ale_frq = pop_ale_frq))
}
#' Compiled version
cmpsimulateHMM <- cmpfun(simulateHMM)

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

#' Standard version
runBPF <- function(s_b, s_c, r, N, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num) {
  # run the BPF
  BPF <- runBPF_arma(s_b, s_c, r, N, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num)

  return(list(lik = BPF$lik,
              wght = BPF$wght,
              pop_frq_pre_resmp = BPF$part_pre_resmp,
              pop_frq_pst_resmp = BPF$part_pst_resmp))
}
#' Compiled version
cmprunBPF <- cmpfun(runBPF)

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

#' Standard version
calculateOptimalParticleNum <- function(s_b, s_c, r, N, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, gap_num) {
  OptNum <- calculateOptimalParticleNum_arma(s_b, s_c, r, N, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, gap_num)

  return(list(opt_pcl_num = as.vector(OptNum$opt_pcl_num),
              log_lik_sdv = as.vector(OptNum$log_lik_sdv)))
}
#' Compiled version
cmpcalculateOptimalParticleNum <- cmpfun(calculateOptimalParticleNum)

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

#' Standard version
runPMMH <- function(s_b, s_c, r, N, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num) {
  # set the recombination rate between ASIP and MC1R to 0.5
  r <- 0.5

  # run the PMMH
  PMMH <- runPMMH_arma(s_b, s_c, r, N, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num)

  return(list(s_b_chn = as.vector(PMMH$s_b_chn),
              s_c_chn = as.vector(PMMH$s_c_chn)))
}
#' Compiled version
cmprunPMMH <- cmpfun(runPMMH)

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

#' Standard version
runBayesianProcedure <- function(s_b, s_c, r, N, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, brn_num, thn_num) {
  # set the recombination rate between ASIP and MC1R to 0.5
  r <- 0.5

  # run the PMMH
  PMMH <- runPMMH_arma(s_b, s_c, r, N, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num)

  # burn-in and thinning
  s_b_chn <- as.vector(PMMH$s_b_chn)
  s_b_chn <- s_b_chn[brn_num:length(s_b_chn)]
  s_b_chn <- s_b_chn[(1:round(length(s_b_chn) / thn_num)) * thn_num]

  s_c_chn <- as.vector(PMMH$s_c_chn)
  s_c_chn <- s_c_chn[brn_num:length(s_c_chn)]
  s_c_chn <- s_c_chn[(1:round(length(s_c_chn) / thn_num)) * thn_num]

  # MMSE estimates for the selection coefficients
  s_b_mmse <- mean(s_b_chn)
  s_c_mmse <- mean(s_c_chn)

  # 95% HPD intervals for the selection coefficients
  s_b_hpd <- HPDinterval(as.mcmc(s_b_chn), prob = 0.95)
  s_c_hpd <- HPDinterval(as.mcmc(s_c_chn), prob = 0.95)

  return(list(s_b_mmse = s_b_mmse,
              s_c_mmse = s_c_mmse,
              s_b_hpd = s_b_hpd,
              s_c_hpd = s_c_hpd,
              s_b_chn = s_b_chn,
              s_c_chn = s_c_chn))
}
#' Compiled version
cmprunBayesianProcedure <- cmpfun(runBayesianProcedure)

################################################################################
