#' @title Inferring natural selection acting on horse coat colours and patterns during the process of domestication from ancient DNA data
#' @author Xiaoyang Dai, Sile Hu, Mark Beaumont, Feng Yu, Zhangyi He

#' version 1.1
#' Horse coat colours (ASIP & MC1R) under constant natural selection and non-constant demographic histories (N/A is not allowed)

#' R functions

#install.packages("purrr")
library("purrr")

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
sourceCpp("./Code/Code v1.0/Code v1.1/CFUN_COL.cpp")

################################################################################

#' Simulate the haplotype frequency trajectories according to the two-locus Wright-Fisher model with selection
#' Parameter setting
#' @param sel_cof the selection coefficients of the black and chestnut phenotypes
#' @param rec_rat the recombination rate between the ASIP and MC1R loci
#' @param pop_siz the size of the horse population (non-constant)
#' @param int_frq the initial haplotype frequencies of the population
#' @param int_gen the generation that the simulated haplotype frequency trajectories started
#' @param lst_gen the generation that the simulated haplotype frequency trajectories ended

#' Standard version
simulateWFM <- function(sel_cof, rec_rat, pop_siz, int_frq, int_gen, lst_gen) {
  fts_mat <- calculateFitnessMat_arma(sel_cof)

  WFM <- simulateWFM_arma(fts_mat, rec_rat, pop_siz, int_frq, int_gen, lst_gen)
  hap_frq_pth <- as.matrix(WFM$hap_frq_pth)
  gen_frq_pth <- as.matrix(WFM$gen_frq_pth)

  return(list(hap_frq_pth = hap_frq_pth,
              gen_frq_pth = gen_frq_pth))
}
#' Compiled version
cmpsimulateWFM <- cmpfun(simulateWFM)

########################################

#' Simulate the haplotype frequency trajectories according to the two-locus Wright-Fisher diffusion with selection using the Euler-Maruyama method
#' Parameter setting
#' @param sel_cof the selection coefficients of the black and chestnut phenotypes
#' @param rec_rat the recombination rate between the ASIP and MC1R loci
#' @param pop_siz the size of the horse population (non-constant)
#' @param ref_siz the reference size of the horse population
#' @param int_frq the initial haplotype frequencies of the population
#' @param int_gen the generation that the simulated haplotype frequency trajectories started
#' @param lst_gen the generation that the simulated haplotype frequency trajectories ended
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param dat_aug = TRUE/FALSE (return the simulated sample trajectory with data augmentation or not)

#' Standard version
simulateWFD <- function(sel_cof, rec_rat, pop_siz, ref_siz, int_frq, int_gen, lst_gen, ptn_num, dat_aug = TRUE) {
  frq_pth <- simulateWFD_arma(sel_cof, rec_rat, pop_siz, ref_siz, int_frq, int_gen, lst_gen, ptn_num)
  frq_pth <- as.matrix(frq_pth)

  if (dat_aug == FALSE) {
    return(frq_pth[, (0:(lst_gen - int_gen)) * ptn_num + 1])
  } else {
    return(frq_pth)
  }
}
#' Compiled version
cmpsimulateWFD <- cmpfun(simulateWFD)

########################################

#' Simulate the hidden Markov model
#' Parameter setting
#' @param model = "WFM"/"WFD" (return the observations from the underlying population evolving according to the WFM or the WFD)
#' @param sel_cof the selection coefficients of the black and chestnut phenotypes
#' @param rec_rat the recombination rate between the ASIP and MC1R loci
#' @param pop_siz the size of the horse population (non-constant)
#' @param int_frq the initial haplotype frequencies of the population
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_siz the count of the horses drawn from the population at all sampling time points
#' @param ref_siz the reference size of the horse population
#' @param ptn_num the number of the subintervals divided per generation in the Euler-Maruyama method for the WFD

#' Standard version
simulateHMM <- function(model, sel_cof, rec_rat, pop_siz, int_frq, smp_gen, smp_siz, ...) {
  int_gen <- min(smp_gen)
  lst_gen <- max(smp_gen)

  # generate the population haplotype and genotype frequency trajectories
  if (model == "WFM") {
    WFM <- cmpsimulateWFM(sel_cof, rec_rat, pop_siz, int_frq, int_gen, lst_gen)
    pop_hap_frq <- as.matrix(WFM$hap_frq_pth)
    pop_gen_frq <- as.matrix(WFM$gen_frq_pth)
  }
  if (model == "WFD") {
    pop_hap_frq <- cmpsimulateWFD(sel_cof, rec_rat, pop_siz, ref_siz, int_frq, int_gen, lst_gen, ptn_num, dat_aug = FALSE)
    pop_hap_frq <- as.matrix(pop_hap_frq)
    pop_gen_frq <- matrix(NA, nrow = 10, ncol = ncol(pop_hap_frq))
    fts_mat <- calculateFitnessMat_arma(sel_cof)
    for (k in 1:ncol(pop_hap_frq)) {
      hap_frq <- pop_hap_frq[, k]
      gen_frq <- fts_mat * (hap_frq %*% t(hap_frq)) / sum(fts_mat * (hap_frq %*% t(hap_frq)))
      gen_frq[lower.tri(gen_frq, diag = FALSE)] <- NA
      pop_gen_frq[, k] <- discard(as.vector(2 * gen_frq - diag(diag(gen_frq), nrow = 4, ncol = 4)), is.na)
    }
    pop_gen_frq <- as.matrix(pop_gen_frq)
  }

  # generate the sample genotype counts at all sampling time points
  # smp_gen_cnt <- matrix(NA, nrow = 9, ncol = length(smp_gen))
  # smp_gen_frq <- matrix(NA, nrow = 9, ncol = length(smp_gen))
  smp_gen_cnt <- matrix(NA, nrow = 10, ncol = length(smp_gen))
  smp_gen_frq <- matrix(NA, nrow = 10, ncol = length(smp_gen))
  for (k in 1:length(smp_gen)) {
    smp_gen_cnt[, k]  <- rmultinom(1, size = smp_siz[k], prob = pop_gen_frq[, smp_gen[k] - int_gen + 1])
    # gen_cnt <- rmultinom(1, size = smp_siz[k], prob = pop_gen_frq[, smp_gen[k] - int_gen + 1])
    # gen_cnt[5] <- gen_cnt[5] + gen_cnt[7]
    # smp_gen_cnt[, k] <- gen_cnt[-7]
    smp_gen_frq[, k] <- smp_gen_cnt[, k] / smp_siz[k]
  }

  return(list(smp_gen = smp_gen,
              smp_siz = smp_siz,
              smp_gen_cnt = smp_gen_cnt,
              smp_gen_frq = smp_gen_frq,
              pop_gen_frq = pop_gen_frq,
              pop_hap_frq = pop_hap_frq))
}
#' Compiled version
cmpsimulateHMM <- cmpfun(simulateHMM)

########################################

#' Run the bootstrap particle filter (BPF) with the two-locus Wright-Fisher diffusion with selection
#' Parameter setting
#' @param sel_cof the selection coefficients of the black and chestnut phenotypes
#' @param rec_rat the recombination rate between the ASIP and MC1R loci
#' @param pop_siz the size of the horse population (non-constant)
#' @param ref_siz the reference size of the horse population
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_siz the count of the horses drawn from the population at all sampling time points
#' @param smp_cnt the count of the genotypes observed in the sample at all sampling time points
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter

#' Standard version
runBPF <- function(sel_cof, rec_rat, pop_siz, ref_siz, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num) {
  # run the BPF
  BPF <- runBPF_arma(sel_cof, rec_rat, pop_siz, ref_siz, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num)

  lik <- BPF$lik
  wght <- BPF$wght
  hap_frq_pre_resmp <- BPF$hap_frq_pre_resmp
  hap_frq_pst_resmp <- BPF$hap_frq_pst_resmp
  gen_frq_pre_resmp <- BPF$gen_frq_pre_resmp
  gen_frq_pre_resmp[5, , ] <- gen_frq_pre_resmp[5, , ] + gen_frq_pre_resmp[7, , ]
  gen_frq_pre_resmp <- gen_frq_pre_resmp[-7, , ]
  gen_frq_pst_resmp <- BPF$gen_frq_pst_resmp
  gen_frq_pst_resmp[5, , ] <- gen_frq_pst_resmp[5, , ] + gen_frq_pst_resmp[7, , ]
  gen_frq_pst_resmp <- gen_frq_pst_resmp[-7, , ]

  return(list(lik = lik,
              wght = wght,
              hap_frq_pre_resmp = hap_frq_pre_resmp,
              hap_frq_pst_resmp = hap_frq_pst_resmp,
              gen_frq_pre_resmp = gen_frq_pre_resmp,
              gen_frq_pst_resmp = gen_frq_pst_resmp))
}
#' Compiled version
cmprunBPF <- cmpfun(runBPF)

########################################

#' Calculate the optimal particle number in the particle marginal Metropolis-Hastings (PMMH)
#' Parameter settings
#' @param sel_cof the selection coefficients of the black and chestnut phenotypes
#' @param rec_rat the recombination rate between the ASIP and MC1R loci
#' @param pop_siz the size of the horse population (non-constant)
#' @param ref_siz the reference size of the horse population
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_siz the count of the horses drawn from the population at all sampling time points
#' @param smp_cnt the count of the genotypes observed in the sample at all sampling time points
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter
#' @param gap_num the number of particles increased or decreased in the optimal particle number search

#' Standard version
calculateOptimalParticleNum <- function(sel_cof, rec_rat, pop_siz, ref_siz, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, gap_num) {
  # calculate the optimal particle number
  OptNum <- calculateOptimalParticleNum_arma(sel_cof, rec_rat, pop_siz, ref_siz, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, gap_num)

  return(list(opt_pcl_num = as.vector(OptNum$opt_pcl_num),
              log_lik_sdv = as.vector(OptNum$log_lik_sdv)))
}
#' Compiled version
cmpcalculateOptimalParticleNum <- cmpfun(calculateOptimalParticleNum)

########################################

#' Run the particle marginal Metropolis-Hastings (PMMH)
#' Parameter settings
#' @param sel_cof the selection coefficients of the black and chestnut phenotypes
#' @param rec_rat the recombination rate between the ASIP and MC1R loci
#' @param pop_siz the size of the horse population (non-constant)
#' @param ref_siz the reference size of the horse population
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_siz the count of the horses drawn from the population at all sampling time points
#' @param smp_cnt the count of the genotypes observed in the sample at all sampling time points
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter
#' @param itn_num the number of the iterations carried out in the particle marginal Metropolis-Hastings

#' Standard version
runPMMH <- function(sel_cof, rec_rat, pop_siz, ref_siz, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num) {
  # run the PMMH
  sel_cof_chn <- runPMMH_arma(sel_cof, rec_rat, pop_siz, ref_siz, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num)
  sel_cof_chn <- as.matrix(sel_cof_chn)

  return(sel_cof_chn)
}
#' Compiled version
cmprunPMMH <- cmpfun(runPMMH)

########################################

#' Run the Bayesian procedure for the inference of natural selection
#' Parameter settings
#' @param sel_cof the selection coefficients of the black and chestnut phenotypes
#' @param rec_rat the recombination rate between the ASIP and MC1R loci
#' @param pop_siz the size of the horse population (non-constant)
#' @param ref_siz the reference size of the horse population
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_siz the count of the chromosomes drawn from the population at all sampling time points
#' @param smp_cnt the count of the mutant alleles observed in the sample at all sampling time points
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter
#' @param itn_num the number of the iterations carried out in the particle marginal Metropolis-Hastings
#' @param brn_num the number of the iterations for burn-in
#' @param thn_num the number of the iterations for thinning

#' Standard version
runBayesianProcedure <- function(sel_cof, rec_rat, pop_siz, ref_siz, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, brn_num, thn_num) {
  # run the PMMH
  sel_cof_chn <- runPMMH_arma(sel_cof, rec_rat, pop_siz, ref_siz, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num)
  sel_cof_chn <- as.matrix(sel_cof_chn)

  # burn-in and thinning
  sel_cof_chn <- sel_cof_chn[, brn_num:dim(sel_cof_chn)[2]]
  sel_cof_chn <- sel_cof_chn[, (1:round(dim(sel_cof_chn)[2] / thn_num)) * thn_num]

  # MMSE estimates for the selection coefficients
  sel_cof_est <- rowMeans(sel_cof_chn)

  # 95% HPD intervals for the selection coefficients
  sel_cof_hpd <- matrix(NA, nrow = 2, ncol = 2)
  sel_cof_hpd[1, ] <- HPDinterval(as.mcmc(sel_cof_chn[1, ]), prob = 0.95)
  sel_cof_hpd[2, ] <- HPDinterval(as.mcmc(sel_cof_chn[2, ]), prob = 0.95)

  return(list(sel_cof_est = sel_cof_est,
              sel_cof_hpd = sel_cof_hpd,
              sel_cof_chn = sel_cof_chn))
}
#' Compiled version
cmprunBayesianProcedure <- cmpfun(runBayesianProcedure)

################################################################################
