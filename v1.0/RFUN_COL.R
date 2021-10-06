#' @title Estimating selection coefficients and testing their changes from ancient DNA data with the flexibility of modelling linkage and epistasis
#' @author Xiaoyang Dai, Wenyang Lyu, Mark Beaumont, Feng Yu, Zhangyi He

#' version 1.0
#' Phenotypes controlled by two gene with genetic linkage
#' Non-constant natural selection and non-constant demographic histories

#' Integrate prior knowledge from modern samples (gene polymorphism)

#' Input: genotype likelihoods
#' Output: posteriors for the selection coefficient, the genotype frequency trajectories of the population and the genotypes of the sample

#' Horse base coat colours (ASIP & MC1R)

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
sourceCpp("./Code/Code v1.0/Code v1.0/CFUN_COL.cpp")

################################################################################

#' Simulate the haplotype frequency trajectories according to the two-locus Wright-Fisher model with selection
#' Parameter setting
#' @param sel_cof the selection coefficients of the black and chestnut against bay
#' @param rec_rat the (artificial) recombination rate between the ASIP and MC1R loci (r = 0.5)
#' @param pop_siz the size of the horse population (non-constant)
#' @param int_frq the initial haplotype frequencies of the population
#' @param evt_gen the generation that the event of interest occurred
#' @param int_gen the generation that the simulated haplotype frequency trajectories started
#' @param lst_gen the generation that the simulated haplotype frequency trajectories ended

#' Standard version
simulateWFM <- function(sel_cof, rec_rat, pop_siz, int_frq, evt_gen, int_gen, lst_gen) {
  if (evt_gen >= lst_gen) {
    fts_mat <- calculateFitnessMat_arma(sel_cof[, 1])
    WFM <- simulateWFM_arma(fts_mat, rec_rat, pop_siz, int_frq, int_gen, lst_gen)
    hap_frq_pth <- as.matrix(WFM$hap_frq_pth)
    gen_frq_pth <- as.matrix(WFM$gen_frq_pth)
  } else if (evt_gen < int_gen) {
    fts_mat <- calculateFitnessMat_arma(sel_cof[, 2])
    WFM <- simulateWFM_arma(fts_mat, rec_rat, pop_siz, int_frq, int_gen, lst_gen)
    hap_frq_pth <- as.matrix(WFM$hap_frq_pth)
    gen_frq_pth <- as.matrix(WFM$gen_frq_pth)
  } else {
    fts_mat <- calculateFitnessMat_arma(sel_cof[, 1])
    WFM <- simulateWFM_arma(fts_mat, rec_rat, pop_siz, int_frq, int_gen, evt_gen)
    hap_frq_pth_pre_evt <- as.matrix(WFM$hap_frq_pth)
    gen_frq_pth_pre_evt <- as.matrix(WFM$gen_frq_pth)

    fts_mat <- calculateFitnessMat_arma(sel_cof[, 2])
    WFM <- simulateWFM_arma(fts_mat, rec_rat, pop_siz, hap_frq_pth_pre_evt[, ncol(hap_frq_pth_pre_evt)], evt_gen, lst_gen)
    hap_frq_pth_pst_evt <- as.matrix(WFM$hap_frq_pth)
    gen_frq_pth_pst_evt <- as.matrix(WFM$gen_frq_pth)

    hap_frq_pth <- cbind(hap_frq_pth_pre_evt, hap_frq_pth_pst_evt[, -1])
    gen_frq_pth <- cbind(gen_frq_pth_pre_evt, gen_frq_pth_pst_evt[, -1])
  }

  return(list(hap_frq_pth = hap_frq_pth,
              gen_frq_pth = gen_frq_pth))
}
#' Compiled version
cmpsimulateWFM <- cmpfun(simulateWFM)

########################################

#' Simulate the haplotype frequency trajectories according to the two-locus Wright-Fisher diffusion with selection using the Euler-Maruyama method
#' Parameter setting
#' @param sel_cof the selection coefficients of the black and chestnut against bay
#' @param rec_rat the (artificial) recombination rate between the ASIP and MC1R loci (r = 0.5)
#' @param pop_siz the size of the horse population (non-constant)
#' @param ref_siz the reference size of the horse population
#' @param int_frq the initial haplotype frequencies of the population
#' @param evt_gen the generation that the event of interest occurred
#' @param int_gen the generation that the simulated haplotype frequency trajectories started
#' @param lst_gen the generation that the simulated haplotype frequency trajectories ended
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param dat_aug = TRUE/FALSE (return the simulated sample trajectory with data augmentation or not)

#' Standard version
simulateWFD <- function(sel_cof, rec_rat, pop_siz, ref_siz, int_frq, evt_gen, int_gen, lst_gen, ptn_num, dat_aug = TRUE) {
  if (evt_gen >= lst_gen) {
    hap_frq_pth <- simulateWFD_arma(sel_cof[, 1], rec_rat, pop_siz, ref_siz, int_frq, int_gen, lst_gen, ptn_num)
    hap_frq_pth <- as.matrix(hap_frq_pth)
  } else if (evt_gen < int_gen) {
    hap_frq_pth <- simulateWFD_arma(sel_cof[, 2], rec_rat, pop_siz, ref_siz, int_frq, int_gen, lst_gen, ptn_num)
    hap_frq_pth <- as.matrix(hap_frq_pth)
  } else {
    hap_frq_pth_pre_evt <- simulateWFD_arma(sel_cof[, 1], rec_rat, pop_siz, ref_siz, int_frq, int_gen, evt_gen, ptn_num)
    hap_frq_pth_pre_evt <- as.matrix(hap_frq_pth_pre_evt)

    hap_frq_pth_pst_evt <- simulateWFD_arma(sel_cof[, 2], rec_rat, pop_siz, ref_siz, hap_frq_pth_pre_evt[, ncol(hap_frq_pth_pre_evt)], evt_gen, lst_gen, ptn_num)
    hap_frq_pth_pst_evt <- as.matrix(hap_frq_pth_pst_evt)

    hap_frq_pth <- cbind(hap_frq_pth_pre_evt, hap_frq_pth_pst_evt[, -1])
  }

  if (dat_aug == FALSE) {
    return(hap_frq_pth[, (0:(lst_gen - int_gen)) * ptn_num + 1])
  } else {
    return(hap_frq_pth)
  }
}
#' Compiled version
cmpsimulateWFD <- cmpfun(simulateWFD)

########################################

#' Simulate the hidden Markov model
#' Parameter setting
#' @param model = "WFM"/"WFD" (return the observations from the underlying population evolving according to the WFM or the WFD)
#' @param sel_cof the selection coefficients of the black and chestnut against bay
#' @param rec_rat the (artificial) recombination rate between the ASIP and MC1R loci (r = 0.5)
#' @param pop_siz the size of the horse population (non-constant)
#' @param int_con the initial haplotype frequencies of the population / the initial mutant allele frequencies and the linkage disequilibrium of the population
#' @param evt_gen the generation that the event of interest occurred
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_siz the count of the horses drawn from the population at all sampling time points
#' @param obs_hap = TRUE/FALSE (return the simulated sample genotypes with haplotype information or not)
#' @param ref_siz the reference size of the horse population
#' @param ptn_num the number of the subintervals divided per generation in the Euler-Maruyama method for the WFD

#' Standard version
simulateHMM <- function(model, sel_cof, rec_rat, pop_siz, int_con, evt_gen, smp_gen, smp_siz, obs_hap = FALSE, ...) {
  int_gen <- min(smp_gen)
  lst_gen <- max(smp_gen)

  # calculate the initial population haplotype frequencies
  int_frq <- rep(0, length.out = 4)
  if (length(int_con) == 3) {
    int_frq[1] <- (1 - int_con[1]) * (1 - int_con[2]) + int_con[3]
    int_frq[2] <- (1 - int_con[1]) * int_con[2] - int_con[3]
    int_frq[3] <- int_con[1] * (1 - int_con[2]) - int_con[3]
    int_frq[4] <- int_con[1] * int_con[2] + int_con[3]
  } else {
    int_frq <- int_con
  }

  # generate the population haplotype and genotype frequency trajectories
  if (model == "WFM") {
    WFM <- cmpsimulateWFM(sel_cof, rec_rat, pop_siz, int_frq, evt_gen, int_gen, lst_gen)
    pop_hap_frq <- WFM$hap_frq_pth
    pop_gen_frq <- WFM$gen_frq_pth
  }
  if (model == "WFD") {
    pop_hap_frq <- cmpsimulateWFD(sel_cof, rec_rat, pop_siz, ref_siz, int_frq, evt_gen, int_gen, lst_gen, ptn_num, dat_aug = FALSE)
    pop_hap_frq <- as.matrix(pop_hap_frq)

    pop_gen_frq <- matrix(NA, nrow = 10, ncol = ncol(pop_hap_frq))
    fts_mat <- calculateFitnessMat_arma(sel_cof[, 1])
    for (k in 1:(evt_gen - int_gen)) {
      hap_frq <- pop_hap_frq[, k]
      gen_frq <- fts_mat * (hap_frq %*% t(hap_frq)) / sum(fts_mat * (hap_frq %*% t(hap_frq)))
      gen_frq[lower.tri(gen_frq, diag = FALSE)] <- NA
      pop_gen_frq[, k] <- discard(as.vector(2 * gen_frq - diag(diag(gen_frq), nrow = 4, ncol = 4)), is.na)
    }
    fts_mat <- calculateFitnessMat_arma(sel_cof[, 2])
    for (k in (evt_gen - int_gen + 1):ncol(pop_hap_frq)) {
      hap_frq <- pop_hap_frq[, k]
      gen_frq <- fts_mat * (hap_frq %*% t(hap_frq)) / sum(fts_mat * (hap_frq %*% t(hap_frq)))
      gen_frq[lower.tri(gen_frq, diag = FALSE)] <- NA
      pop_gen_frq[, k] <- discard(as.vector(2 * gen_frq - diag(diag(gen_frq), nrow = 4, ncol = 4)), is.na)
    }
  }

  # generate the sample genotype counts at all sampling time points
  if (obs_hap == TRUE) {
    smp_gen_cnt <- matrix(NA, nrow = 10, ncol = length(smp_gen))
    smp_gen_frq <- matrix(NA, nrow = 10, ncol = length(smp_gen))
    for (k in 1:length(smp_gen)) {
      smp_gen_cnt[, k]  <- rmultinom(1, size = smp_siz[k], prob = pop_gen_frq[, smp_gen[k] - int_gen + 1])
      smp_gen_frq[, k] <- smp_gen_cnt[, k] / smp_siz[k]
    }
  } else {
    smp_gen_cnt <- matrix(NA, nrow = 9, ncol = length(smp_gen))
    smp_gen_frq <- matrix(NA, nrow = 9, ncol = length(smp_gen))
    for (k in 1:length(smp_gen)) {
      gen_cnt <- rmultinom(1, size = smp_siz[k], prob = pop_gen_frq[, smp_gen[k] - int_gen + 1])
      gen_cnt[5] <- gen_cnt[5] + gen_cnt[7]
      smp_gen_cnt[, k] <- gen_cnt[-7]
      smp_gen_frq[, k] <- smp_gen_cnt[, k] / smp_siz[k]
    }
    pop_gen_frq[5, ] <- pop_gen_frq[5, ] + pop_gen_frq[7, ]
    pop_gen_frq <- pop_gen_frq[-7, ]
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
#' @param sel_cof the selection coefficients of the black and chestnut against bay
#' @param rec_rat the (artificial) recombination rate between the ASIP and MC1R loci (r = 0.5)
#' @param pop_siz the size of the horse population (non-constant)
#' @param ref_siz the reference size of the horse population
#' @param evt_gen the generation that the event of interest occurred
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_siz the count of the horses drawn from the population at all sampling time points
#' @param smp_cnt the count of the genotypes observed in the sample at all sampling time points
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter

#' Standard version
runBPF <- function(sel_cof, rec_rat, pop_siz, ref_siz, evt_gen, frq_pth, raw_smp, ptn_num, pcl_num) {
  # preprocess the raw sample
  if (ncol(raw_smp) == 7) {
    rownames(raw_smp) <- NULL
    colnames(raw_smp) <- c("generation", "A0A0", "A0A1", "A1A1", "B0B0", "B0B1", "B1B1")

    raw_smp <- as.data.frame(cbind(raw_smp[, 1], 
      raw_smp[, 2] * raw_smp[, 5], 
      raw_smp[, 2] * raw_smp[, 6], 
      raw_smp[, 2] * raw_smp[, 7], 
      raw_smp[, 3] * raw_smp[, 5], 
      raw_smp[, 3] * raw_smp[, 6], 
      raw_smp[, 4] * raw_smp[, 5], 
      raw_smp[, 3] * raw_smp[, 7], 
      raw_smp[, 4] * raw_smp[, 6], 
      raw_smp[, 4] * raw_smp[, 7]))
    colnames(raw_smp) <- c("generation", "A0A0/B0B0", "A0A0/B0B1", "A0A0/B1B1", "A0A1/B0B0", "A0A1/B0B1", "A1A1/B0B0", "A0A1/B1B1", "A1A1/B0B1", "A1A1/B1B1")

    raw_smp <- raw_smp[order(raw_smp$generation), ]
    evt_gen <- evt_gen - min(raw_smp$generation)
    raw_smp$generation <- raw_smp$generation - min(raw_smp$generation)
    rownames(raw_smp) <- NULL
    colnames(raw_smp) <- NULL
    raw_smp <- t(as.matrix(raw_smp))
  } else {
    raw_smp <- raw_smp[, -(2:3)]
    rownames(raw_smp) <- NULL
    colnames(raw_smp) <- c("generation", "A0A0", "A0A1", "A1A1", "B0B0", "B0B1", "B1B1")

    raw_smp <- as.data.frame(cbind(raw_smp[, 1], 
      raw_smp[, 2] * raw_smp[, 5], 
      raw_smp[, 2] * raw_smp[, 6], 
      raw_smp[, 2] * raw_smp[, 7], 
      raw_smp[, 3] * raw_smp[, 5], 
      raw_smp[, 3] * raw_smp[, 6], 
      raw_smp[, 4] * raw_smp[, 5], 
      raw_smp[, 3] * raw_smp[, 7], 
      raw_smp[, 4] * raw_smp[, 6], 
      raw_smp[, 4] * raw_smp[, 7]))
    colnames(raw_smp) <- c("generation", "A0A0/B0B0", "A0A0/B0B1", "A0A0/B1B1", "A0A1/B0B0", "A0A1/B0B1", "A1A1/B0B0", "A0A1/B1B1", "A1A1/B0B1", "A1A1/B1B1")

    raw_smp <- raw_smp[order(raw_smp$generation), ]
    evt_gen <- evt_gen - min(raw_smp$generation)
    raw_smp$generation <- raw_smp$generation - min(raw_smp$generation)
    rownames(raw_smp) <- NULL
    colnames(raw_smp) <- NULL
    raw_smp <- t(as.matrix(raw_smp))
  }

  # run the BPF
  BPF <- runBPF_arma(sel_cof, rec_rat, pop_siz, ref_siz, evt_gen, frq_pth, raw_smp, ptn_num, pcl_num)

  lik <- BPF$lik
  wght <- BPF$wght
  hap_frq_pth <- BPF$hap_frq_pth
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
              hap_frq_pth = hap_frq_pth,
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
#' @param sel_cof the selection coefficients of the black and chestnut against bay
#' @param rec_rat the (artificial) recombination rate between the ASIP and MC1R loci (r = 0.5)
#' @param pop_siz the size of the horse population (non-constant)
#' @param ref_siz the reference size of the horse population
#' @param evt_gen the generation that the event of interest occurred
#' @param frq_pth the mutant allele frequency of the population
#' @param raw_smp the sample of ancient horses (sampling times and individual genotypes)
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter
#' @param gap_num the number of particles increased or decreased in the optimal particle number search

#' Standard version
calculateOptimalParticleNum <- function(sel_cof, rec_rat, pop_siz, ref_siz, evt_gen, frq_pth, raw_smp, ptn_num, pcl_num, gap_num) {
  # preprocess the raw sample
  if (ncol(raw_smp) == 7) {
    rownames(raw_smp) <- NULL
    colnames(raw_smp) <- c("generation", "A0A0", "A0A1", "A1A1", "B0B0", "B0B1", "B1B1")

    raw_smp <- as.data.frame(cbind(raw_smp[, 1], 
      raw_smp[, 2] * raw_smp[, 5], 
      raw_smp[, 2] * raw_smp[, 6], 
      raw_smp[, 2] * raw_smp[, 7], 
      raw_smp[, 3] * raw_smp[, 5], 
      raw_smp[, 3] * raw_smp[, 6], 
      raw_smp[, 4] * raw_smp[, 5], 
      raw_smp[, 3] * raw_smp[, 7], 
      raw_smp[, 4] * raw_smp[, 6], 
      raw_smp[, 4] * raw_smp[, 7]))
    colnames(raw_smp) <- c("generation", "A0A0/B0B0", "A0A0/B0B1", "A0A0/B1B1", "A0A1/B0B0", "A0A1/B0B1", "A1A1/B0B0", "A0A1/B1B1", "A1A1/B0B1", "A1A1/B1B1")

    raw_smp <- raw_smp[order(raw_smp$generation), ]
    evt_gen <- evt_gen - min(raw_smp$generation)
    raw_smp$generation <- raw_smp$generation - min(raw_smp$generation)
    rownames(raw_smp) <- NULL
    colnames(raw_smp) <- NULL
    raw_smp <- t(as.matrix(raw_smp))
  } else {
    raw_smp <- raw_smp[, -(2:3)]
    rownames(raw_smp) <- NULL
    colnames(raw_smp) <- c("generation", "A0A0", "A0A1", "A1A1", "B0B0", "B0B1", "B1B1")

    raw_smp <- as.data.frame(cbind(raw_smp[, 1], 
      raw_smp[, 2] * raw_smp[, 5], 
      raw_smp[, 2] * raw_smp[, 6], 
      raw_smp[, 2] * raw_smp[, 7], 
      raw_smp[, 3] * raw_smp[, 5], 
      raw_smp[, 3] * raw_smp[, 6], 
      raw_smp[, 4] * raw_smp[, 5], 
      raw_smp[, 3] * raw_smp[, 7], 
      raw_smp[, 4] * raw_smp[, 6], 
      raw_smp[, 4] * raw_smp[, 7]))
    colnames(raw_smp) <- c("generation", "A0A0/B0B0", "A0A0/B0B1", "A0A0/B1B1", "A0A1/B0B0", "A0A1/B0B1", "A1A1/B0B0", "A0A1/B1B1", "A1A1/B0B1", "A1A1/B1B1")

    raw_smp <- raw_smp[order(raw_smp$generation), ]
    evt_gen <- evt_gen - min(raw_smp$generation)
    raw_smp$generation <- raw_smp$generation - min(raw_smp$generation)
    rownames(raw_smp) <- NULL
    colnames(raw_smp) <- NULL
    raw_smp <- t(as.matrix(raw_smp))
  }

  # calculate the optimal particle number
  OptNum <- calculateOptimalParticleNum_arma(sel_cof, rec_rat, pop_siz, ref_siz, evt_gen, frq_pth, raw_smp, ptn_num, pcl_num, gap_num)

  return(list(opt_pcl_num = as.vector(OptNum$opt_pcl_num),
              log_lik_sdv = as.vector(OptNum$log_lik_sdv)))
}
#' Compiled version
cmpcalculateOptimalParticleNum <- cmpfun(calculateOptimalParticleNum)

########################################

#' Run the particle marginal Metropolis-Hastings (PMMH)
#' Parameter settings
#' @param sel_cof the selection coefficients of the black and chestnut against bay
#' @param rec_rat the (artificial) recombination rate between the ASIP and MC1R loci (r = 0.5)
#' @param pop_siz the size of the horse population (non-constant)
#' @param ref_siz the reference size of the horse population
#' @param evt_gen the generation that the event of interest occurred
#' @param raw_smp the sample of ancient horses (sampling times and individual genotypes)
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter
#' @param itn_num the number of the iterations carried out in the PMMH

#' Standard version
runPMMH <- function(sel_cof, rec_rat, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num) {
  # preprocess the raw sample
  if (ncol(raw_smp) == 7) {
    rownames(raw_smp) <- NULL
    colnames(raw_smp) <- c("generation", "A0A0", "A0A1", "A1A1", "B0B0", "B0B1", "B1B1")

    raw_smp <- as.data.frame(cbind(raw_smp[, 1], 
      raw_smp[, 2] * raw_smp[, 5], 
      raw_smp[, 2] * raw_smp[, 6], 
      raw_smp[, 2] * raw_smp[, 7], 
      raw_smp[, 3] * raw_smp[, 5], 
      raw_smp[, 3] * raw_smp[, 6], 
      raw_smp[, 4] * raw_smp[, 5], 
      raw_smp[, 3] * raw_smp[, 7], 
      raw_smp[, 4] * raw_smp[, 6], 
      raw_smp[, 4] * raw_smp[, 7]))
    colnames(raw_smp) <- c("generation", "A0A0/B0B0", "A0A0/B0B1", "A0A0/B1B1", "A0A1/B0B0", "A0A1/B0B1", "A1A1/B0B0", "A0A1/B1B1", "A1A1/B0B1", "A1A1/B1B1")

    raw_smp <- raw_smp[order(raw_smp$generation), ]
    evt_gen <- evt_gen - min(raw_smp$generation)
    raw_smp$generation <- raw_smp$generation - min(raw_smp$generation)
    rownames(raw_smp) <- NULL
    colnames(raw_smp) <- NULL
    raw_smp <- t(as.matrix(raw_smp))
  } else {
    raw_smp <- raw_smp[, -(2:3)]
    rownames(raw_smp) <- NULL
    colnames(raw_smp) <- c("generation", "A0A0", "A0A1", "A1A1", "B0B0", "B0B1", "B1B1")

    raw_smp <- as.data.frame(cbind(raw_smp[, 1], 
      raw_smp[, 2] * raw_smp[, 5], 
      raw_smp[, 2] * raw_smp[, 6], 
      raw_smp[, 2] * raw_smp[, 7], 
      raw_smp[, 3] * raw_smp[, 5], 
      raw_smp[, 3] * raw_smp[, 6], 
      raw_smp[, 4] * raw_smp[, 5], 
      raw_smp[, 3] * raw_smp[, 7], 
      raw_smp[, 4] * raw_smp[, 6], 
      raw_smp[, 4] * raw_smp[, 7]))
    colnames(raw_smp) <- c("generation", "A0A0/B0B0", "A0A0/B0B1", "A0A0/B1B1", "A0A1/B0B0", "A0A1/B0B1", "A1A1/B0B0", "A0A1/B1B1", "A1A1/B0B1", "A1A1/B1B1")

    raw_smp <- raw_smp[order(raw_smp$generation), ]
    evt_gen <- evt_gen - min(raw_smp$generation)
    raw_smp$generation <- raw_smp$generation - min(raw_smp$generation)
    rownames(raw_smp) <- NULL
    colnames(raw_smp) <- NULL
    raw_smp <- t(as.matrix(raw_smp))
  }

  # run the PMMH
  PMMH <- runPMMH_arma(sel_cof, rec_rat, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num)

  return(list(sel_cof_chn = as.array(PMMH$sel_cof_chn),
              frq_pth_chn = as.array(PMMH$frq_pth_chn),
              imp_smp_chn = as.array(PMMH$imp_smp_chn)))
}
#' Compiled version
cmprunPMMH <- cmpfun(runPMMH)

########################################

#' Run the adaptive particle marginal Metropolis-Hastings (AdaptPMMH)
#' Parameter settings
#' @param sel_cof the selection coefficients of the black and chestnut against bay
#' @param rec_rat the (artificial) recombination rate between the ASIP and MC1R loci (r = 0.5)
#' @param pop_siz the size of the horse population (non-constant)
#' @param ref_siz the reference size of the horse population
#' @param evt_gen the generation that the event of interest occurred
#' @param raw_smp the sample of ancient horses (sampling times and individual genotypes)
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter
#' @param itn_num the number of the iterations carried out in the PMMH
#' @param stp_siz the step size sequence in the adaptive setting (decaying to zero)
#' @param apt_rto the target mean acceptance probability of the adaptive setting

#' Standard version
runAdaptPMMH <- function(sel_cof, rec_rat, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, stp_siz, apt_rto) {
  # preprocess the raw sample
  if (ncol(raw_smp) == 7) {
    rownames(raw_smp) <- NULL
    colnames(raw_smp) <- c("generation", "A0A0", "A0A1", "A1A1", "B0B0", "B0B1", "B1B1")

    raw_smp <- as.data.frame(cbind(raw_smp[, 1], 
      raw_smp[, 2] * raw_smp[, 5], 
      raw_smp[, 2] * raw_smp[, 6], 
      raw_smp[, 2] * raw_smp[, 7], 
      raw_smp[, 3] * raw_smp[, 5], 
      raw_smp[, 3] * raw_smp[, 6], 
      raw_smp[, 4] * raw_smp[, 5], 
      raw_smp[, 3] * raw_smp[, 7], 
      raw_smp[, 4] * raw_smp[, 6], 
      raw_smp[, 4] * raw_smp[, 7]))
    colnames(raw_smp) <- c("generation", "A0A0/B0B0", "A0A0/B0B1", "A0A0/B1B1", "A0A1/B0B0", "A0A1/B0B1", "A1A1/B0B0", "A0A1/B1B1", "A1A1/B0B1", "A1A1/B1B1")

    raw_smp <- raw_smp[order(raw_smp$generation), ]
    evt_gen <- evt_gen - min(raw_smp$generation)
    raw_smp$generation <- raw_smp$generation - min(raw_smp$generation)
    rownames(raw_smp) <- NULL
    colnames(raw_smp) <- NULL
    raw_smp <- t(as.matrix(raw_smp))
  } else {
    raw_smp <- raw_smp[, -(2:3)]
    rownames(raw_smp) <- NULL
    colnames(raw_smp) <- c("generation", "A0A0", "A0A1", "A1A1", "B0B0", "B0B1", "B1B1")

    raw_smp <- as.data.frame(cbind(raw_smp[, 1], 
      raw_smp[, 2] * raw_smp[, 5], 
      raw_smp[, 2] * raw_smp[, 6], 
      raw_smp[, 2] * raw_smp[, 7], 
      raw_smp[, 3] * raw_smp[, 5], 
      raw_smp[, 3] * raw_smp[, 6], 
      raw_smp[, 4] * raw_smp[, 5], 
      raw_smp[, 3] * raw_smp[, 7], 
      raw_smp[, 4] * raw_smp[, 6], 
      raw_smp[, 4] * raw_smp[, 7]))
    colnames(raw_smp) <- c("generation", "A0A0/B0B0", "A0A0/B0B1", "A0A0/B1B1", "A0A1/B0B0", "A0A1/B0B1", "A1A1/B0B0", "A0A1/B1B1", "A1A1/B0B1", "A1A1/B1B1")

    raw_smp <- raw_smp[order(raw_smp$generation), ]
    evt_gen <- evt_gen - min(raw_smp$generation)
    raw_smp$generation <- raw_smp$generation - min(raw_smp$generation)
    rownames(raw_smp) <- NULL
    colnames(raw_smp) <- NULL
    raw_smp <- t(as.matrix(raw_smp))
  }

  # run the adaptive PMMH
  PMMH <- runAdaptPMMH_arma(sel_cof, rec_rat, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, stp_siz, apt_rto)

  return(list(sel_cof_chn = as.array(PMMH$sel_cof_chn),
              frq_pth_chn = as.array(PMMH$frq_pth_chn),
              imp_smp_chn = as.array(PMMH$imp_smp_chn)))
}
#' Compiled version
cmprunAdaptPMMH <- cmpfun(runAdaptPMMH)

########################################

#' Run the Bayesian procedure for the inference of natural selection
#' Parameter settings
#' @param sel_cof the selection coefficients of the black and chestnut against bay
#' @param rec_rat the (artificial) recombination rate between the ASIP and MC1R loci (r = 0.5)
#' @param pop_siz the size of the horse population (non-constant)
#' @param ref_siz the reference size of the horse population
#' @param evt_gen the generation that the event of interest occurred
#' @param raw_smp the sample of ancient horses (sampling times and individual genotypes)
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter
#' @param itn_num the number of the iterations carried out in the PMMH
#' @param brn_num the number of the iterations for burn-in
#' @param thn_num the number of the iterations for thinning
#' @param adp_set = TRUE/FALSE (return the result with the adaptive setting or not)
#' @param stp_siz the step size sequence in the adaptive setting (decaying to zero)
#' @param apt_rto the target mean acceptance probability of the adaptive setting

#' Standard version
runBayesianProcedure <- function(sel_cof, rec_rat, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, brn_num, thn_num, adp_set, ...) {
  # preprocess the raw sample
  if (ncol(raw_smp) == 7) {
    rownames(raw_smp) <- NULL
    colnames(raw_smp) <- c("generation", "A0A0", "A0A1", "A1A1", "B0B0", "B0B1", "B1B1")

    raw_smp <- as.data.frame(cbind(raw_smp[, 1], 
      raw_smp[, 2] * raw_smp[, 5], 
      raw_smp[, 2] * raw_smp[, 6], 
      raw_smp[, 2] * raw_smp[, 7], 
      raw_smp[, 3] * raw_smp[, 5], 
      raw_smp[, 3] * raw_smp[, 6], 
      raw_smp[, 4] * raw_smp[, 5], 
      raw_smp[, 3] * raw_smp[, 7], 
      raw_smp[, 4] * raw_smp[, 6], 
      raw_smp[, 4] * raw_smp[, 7]))
    colnames(raw_smp) <- c("generation", "A0A0/B0B0", "A0A0/B0B1", "A0A0/B1B1", "A0A1/B0B0", "A0A1/B0B1", "A1A1/B0B0", "A0A1/B1B1", "A1A1/B0B1", "A1A1/B1B1")

    raw_smp <- raw_smp[order(raw_smp$generation), ]
    evt_gen <- evt_gen - min(raw_smp$generation)
    raw_smp$generation <- raw_smp$generation - min(raw_smp$generation)
    rownames(raw_smp) <- NULL
    colnames(raw_smp) <- NULL
    raw_smp <- t(as.matrix(raw_smp))
  } else {
    raw_smp <- raw_smp[, -(2:3)]
    rownames(raw_smp) <- NULL
    colnames(raw_smp) <- c("generation", "A0A0", "A0A1", "A1A1", "B0B0", "B0B1", "B1B1")

    raw_smp <- as.data.frame(cbind(raw_smp[, 1], 
      raw_smp[, 2] * raw_smp[, 5], 
      raw_smp[, 2] * raw_smp[, 6], 
      raw_smp[, 2] * raw_smp[, 7], 
      raw_smp[, 3] * raw_smp[, 5], 
      raw_smp[, 3] * raw_smp[, 6], 
      raw_smp[, 4] * raw_smp[, 5], 
      raw_smp[, 3] * raw_smp[, 7], 
      raw_smp[, 4] * raw_smp[, 6], 
      raw_smp[, 4] * raw_smp[, 7]))
    colnames(raw_smp) <- c("generation", "A0A0/B0B0", "A0A0/B0B1", "A0A0/B1B1", "A0A1/B0B0", "A0A1/B0B1", "A1A1/B0B0", "A0A1/B1B1", "A1A1/B0B1", "A1A1/B1B1")

    raw_smp <- raw_smp[order(raw_smp$generation), ]
    evt_gen <- evt_gen - min(raw_smp$generation)
    raw_smp$generation <- raw_smp$generation - min(raw_smp$generation)
    rownames(raw_smp) <- NULL
    colnames(raw_smp) <- NULL
    raw_smp <- t(as.matrix(raw_smp))
  }

  if (adp_set == TRUE) {
    # run the adaptive PMMH
    PMMH <- runAdaptPMMH_arma(sel_cof, rec_rat, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, stp_siz, apt_rto)
  } else {
    # run the PMMH
    PMMH <- runPMMH_arma(sel_cof, rec_rat, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num)
  }
  sel_cof_chn <- as.array(PMMH$sel_cof_chn)
  frq_pth_chn <- as.array(PMMH$frq_pth_chn)
  imp_smp_chn <- as.array(PMMH$imp_smp_chn)

  # burn-in and thinning
  sel_cof_chn <- sel_cof_chn[, , brn_num:dim(sel_cof_chn)[3]]
  sel_cof_chn <- sel_cof_chn[, , (1:round(dim(sel_cof_chn)[3] / thn_num)) * thn_num]
  frq_pth_chn <- frq_pth_chn[, , brn_num:dim(frq_pth_chn)[3]]
  frq_pth_chn <- frq_pth_chn[, , (1:round(dim(frq_pth_chn)[3] / thn_num)) * thn_num]
  smp_dat_chn <- smp_dat_chn[, , brn_num:dim(smp_dat_chn)[3]]
  smp_dat_chn <- smp_dat_chn[, , (1:round(dim(smp_dat_chn)[3] / thn_num)) * thn_num]
  
  # MMSE estimates for selection coefficients and haplotype frequencies
  sel_cof_est <- matrix(NA, nrow = 3, ncol = 2)
  sel_cof_est[, 1] <- rowMeans(sel_cof_chn[, 1, ])
  sel_cof_est[, 2] <- rowMeans(sel_cof_chn[, 2, ])
  frq_pth_est <- matrix(NA, nrow = 4, ncol = dim(frq_pth_chn)[2])
  frq_pth_est[1, ] <- rowMeans(frq_pth_chn[1, , ])
  frq_pth_est[2, ] <- rowMeans(frq_pth_chn[2, , ])
  frq_pth_est[3, ] <- rowMeans(frq_pth_chn[3, , ])
  frq_pth_est[4, ] <- rowMeans(frq_pth_chn[4, , ])

  # 95% HPD intervals for selection coefficients and haplotype frequencies
  sel_cof_hpd <- array(NA, dim = c(3, 2, 2))
  sel_cof_hpd[1, , 1] <- HPDinterval(as.mcmc(sel_cof_chn[1, 1, ]), prob = 0.95)
  sel_cof_hpd[2, , 1] <- HPDinterval(as.mcmc(sel_cof_chn[2, 1, ]), prob = 0.95)
  sel_cof_hpd[3, , 1] <- HPDinterval(as.mcmc(sel_cof_chn[3, 1, ]), prob = 0.95)
  sel_cof_hpd[1, , 2] <- HPDinterval(as.mcmc(sel_cof_chn[1, 2, ]), prob = 0.95)
  sel_cof_hpd[2, , 2] <- HPDinterval(as.mcmc(sel_cof_chn[2, 2, ]), prob = 0.95)
  sel_cof_hpd[3, , 2] <- HPDinterval(as.mcmc(sel_cof_chn[3, 2, ]), prob = 0.95)
  frq_pth_hpd <- array(NA, dim = c(4, 2, dim(frq_pth_chn)[2]))
  for (i in 1:dim(frq_pth_chn)[2]) {
    frq_pth_hpd[1, , i] <- HPDinterval(as.mcmc(frq_pth_chn[1, i, ]), prob = 0.95)
    frq_pth_hpd[2, , i] <- HPDinterval(as.mcmc(frq_pth_chn[2, i, ]), prob = 0.95)
    frq_pth_hpd[3, , i] <- HPDinterval(as.mcmc(frq_pth_chn[3, i, ]), prob = 0.95)
    frq_pth_hpd[4, , i] <- HPDinterval(as.mcmc(frq_pth_chn[4, i, ]), prob = 0.95)
  }

  # calculate the changes in the selection coefficients before and after the event
  dif_sel_chn <- sel_cof_chn[, 2, ] - sel_cof_chn[, 1, ]

  # MMSE estimates for the changes in the selection coefficients before and after the event
  dif_sel_est <- rowMeans(dif_sel_chn)

  # 95% HPD intervals for the changes in the selection coefficients before and after the event
  dif_sel_hpd <- matrix(NA, nrow = 3, ncol = 2)
  dif_sel_hpd[1, ] <- HPDinterval(as.mcmc(dif_sel_chn[1, ]), prob = 0.95)
  dif_sel_hpd[2, ] <- HPDinterval(as.mcmc(dif_sel_chn[2, ]), prob = 0.95)
  dif_sel_hpd[3, ] <- HPDinterval(as.mcmc(dif_sel_chn[3, ]), prob = 0.95)

  # posterior probability for genotypes
  imp_smp_est <- matrix(NA, nrow = ncol(raw_smp), ncol = 9)
  for (i in 1:nrow(imp_smp_est)) {
    imp_smp_est[i, ] <- rowSums(imp_smp_chn[, i, ]) / dim(imp_smp_chn)[3]
  }
  imp_smp_est <- cbind(raw_smp[1, ], imp_smp_est)
  imp_smp_est <- as.data.frame(imp_smp_est)
  rownames(imp_smp_est) <- NULL
  colnames(imp_smp_est) <- c("generation", "A0A0/B0B0", "A0A0/B0B1", "A0A0/B1B1", "A0A1/B0B0", "A0A1/B0B1", "A1A1/B0B0", "A0A1/B1B1", "A1A1/B0B1", "A1A1/B1B1")

  imp_smp_est <- as.data.frame(cbind(imp_smp_est[, 1], 
    imp_smp_est[, 2] + imp_smp_est[, 3] + imp_smp_est[, 4], 
    imp_smp_est[, 5] + imp_smp_est[, 6] + imp_smp_est[, 8], 
    imp_smp_est[, 7] + imp_smp_est[, 9] + imp_smp_est[, 10], 
    imp_smp_est[, 2] + imp_smp_est[, 5] + imp_smp_est[, 7], 
    imp_smp_est[, 3] + imp_smp_est[, 6] + imp_smp_est[, 9], 
    imp_smp_est[, 4] + imp_smp_est[, 8] + imp_smp_est[, 10]))
  colnames(raw_smp) <- c("generation", "A0A0", "A0A1", "A1A1", "B0B0", "B0B1", "B1B1")

  return(list(sel_cof_est = sel_cof_est,
              sel_cof_hpd = sel_cof_hpd,
              sel_cof_chn = sel_cof_chn,
              dif_sel_est = dif_sel_est,
              dif_sel_hpd = dif_sel_hpd,
              dif_sel_chn = dif_sel_chn,
              frq_pth_est = frq_pth_est,
              frq_pth_hpd = frq_pth_hpd,
              frq_pth_chn = frq_pth_chn,
              imp_smp_est = imp_smp_est,
              imp_smp_chn = imp_smp_chn))
}
#' Compiled version
cmprunBayesianProcedure <- cmpfun(runBayesianProcedure)

################################################################################
