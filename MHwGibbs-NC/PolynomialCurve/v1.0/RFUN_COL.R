#' @title A Bayesian approach for estimating time-varying selection coefficients from ancient DNA data
#' @author Xiaoyang Dai, Mark Beaumont, Feng Yu, Ludovic Orlando, Zhangyi He

#' version 1.0
#' Polynomial curves of the selection coefficient
#' Horse coat colours (ASIP & MC1R) under constant demographic histories (N/A is not allowed)

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
#' @param sel_cof the selection coefficients of the black and chestnut phenotypes
#' @param rec_rat the recombination rate between the ASIP and MC1R loci
#' @param pop_siz the size of the horse population (constant)
#' @param int_frq the initial haplotype frequencies of the population
#' @param int_gen the generation that the simulated haplotype frequency trajectories started
#' @param lst_gen the generation that the simulated haplotype frequency trajectories ended

#' Standard version
simulateWFM <- function(sel_cof, rec_rat, pop_siz, int_frq, int_gen, lst_gen) {
  WFM <- simulateWFM_arma(sel_cof, rec_rat, pop_siz, int_frq, int_gen, lst_gen)
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
#' @param pop_siz the size of the horse population (constant)
#' @param int_frq the initial haplotype frequencies of the population
#' @param int_gen the generation that the simulated haplotype frequency trajectories started
#' @param lst_gen the generation that the simulated haplotype frequency trajectories ended
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param dat_aug = TRUE/FALSE (return the simulated sample trajectory with data augmentation or not)

#' Standard version
simulateWFD <- function(sel_cof, rec_rat, pop_siz, int_frq, int_gen, lst_gen, ptn_num, dat_aug = TRUE) {
  sel_cof_aug <- matrix(NA, nrow = nrow(sel_cof), ncol = ncol(sel_cof) * ptn_num)
  sel_cof_aug[1, ] <- rep(sel_cof[1, ], each = ptn_num)
  sel_cof_aug[2, ] <- rep(sel_cof[2, ], each = ptn_num)
  sel_cof_aug <- sel_cof_aug[, 1:((ncol(sel_cof) - 1) * ptn_num + 1)]

  frq_pth <- simulateWFD_arma(sel_cof_aug, rec_rat, pop_siz, int_frq, int_gen, lst_gen, ptn_num)
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
#' @param pop_siz the size of the horse population (constant)
#' @param int_con the initial haplotype frequencies of the population / the initial mutant allele frequencies and the linkage disequilibrium of the population
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_siz the count of the horses drawn from the population at all sampling time points
#' @param obs_hap = TRUE/FALSE (return the simulated sample genotypes with haplotype information or not)
#' @param ptn_num the number of the subintervals divided per generation in the Euler-Maruyama method for the WFD

#' Standard version
simulateHMM <- function(model, sel_cof, rec_rat, pop_siz, int_con, smp_gen, smp_siz, obs_hap = FALSE, ...) {
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
    WFM <- cmpsimulateWFM(sel_cof, rec_rat, pop_siz, int_frq, int_gen, lst_gen)
    pop_hap_frq <- as.matrix(WFM$hap_frq_pth)
    pop_gen_frq <- as.matrix(WFM$gen_frq_pth)
  }
  if (model == "WFD") {
    sel_cof_aug <- matrix(NA, nrow = nrow(sel_cof), ncol = ncol(sel_cof) * ptn_num)
    sel_cof_aug[1, ] <- rep(sel_cof[1, ], each = ptn_num)
    sel_cof_aug[2, ] <- rep(sel_cof[2, ], each = ptn_num)
    sel_cof_aug <- sel_cof_aug[, 1:((ncol(sel_cof) - 1) * ptn_num + 1)]

    pop_hap_frq <- cmpsimulateWFD(sel_cof_aug, rec_rat, pop_siz, int_frq, int_gen, lst_gen, ptn_num, dat_aug = FALSE)
    pop_hap_frq <- as.matrix(pop_hap_frq)
    pop_gen_frq <- matrix(NA, nrow = 10, ncol = ncol(pop_hap_frq))
    fts_mat <- calculateFitnessMat_arma(sel_cof_aug)
    for (k in 1:ncol(pop_hap_frq)) {
      hap_frq <- pop_hap_frq[, k]
      gen_frq <- fts_mat * (hap_frq %*% t(hap_frq)) / sum(fts_mat * (hap_frq %*% t(hap_frq)))
      gen_frq[lower.tri(gen_frq, diag = FALSE)] <- NA
      pop_gen_frq[, k] <- discard(as.vector(2 * gen_frq - diag(diag(gen_frq), nrow = 4, ncol = 4)), is.na)
    }
    pop_gen_frq <- as.matrix(pop_gen_frq)
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
#' @param sel_cof the selection coefficients of the black and chestnut phenotypes
#' @param rec_rat the recombination rate between the ASIP and MC1R loci
#' @param pop_siz the size of the horse population (constant)
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_siz the count of the horses drawn from the population at all sampling time points
#' @param smp_cnt the count of the genotypes observed in the sample at all sampling time points
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter

#' Standard version
runBPF <- function(sel_cof, rec_rat, pop_siz, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num) {
  smp_gen <- smp_gen - min(smp_gen)

  sel_cof_aug <- matrix(NA, nrow = nrow(sel_cof), ncol = ncol(sel_cof) * ptn_num)
  sel_cof_aug[1, ] <- rep(sel_cof[1, ], each = ptn_num)
  sel_cof_aug[2, ] <- rep(sel_cof[2, ], each = ptn_num)
  sel_cof_aug <- sel_cof_aug[, 1:((ncol(sel_cof) - 1) * ptn_num + 1)]

  # run the BPF
  BPF <- runBPF_arma(sel_cof_aug, rec_rat, pop_siz, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num)

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
#' @param pop_siz the size of the horse population (constant)
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_siz the count of the horses drawn from the population at all sampling time points
#' @param smp_cnt the count of the genotypes observed in the sample at all sampling time points
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter
#' @param gap_num the number of particles increased or decreased in the optimal particle number search

#' Standard version
calculateOptimalParticleNum <- function(sel_cof, rec_rat, pop_siz, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, gap_num) {
  smp_gen <- smp_gen - min(smp_gen)

  sel_cof_aug <- matrix(NA, nrow = nrow(sel_cof), ncol = ncol(sel_cof) * ptn_num)
  sel_cof_aug[1, ] <- rep(sel_cof[1, ], each = ptn_num)
  sel_cof_aug[2, ] <- rep(sel_cof[2, ], each = ptn_num)
  sel_cof_aug <- sel_cof_aug[, 1:((ncol(sel_cof) - 1) * ptn_num + 1)]

  # calculate the optimal particle number
  OptNum <- calculateOptimalParticleNum_arma(sel_cof_aug, rec_rat, pop_siz, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, gap_num)

  return(list(opt_pcl_num = as.vector(OptNum$opt_pcl_num),
              log_lik_sdv = as.vector(OptNum$log_lik_sdv)))
}
#' Compiled version
cmpcalculateOptimalParticleNum <- cmpfun(calculateOptimalParticleNum)

########################################

#' Run the particle marginal Metropolis-Hastings (PMMH)
#' Parameter settings
#' @param reg_cof
#' @param rec_rat the recombination rate between the ASIP and MC1R loci
#' @param pop_siz the size of the horse population (constant)
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_siz the count of the horses drawn from the population at all sampling time points
#' @param smp_cnt the count of the genotypes observed in the sample at all sampling time points
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter
#' @param itn_num the number of the iterations carried out in the PMMH

#' Standard version
runPMMH <- function(reg_cof, rec_rat, pop_siz, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num) {
  reg_var <- matrix(NA, nrow = 4, ncol = (max(smp_gen) - min(smp_gen)) * ptn_num + 1)
  reg_var[1, ] <- rep(1, length.out = (max(smp_gen) - min(smp_gen)) * ptn_num + 1)
  reg_var[2, ] <- (min(smp_gen) + 0:((max(smp_gen) - min(smp_gen)) * ptn_num) / ptn_num) / (max(smp_gen) - min(smp_gen))
  reg_var[3, ] <- reg_var[2, ]^2
  reg_var[4, ] <- reg_var[2, ]^3

  # run the PMMH
  reg_cof_chn <- runPMMH_arma(reg_cof, reg_var, rec_rat, pop_siz, smp_gen - min(smp_gen), smp_siz, smp_cnt, ptn_num, pcl_num, itn_num)
  reg_cof_chn <- as.array(reg_cof_chn)

  # reg_cof_chn[, 1, ]
  reg_cof_chn[, 2, ] <- reg_cof_chn[, 2, ] / (max(smp_gen) - min(smp_gen))
  reg_cof_chn[, 3, ] <- reg_cof_chn[, 3, ] / (max(smp_gen) - min(smp_gen))^2
  reg_cof_chn[, 4, ] <- reg_cof_chn[, 4, ] / (max(smp_gen) - min(smp_gen))^3

  reg_var <- matrix(NA, nrow = max(smp_gen) - min(smp_gen) + 1, ncol = 4)
  reg_var[, 1] <- rep(1, length.out = max(smp_gen) - min(smp_gen) + 1)
  reg_var[, 2] <- min(smp_gen):max(smp_gen)
  reg_var[, 3] <- (min(smp_gen):max(smp_gen))^2
  reg_var[, 4] <- (min(smp_gen):max(smp_gen))^3

  sel_cof_chn <- array(NA, dim = c(2, max(smp_gen) - min(smp_gen) + 1, itn_num))
  sel_cof_chn[1, , ] <- reg_var %*% reg_cof_chn[1, , ]
  sel_cof_chn[2, , ] <- reg_var %*% reg_cof_chn[2, , ]

  return(list(reg_cof_chn = reg_cof_chn,
              sel_cof_chn = sel_cof_chn))
}
#' Compiled version
cmprunPMMH <- cmpfun(runPMMH)

########################################

#' Run the adaptive particle marginal Metropolis-Hastings (AdaptPMMH)
#' Parameter settings
#' @param reg_cof
#' @param rec_rat the recombination rate between the ASIP and MC1R loci
#' @param pop_siz the size of the horse population (constant)
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_siz the count of the horses drawn from the population at all sampling time points
#' @param smp_cnt the count of the genotypes observed in the sample at all sampling time points
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter
#' @param itn_num the number of the iterations carried out in the PMMH
#' @param stp_siz the step size sequence in the adaptive setting (decaying to zero)
#' @param apt_rto the target mean acceptance probability of the adaptive setting

#' Standard version
runAdaptPMMH <- function(reg_cof, rec_rat, pop_siz, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, stp_siz, apt_rto) {
  reg_var <- matrix(NA, nrow = 4, ncol = (max(smp_gen) - min(smp_gen)) * ptn_num + 1)
  reg_var[1, ] <- rep(1, length.out = (max(smp_gen) - min(smp_gen)) * ptn_num + 1)
  reg_var[2, ] <- (min(smp_gen) + 0:((max(smp_gen) - min(smp_gen)) * ptn_num) / ptn_num) / (max(smp_gen) - min(smp_gen))
  reg_var[3, ] <- reg_var[2, ]^2
  reg_var[4, ] <- reg_var[2, ]^3

  # run the adaptive PMMH
  reg_cof_chn <- runAdaptPMMH_arma(reg_cof, reg_var, rec_rat, pop_siz, smp_gen - min(smp_gen), smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, stp_siz, apt_rto)
  reg_cof_chn <- as.array(reg_cof_chn)

  # reg_cof_chn[, 1, ]
  reg_cof_chn[, 2, ] <- reg_cof_chn[, 2, ] / (max(smp_gen) - min(smp_gen))
  reg_cof_chn[, 3, ] <- reg_cof_chn[, 3, ] / (max(smp_gen) - min(smp_gen))^2
  reg_cof_chn[, 4, ] <- reg_cof_chn[, 4, ] / (max(smp_gen) - min(smp_gen))^3

  reg_var <- matrix(NA, nrow = max(smp_gen) - min(smp_gen) + 1, ncol = 4)
  reg_var[, 1] <- rep(1, length.out = max(smp_gen) - min(smp_gen) + 1)
  reg_var[, 2] <- min(smp_gen):max(smp_gen)
  reg_var[, 3] <- reg_var[, 2]^2
  reg_var[, 4] <- reg_var[, 2]^3

  sel_cof_chn <- array(NA, dim = c(2, max(smp_gen) - min(smp_gen) + 1, itn_num))
  sel_cof_chn[1, , ] <- reg_var %*% reg_cof_chn[1, , ]
  sel_cof_chn[2, , ] <- reg_var %*% reg_cof_chn[2, , ]

  return(list(reg_cof_chn = reg_cof_chn,
              sel_cof_chn = sel_cof_chn))
}
#' Compiled version
cmprunAdaptPMMH <- cmpfun(runAdaptPMMH)

########################################

#' Run the Bayesian procedure for the inference of natural selection
#' Parameter settings
#' @param reg_cof
#' @param rec_rat the recombination rate between the ASIP and MC1R loci
#' @param pop_siz the size of the horse population (constant)
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_siz the count of the horses drawn from the population at all sampling time points
#' @param smp_cnt the count of the genotypes observed in the sample at all sampling time points
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter
#' @param itn_num the number of the iterations carried out in the PMMH
#' @param brn_num the number of the iterations for burn-in
#' @param thn_num the number of the iterations for thinning
#' @param adp_set = TRUE/FALSE (return the result with the adaptive setting or not)
#' @param stp_siz the step size sequence in the adaptive setting (decaying to zero)
#' @param apt_rto the target mean acceptance probability of the adaptive setting

#' Standard version
runBayesianProcedure <- function(reg_cof, rec_rat, pop_siz, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, brn_num, thn_num, adp_set, ...) {
  reg_var <- matrix(NA, nrow = 4, ncol = (max(smp_gen) - min(smp_gen)) * ptn_num + 1)
  reg_var[1, ] <- rep(1, length.out = (max(smp_gen) - min(smp_gen)) * ptn_num + 1)
  reg_var[2, ] <- (min(smp_gen) + 0:((max(smp_gen) - min(smp_gen)) * ptn_num) / ptn_num) / (max(smp_gen) - min(smp_gen))
  reg_var[3, ] <- reg_var[2, ]^2
  reg_var[4, ] <- reg_var[2, ]^3

  if (adp_set == TRUE) {
    # run the adaptive PMMH
    reg_cof_chn <- runAdaptPMMH_arma(reg_cof, reg_var, rec_rat, pop_siz, smp_gen - min(smp_gen), smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, stp_siz, apt_rto)
  } else {
    # run the PMMH
    reg_cof_chn <- runPMMH_arma(reg_cof, reg_var, rec_rat, pop_siz, smp_gen - min(smp_gen), smp_siz, smp_cnt, ptn_num, pcl_num, itn_num)
  }
  reg_cof_chn <- as.array(reg_cof_chn)

  # reg_cof_chn[, 1, ]
  reg_cof_chn[, 2, ] <- reg_cof_chn[, 2, ] / (max(smp_gen) - min(smp_gen))
  reg_cof_chn[, 3, ] <- reg_cof_chn[, 3, ] / (max(smp_gen) - min(smp_gen))^2
  reg_cof_chn[, 4, ] <- reg_cof_chn[, 4, ] / (max(smp_gen) - min(smp_gen))^3

  reg_var <- matrix(NA, nrow = max(smp_gen) - min(smp_gen) + 1, ncol = 4)
  reg_var[, 1] <- rep(1, length.out = max(smp_gen) - min(smp_gen) + 1)
  reg_var[, 2] <- min(smp_gen):max(smp_gen)
  reg_var[, 3] <- reg_var[, 2]^2
  reg_var[, 4] <- reg_var[, 2]^3

  # selection coefficients
  sel_cof_chn <- array(NA, dim = c(2, max(smp_gen) - min(smp_gen) + 1, itn_num))
  sel_cof_chn[1, , ] <- reg_var %*% reg_cof_chn[1, , ]
  sel_cof_chn[2, , ] <- reg_var %*% reg_cof_chn[2, , ]

  sel_cof_chn <- sel_cof_chn[, , brn_num:dim(sel_cof_chn)[3]]
  sel_cof_chn <- sel_cof_chn[, , (1:round(dim(sel_cof_chn)[3] / thn_num)) * thn_num]

  sel_cof_est <- matrix(NA, nrow = 2, ncol = max(smp_gen) - min(smp_gen) + 1)
  for (i in 1:2) {
    sel_cof_est[i, ] <- rowMeans(sel_cof_chn[i, , ])
  }
  sel_cof_hpd <- array(NA, dim = c(2, max(smp_gen) - min(smp_gen) + 1, 2))
  for (i in 1:2) {
    for (k in 1:(max(smp_gen) - min(smp_gen) + 1)) {
      sel_cof_hpd[i, k, ] <- HPDinterval(as.mcmc(sel_cof_chn[i, k, ]), prob = 0.95)
    }
  }

  # Taylor series coefficients
  reg_cof_chn <- reg_cof_chn[, , brn_num:dim(reg_cof_chn)[3]]
  reg_cof_chn <- reg_cof_chn[, , (1:round(dim(reg_cof_chn)[3] / thn_num)) * thn_num]

  reg_cof_est <- matrix(NA, nrow = 2, ncol = 4)
  for (i in 1:2) {
    reg_cof_est[i, ] <- rowMeans(reg_cof_chn[i, , ])
  }

  reg_cof_hpd <- array(NA, dim = c(2, 2, 4))
  for (i in 1:2) {
    for (j in 1:4) {
      reg_cof_hpd[i, j, ] <- HPDinterval(as.mcmc(reg_cof_chn[i, j, ]), prob = 0.95)
    }
  }

  return(list(reg_cof_est = reg_cof_est,
              reg_cof_hpd = reg_cof_hpd,
              reg_cof_chn = reg_cof_chn,
              sel_cof_est = sel_cof_est,
              sel_cof_hpd = sel_cof_hpd,
              sel_cof_chn = sel_cof_chn))
}
#' Compiled version
cmprunBayesianProcedure <- cmpfun(runBayesianProcedure)

################################################################################
