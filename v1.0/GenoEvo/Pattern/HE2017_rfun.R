#' @title Inferring fluctuating selection acting on horse coat colours and patterns from ancient DNA data
#' @author Zhangyi He, Xiaoyang Dai, Mark Beaumont and Feng Yu

#' Genotype evolution (Wright-Fisher diffusion)

#' Horse base coat colour: KIT13 and KIT16 (constant selection)

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
sourceCpp("./Code/Code v1.0/GenoEvo/Horse coat pattern/HE2017_cfun.cpp")

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

#' Standard version
simulateTLWFMS <- function(type, phased, sel_cof, rec_rat, pop_siz, int_frq, int_gen, lst_gen) {
  if (type == "allele") {
    frq_pth <- simulateTLWFMS_arma(sel_cof, rec_rat, pop_siz, int_frq, int_gen, lst_gen)
    frq_pth <- convertPhasedGenoPath2AllelePath_arma(frq_pth)
  } else if (type == "haplotype") {
    frq_pth <- simulateTLWFMS_arma(sel_cof, rec_rat, pop_siz, int_frq, int_gen, lst_gen)
    frq_pth <- convertPhasedGenoPath2HaploPath_arma(frq_pth)
  } else {
    if (phased == FALSE) {
      frq_pth <- simulateTLWFMS_arma(sel_cof, rec_rat, pop_siz, int_frq, int_gen, lst_gen)
      frq_pth <- convertPhasedGenoPath2UnphasedGenoPath_arma(frq_pth)
    } else {
      frq_pth <- simulateTLWFMS_arma(sel_cof, rec_rat, pop_siz, int_frq, int_gen, lst_gen)
    }
  }
  
  return(frq_pth)
}
#' Compiled version
cmpsimulateTLWFMS <- cmpfun(simulateTLWFMS)

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

#' Standard version
simulateTLWFDS <- function(type, phased, sel_cof, rec_rat, pop_siz, int_frq, int_gen, lst_gen, ptn_num) {
  if (type == "allele") {
    frq_pth <- simulateTLWFDS_arma(sel_cof, rec_rat, pop_siz, int_frq, int_gen, lst_gen, ptn_num)
    frq_pth <- convertPhasedGenoPath2AllelePath_arma(frq_pth)
  } else if (type == "haplotype") {
    frq_pth <- simulateTLWFDS_arma(sel_cof, rec_rat, pop_siz, int_frq, int_gen, lst_gen, ptn_num)
    frq_pth <- convertPhasedGenoPath2HaploPath_arma(frq_pth)
  } else {
    if (phased == FALSE) {
      frq_pth <- simulateTLWFDS_arma(sel_cof, rec_rat, pop_siz, int_frq, int_gen, lst_gen, ptn_num)
      frq_pth <- convertPhasedGenoPath2UnphasedGenoPath_arma(frq_pth)
    } else {
      frq_pth <- simulateTLWFDS_arma(sel_cof, rec_rat, pop_siz, int_frq, int_gen, lst_gen, ptn_num)
    }
  }
  
  return(frq_pth)
}
#' Compiled version
cmpsimulateTLWFDS <- cmpfun(simulateTLWFDS)

########################################

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

#' Standard version
simulateHMM <- function(model, missing, sel_cof, rec_rat, pop_siz, int_frq, smp_gen, smp_siz, ...) {
  int_gen <- min(smp_gen)
  lst_gen <- max(smp_gen)
  
  # generate the population genotype frequency trajectories
  if (model == "WFM") {
    pop_frq <- cmpsimulateTLWFMS(type = "genotype", phased = TRUE, sel_cof, rec_rat, pop_siz, int_frq, int_gen, lst_gen)
  }
  if (model == "WFD") {
    pop_frq <- cmpsimulateTLWFDS(type = "genotype", phased = TRUE, sel_cof, rec_rat, pop_siz, int_frq, int_gen, lst_gen, ptn_num)
  }
  
  smp_cnt <- matrix(0, nrow = 9, ncol = length(smp_gen))
  mis_cnt <- matrix(0, nrow = 27, ncol = length(smp_gen))
  if (missing == TRUE) {
    for (k in 1:length(smp_gen)) {
      smp_gen_cnt <- rmultinom(1, size = smp_siz[k], prob = pop_frq[, smp_gen[k] - int_gen + 1])
      
      mis_gen_cnt <- sample(0:round(smp_gen_cnt[1] * 0.1 / k), size = 1)
      smp_cnt[1, k] <- smp_gen_cnt[1] - mis_gen_cnt
      mis_cnt[, k] <- mis_cnt[, k] + rmultinom(1, size = mis_gen_cnt, prob = c(0.4, 0, 0, 0, 0, 0, 0.4, 0, 0, 0, 0, 0, 0.05, 0, 0, 0.05, 0, 0, 0.05, 0, 0, 0, 0.02, 0, 0.02, 0, 0.01))
      
      mis_gen_cnt <- sample(0:round(smp_gen_cnt[2] * 0.1 / k), size = 1)
      smp_cnt[2, k] <- smp_gen_cnt[2] - mis_gen_cnt
      mis_cnt[, k] <- mis_cnt[, k] + rmultinom(1, size = mis_gen_cnt, prob = c(0, 0, 0.23, 0, 0, 0, 0.23, 0.23, 0, 0, 0, 0, 0, 0.06, 0, 0.06, 0, 0, 0.06, 0, 0.06, 0, 0.02, 0.02, 0.02, 0, 0.01))
      
      mis_gen_cnt <- sample(0:round(smp_gen_cnt[3] * 0.1 / k), size = 1)
      smp_cnt[3, k] <- smp_gen_cnt[3] - mis_gen_cnt
      mis_cnt[, k] <- mis_cnt[, k] + rmultinom(1, size = mis_gen_cnt, prob = c(0, 0, 0, 0, 0.4, 0, 0, 0.4, 0, 0, 0, 0, 0, 0, 0.05, 0.05, 0, 0, 0, 0, 0.05, 0, 0, 0.02, 0.02, 0, 0.01))
      
      mis_gen_cnt <- sample(0:round(smp_gen_cnt[4] * 0.1 / k), size = 1)
      smp_cnt[4, k] <- smp_gen_cnt[4] - mis_gen_cnt
      mis_cnt[, k] <- mis_cnt[, k] + rmultinom(1, size = mis_gen_cnt, prob = c(0.23, 0.23, 0, 0, 0, 0, 0, 0, 0.23, 0, 0, 0, 0.06, 0, 0, 0.06, 0, 0, 0.06, 0.06, 0, 0, 0.02, 0, 0.02, 0.02, 0.01))
      
      mis_gen_cnt <- sample(0:round((smp_gen_cnt[5] + smp_gen_cnt[10]) * 0.1 / k), size = 1)
      smp_cnt[5, k] <- (smp_gen_cnt[5] + smp_gen_cnt[10]) - mis_gen_cnt
      mis_cnt[, k] <- mis_cnt[, k] + rmultinom(1, size = mis_gen_cnt, prob = c(0, 0, 0.1825, 0.1825, 0, 0, 0, 0, 0.1825, 0.1825, 0, 0, 0, 0.03, 0, 0, 0.03, 0, 0.03, 0.03, 0.03, 0.03, 0.02, 0.02, 0.02, 0.02, 0.01))
      
      mis_gen_cnt <- sample(0:round(smp_gen_cnt[6] * 0.1 / k), size = 1)
      smp_cnt[6, k] <- smp_gen_cnt[6] - mis_gen_cnt
      mis_cnt[, k] <- mis_cnt[, k] + rmultinom(1, size = mis_gen_cnt, prob = c(0, 0, 0, 0, 0.23, 0.23, 0, 0, 0, 0.23, 0, 0, 0, 0, 0.06, 0, 0.06, 0, 0, 0, 0.06, 0.06, 0, 0.02, 0.02, 0.02, 0.01))
      
      mis_gen_cnt <- sample(0:round(smp_gen_cnt[7] * 0.1 / k), size = 1)
      smp_cnt[7, k] <- smp_gen_cnt[7] - mis_gen_cnt
      mis_cnt[, k] <- mis_cnt[, k] + rmultinom(1, size = mis_gen_cnt, prob = c(0, 0.4, 0, 0, 0, 0, 0, 0, 0, 0, 0.4, 0, 0.05, 0, 0, 0, 0, 0.05, 0, 0.05, 0, 0, 0.02, 0, 0, 0.02, 0.01))
      
      mis_gen_cnt <- sample(0:round(smp_gen_cnt[8] * 0.1 / k), size = 1)
      smp_cnt[8, k] <- smp_gen_cnt[8] - mis_gen_cnt
      mis_cnt[, k] <- mis_cnt[, k] + rmultinom(1, size = mis_gen_cnt, prob = c(0, 0, 0, 0.23, 0, 0, 0, 0, 0, 0, 0.23, 0.23, 0, 0.06, 0, 0, 0, 0.06, 0, 0.06, 0, 0.06, 0.02, 0.02, 0, 0.02, 0.01))
      
      mis_gen_cnt <- sample(0:round(smp_gen_cnt[9] * 0.1 / k), size = 1)
      smp_cnt[9, k] <- smp_gen_cnt[9] - mis_gen_cnt
      mis_cnt[, k] <- mis_cnt[, k] + rmultinom(1, size = mis_gen_cnt, prob = c(0, 0, 0, 0, 0, 0.4, 0, 0, 0, 0, 0, 0.4, 0, 0, 0.05, 0, 0, 0.05, 0, 0, 0, 0.05, 0, 0.02, 0, 0.02, 0.01))
    }
  } else {
    for (k in 1:length(smp_gen)) {
      smp_gen_cnt <- rmultinom(1, size = smp_siz[k], prob = pop_frq[, smp_gen[k] - int_gen + 1])
      
      smp_cnt[1, k] <- smp_gen_cnt[1]
      smp_cnt[2, k] <- smp_gen_cnt[2]
      smp_cnt[3, k] <- smp_gen_cnt[3]
      smp_cnt[4, k] <- smp_gen_cnt[4]
      smp_cnt[5, k] <- smp_gen_cnt[5] + smp_gen_cnt[10]
      smp_cnt[6, k] <- smp_gen_cnt[6]
      smp_cnt[7, k] <- smp_gen_cnt[7]
      smp_cnt[8, k] <- smp_gen_cnt[8]
      smp_cnt[9, k] <- smp_gen_cnt[9]
    }
  }
  
  return(list(smp_gen = smp_gen, 
              smp_siz = smp_siz, 
              smp_cnt = smp_cnt, 
              mis_cnt = mis_cnt, 
              pop_frq = pop_frq))
}
#' Compiled version
cmpsimulateHMM <- cmpfun(simulateHMM)

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

#' Standard version
runBPF <- function(sel_cof, rec_rat, pop_siz, smp_gen, smp_siz, smp_cnt, mis_cnt, ptn_num, pcl_num) {
  # run the BPF
  BPF <- runBPF_arma(sel_cof, rec_rat, pop_siz, smp_gen, smp_siz, smp_cnt, mis_cnt, ptn_num, pcl_num)
  
  lik = BPF$lik
  wght = BPF$wght
  pop_gen_frq_pre_resmp = BPF$part_pre_resmp
  pop_gen_frq_pst_resmp = BPF$part_pst_resmp
  
  # convert haplotype particles to genotype particles
  pop_frq_pre_resmp <- array(NA, dim = c(9, pcl_num, dim(pop_gen_frq_pre_resmp)[3]))
  pop_frq_pst_resmp <- array(NA, dim = c(9, pcl_num, dim(pop_gen_frq_pst_resmp)[3]))
  for (i in 1:pcl_num) {
    pop_frq_pre_resmp[, i, ] <- convertPhasedGenoPath2UnphasedGenoPath_arma(pop_gen_frq_pre_resmp[, i, ])
    pop_frq_pst_resmp[, i, ] <- convertPhasedGenoPath2UnphasedGenoPath_arma(pop_gen_frq_pst_resmp[, i, ])
  }
  
  return(list(lik = lik, 
              wght = wght, 
              pop_frq_pre_resmp = pop_frq_pre_resmp, 
              pop_frq_pst_resmp = pop_frq_pst_resmp))
}
#' Compiled version
cmprunBPF <- cmpfun(runBPF)

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

#' Standard version
calculateOptimalParticleNum <- function(sel_cof, rec_rat, pop_siz, smp_gen, smp_siz, smp_cnt, mis_cnt, ptn_num, pcl_num, gap_num) {
  OptNum <- calculateOptimalParticleNum_arma(sel_cof, rec_rat, pop_siz, smp_gen, smp_siz, smp_cnt, mis_cnt, ptn_num, pcl_num, gap_num)
  
  return(list(opt_pcl_num = as.vector(OptNum$opt_pcl_num), 
              log_lik_sdv = as.vector(OptNum$log_lik_sdv)))
}
#' Compiled version
cmpcalculateOptimalParticleNum <- cmpfun(calculateOptimalParticleNum)

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

#' Standard version
runPMMH <- function(sel_cof, rec_rat, pop_siz, smp_gen, smp_siz, smp_cnt, mis_cnt, ptn_num, pcl_num, itn_num, ER = TRUE, PA = TRUE, ...) {
  # run the PMMH
  if (ER == TRUE) {
    if (PA == TRUE) {
      PMMH <- runERPAPMMH_arma(sel_cof, rec_rat, pop_siz, smp_gen, smp_siz, smp_cnt, mis_cnt, ptn_num, pcl_num, itn_num, nap_num)
    } else {
      PMMH <- runERPMMH_arma(sel_cof, rec_rat, pop_siz, smp_gen, smp_siz, smp_cnt, mis_cnt, ptn_num, pcl_num, itn_num)
    }
  } else {
    if (PA == TRUE) {
      PMMH <- runPAPMMH_arma(sel_cof, rec_rat, pop_siz, smp_gen, smp_siz, smp_cnt, mis_cnt, ptn_num, pcl_num, itn_num, nap_num)
    } else {
      PMMH <- runPMMH_arma(sel_cof, rec_rat, pop_siz, smp_gen, smp_siz, smp_cnt, mis_cnt, ptn_num, pcl_num, itn_num)
    }
  }
  
  return(list(sel_cof_to_chn = as.vector(PMMH[1, ]), 
              sel_cof_sb_chn = as.vector(PMMH[2, ]), 
              sel_cof_ts_chn = as.vector(PMMH[3, ])))
}
#' Compiled version
cmprunPMMH <- cmpfun(runPMMH)

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

#' Standard version
runBayesianProcedure <- function(sel_cof, rec_rat, pop_siz, smp_gen, smp_siz, smp_cnt, mis_cnt, ptn_num, pcl_num, itn_num, ER = TRUE, PA = TRUE, brn_num, thn_num, est = "MMSE", ...) {
  # run the PMMH
  if (ER == TRUE) {
    if (PA == TRUE) {
      PMMH <- runERPAPMMH_arma(sel_cof, rec_rat, pop_siz, smp_gen, smp_siz, smp_cnt, mis_cnt, ptn_num, pcl_num, itn_num, nap_num)
    } else {
      PMMH <- runERPMMH_arma(sel_cof, rec_rat, pop_siz, smp_gen, smp_siz, smp_cnt, mis_cnt, ptn_num, pcl_num, itn_num)
    }
  } else {
    if (PA == TRUE) {
      PMMH <- runPAPMMH_arma(sel_cof, rec_rat, pop_siz, smp_gen, smp_siz, smp_cnt, mis_cnt, ptn_num, pcl_num, itn_num, nap_num)
    } else {
      PMMH <- runPMMH_arma(sel_cof, rec_rat, pop_siz, smp_gen, smp_siz, smp_cnt, mis_cnt, ptn_num, pcl_num, itn_num)
    }
  }
  
  # burn-in and thinning
  sel_cof_to_chn <- as.vector(PMMH[1, ])
  sel_cof_to_chn <- sel_cof_to_chn[brn_num:length(sel_cof_to_chn)]
  sel_cof_to_chn <- sel_cof_to_chn[(1:round(length(sel_cof_to_chn) / thn_num)) * thn_num]
  sel_cof_sb_chn <- as.vector(PMMH[2, ])
  sel_cof_sb_chn <- sel_cof_sb_chn[brn_num:length(sel_cof_sb_chn)]
  sel_cof_sb_chn <- sel_cof_sb_chn[(1:round(length(sel_cof_sb_chn) / thn_num)) * thn_num]
  sel_cof_ts_chn <- as.vector(PMMH[3, ])
  sel_cof_ts_chn <- sel_cof_ts_chn[brn_num:length(sel_cof_ts_chn)]
  sel_cof_ts_chn <- sel_cof_ts_chn[(1:round(length(sel_cof_ts_chn) / thn_num)) * thn_num]
  
  if (est == "MAP") {
    # MAP estimates for the selection coefficients
    if (max(length(sel_cof_to_chn), length(sel_cof_sb_chn), length(sel_cof_ts_chn)) < 1e+05) {
      sel_cof_pdf <- kde3d(sel_cof_to_chn, sel_cof_sb_chn, sel_cof_ts_chn, n = grd_num)
      sel_cof_to_grd <- sel_cof_pdf$x
      sel_cof_sb_grd <- sel_cof_pdf$y
      sel_cof_ts_grd <- sel_cof_pdf$z
    } else {
      sel_cof_pdf <- kde3d(tail(sel_cof_to_chn, 1e+05), tail(sel_cof_sb_chn, 1e+05), tail(sel_cof_ts_chn, 1e+05), n = grd_num)
      sel_cof_to_grd <- sel_cof_pdf$x
      sel_cof_sb_grd <- sel_cof_pdf$y
      sel_cof_ts_grd <- sel_cof_pdf$z
    }
    sel_cof_to_map <- sel_cof_to_grd[which(sel_cof_pdf$d == max(sel_cof_pdf$d), arr.ind = TRUE)[1]]
    sel_cof_sb_map <- sel_cof_sb_grd[which(sel_cof_pdf$d == max(sel_cof_pdf$d), arr.ind = TRUE)[2]]
    sel_cof_ts_map <- sel_cof_ts_grd[which(sel_cof_pdf$d == max(sel_cof_pdf$d), arr.ind = TRUE)[3]]
    
    # 95% HPD intervals for the selection coefficients
    sel_cof_to_hpd <- HPDinterval(as.mcmc(sel_cof_to_chn), prob = 0.95)
    sel_cof_sb_hpd <- HPDinterval(as.mcmc(sel_cof_sb_chn), prob = 0.95)
    sel_cof_ts_hpd <- HPDinterval(as.mcmc(sel_cof_ts_chn), prob = 0.95)
    
    return(list(sel_cof_to_map = sel_cof_to_map, 
                sel_cof_sb_map = sel_cof_sb_map, 
                sel_cof_ts_map = sel_cof_ts_map, 
                sel_cof_to_hpd = sel_cof_to_hpd, 
                sel_cof_sb_hpd = sel_cof_sb_hpd, 
                sel_cof_ts_hpd = sel_cof_ts_hpd, 
                sel_cof_to_chn = sel_cof_to_chn, 
                sel_cof_sb_chn = sel_cof_sb_chn, 
                sel_cof_ts_chn = sel_cof_ts_chn))
  }
  
  if (est == "MMSE") {
    # MMSE estimates for the selection coefficients
    sel_cof_to_mmse <- mean(sel_cof_to_chn)
    sel_cof_sb_mmse <- mean(sel_cof_sb_chn)
    sel_cof_ts_mmse <- mean(sel_cof_ts_chn)
    
    # 95% HPD intervals for the selection coefficients
    sel_cof_to_hpd <- HPDinterval(as.mcmc(sel_cof_to_chn), prob = 0.95)
    sel_cof_sb_hpd <- HPDinterval(as.mcmc(sel_cof_sb_chn), prob = 0.95)
    sel_cof_ts_hpd <- HPDinterval(as.mcmc(sel_cof_ts_chn), prob = 0.95)
    
    return(list(sel_cof_to_mmse = sel_cof_to_mmse, 
                sel_cof_sb_mmse = sel_cof_sb_mmse, 
                sel_cof_ts_mmse = sel_cof_ts_mmse, 
                sel_cof_to_hpd = sel_cof_to_hpd, 
                sel_cof_sb_hpd = sel_cof_sb_hpd, 
                sel_cof_ts_hpd = sel_cof_ts_hpd, 
                sel_cof_to_chn = sel_cof_to_chn, 
                sel_cof_sb_chn = sel_cof_sb_chn, 
                sel_cof_ts_chn = sel_cof_ts_chn))
  }
  
  if (est == "both") {
    # MAP estimates for the selection coefficients
    if (max(length(sel_cof_to_chn), length(sel_cof_sb_chn), length(sel_cof_ts_chn)) < 1e+05) {
      sel_cof_pdf <- kde3d(sel_cof_to_chn, sel_cof_sb_chn, sel_cof_ts_chn, n = grd_num)
      sel_cof_to_grd <- sel_cof_pdf$x
      sel_cof_sb_grd <- sel_cof_pdf$y
      sel_cof_ts_grd <- sel_cof_pdf$z
    } else {
      sel_cof_pdf <- kde3d(tail(sel_cof_to_chn, 1e+05), tail(sel_cof_sb_chn, 1e+05), tail(sel_cof_ts_chn, 1e+05), n = grd_num)
      sel_cof_to_grd <- sel_cof_pdf$x
      sel_cof_sb_grd <- sel_cof_pdf$y
      sel_cof_ts_grd <- sel_cof_pdf$z
    }
    sel_cof_to_map <- sel_cof_to_grd[which(sel_cof_pdf$d == max(sel_cof_pdf$d), arr.ind = TRUE)[1]]
    sel_cof_sb_map <- sel_cof_sb_grd[which(sel_cof_pdf$d == max(sel_cof_pdf$d), arr.ind = TRUE)[2]]
    sel_cof_ts_map <- sel_cof_ts_grd[which(sel_cof_pdf$d == max(sel_cof_pdf$d), arr.ind = TRUE)[3]]
    
    # MMSE estimates for the selection coefficients
    sel_cof_to_mmse <- mean(sel_cof_to_chn)
    sel_cof_sb_mmse <- mean(sel_cof_sb_chn)
    sel_cof_ts_mmse <- mean(sel_cof_ts_chn)
    
    # 95% HPD intervals for the selection coefficients
    sel_cof_to_hpd <- HPDinterval(as.mcmc(sel_cof_to_chn), prob = 0.95)
    sel_cof_sb_hpd <- HPDinterval(as.mcmc(sel_cof_sb_chn), prob = 0.95)
    sel_cof_ts_hpd <- HPDinterval(as.mcmc(sel_cof_ts_chn), prob = 0.95)
    
    return(list(sel_cof_to_map = sel_cof_to_map, 
                sel_cof_sb_map = sel_cof_sb_map, 
                sel_cof_ts_map = sel_cof_ts_map, 
                sel_cof_to_mmse = sel_cof_to_mmse, 
                sel_cof_sb_mmse = sel_cof_sb_mmse, 
                sel_cof_ts_mmse = sel_cof_ts_mmse, 
                sel_cof_to_hpd = sel_cof_to_hpd, 
                sel_cof_sb_hpd = sel_cof_sb_hpd, 
                sel_cof_ts_hpd = sel_cof_ts_hpd, 
                sel_cof_to_chn = sel_cof_to_chn, 
                sel_cof_sb_chn = sel_cof_sb_chn, 
                sel_cof_ts_chn = sel_cof_ts_chn))
  }
}
#' Compiled version
cmprunBayesianProcedure <- cmpfun(runBayesianProcedure)

################################################################################
