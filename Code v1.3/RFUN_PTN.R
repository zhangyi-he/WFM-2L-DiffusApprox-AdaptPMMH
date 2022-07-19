#' @title Estimating temporally variable selection intensity from ancient DNA data with the flexibility of modelling linkage and epistasis
#' @author Zhangyi He, Xiaoyang Dai, Wenyang Lyu, Mark Beaumont, Feng Yu

#' version 1.3
#' Phenotypes controlled by two genes with genetic linkage
#' Non-constant natural selection and non-constant demographic histories

#' Use the flat Dirichlet prior for the starting haplotype frequencies of the underlying population

#' Input: genotype likelihoods
#' Output: posteriors for the selection coefficient and the genotype frequency trajectories of the population

#' Horse pinto coat patterns (KIT13 & KIT116)

#' R functions

#install.packages("gtools")
library("gtools")

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
sourceCpp("./CFUN_PTN.cpp")

################################################################################

#' Simulate the haplotype frequency trajectories according to the two-locus Wright-Fisher model with selection
#' Parameter setting
#' @param sel_cof the selection coefficients of the tobiano, sabino and mixed against solid
#' @param rec_rat the recombination rate between the KIT13 and KIT16 loci
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
    WFM <- simulateWFM_arma(fts_mat, rec_rat, pop_siz[1:(evt_gen - int_gen + 1)], int_frq, int_gen, evt_gen)
    hap_frq_pth_pre_evt <- as.matrix(WFM$hap_frq_pth)
    gen_frq_pth_pre_evt <- as.matrix(WFM$gen_frq_pth)

    fts_mat <- calculateFitnessMat_arma(sel_cof[, 2])
    WFM <- simulateWFM_arma(fts_mat, rec_rat, pop_siz[(evt_gen - int_gen + 1):(lst_gen - int_gen + 1)], hap_frq_pth_pre_evt[, ncol(hap_frq_pth_pre_evt)], evt_gen, lst_gen)
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
#' @param sel_cof the selection coefficients of the tobiano, sabino and mixed against solid
#' @param rec_rat the recombination rate between the KIT13 and KIT16 loci
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
    hap_frq_pth_pre_evt <- simulateWFD_arma(sel_cof[, 1], rec_rat, pop_siz[1:(evt_gen - int_gen + 1)], ref_siz, int_frq, int_gen, evt_gen, ptn_num)
    hap_frq_pth_pre_evt <- as.matrix(hap_frq_pth_pre_evt)

    hap_frq_pth_pst_evt <- simulateWFD_arma(sel_cof[, 2], rec_rat, pop_siz[(evt_gen - int_gen + 1):(lst_gen - int_gen + 1)], ref_siz, hap_frq_pth_pre_evt[, ncol(hap_frq_pth_pre_evt)], evt_gen, lst_gen, ptn_num)
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
#' @param sel_cof the selection coefficients of the tobiano, sabino and mixed against solid
#' @param rec_rat the recombination rate between the KIT13 and KIT16 loci
#' @param pop_siz the size of the horse population (non-constant)
#' @param int_con the initial haplotype frequencies of the population / the initial mutant allele frequencies and the linkage disequilibrium of the population
#' @param evt_gen the generation that the event of interest occurred
#' @param smp_lab the identifier of the sample assigned
#' @param smp_gen the generation of the sample drawn
#' @param smp_qua the quality of the sample tested
#' @param thr_val the threshold for genotype calling
#' @param ref_siz the reference size of the horse population
#' @param ptn_num the number of the subintervals divided per generation in the Euler-Maruyama method for the WFD

#' Standard version
simulateHMM <- function(model, sel_cof, rec_rat, pop_siz, int_con, evt_gen, smp_lab, smp_gen, smp_qua, thr_val, ...) {
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
    hap_frq <- WFM$hap_frq_pth
    gen_frq <- WFM$gen_frq_pth
  }
  if (model == "WFD") {
    hap_frq <- cmpsimulateWFD(sel_cof, rec_rat, pop_siz, ref_siz, int_frq, evt_gen, int_gen, lst_gen, ptn_num, dat_aug = FALSE)
    hap_frq <- as.matrix(hap_frq)
    gen_frq <- matrix(NA, nrow = 10, ncol = ncol(hap_frq))
    if (evt_gen >= lst_gen) {
      fts_mat <- calculateFitnessMat_arma(sel_cof[, 1])
      for (k in 1:ncol(hap_frq)) {
        pop_hap_frq <- hap_frq[, k]
        pop_gen_frq <- fts_mat * (pop_hap_frq %*% t(pop_hap_frq)) / sum(fts_mat * (pop_hap_frq %*% t(pop_hap_frq)))
        pop_gen_frq[lower.tri(pop_gen_frq, diag = FALSE)] <- NA
        gen_frq[, k] <- discard(as.vector(2 * pop_gen_frq - diag(diag(pop_gen_frq), nrow = 4, ncol = 4)), is.na)
      }
    } else if (evt_gen < int_gen) {
      fts_mat <- calculateFitnessMat_arma(sel_cof[, 2])
      for (k in 1:ncol(hap_frq)) {
        pop_hap_frq <- hap_frq[, k]
        pop_gen_frq <- fts_mat * (pop_hap_frq %*% t(pop_hap_frq)) / sum(fts_mat * (pop_hap_frq %*% t(pop_hap_frq)))
        pop_gen_frq[lower.tri(pop_gen_frq, diag = FALSE)] <- NA
        gen_frq[, k] <- discard(as.vector(2 * pop_gen_frq - diag(diag(pop_gen_frq), nrow = 4, ncol = 4)), is.na)
      }
    } else {
      fts_mat <- calculateFitnessMat_arma(sel_cof[, 1])
      for (k in 1:(evt_gen - int_gen + 1)) {
        pop_hap_frq <- hap_frq[, k]
        pop_gen_frq <- fts_mat * (pop_hap_frq %*% t(pop_hap_frq)) / sum(fts_mat * (pop_hap_frq %*% t(pop_hap_frq)))
        pop_gen_frq[lower.tri(pop_gen_frq, diag = FALSE)] <- NA
        gen_frq[, k] <- discard(as.vector(2 * pop_gen_frq - diag(diag(pop_gen_frq), nrow = 4, ncol = 4)), is.na)
      }

      fts_mat <- calculateFitnessMat_arma(sel_cof[, 2])
      for (k in (evt_gen - int_gen + 2):ncol(hap_frq)) {
        pop_hap_frq <- hap_frq[, k]
        pop_gen_frq <- fts_mat * (pop_hap_frq %*% t(pop_hap_frq)) / sum(fts_mat * (pop_hap_frq %*% t(pop_hap_frq)))
        pop_gen_frq[lower.tri(pop_gen_frq, diag = FALSE)] <- NA
        gen_frq[, k] <- discard(as.vector(2 * pop_gen_frq - diag(diag(pop_gen_frq), nrow = 4, ncol = 4)), is.na)
      }
    }
    gen_frq <- as.matrix(gen_frq)
  }

  # generate the sample individual genotypes at all sampling time points
  tru_smp <- NULL
  for (k in 1:length(smp_gen)) {
    tru_smp <- cbind(tru_smp, rbind(rep(smp_gen[k], times = 1), rmultinom(1, size = 1, prob = gen_frq[, smp_gen[k] - int_gen + 1])))
  }
  tru_smp[5, ] <- tru_smp[5, ] + tru_smp[7, ]
  tru_smp <- tru_smp[-7, ]
  tru_smp <- as.data.frame(t(tru_smp))
  tru_smp <- as.data.frame(cbind(tru_smp[, 1],
    tru_smp[, 2] + tru_smp[, 3] + tru_smp[, 4],
    tru_smp[, 5] + tru_smp[, 6] + tru_smp[, 8],
    tru_smp[, 7] + tru_smp[, 9] + tru_smp[, 10],
    tru_smp[, 2] + tru_smp[, 5] + tru_smp[, 7],
    tru_smp[, 3] + tru_smp[, 6] + tru_smp[, 9],
    tru_smp[, 4] + tru_smp[, 8] + tru_smp[, 10]))
  rownames(tru_smp) <- smp_lab
  colnames(tru_smp) <- c("generation", "KM0KM0", "KM0KM1", "KM1KM1", "sb1sb1", "sb1SB1", "SB1SB1")
  tru_smp <- tru_smp[order(tru_smp$generation), ]

  raw_smp <- tru_smp
  sel_idx <- which(tru_smp$'KM0KM0' == 1)
  raw_smp[sel_idx, 2:4] <- rdirichlet(length(sel_idx), alpha = c(smp_qua[1], (1 - smp_qua[1]) / 2, (1 - smp_qua[1]) / 2) * smp_qua[2])
  sel_idx <- which(tru_smp$'KM0KM1' == 1)
  raw_smp[sel_idx, 2:4] <- rdirichlet(length(sel_idx), alpha = c((1 - smp_qua[1]) / 2, smp_qua[1], (1 - smp_qua[1]) / 2) * smp_qua[2])
  sel_idx <- which(tru_smp$'KM1KM1' == 1)
  raw_smp[sel_idx, 2:4] <- rdirichlet(length(sel_idx), alpha = c((1 - smp_qua[1]) / 2, (1 - smp_qua[1]) / 2, smp_qua[1]) * smp_qua[2])
  sel_idx <- which(tru_smp$'sb1sb1' == 1)
  raw_smp[sel_idx, 5:7] <- rdirichlet(length(sel_idx), alpha = c(smp_qua[1], (1 - smp_qua[1]) / 2, (1 - smp_qua[1]) / 2) * smp_qua[2])
  sel_idx <- which(tru_smp$'sb1SB1' == 1)
  raw_smp[sel_idx, 5:7] <- rdirichlet(length(sel_idx), alpha = c((1 - smp_qua[1]) / 2, smp_qua[1], (1 - smp_qua[1]) / 2) * smp_qua[2])
  sel_idx <- which(tru_smp$'SB1SB1' == 1)
  raw_smp[sel_idx, 5:7] <- rdirichlet(length(sel_idx), alpha = c((1 - smp_qua[1]) / 2, (1 - smp_qua[1]) / 2, smp_qua[1]) * smp_qua[2])

  cal_smp <- raw_smp
  mis_idx <- 1:nrow(raw_smp)
  sel_idx <- which(raw_smp$'KM0KM0' / raw_smp$'KM0KM1' >= thr_val & raw_smp$'KM0KM0' / raw_smp$'KM1KM1' >= thr_val)
  if (length(sel_idx) > 1) {
    cal_smp[sel_idx, 2:4] <- matrix(rep(c(1, 0, 0), times = length(sel_idx)), nrow = length(sel_idx), ncol = 3, byrow = TRUE)
    mis_idx <- mis_idx[-which(mis_idx %in% sel_idx)]
  }
  sel_idx <- which(raw_smp$'KM0KM1' / raw_smp$'KM0KM0' >= thr_val & raw_smp$'KM0KM1' / raw_smp$'KM1KM1' >= thr_val)
  if (length(sel_idx) > 1) {
    cal_smp[sel_idx, 2:4] <- matrix(rep(c(0, 1, 0), times = length(sel_idx)), nrow = length(sel_idx), ncol = 3, byrow = TRUE)
    mis_idx <- mis_idx[-which(mis_idx %in% sel_idx)]
  }
  sel_idx <- which(raw_smp$'KM1KM1' / raw_smp$'KM0KM0' >= thr_val & raw_smp$'KM1KM1' / raw_smp$'KM0KM1' >= thr_val)
  if (length(sel_idx) > 1) {
    cal_smp[sel_idx, 2:4] <- matrix(rep(c(0, 0, 1), times = length(sel_idx)), nrow = length(sel_idx), ncol = 3, byrow = TRUE)
    mis_idx <- mis_idx[-which(mis_idx %in% sel_idx)]
  }
  if (length(mis_idx) > 1) {
    cal_smp[mis_idx, 2:4] <- matrix(rep(c(NA, NA, NA), times = length(mis_idx)), nrow = length(mis_idx), ncol = 3, byrow = TRUE)
  }
  mis_idx <- 1:nrow(raw_smp)
  sel_idx <- which(raw_smp$'sb1sb1' / raw_smp$'sb1SB1' >= thr_val & raw_smp$'sb1sb1' / raw_smp$'SB1SB1' >= thr_val)
  if (length(sel_idx) > 1) {
    cal_smp[sel_idx, 5:7] <- matrix(rep(c(1, 0, 0), times = length(sel_idx)), nrow = length(sel_idx), ncol = 3, byrow = TRUE)
    mis_idx <- mis_idx[-which(mis_idx %in% sel_idx)]
  }
  sel_idx <- which(raw_smp$'sb1SB1' / raw_smp$'sb1sb1' >= thr_val & raw_smp$'sb1SB1' / raw_smp$'SB1SB1' >= thr_val)
  if (length(sel_idx) > 1) {
    cal_smp[sel_idx, 5:7] <- matrix(rep(c(0, 1, 0), times = length(sel_idx)), nrow = length(sel_idx), ncol = 3, byrow = TRUE)
    mis_idx <- mis_idx[-which(mis_idx %in% sel_idx)]
  }
  sel_idx <- which(raw_smp$'SB1SB1' / raw_smp$'sb1sb1' >= thr_val & raw_smp$'SB1SB1' / raw_smp$'sb1SB1' >= thr_val)
  if (length(sel_idx) > 1) {
    cal_smp[sel_idx, 5:7] <- matrix(rep(c(0, 0, 1), times = length(sel_idx)), nrow = length(sel_idx), ncol = 3, byrow = TRUE)
    mis_idx <- mis_idx[-which(mis_idx %in% sel_idx)]
  }
  if (length(mis_idx) > 1) {
    cal_smp[mis_idx, 5:7] <- matrix(rep(c(NA, NA, NA), times = length(mis_idx)), nrow = length(mis_idx), ncol = 3, byrow = TRUE)
  }

  fil_smp <- cal_smp[complete.cases(cal_smp), ]

  cal_rat <- nrow(fil_smp) / nrow(cal_smp)
  # cal_rat
  err_rat <- sum(rowSums(cal_smp[, -1] == tru_smp[, -1]) == 1, na.rm = TRUE) / nrow(fil_smp)
  # err_rat

  return(list(tru_smp = tru_smp,
              raw_smp = raw_smp,
              cal_smp = cal_smp,
              fil_smp = fil_smp,
              cal_rat = cal_rat,
              err_rat = err_rat,
              gen_frq = gen_frq,
              hap_frq = hap_frq))
}
#' Compiled version
cmpsimulateHMM <- cmpfun(simulateHMM)
########################################

#' Run the bootstrap particle filter (BPF) with the two-locus Wright-Fisher diffusion with selection
#' Parameter setting
#' @param sel_cof the selection coefficients of the tobiano, sabino and mixed against solid
#' @param rec_rat the recombination rate between the KIT13 and KIT16 loci
#' @param pop_siz the size of the horse population (non-constant)
#' @param ref_siz the reference size of the horse population
#' @param evt_gen the generation that the event of interest occurred
#' @param raw_smp the sample of ancient horses (sampling times and individual genotypes)
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter

#' Standard version
runBPF <- function(sel_cof, rec_rat, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num) {
  # preprocess the raw sample
  raw_smp <- as.data.frame(raw_smp)
  # rownames(raw_smp) <- NULL
  # colnames(raw_smp) <- c("generation", "KM0KM0", "KM0KM1", "KM1KM1", "sb1sb1", "sb1SB1", "SB1SB1")
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
  colnames(raw_smp) <- c("generation", "KM0KM0/sb1sb1", "KM0KM0/sb1SB1", "KM0KM0/SB1SB1", "KM0KM1/sb1sb1", "KM0KM1/sb1SB1", "KM1KM1/sb1sb1", "KM0KM1/SB1SB1", "KM1KM1/sb1SB1", "KM1KM1/SB1SB1")
  raw_smp <- raw_smp[order(raw_smp$generation), ]
  evt_gen <- evt_gen - min(raw_smp$generation)
  raw_smp$generation <- raw_smp$generation - min(raw_smp$generation)

  rownames(raw_smp) <- NULL
  colnames(raw_smp) <- NULL
  raw_smp <- t(as.matrix(raw_smp))

  # run the BPF
  BPF <- runBPF_arma(sel_cof, rec_rat, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num)

  return(list(lik = BPF$lik,
              wght = BPF$wght,
              hap_frq_pth = BPF$hap_frq_pth,
              hap_frq_pre_resmp = BPF$hap_frq_pre_resmp,
              hap_frq_pst_resmp = BPF$hap_frq_pst_resmp,
              gen_frq_pre_resmp = BPF$gen_frq_pre_resmp,
              gen_frq_pst_resmp = BPF$gen_frq_pst_resmp))
}
#' Compiled version
cmprunBPF <- cmpfun(runBPF)

########################################

#' Calculate the optimal particle number in the particle marginal Metropolis-Hastings (PMMH)
#' Parameter settings
#' @param sel_cof the selection coefficients of the tobiano, sabino and mixed against solid
#' @param rec_rat the recombination rate between the KIT13 and KIT16 loci
#' @param pop_siz the size of the horse population (non-constant)
#' @param ref_siz the reference size of the horse population
#' @param evt_gen the generation that the event of interest occurred
#' @param raw_smp the sample of ancient horses (sampling times and individual genotypes)
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter
#' @param gap_num the number of particles increased or decreased in the optimal particle number search

#' Standard version
calculateOptimalParticleNum <- function(sel_cof, rec_rat, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, gap_num) {
  # preprocess the raw sample
  raw_smp <- as.data.frame(raw_smp)
  # rownames(raw_smp) <- NULL
  # colnames(raw_smp) <- c("generation", "KM0KM0", "KM0KM1", "KM1KM1", "sb1sb1", "sb1SB1", "SB1SB1")
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
  colnames(raw_smp) <- c("generation", "KM0KM0/sb1sb1", "KM0KM0/sb1SB1", "KM0KM0/SB1SB1", "KM0KM1/sb1sb1", "KM0KM1/sb1SB1", "KM1KM1/sb1sb1", "KM0KM1/SB1SB1", "KM1KM1/sb1SB1", "KM1KM1/SB1SB1")
  raw_smp <- raw_smp[order(raw_smp$generation), ]
  evt_gen <- evt_gen - min(raw_smp$generation)
  raw_smp$generation <- raw_smp$generation - min(raw_smp$generation)

  rownames(raw_smp) <- NULL
  colnames(raw_smp) <- NULL
  raw_smp <- t(as.matrix(raw_smp))

  # calculate the optimal particle number
  OptNum <- calculateOptimalParticleNum_arma(sel_cof, rec_rat, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, gap_num)

  return(list(opt_pcl_num = as.vector(OptNum$opt_pcl_num),
              log_lik_sdv = as.vector(OptNum$log_lik_sdv)))
}
#' Compiled version
cmpcalculateOptimalParticleNum <- cmpfun(calculateOptimalParticleNum)

########################################

#' Run the particle marginal Metropolis-Hastings (PMMH)
#' Parameter settings
#' @param sel_cof the selection coefficients of the tobiano, sabino and mixed against solid
#' @param rec_rat the recombination rate between the KIT13 and KIT16 loci
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
  raw_smp <- as.data.frame(raw_smp)
  # rownames(raw_smp) <- NULL
  # colnames(raw_smp) <- c("generation", "KM0KM0", "KM0KM1", "KM1KM1", "sb1sb1", "sb1SB1", "SB1SB1")
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
  colnames(raw_smp) <- c("generation", "KM0KM0/sb1sb1", "KM0KM0/sb1SB1", "KM0KM0/SB1SB1", "KM0KM1/sb1sb1", "KM0KM1/sb1SB1", "KM1KM1/sb1sb1", "KM0KM1/SB1SB1", "KM1KM1/sb1SB1", "KM1KM1/SB1SB1")
  raw_smp <- raw_smp[order(raw_smp$generation), ]
  evt_gen <- evt_gen - min(raw_smp$generation)
  raw_smp$generation <- raw_smp$generation - min(raw_smp$generation)

  rownames(raw_smp) <- NULL
  colnames(raw_smp) <- NULL
  raw_smp <- t(as.matrix(raw_smp))

  # run the PMMH
  PMMH <- runPMMH_arma(sel_cof, rec_rat, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num)

  return(list(sel_cof_chn = as.array(PMMH$sel_cof_chn),
              frq_pth_chn = as.array(PMMH$frq_pth_chn)))
}
#' Compiled version
cmprunPMMH <- cmpfun(runPMMH)

########################################

#' Run the adaptive particle marginal Metropolis-Hastings (AdaptPMMH)
#' Parameter settings
#' @param sel_cof the selection coefficients of the tobiano, sabino and mixed against solid
#' @param rec_rat the recombination rate between the KIT13 and KIT16 loci
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
  raw_smp <- as.data.frame(raw_smp)
  # rownames(raw_smp) <- NULL
  # colnames(raw_smp) <- c("generation", "KM0KM0", "KM0KM1", "KM1KM1", "sb1sb1", "sb1SB1", "SB1SB1")
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
  colnames(raw_smp) <- c("generation", "KM0KM0/sb1sb1", "KM0KM0/sb1SB1", "KM0KM0/SB1SB1", "KM0KM1/sb1sb1", "KM0KM1/sb1SB1", "KM1KM1/sb1sb1", "KM0KM1/SB1SB1", "KM1KM1/sb1SB1", "KM1KM1/SB1SB1")
  raw_smp <- raw_smp[order(raw_smp$generation), ]
  evt_gen <- evt_gen - min(raw_smp$generation)
  raw_smp$generation <- raw_smp$generation - min(raw_smp$generation)

  rownames(raw_smp) <- NULL
  colnames(raw_smp) <- NULL
  raw_smp <- t(as.matrix(raw_smp))

  # run the adaptive PMMH
  PMMH <- runAdaptPMMH_arma(sel_cof, rec_rat, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, stp_siz, apt_rto)

  return(list(sel_cof_chn = as.array(PMMH$sel_cof_chn),
              frq_pth_chn = as.array(PMMH$frq_pth_chn)))
}
#' Compiled version
cmprunAdaptPMMH <- cmpfun(runAdaptPMMH)

########################################

#' Run the Bayesian procedure for the inference of selection
#' Parameter settings
#' @param sel_cof the selection coefficients of the tobiano, sabino and mixed against solid
#' @param rec_rat the recombination rate between the KIT13 and KIT16 loci
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
  raw_smp <- as.data.frame(raw_smp)
  # rownames(raw_smp) <- NULL
  # colnames(raw_smp) <- c("generation", "KM0KM0", "KM0KM1", "KM1KM1", "sb1sb1", "sb1SB1", "SB1SB1")
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
  colnames(raw_smp) <- c("generation", "KM0KM0/sb1sb1", "KM0KM0/sb1SB1", "KM0KM0/SB1SB1", "KM0KM1/sb1sb1", "KM0KM1/sb1SB1", "KM1KM1/sb1sb1", "KM0KM1/SB1SB1", "KM1KM1/sb1SB1", "KM1KM1/SB1SB1")
  raw_smp <- raw_smp[order(raw_smp$generation), ]
  evt_gen <- evt_gen - min(raw_smp$generation)
  raw_smp$generation <- raw_smp$generation - min(raw_smp$generation)

  rownames(raw_smp) <- NULL
  colnames(raw_smp) <- NULL
  raw_smp <- t(as.matrix(raw_smp))

  if (adp_set == TRUE) {
    # run the adaptive PMMH
    PMMH <- runAdaptPMMH_arma(sel_cof, rec_rat, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, stp_siz, apt_rto)
  } else {
    # run the PMMH
    PMMH <- runPMMH_arma(sel_cof, rec_rat, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num)
  }
  sel_cof_chn <- as.array(PMMH$sel_cof_chn)
  frq_pth_chn <- as.array(PMMH$frq_pth_chn)
  cal_smp_chn <- as.array(PMMH$cal_smp_chn)

  # burn-in and thinning
  sel_cof_chn <- sel_cof_chn[, , brn_num:dim(sel_cof_chn)[3]]
  sel_cof_chn <- sel_cof_chn[, , (1:round(dim(sel_cof_chn)[3] / thn_num)) * thn_num]

  dif_sel_chn <- sel_cof_chn[, 2, ] - sel_cof_chn[, 1, ]

  frq_pth_chn <- frq_pth_chn[, , brn_num:dim(frq_pth_chn)[3]]
  frq_pth_chn <- frq_pth_chn[, , (1:round(dim(frq_pth_chn)[3] / thn_num)) * thn_num]

  cal_smp_chn <- cal_smp_chn[, , brn_num:dim(cal_smp_chn)[3]]
  cal_smp_chn <- cal_smp_chn[, , (1:round(dim(cal_smp_chn)[3] / thn_num)) * thn_num]

  # MMSE estimates for selection coefficients with their 95% HPD intervals
  sel_cof_est <- matrix(NA, nrow = 3, ncol = 2)
  sel_cof_est[, 1] <- rowMeans(sel_cof_chn[, 1, ])
  sel_cof_est[, 2] <- rowMeans(sel_cof_chn[, 2, ])
  sel_cof_est[, 3] <- rowMeans(sel_cof_chn[, 3, ])

  sel_cof_hpd <- array(NA, dim = c(3, 2, 2))
  sel_cof_hpd[1, , 1] <- HPDinterval(as.mcmc(sel_cof_chn[1, 1, ]), prob = 0.95)
  sel_cof_hpd[2, , 1] <- HPDinterval(as.mcmc(sel_cof_chn[2, 1, ]), prob = 0.95)
  sel_cof_hpd[3, , 1] <- HPDinterval(as.mcmc(sel_cof_chn[3, 1, ]), prob = 0.95)
  sel_cof_hpd[1, , 2] <- HPDinterval(as.mcmc(sel_cof_chn[1, 2, ]), prob = 0.95)
  sel_cof_hpd[2, , 2] <- HPDinterval(as.mcmc(sel_cof_chn[2, 2, ]), prob = 0.95)
  sel_cof_hpd[3, , 2] <- HPDinterval(as.mcmc(sel_cof_chn[3, 2, ]), prob = 0.95)

  # MMSE estimates for selection changes with their 95% HPD intervals
  dif_sel_est <- rowMeans(dif_sel_chn)

  dif_sel_hpd <- matrix(NA, nrow = 3, ncol = 2)
  dif_sel_hpd[1, ] <- HPDinterval(as.mcmc(dif_sel_chn[1, ]), prob = 0.95)
  dif_sel_hpd[2, ] <- HPDinterval(as.mcmc(dif_sel_chn[2, ]), prob = 0.95)
  dif_sel_hpd[3, ] <- HPDinterval(as.mcmc(dif_sel_chn[3, ]), prob = 0.95)

  # MMSE estimates for population haplotype frequency trajectories with their 95% HPD intervals
  frq_pth_est <- matrix(NA, nrow = 4, ncol = dim(frq_pth_chn)[2])
  frq_pth_est[1, ] <- rowMeans(frq_pth_chn[1, , ])
  frq_pth_est[2, ] <- rowMeans(frq_pth_chn[2, , ])
  frq_pth_est[3, ] <- rowMeans(frq_pth_chn[3, , ])
  frq_pth_est[4, ] <- rowMeans(frq_pth_chn[4, , ])

  frq_pth_hpd <- array(NA, dim = c(4, 2, dim(frq_pth_chn)[2]))
  for (i in 1:dim(frq_pth_chn)[2]) {
    frq_pth_hpd[1, , i] <- HPDinterval(as.mcmc(frq_pth_chn[1, i, ]), prob = 0.95)
    frq_pth_hpd[2, , i] <- HPDinterval(as.mcmc(frq_pth_chn[2, i, ]), prob = 0.95)
    frq_pth_hpd[3, , i] <- HPDinterval(as.mcmc(frq_pth_chn[3, i, ]), prob = 0.95)
    frq_pth_hpd[4, , i] <- HPDinterval(as.mcmc(frq_pth_chn[4, i, ]), prob = 0.95)
  }

  return(list(sel_cof_est = sel_cof_est,
              sel_cof_hpd = sel_cof_hpd,
              sel_cof_chn = sel_cof_chn,
              dif_sel_est = dif_sel_est,
              dif_sel_hpd = dif_sel_hpd,
              dif_sel_chn = dif_sel_chn,
              frq_pth_est = frq_pth_est,
              frq_pth_hpd = frq_pth_hpd,
              frq_pth_chn = frq_pth_chn))
}
#' Compiled version
cmprunBayesianProcedure <- cmpfun(runBayesianProcedure)

########################################

#' Convert the haplotype frequency trajectories to the genotype/phenotype frequency trajectories
#' Parameter settings
#' @param hap_pth the haplotype frequency trajectories
#' @param sel_cof the selection coefficients of the tobiano, sabino and mixed against solid
#' @param evt_gen the generation that the event of interest occurred
#' @param int_gen the generation that the haplotype frequency trajectories started
#' @param lst_gen the generation that the haplotype frequency trajectories ended

#' Standard version
convertHaploFreq <- function(hap_pth, sel_cof, evt_gen, int_gen, lst_gen) {
  frq_pth <- convertHaploFreq_arma(hap_pth, sel_cof, evt_gen, int_gen, lst_gen)

  return(list(gen_pth = as.array(frq_pth$gen_pth),
              phe_pth = as.array(frq_pth$phe_pth)))
}
#' Compiled version
cmpconvertHaploFreq <- cmpfun(convertHaploFreq)

################################################################################
