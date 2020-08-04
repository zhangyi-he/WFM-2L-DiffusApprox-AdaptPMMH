// Inferring natural selection acting on horse coat colours and patterns during the process of domestication
// Zhangyi He, Xiaoyang Dai, Mark Beaumont and Feng Yu

// version 1.0

// C functions

#define ARMA_64BIT_WORD 1

#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::plugins(cpp11)]]
#include <math.h>

using namespace Rcpp;
using namespace std;
// using namespace arma;

/********** WFM **********/
// Calculate the fitness matrix for the Wright-Fisher model
// [[Rcpp::export]]
arma::dmat calculateFitnessMat_arma(const double& s_b, const double& s_c) {
  // ensure RNG gets set/reset
  RNGScope scope;

  // declare the fitness
  arma::dmat fts_mat = arma::ones<arma::dmat>(4, 4);
  fts_mat(1, 1) = 1 + s_c;
  fts_mat(3, 1) = 1 + s_c;
  fts_mat(2, 2) = 1 + s_b;
  fts_mat(3, 2) = 1 + s_b;
  fts_mat(1, 3) = 1 + s_c;
  fts_mat(2, 3) = 1 + s_b;
  fts_mat(3, 3) = 1 + s_c;

  return fts_mat;
}

// Simulate the haplotype frequency trajectories according to the Wright-Fisher model
// [[Rcpp::export]]
arma::dmat simulateTLWFMS_arma(const arma::dmat& fts_mat, const double& r, const int& N, const arma::dcolvec& int_frq, const int& int_gen, const int& lst_gen) {
  // ensure RNG gets set/reset
  RNGScope scope;

  // declare the haplotype frequency trajectories
  arma::dmat frq_pth = arma::zeros<arma::dmat>(4, arma::uword(lst_gen - int_gen) + 1);

  // initialise the haplotype frequencies in generation 0
  frq_pth.col(0) = int_frq;

  // declare eta
  arma::dcolvec eta = arma::ones<arma::dcolvec>(4);
  eta(0) = -1;
  // eta(1) = 1;
  // eta(2) = 1;
  eta(3) = -1;

  // declare and initialise the sampling probabilities
  arma::dcolvec prob = int_frq;

  for(arma::uword t = 1; t < arma::uword(lst_gen - int_gen) + 1; t++) {
    // calculate the haplotype frequencies after natural selection
    prob = prob % (fts_mat * prob) / arma::as_scalar(prob.t() * fts_mat * prob);

    // calculate the haplotype frequencies after genetic recombination
    prob = prob + eta * r * (prob(0) * prob(3) - prob(1) * prob(2));

    // do the Wright-Fisher sampling (reproduction)
    IntegerVector hap_cnt(4);
    R::rmultinom(2 * N, prob.begin(), 4, hap_cnt.begin());
    prob = as<arma::dcolvec>(hap_cnt) / 2 / N;

    frq_pth.col(t) = prob;
  }

  return frq_pth;
}
/*************************/

/********** WFD **********/
// Simulate the haplotype frequency trajectories according to the Wright-Fisher diffusion using the Euler-Maruyama method
// [[Rcpp::export]]
arma::dmat simulateTLWFDS_arma(const double& alpha_b, const double& alpha_c, const double& rho, const int& N, const arma::dcolvec& int_frq, const int& int_gen, const int& lst_gen, const arma::uword& ptn_num) {
  // ensure RNG gets set/reset
  RNGScope scope;

  // declare delta t
  double dt = 1.0 / (2 * N) / ptn_num;
  // declare delta W
  arma::dmat dW = pow(dt, 0.5) * arma::randn<arma::dmat>(6, arma::uword(lst_gen - int_gen) * ptn_num);

  // declare the haplotype frequency trajectories
  arma::dmat frq_pth = arma::zeros<arma::dmat>(4, arma::uword(lst_gen - int_gen) * ptn_num + 1);

  // initialise the haplotype frequencies in generation 0
  frq_pth.col(0) = int_frq;

  for (arma::uword t = 1; t < arma::uword(lst_gen - int_gen) * ptn_num + 1; t++) {
    // calculate the drift coefficient vector
    arma::dcolvec mu = arma::zeros<arma::dcolvec>(4);
    double LD = frq_pth(0, t - 1) * frq_pth(3, t - 1) - frq_pth(1, t - 1) * frq_pth(2, t - 1);
    mu(0) = -alpha_b * frq_pth(2, t - 1) * (frq_pth(0, t - 1) * frq_pth(3, t - 1) + frq_pth(0, t - 1) * (frq_pth(2, t - 1) + frq_pth(3, t - 1))) -
      alpha_c * frq_pth(0, t - 1) * (frq_pth(1, t - 1) + frq_pth(3, t - 1)) * (frq_pth(1, t - 1) + frq_pth(3, t - 1)) -
      0.5 * rho * LD;
    mu(1) = -alpha_b * frq_pth(2, t - 1) * (frq_pth(1, t - 1) * frq_pth(3, t - 1) + frq_pth(1, t - 1) * (frq_pth(2, t - 1) + frq_pth(3, t - 1))) +
      alpha_c * frq_pth(1, t - 1) * (frq_pth(0, t - 1) + frq_pth(2, t - 1)) * (frq_pth(1, t - 1) + frq_pth(3, t - 1)) +
      0.5 * rho * LD;
    mu(2) = -alpha_b * frq_pth(2, t - 1) * (frq_pth(2, t - 1) * frq_pth(3, t - 1) + frq_pth(2, t - 1) * (frq_pth(2, t - 1) + frq_pth(3, t - 1)) - (frq_pth(2, t - 1) + frq_pth(3, t - 1))) -
      alpha_c * frq_pth(2, t - 1) * (frq_pth(1, t - 1) + frq_pth(3, t - 1)) * (frq_pth(1, t - 1) + frq_pth(3, t - 1)) +
      0.5 * rho * LD;
    mu(3) = -alpha_b * frq_pth(2, t - 1) * (frq_pth(3, t - 1) * frq_pth(3, t - 1) + frq_pth(3, t - 1) * (frq_pth(2, t - 1) + frq_pth(3, t - 1)) - frq_pth(3, t - 1)) +
      alpha_c * frq_pth(3, t - 1) * (frq_pth(0, t - 1) + frq_pth(2, t - 1)) * (frq_pth(1, t - 1) + frq_pth(3, t - 1)) -
      0.5 * rho * LD;

    // calculate the diffusion coefficient matrix
    arma::dmat sigma = arma::zeros<arma::dmat>(4, 6);
    sigma(0, 0) = pow(frq_pth(0, t - 1) * frq_pth(1, t - 1), 0.5);
    sigma(0, 1) = pow(frq_pth(0, t - 1) * frq_pth(2, t - 1), 0.5);
    sigma(0, 2) = pow(frq_pth(0, t - 1) * frq_pth(3, t - 1), 0.5);
    // sigma(0, 3) = 0;
    // sigma(0, 4) = 0;
    // sigma(0, 5) = 0;
    sigma(1, 0) = -pow(frq_pth(1, t - 1) * frq_pth(0, t - 1), 0.5);
    // sigma(1, 1) = 0;
    // sigma(1, 2) = 0;
    sigma(1, 3) = pow(frq_pth(1, t - 1) * frq_pth(2, t - 1), 0.5);
    sigma(1, 4) = pow(frq_pth(1, t - 1) * frq_pth(3, t - 1), 0.5);
    // sigma(1, 5) = 0;
    // sigma(2, 0) = 0;
    sigma(2, 1) = -pow(frq_pth(2, t - 1) * frq_pth(0, t - 1), 0.5);
    // sigma(2, 2) = 0;
    sigma(2, 3) = -pow(frq_pth(2, t - 1) * frq_pth(1, t - 1), 0.5);
    // sigma(2, 4) = 0;
    sigma(2, 5) = pow(frq_pth(2, t - 1) * frq_pth(3, t - 1), 0.5);
    // sigma(3, 0) = 0;
    // sigma(3, 1) = 0;
    sigma(3, 2) = -pow(frq_pth(3, t - 1) * frq_pth(0, t - 1), 0.5);
    // sigma(3, 3) = 0;
    sigma(3, 4) = -pow(frq_pth(3, t - 1) * frq_pth(1, t - 1), 0.5);
    sigma(3, 5) = -pow(frq_pth(3, t - 1) * frq_pth(2, t - 1), 0.5);

    // proceed the Euler-Maruyama scheme
    frq_pth.col(t) = frq_pth.col(t - 1) + mu * dt + sigma * dW.col(t - 1);

    // remove the noise from the numerical techniques
    // for(arma::uword i = 0; i < 4; i++) {
    //   if(frq_pth(i, t) < 0) {
    //     frq_pth(i, t) = 0;
    //   }
    //   if(frq_pth(i, t) > 1) {
    //     frq_pth(i, t) = 1;
    //   }
    // }
    // frq_pth.col(t) = frq_pth.col(t) / sum(frq_pth.col(t));
  }

  return frq_pth;
}
/*************************/

/********** BPF **********/
// Calculate the possible haplotype counts in the sample of genotypes
// [[Rcpp::export]]
arma::imat calculateHaploCnt_arma(const int& smp_siz, const arma::icolvec& smp_cnt) {
  // ensure RNG gets set/reset
  RNGScope scope;

  if (smp_cnt.n_elem == 2) {
    arma::imat smp_hap_cnt = arma::zeros<arma::imat>(4, 1);

    for (int i = 0; i <= min(smp_cnt(0), smp_cnt(1)); i++) {
      int j = smp_cnt(0) - i;
      int k = smp_cnt(1) - i;
      if (i + j + k <= smp_siz) {
        smp_hap_cnt(0, 0) = i;
        smp_hap_cnt(1, 0) = j;
        smp_hap_cnt(2, 0) = k;
        smp_hap_cnt(3, 0) = smp_siz - i - j - k;
        smp_hap_cnt.insert_cols(0, 1);
      }
    }
    smp_hap_cnt.shed_cols(0, 0);

    return smp_hap_cnt;
  } else {
    arma::imat smp_hap_cnt = arma::zeros<arma::imat>(4, 1);
    smp_hap_cnt.col(0) = smp_cnt;

    return smp_hap_cnt;
  }
}

// Calculate the multinomial probabilities
// [[Rcpp::export]]
double calculateMultinomProb_arma(const arma::icolvec& smp_cnt, const int& smp_siz, const arma::dcolvec& pop_frq) {
  // ensure RNG gets set/reset
  RNGScope scope;

  if (arma::any(pop_frq == 0)) {
    if (arma::any(smp_cnt.elem(arma::find(pop_frq == 0)) != 0)) {
      return 0;
    }

    arma::icolvec smp_cnt_nonzero = smp_cnt.elem(arma::find(pop_frq != 0));
    arma::dcolvec pop_frq_nonzero = pop_frq.elem(arma::find(pop_frq != 0));

    return exp(lgamma(smp_siz + 1) + sum(arma::conv_to<arma::dcolvec>::from(smp_cnt_nonzero) % log(pop_frq_nonzero) - lgamma(arma::conv_to<arma::dcolvec>::from(smp_cnt_nonzero) + 1)));
  } else {
    return exp(lgamma(smp_siz + 1) + sum(arma::conv_to<arma::dcolvec>::from(smp_cnt) % log(pop_frq) - lgamma(arma::conv_to<arma::dcolvec>::from(smp_cnt) + 1)));
  }
}

// Initialise the particles in the particle filter (uniform generation from the flat Dirichlet distribution)
// [[Rcpp::export]]
arma::dmat initialiseParticle(const arma::uword& pcl_num){
  // ensure RNG gets set/reset
  RNGScope scope;

  NumericMatrix part(pcl_num, 4);
  for(int j = 0; j < 4; j++){
    part(_, j) = rgamma(pcl_num, 1.0, 1.0);
  }
  for(int i = 0; i < pcl_num; i++){
    part(i, _) = part(i, _) / sum(part(i, _));
  }

  return as<arma::dmat>(transpose(part));
}

// Run the bootstrap particle filter
// [[Rcpp::export]]
List runBPF_arma(const double& s_b, const double& s_c, const double& r, const int& N, const arma::irowvec& smp_gen, const arma::irowvec& smp_siz, const arma::imat& smp_cnt, const arma::uword& ptn_num, const arma::uword& pcl_num) {
  // ensure RNG gets set/reset
  RNGScope scope;

  // rescale the selection coefficients and the recombination rate
  double alpha_b = 2 * N * s_b;
  double alpha_c = 2 * N * s_c;
  double rho = 4 * N * r;

  double lik = 1;

  arma::dmat wght = arma::zeros<arma::dmat>(pcl_num, smp_gen.n_elem);
  arma::dcube part_pre = arma::zeros<arma::dcube>(4, pcl_num, smp_gen.n_elem);
  arma::dcube part_pst = arma::zeros<arma::dcube>(4, pcl_num, smp_gen.n_elem);

  arma::dcolvec wght_tmp = arma::zeros<arma::dcolvec>(pcl_num);
  arma::dmat part_tmp = arma::zeros<arma::dmat>(4, pcl_num);

  // initialise the particles
  cout << "generation: " << smp_gen(0) << endl;
  part_tmp = initialiseParticle(pcl_num);
  arma::imat smp_hap_cnt = calculateHaploCnt_arma(smp_siz(0), smp_cnt.col(0));
  for (arma::uword i = 0; i < pcl_num; i++) {
    for (arma::uword j = 0; j < smp_hap_cnt.n_cols; j++) {
      wght_tmp(i) = wght_tmp(i) + calculateMultinomProb_arma(smp_hap_cnt.col(j), smp_siz(0), part_tmp.col(i));
    }
  }

  if (arma::sum(wght_tmp) > 0) {
    arma::dcolvec prob = arma::normalise(wght_tmp, 1);
    arma::ucolvec elem = arma::linspace<arma::ucolvec>(0, pcl_num - 1, pcl_num);
    arma::ucolvec indx = RcppArmadillo::sample(elem, pcl_num, TRUE, prob);

    lik = lik * arma::mean(wght_tmp);
    wght.col(0) = wght_tmp;
    part_pre.slice(0) = part_tmp;
    part_pst.slice(0) = part_tmp.cols(indx);
  } else {
    lik = 0;
    wght.shed_cols(0, smp_gen.n_elem - 1);
    part_pre.shed_slices(0, smp_gen.n_elem - 1);
    part_pst.shed_slices(0, smp_gen.n_elem - 1);

    return List::create(Named("lik", lik),
                        Named("wght", wght),
                        Named("part_pre_resmp", part_pre),
                        Named("part_pst_resmp", part_pst));
  }

  // run the bootstrap particle filter
  for (arma::uword k = 1; k < smp_gen.n_elem; k++) {
    cout << "generation: " << smp_gen(k) << endl;
    wght_tmp = arma::zeros<arma::dcolvec>(pcl_num);
    part_tmp = part_pst.slice(k - 1);
    arma::imat smp_hap_cnt = calculateHaploCnt_arma(smp_siz(k), smp_cnt.col(k));
    for (arma::uword i = 0; i < pcl_num; i++) {
      arma::dmat path = simulateTLWFDS_arma(alpha_b, alpha_c, rho, N, part_tmp.col(i), smp_gen(k - 1), smp_gen(k), ptn_num);
      part_tmp.col(i) = arma::vectorise(path.tail_cols(1), 0);
      for (arma::uword j = 0; j < smp_hap_cnt.n_cols; j++) {
        wght_tmp(i) = wght_tmp(i) + calculateMultinomProb_arma(smp_hap_cnt.col(j), smp_siz(k), part_tmp.col(i));
      }
    }

    if (arma::sum(wght_tmp) > 0) {
      arma::dcolvec prob = arma::normalise(wght_tmp, 1);
      arma::ucolvec elem = arma::linspace<arma::ucolvec>(0, pcl_num - 1, pcl_num);
      arma::ucolvec indx = RcppArmadillo::sample(elem, pcl_num, TRUE, prob);

      lik = lik * arma::mean(wght_tmp);
      wght.col(k) = wght_tmp;
      part_pre.slice(k) = part_tmp;
      part_pst.slice(k) = part_tmp.cols(indx);
    } else {
      lik = 0;
      wght.shed_cols(k, smp_gen.n_elem - 1);
      part_pre.shed_slices(k, smp_gen.n_elem - 1);
      part_pst.shed_slices(k, smp_gen.n_elem - 1);
      break;
    }
  }

  return List::create(Named("lik", lik),
                      Named("wght", wght),
                      Named("part_pre_resmp", part_pre),
                      Named("part_pst_resmp", part_pst));
}
/*************************/

/********** PMMH **********/
// Calculate the log-likelihood using the bootstrap particle filter
// [[Rcpp::export]]
double calculateLogLikelihood_arma(const double& s_b, const double& s_c, const double& r, const int& N, const arma::irowvec& smp_gen, const arma::irowvec& smp_siz, const arma::field<arma::imat>& ptl_hap_cnt, const arma::uword& ptn_num, const arma::uword& pcl_num) {
  // ensure RNG gets set/reset
  RNGScope scope;

  // rescale the selection coefficients and the recombination rate
  double alpha_b = 2 * N * s_b;
  double alpha_c = 2 * N * s_c;
  double rho = 4 * N * r;

  double log_lik = 0;

  arma::dcolvec wght = arma::zeros<arma::dcolvec>(pcl_num);
  arma::dmat part_pre = arma::zeros<arma::dmat>(4, pcl_num);
  arma::dmat part_pst = arma::zeros<arma::dmat>(4, pcl_num);

  // initialise the particles
  part_pre = initialiseParticle(pcl_num);
  arma::imat smp_hap_cnt = ptl_hap_cnt(0);
  for (arma::uword i = 0; i < pcl_num; i++) {
    for (arma::uword j = 0; j < smp_hap_cnt.n_cols; j++) {
      wght(i) = wght(i) + calculateMultinomProb_arma(smp_hap_cnt.col(j), smp_siz(0), part_pre.col(i));
    }
  }

  if (arma::mean(wght) > 0) {
    log_lik = log_lik + log(arma::mean(wght));
    arma::dcolvec prob = arma::normalise(wght, 1);
    arma::ucolvec elem = arma::linspace<arma::ucolvec>(0, pcl_num - 1, pcl_num);
    arma::ucolvec indx = RcppArmadillo::sample(elem, pcl_num, TRUE, prob);
    part_pst = part_pre.cols(indx);
  } else {
    log_lik = -(arma::datum::inf);
    return log_lik;
  }

  // run the bootstrap particle filter
  for (arma::uword k = 1; k < smp_gen.n_elem; k++) {
    wght = arma::zeros<arma::dcolvec>(pcl_num);
    arma::imat smp_hap_cnt = ptl_hap_cnt(k);
    for (arma::uword i = 0; i < pcl_num; i++) {
      arma::dmat path = simulateTLWFDS_arma(alpha_b, alpha_c, rho, N, part_pst.col(i), smp_gen(k - 1), smp_gen(k), ptn_num);
      part_pre.col(i) = arma::vectorise(path.tail_cols(1), 0);
      for (arma::uword j = 0; j < smp_hap_cnt.n_cols; j++) {
        wght(i) = wght(i) + calculateMultinomProb_arma(smp_hap_cnt.col(j), smp_siz(k), part_pre.col(i));
      }
    }

    if (arma::mean(wght) > 0) {
      log_lik = log_lik + log(arma::mean(wght));
      arma::dcolvec prob = arma::normalise(wght, 1);
      arma::ucolvec elem = arma::linspace<arma::ucolvec>(0, pcl_num - 1, pcl_num);
      arma::ucolvec indx = RcppArmadillo::sample(elem, pcl_num, TRUE, prob);
      part_pst = part_pre.cols(indx);
    } else {
      log_lik = -(arma::datum::inf);
      return log_lik;
    }
  }

  return log_lik;
}

// Calculate the optimal particle number in the particle marginal Metropolis-Hastings
// [[Rcpp::export]]
List calculateOptimalParticleNum_arma(const double& s_b, const double& s_c, const double& r, const int& N, const arma::irowvec& smp_gen, const arma::irowvec& smp_siz, const arma::imat& smp_cnt, const arma::uword& ptn_num, const arma::uword& pcl_num, const arma::uword& gap_num) {
  // ensure RNG gets set/reset
  RNGScope scope;

  arma::field<arma::imat> ptl_hap_cnt(smp_gen.n_elem);
  for (arma::uword k = 0; k < smp_gen.n_elem; k++) {
    ptl_hap_cnt(k) = calculateHaploCnt_arma(smp_siz(k), smp_cnt.col(k));
  }

  arma::drowvec log_lik(300);
  for (arma::uword i = 0; i < 300; i++) {
    log_lik(i) = calculateLogLikelihood_arma(s_b, s_c, r, N, smp_gen, smp_siz, ptl_hap_cnt, ptn_num, pcl_num);
  }

  arma::drowvec log_lik_sdv(1);
  log_lik_sdv(0) = arma::stddev(log_lik);
  log_lik_sdv.print();
  arma::urowvec opt_pcl_num(1);
  opt_pcl_num(0) = pcl_num;

  if (log_lik_sdv(0) > 1.7) {
    while (log_lik_sdv(0) > 1.0) {
      opt_pcl_num.insert_cols(0, 1);
      opt_pcl_num(0) = opt_pcl_num(1) + gap_num;
      for (arma::uword i = 0; i < 300; i++) {
        log_lik(i) = calculateLogLikelihood_arma(s_b, s_c, r, N, smp_gen, smp_siz, ptl_hap_cnt, ptn_num, opt_pcl_num(0));
      }
      log_lik_sdv.insert_cols(0, 1);
      log_lik_sdv(0) = arma::stddev(log_lik);
      log_lik_sdv.print();
    }
    opt_pcl_num = arma::reverse(opt_pcl_num);
    log_lik_sdv = arma::reverse(log_lik_sdv);
  } else if (log_lik_sdv(0) < 1.0) {
    while (log_lik_sdv(0) < 1.7 && opt_pcl_num(0) > gap_num) {
      opt_pcl_num.insert_cols(0, 1);
      opt_pcl_num(0) = opt_pcl_num(1) - gap_num;
      for (arma::uword i = 0; i < 300; i++) {
        log_lik(i) = calculateLogLikelihood_arma(s_b, s_c, r, N, smp_gen, smp_siz, ptl_hap_cnt, ptn_num, opt_pcl_num(0));
      }
      log_lik_sdv.insert_cols(0, 1);
      log_lik_sdv(0) = arma::stddev(log_lik);
      log_lik_sdv.print();
    }
  } else {
    while (log_lik_sdv(0) > 1.0) {
      opt_pcl_num.insert_cols(0, 1);
      opt_pcl_num(0) = opt_pcl_num(1) + gap_num;
      for (arma::uword i = 0; i < 300; i++) {
        log_lik(i) = calculateLogLikelihood_arma(s_b, s_c, r, N, smp_gen, smp_siz, ptl_hap_cnt, ptn_num, opt_pcl_num(0));
      }
      log_lik_sdv.insert_cols(0, 1);
      log_lik_sdv(0) = arma::stddev(log_lik);
      log_lik_sdv.print();
    }
    opt_pcl_num = arma::reverse(opt_pcl_num);
    log_lik_sdv = arma::reverse(log_lik_sdv);

    while (log_lik_sdv(0) < 1.7 && opt_pcl_num(0) > gap_num) {
      opt_pcl_num.insert_cols(0, 1);
      opt_pcl_num(0) = opt_pcl_num(1) - gap_num;
      for (arma::uword i = 0; i < 300; i++) {
        log_lik(i) = calculateLogLikelihood_arma(s_b, s_c, r, N, smp_gen, smp_siz, ptl_hap_cnt, ptn_num, opt_pcl_num(0));
      }
      log_lik_sdv.insert_cols(0, 1);
      log_lik_sdv(0) = arma::stddev(log_lik);
      log_lik_sdv.print();
    }
  }

  return List::create(Named("opt_pcl_num", opt_pcl_num),
                      Named("log_lik_sdv", log_lik_sdv));
}

// Run the particle marginal Metropolis-Hastings
//[[Rcpp::export]]
List runPMMH_arma(const double& s_b, const double& s_c, const double& r, const int& N, const arma::irowvec& smp_gen, const arma::irowvec& smp_siz, const arma::imat& smp_cnt, const arma::uword& ptn_num, const arma::uword& pcl_num, const arma::uword& itn_num) {
  // ensure RNG gets set/reset
  RNGScope scope;

  arma::field<arma::imat> ptl_hap_cnt(smp_gen.n_elem);
  for (arma::uword k = 0; k < smp_gen.n_elem; k++) {
    ptl_hap_cnt(k) = calculateHaploCnt_arma(smp_siz(k), smp_cnt.col(k));
  }

  arma::drowvec s_b_chn = arma::zeros<arma::drowvec>(itn_num);
  arma::drowvec s_c_chn = arma::zeros<arma::drowvec>(itn_num);

  //arma::drowvec log_pri_chn = arma::zeros<arma::drowvec>(itn_num);
  arma::drowvec log_lik_chn = arma::zeros<arma::drowvec>(itn_num);

  double sel_cof_sd = 1e-02;

  // initialise the population genetic parameters
  cout << "iteration: " << 1 << endl;
  // take the uniform prior and fix the initial value of the selection coefficient to zero
  // or take the beta prior with alpha = 1 and beta = 3
  s_b_chn(0) = s_b;
  s_c_chn(0) = s_c;

  log_lik_chn(0) = calculateLogLikelihood_arma(s_b_chn(0), s_c_chn(0), r, N, smp_gen, smp_siz, ptl_hap_cnt, ptn_num, pcl_num);

  double apt_rto = 0;
  for (arma::uword i = 1; i < itn_num; i++) {
    cout << "iteration: " << i + 1 << endl;

    // draw the candidates of the selection coefficients from the random walk proposal
    s_b_chn(i) = s_b_chn(i - 1) + sel_cof_sd * arma::randn();
    s_c_chn(i) = s_c_chn(i - 1) + sel_cof_sd * arma::randn();

    if (s_b_chn(i) > 1 || s_c_chn(i) > 1) {
      s_b_chn(i) = s_b_chn(i - 1);
      s_c_chn(i) = s_c_chn(i - 1);
      log_lik_chn(i) = log_lik_chn(i - 1);
    } else {
      // calculate the proposal
      //double log_psl_old_new = log(arma::normpdf(s_b_chn(i - 1), s_b_chn(i), sel_cof_sd)) + log(arma::normpdf(s_c_chn(i - 1), s_c_chn(i), sel_cof_sd));
      //double log_psl_new_old = log(arma::normpdf(s_b_chn(i), s_b_chn(i - 1), sel_cof_sd)) + log(arma::normpdf(s_c_chn(i), s_c_chn(i - 1), sel_cof_sd));

      // calculate the likelihood
      log_lik_chn(i) = calculateLogLikelihood_arma(s_b_chn(i), s_c_chn(i), r, N, smp_gen, smp_siz, ptl_hap_cnt, ptn_num, pcl_num);

      // calculate the acceptance ratio
      apt_rto = exp(log_lik_chn(i) - log_lik_chn(i - 1));
      //apt_rto = exp((log_pri_chn(i) + log_lik_chn(i) + log_psl_old_new) - (log_pri_chn(i - 1) + log_lik_chn(i - 1) + log_psl_new_old));

      if (arma::randu() > apt_rto) {
        s_b_chn(i) = s_b_chn(i - 1);
        s_c_chn(i) = s_c_chn(i - 1);
        log_lik_chn(i) = log_lik_chn(i - 1);
      }
    }
  }

  return List::create(Named("s_b_chn", s_b_chn),
                      Named("s_c_chn", s_c_chn));
}
/*************************/
