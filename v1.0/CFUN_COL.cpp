// Inferring natural selection acting on horse coat colours and patterns during the process of domestication from ancient DNA data
// Xiaoyang Dai, Sile Hu, Mark Beaumont, Feng Yu, Zhangyi He

// version 1.0
// Horse coat colours (ASIP & MC1R) under constant natural selection and constant demographic histories (N/A is not allowed)

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
// Calculate the fitness matrix for the two-locus Wright-Fisher model with selection
// [[Rcpp::export]]
arma::dmat calculateFitnessMat_arma(const arma::dcolvec& sel_cof) {
  // ensure RNG gets set/reset
  RNGScope scope;

  // calculate the fitness matrix
  arma::dmat fts_mat = arma::ones<arma::dmat>(4, 4);
  // fts_mat(0, 0) = 1;
  // fts_mat(1, 0) = 1;
  // fts_mat(2, 0) = 1;
  // fts_mat(3, 0) = 1;
  // fts_mat(0, 1) = 1;
  fts_mat(1, 1) = 1 + sel_cof(1);
  // fts_mat(2, 1) = 1;
  fts_mat(3, 1) = 1 + sel_cof(1);
  // fts_mat(0, 2) = 1;
  // fts_mat(1, 2) = 1;
  fts_mat(2, 2) = 1 + sel_cof(0);
  fts_mat(3, 2) = 1 + sel_cof(0);
  // fts_mat(0, 3) = 1;
  fts_mat(1, 3) = 1 + sel_cof(1);
  fts_mat(2, 3) = 1 + sel_cof(0);
  fts_mat(3, 3) = 1 + sel_cof(1);

  // return the fitness matrix for the Wright-Fisher model
  return fts_mat;
}

// Simulate the haplotype frequency trajectories according to the two-locus Wright-Fisher model with selection
// [[Rcpp::export]]
List simulateWFM_arma(const arma::dmat& fts_mat, const double& rec_rat, const int& pop_siz, const arma::dcolvec& int_frq, const int& int_gen, const int& lst_gen) {
  // ensure RNG gets set/reset
  RNGScope scope;

  // declare the haplotype and genotype frequency trajectories
  arma::dmat hap_frq_pth = arma::zeros<arma::dmat>(4, arma::uword(lst_gen - int_gen) + 1);
  arma::dmat gen_frq_pth = arma::zeros<arma::dmat>(10, arma::uword(lst_gen - int_gen) + 1);

  // initialise the haplotype and genotype frequencies in generation 0
  arma::dcolvec hap_frq = int_frq;
  hap_frq_pth.col(0) = hap_frq;
  arma::dmat gen_frq = fts_mat % (hap_frq * hap_frq.t()) / arma::as_scalar(hap_frq.t() * fts_mat * hap_frq);
  gen_frq = arma::diagmat(gen_frq) + 2 * arma::trimatu(gen_frq, 1);
  gen_frq_pth.col(0) = gen_frq(arma::trimatu_ind(arma::size(gen_frq)));

  // declare eta
  arma::dcolvec eta = {-1.0, 1.0, 1.0, -1.0};

  // simulate the haplotype and genotype frequency trajectories
  for (arma::uword k = 1; k < arma::uword(lst_gen - int_gen) + 1; k++) {
    // calculate the sampling probabilities
    arma::dcolvec prob = hap_frq;
    prob = hap_frq % (fts_mat * hap_frq) / arma::as_scalar(hap_frq.t() * fts_mat * hap_frq);
    prob = prob + eta * rec_rat * (prob(0) * prob(3) - prob(1) * prob(2));

    // proceed the Wright-Fisher sampling
    IntegerVector hap_cnt(4);
    R::rmultinom(2 * pop_siz, prob.begin(), 4, hap_cnt.begin());
    hap_frq = as<arma::dcolvec>(hap_cnt) / 2 / pop_siz;
    hap_frq_pth.col(k) = hap_frq;

    gen_frq = fts_mat % (hap_frq * hap_frq.t()) / arma::as_scalar(hap_frq.t() * fts_mat * hap_frq);
    gen_frq = arma::diagmat(gen_frq) + 2 * arma::trimatu(gen_frq, 1);
    gen_frq_pth.col(k) = gen_frq(arma::trimatu_ind(arma::size(gen_frq)));
  }

  // return the haplotype and genotype frequency trajectories under the Wright-Fisher model
  return List::create(Named("hap_frq_pth", hap_frq_pth),
                      Named("gen_frq_pth", gen_frq_pth));
}
/*************************/


/********** WFD **********/
// Simulate the haplotype frequency trajectories according to the two-locus Wright-Fisher diffusion with selection using the Euler-Maruyama method
// [[Rcpp::export]]
arma::dmat simulateWFD_arma(const arma::dcolvec& sel_cof, const double& rec_rat, const int& pop_siz, const arma::dcolvec& int_frq, const int& int_gen, const int& lst_gen, const arma::uword& ptn_num) {
  // ensure RNG gets set/reset
  RNGScope scope;

  // rescale the selection coefficients and the recombination rate
  arma::dcolvec scl_sel_cof = 2 * pop_siz * sel_cof;
  double scl_rec_rat = 4 * pop_siz * rec_rat;

  // calculate delta t
  double dt = 1.0 / (2 * pop_siz) / ptn_num;
  // generate delta W
  arma::dmat dW = pow(dt, 0.5) * arma::randn<arma::dmat>(6, arma::uword(lst_gen - int_gen) * ptn_num);

  // declare the haplotype frequency trajectories
  arma::dmat frq_pth = arma::zeros<arma::dmat>(4, arma::uword(lst_gen - int_gen) * ptn_num + 1);

  // initialise the haplotype frequencies in generation 0
  frq_pth.col(0) = int_frq;

  // simulate the haplotype frequency trajectories
  for (arma::uword t = 1; t < arma::uword(lst_gen - int_gen) * ptn_num + 1; t++) {
    // calculate the drift coefficient vector
    arma::dcolvec mu = arma::zeros<arma::dcolvec>(4);
    double LD = frq_pth(0, t - 1) * frq_pth(3, t - 1) - frq_pth(1, t - 1) * frq_pth(2, t - 1);
    mu(0) = -scl_sel_cof(0) * frq_pth(2, t - 1) * (frq_pth(0, t - 1) * frq_pth(3, t - 1) + frq_pth(0, t - 1) * (frq_pth(2, t - 1) + frq_pth(3, t - 1))) -
      scl_sel_cof(1) * frq_pth(0, t - 1) * (frq_pth(1, t - 1) + frq_pth(3, t - 1)) * (frq_pth(1, t - 1) + frq_pth(3, t - 1)) -
      0.5 * scl_rec_rat * LD;
    mu(1) = -scl_sel_cof(0) * frq_pth(2, t - 1) * (frq_pth(1, t - 1) * frq_pth(3, t - 1) + frq_pth(1, t - 1) * (frq_pth(2, t - 1) + frq_pth(3, t - 1))) +
      scl_sel_cof(1) * frq_pth(1, t - 1) * (frq_pth(0, t - 1) + frq_pth(2, t - 1)) * (frq_pth(1, t - 1) + frq_pth(3, t - 1)) +
      0.5 * scl_rec_rat * LD;
    mu(2) = -scl_sel_cof(0) * frq_pth(2, t - 1) * (frq_pth(2, t - 1) * frq_pth(3, t - 1) + frq_pth(2, t - 1) * (frq_pth(2, t - 1) + frq_pth(3, t - 1)) - (frq_pth(2, t - 1) + frq_pth(3, t - 1))) -
      scl_sel_cof(1) * frq_pth(2, t - 1) * (frq_pth(1, t - 1) + frq_pth(3, t - 1)) * (frq_pth(1, t - 1) + frq_pth(3, t - 1)) +
      0.5 * scl_rec_rat * LD;
    mu(3) = -scl_sel_cof(0) * frq_pth(2, t - 1) * (frq_pth(3, t - 1) * frq_pth(3, t - 1) + frq_pth(3, t - 1) * (frq_pth(2, t - 1) + frq_pth(3, t - 1)) - frq_pth(3, t - 1)) +
      scl_sel_cof(1) * frq_pth(3, t - 1) * (frq_pth(0, t - 1) + frq_pth(2, t - 1)) * (frq_pth(1, t - 1) + frq_pth(3, t - 1)) -
      0.5 * scl_rec_rat * LD;

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
    sigma(2, 4) = 0;
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
    for (arma::uword i = 0; i < 4; i++) {
      if (frq_pth(i, t) < 0) {
        frq_pth(i, t) = 0;
      }
      if (frq_pth(i, t) > 1) {
        frq_pth(i, t) = 1;
      }
    }
    frq_pth.col(t) = frq_pth.col(t) / sum(frq_pth.col(t));
  }

  // return the haplotype frequency trajectories under the Wright-Fisher diffusion
  return frq_pth;
}
/*************************/


/********** BPF **********/
// Calculate the possible genotype counts in the sample
// [[Rcpp::export]]
arma::imat calculateGenoCnt_arma(const arma::icolvec& smp_cnt) {
  // ensure RNG gets set/reset
  RNGScope scope;

  if (smp_cnt.n_elem == 9) {
    arma::imat gen_cnt = arma::zeros<arma::imat>(10, smp_cnt(4) + 1);

    for (int j = 0; j <= smp_cnt(4); j++) {
      gen_cnt(0, j) = smp_cnt(0);
      gen_cnt(1, j) = smp_cnt(1);
      gen_cnt(2, j) = smp_cnt(2);
      gen_cnt(3, j) = smp_cnt(3);
      gen_cnt(4, j) = j;
      gen_cnt(5, j) = smp_cnt(5);
      gen_cnt(6, j) = smp_cnt(4) - j;
      gen_cnt(7, j) = smp_cnt(6);
      gen_cnt(8, j) = smp_cnt(7);
      gen_cnt(9, j) = smp_cnt(8);
    }

    return gen_cnt;
  } else {
    arma::imat gen_cnt = arma::zeros<arma::imat>(10, 1);
    gen_cnt.col(0) = smp_cnt;

    return gen_cnt;
  }
}

// Calculate the genotype frequencies in the adult with the haplotype frequencies
// [[Rcpp::export]]
arma::dcolvec calculateGenoFrq_arma(const arma::dmat& fts_mat, const arma::dcolvec& hap_frq) {
  // ensure RNG gets set/reset
  RNGScope scope;

  arma::dmat gen_frq = fts_mat % (hap_frq * hap_frq.t()) / arma::as_scalar(hap_frq.t() * fts_mat * hap_frq);
  gen_frq = arma::diagmat(gen_frq) + 2 * arma::trimatu(gen_frq, 1);
  arma::dcolvec pop_frq = gen_frq(arma::trimatu_ind(arma::size(gen_frq)));

  return pop_frq;
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

// Calculate the emission probabilities
// [[Rcpp::export]]
double calculateEmissionProb_arma(const arma::icolvec& smp_cnt, const int& smp_siz, const arma::dmat& fts_mat, const arma::dcolvec& hap_frq) {
  // ensure RNG gets set/reset
  RNGScope scope;

  arma::dcolvec pop_frq = calculateGenoFrq_arma(fts_mat, hap_frq);
  double prob = calculateMultinomProb_arma(smp_cnt, smp_siz, pop_frq);

  return prob;
}

// Initialise the particles in the particle filter (uniform generation from the flat Dirichlet distribution)
// [[Rcpp::export]]
arma::dmat initialiseParticle(const arma::uword& pcl_num) {
  // ensure RNG gets set/reset
  RNGScope scope;

  NumericMatrix part(4, pcl_num);
  for (int i = 0; i < 4; i++) {
    part(i, _) = rgamma(pcl_num, 1.0, 1.0);
  }
  for (int j = 0; j < pcl_num; j++) {
    part(_, j) = part(_, j) / sum(part(_, j));
  }

  return as<arma::dmat>(part);
}

// Run the bootstrap particle filter
// [[Rcpp::export]]
List runBPF_arma(const arma::dcolvec& sel_cof, const double& rec_rat, const int& pop_siz, const arma::irowvec& smp_gen, const arma::irowvec& smp_siz, const arma::imat& smp_cnt, const arma::uword& ptn_num, const arma::uword& pcl_num) {
  // ensure RNG gets set/reset
  RNGScope scope;

  double lik = 1;

  arma::dmat wght = arma::zeros<arma::dmat>(pcl_num, smp_gen.n_elem);
  arma::dcube hap_frq_pre = arma::zeros<arma::dcube>(4, pcl_num, smp_gen.n_elem);
  arma::dcube hap_frq_pst = arma::zeros<arma::dcube>(4, pcl_num, smp_gen.n_elem);
  arma::dcube gen_frq_pre = arma::zeros<arma::dcube>(10, pcl_num, smp_gen.n_elem);
  arma::dcube gen_frq_pst = arma::zeros<arma::dcube>(10, pcl_num, smp_gen.n_elem);

  arma::dcolvec wght_tmp = arma::zeros<arma::dcolvec>(pcl_num);
  arma::dmat hap_frq_tmp = arma::zeros<arma::dmat>(4, pcl_num);
  arma::dmat gen_frq_tmp = arma::zeros<arma::dmat>(10, pcl_num);

  arma::dmat fts_mat = calculateFitnessMat_arma(sel_cof);

  // initialise the particles
  cout << "generation: " << smp_gen(0) << endl;
  hap_frq_tmp = initialiseParticle(pcl_num);
  arma::imat gen_cnt = calculateGenoCnt_arma(smp_cnt.col(0));
  for (arma::uword i = 0; i < pcl_num; i++) {
    gen_frq_tmp.col(i) = calculateGenoFrq_arma(fts_mat, hap_frq_tmp.col(i));
    for (arma::uword j = 0; j < gen_cnt.n_cols; j++) {
      wght_tmp(i) = wght_tmp(i) + calculateMultinomProb_arma(gen_cnt.col(j), smp_siz(0), gen_frq_tmp.col(i));
    }
  }

  if (arma::sum(wght_tmp) > 0) {
    arma::dcolvec prob = arma::normalise(wght_tmp, 1);
    arma::ucolvec elem = arma::linspace<arma::ucolvec>(0, pcl_num - 1, pcl_num);
    arma::ucolvec indx = RcppArmadillo::sample(elem, pcl_num, TRUE, prob);

    lik = lik * arma::mean(wght_tmp);
    wght.col(0) = wght_tmp;
    hap_frq_pre.slice(0) = hap_frq_tmp;
    hap_frq_pst.slice(0) = hap_frq_tmp.cols(indx);
    gen_frq_pre.slice(0) = gen_frq_tmp;
    gen_frq_pst.slice(0) = gen_frq_tmp.cols(indx);
  } else {
    lik = 0;
    wght.shed_cols(0, smp_gen.n_elem - 1);
    hap_frq_pre.shed_slices(0, smp_gen.n_elem - 1);
    hap_frq_pst.shed_slices(0, smp_gen.n_elem - 1);
    gen_frq_pre.shed_slices(0, smp_gen.n_elem - 1);
    gen_frq_pst.shed_slices(0, smp_gen.n_elem - 1);

    return List::create(Named("lik", lik),
                        Named("wght", wght),
                        Named("hap_frq_pre_resmp", hap_frq_pre),
                        Named("hap_frq_pst_resmp", hap_frq_pst),
                        Named("gen_frq_pre_resmp", gen_frq_pre),
                        Named("gen_frq_pst_resmp", gen_frq_pst));
  }

  // run the bootstrap particle filter
  for (arma::uword k = 1; k < smp_gen.n_elem; k++) {
    cout << "generation: " << smp_gen(k) << endl;
    wght_tmp = arma::zeros<arma::dcolvec>(pcl_num);
    hap_frq_tmp = hap_frq_pst.slice(k - 1);
    arma::imat gen_cnt = calculateGenoCnt_arma(smp_cnt.col(k));
    for (arma::uword i = 0; i < pcl_num; i++) {
      arma::dmat path = simulateWFD_arma(sel_cof, rec_rat, pop_siz, hap_frq_tmp.col(i), smp_gen(k - 1), smp_gen(k), ptn_num);
      hap_frq_tmp.col(i) = arma::vectorise(path.tail_cols(1), 0);
      gen_frq_tmp.col(i) = calculateGenoFrq_arma(fts_mat, hap_frq_tmp.col(i));
      for (arma::uword j = 0; j < gen_cnt.n_cols; j++) {
        wght_tmp(i) = wght_tmp(i) + calculateMultinomProb_arma(gen_cnt.col(j), smp_siz(k), gen_frq_tmp.col(i));
      }
    }

    if (arma::sum(wght_tmp) > 0) {
      arma::dcolvec prob = arma::normalise(wght_tmp, 1);
      arma::ucolvec elem = arma::linspace<arma::ucolvec>(0, pcl_num - 1, pcl_num);
      arma::ucolvec indx = RcppArmadillo::sample(elem, pcl_num, TRUE, prob);

      lik = lik * arma::mean(wght_tmp);
      wght.col(k) = wght_tmp;
      hap_frq_pre.slice(k) = hap_frq_tmp;
      hap_frq_pst.slice(k) = hap_frq_tmp.cols(indx);
      gen_frq_pre.slice(k) = gen_frq_tmp;
      gen_frq_pst.slice(k) = gen_frq_tmp.cols(indx);
    } else {
      lik = 0;
      wght.shed_cols(k, smp_gen.n_elem - 1);
      hap_frq_pre.shed_slices(k, smp_gen.n_elem - 1);
      hap_frq_pst.shed_slices(k, smp_gen.n_elem - 1);
      gen_frq_pre.shed_slices(k, smp_gen.n_elem - 1);
      gen_frq_pst.shed_slices(k, smp_gen.n_elem - 1);

      break;
    }
  }

  return List::create(Named("lik", lik),
                      Named("wght", wght),
                      Named("hap_frq_pre_resmp", hap_frq_pre),
                      Named("hap_frq_pst_resmp", hap_frq_pst),
                      Named("gen_frq_pre_resmp", gen_frq_pre),
                      Named("gen_frq_pst_resmp", gen_frq_pst));
}
/*************************/


/********** PMMH **********/
// Calculate the log-likelihood using the bootstrap particle filter
// [[Rcpp::export]]
double calculateLogLikelihood_arma(const arma::dcolvec& sel_cof, const double& rec_rat, const int& pop_siz, const arma::irowvec& smp_gen, const arma::irowvec& smp_siz, const arma::field<arma::imat>& ptl_cnt, const arma::uword& ptn_num, const arma::uword& pcl_num) {
  // ensure RNG gets set/reset
  RNGScope scope;

  double log_lik = 0;

  arma::dcolvec wght = arma::zeros<arma::dcolvec>(pcl_num);
  arma::dmat hap_frq_pre = arma::zeros<arma::dmat>(4, pcl_num);
  arma::dmat hap_frq_pst = arma::zeros<arma::dmat>(4, pcl_num);

  arma::dmat fts_mat = calculateFitnessMat_arma(sel_cof);

  // initialise the particles
  hap_frq_pre = initialiseParticle(pcl_num);
  arma::imat gen_cnt = ptl_cnt(0);
  for (arma::uword i = 0; i < pcl_num; i++) {
    for (arma::uword j = 0; j < gen_cnt.n_cols; j++) {
      wght(i) = wght(i) + calculateEmissionProb_arma(gen_cnt.col(j), smp_siz(0), fts_mat, hap_frq_pre.col(i));
    }
  }

  if (arma::mean(wght) > 0) {
    log_lik = log_lik + log(arma::mean(wght));
    arma::dcolvec prob = arma::normalise(wght, 1);
    arma::ucolvec elem = arma::linspace<arma::ucolvec>(0, pcl_num - 1, pcl_num);
    arma::ucolvec indx = RcppArmadillo::sample(elem, pcl_num, TRUE, prob);
    hap_frq_pst = hap_frq_pre.cols(indx);
  } else {
    log_lik = -(arma::datum::inf);

    return log_lik;
  }

  // run the bootstrap particle filter
  for (arma::uword k = 1; k < smp_gen.n_elem; k++) {
    wght = arma::zeros<arma::dcolvec>(pcl_num);
    arma::imat gen_cnt = ptl_cnt(k);
    for (arma::uword i = 0; i < pcl_num; i++) {
      arma::dmat path = simulateWFD_arma(sel_cof, rec_rat, pop_siz, hap_frq_pst.col(i), smp_gen(k - 1), smp_gen(k), ptn_num);
      hap_frq_pre.col(i) = arma::vectorise(path.tail_cols(1), 0);
      for (arma::uword j = 0; j < gen_cnt.n_cols; j++) {
        wght(i) = wght(i) + calculateEmissionProb_arma(gen_cnt.col(j), smp_siz(k), fts_mat, hap_frq_pre.col(i));
      }
    }

    if (arma::mean(wght) > 0) {
      log_lik = log_lik + log(arma::mean(wght));
      arma::dcolvec prob = arma::normalise(wght, 1);
      arma::ucolvec elem = arma::linspace<arma::ucolvec>(0, pcl_num - 1, pcl_num);
      arma::ucolvec indx = RcppArmadillo::sample(elem, pcl_num, TRUE, prob);
      hap_frq_pst = hap_frq_pre.cols(indx);
    } else {
      log_lik = -(arma::datum::inf);

      break;
    }
  }

  return log_lik;
}

// Calculate the optimal particle number in the particle marginal Metropolis-Hastings
// [[Rcpp::export]]
List calculateOptimalParticleNum_arma(const arma::dcolvec& sel_cof, const double& rec_rat, const int& pop_siz, const arma::irowvec& smp_gen, const arma::irowvec& smp_siz, const arma::imat& smp_cnt, const arma::uword& ptn_num, const arma::uword& pcl_num, const arma::uword& gap_num) {
  // ensure RNG gets set/reset
  RNGScope scope;

  arma::field<arma::imat> ptl_cnt(smp_gen.n_elem);
  for (arma::uword k = 0; k < smp_gen.n_elem; k++) {
    ptl_cnt(k) = calculateGenoCnt_arma(smp_cnt.col(k));
  }

  arma::drowvec log_lik(300);
  for (arma::uword i = 0; i < 300; i++) {
    log_lik(i) = calculateLogLikelihood_arma(sel_cof, rec_rat, pop_siz, smp_gen, smp_siz, ptl_cnt, ptn_num, pcl_num);
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
        log_lik(i) = calculateLogLikelihood_arma(sel_cof, rec_rat, pop_siz, smp_gen, smp_siz, ptl_cnt, ptn_num, opt_pcl_num(0));
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
        log_lik(i) = calculateLogLikelihood_arma(sel_cof, rec_rat, pop_siz, smp_gen, smp_siz, ptl_cnt, ptn_num, opt_pcl_num(0));
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
        log_lik(i) = calculateLogLikelihood_arma(sel_cof, rec_rat, pop_siz, smp_gen, smp_siz, ptl_cnt, ptn_num, opt_pcl_num(0));
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
        log_lik(i) = calculateLogLikelihood_arma(sel_cof, rec_rat, pop_siz, smp_gen, smp_siz, ptl_cnt, ptn_num, opt_pcl_num(0));
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
arma::dmat runPMMH_arma(const arma::dcolvec& sel_cof, const double& rec_rat, const int& pop_siz, const arma::irowvec& smp_gen, const arma::irowvec& smp_siz, const arma::imat& smp_cnt, const arma::uword& ptn_num, const arma::uword& pcl_num, const arma::uword& itn_num) {
  // ensure RNG gets set/reset
  RNGScope scope;

  arma::field<arma::imat> ptl_cnt(smp_gen.n_elem);
  for (arma::uword k = 0; k < smp_gen.n_elem; k++) {
    ptl_cnt(k) = calculateGenoCnt_arma(smp_cnt.col(k));
  }

  arma::dmat sel_cof_chn = arma::zeros<arma::dmat>(2, itn_num);

  // arma::drowvec log_pri = arma::zeros<arma::drowvec>(2);
  arma::drowvec log_lik = arma::zeros<arma::drowvec>(2);

  arma::dcolvec sel_cof_sd = {5e-03, 5e-03};

  // initialise the population genetic parameters
  cout << "iteration: " << 1 << endl;
  // take the uniform prior and fix the initial value of the selection coefficient to zero
  // or take the beta prior with alpha = 1 and beta = 3
  sel_cof_chn.col(0) = sel_cof;

  log_lik(0) = calculateLogLikelihood_arma(sel_cof_chn.col(0), rec_rat, pop_siz, smp_gen, smp_siz, ptl_cnt, ptn_num, pcl_num);

  double apt_cnt = 0;
  for (arma::uword i = 1; i < itn_num; i++) {
    cout << "iteration: " << i + 1 << endl;

    // draw the candidates of the selection coefficients from the random walk proposal
    sel_cof_chn.col(i) = sel_cof_chn.col(i - 1) + sel_cof_sd % arma::randn<arma::dcolvec>(2);

    if (arma::any(sel_cof_chn.col(i) < -1)) {
      sel_cof_chn.col(i) = sel_cof_chn.col(i - 1);
      log_lik(1) = log_lik(0);
      // apt_cnt = apt_cnt + 0;
      cout << "acceptance: " << apt_cnt / i << endl;
    } else {
      // calculate the proposal
      // arma::drowvec log_psl = arma::zeros<arma::drowvec>(2);

      // calculate the likelihood
      log_lik(1) = calculateLogLikelihood_arma(sel_cof_chn.col(i), rec_rat, pop_siz, smp_gen, smp_siz, ptl_cnt, ptn_num, pcl_num);

      // calculate the acceptance ratio
      // double apt_rto = exp(log_lik(1) - log_lik(0));
      // double apt_rto = exp((log_pri(1) + log_lik(1) + log_psl(1)) - (log_pri(0) + log_lik(0) + log_psl(0)));

      if (arma::randu() > exp(log_lik(1) - log_lik(0))) {
        sel_cof_chn.col(i) = sel_cof_chn.col(i - 1);
        log_lik(1) = log_lik(0);
        // apt_cnt = apt_cnt + 0;
        cout << "acceptance: " << apt_cnt / i << endl;
      } else {
        log_lik(0) = log_lik(1);
        apt_cnt = apt_cnt + 1;
        cout << "acceptance: " << apt_cnt / i << endl;
      }
    }
  }

  return sel_cof_chn;
}

// Run the adaptive particle marginal Metropolis-Hastings
//[[Rcpp::export]]
arma::dmat runAdaptivePMMH_arma(const arma::dcolvec& sel_cof, const double& rec_rat, const int& pop_siz, const arma::irowvec& smp_gen, const arma::irowvec& smp_siz, const arma::imat& smp_cnt, const arma::uword& ptn_num, const arma::uword& pcl_num, const arma::uword& itn_num, const arma::drowvec& stp_siz, double& apt_rto) {
  // ensure RNG gets set/reset
  RNGScope scope;

  arma::field<arma::imat> ptl_cnt(smp_gen.n_elem);
  for (arma::uword k = 0; k < smp_gen.n_elem; k++) {
    ptl_cnt(k) = calculateGenoCnt_arma(smp_cnt.col(k));
  }

  arma::dmat sel_cof_chn = arma::zeros<arma::dmat>(2, itn_num);

  // arma::drowvec log_pri = arma::zeros<arma::drowvec>(2);
  arma::drowvec log_lik = arma::zeros<arma::drowvec>(2);

  arma::dcolvec U = arma::zeros<arma::dcolvec>(2);
  arma::dmat S = {{5e-03, 0e-00}, {0e-00, 5e-03}};
  arma::dmat M = arma::zeros<arma::dmat>(2, 2);
  arma::dmat I = arma::eye<arma::dmat>(2, 2);

  // initialise the population genetic parameters
  cout << "iteration: " << 1 << endl;
  // take the uniform prior and fix the initial value of the selection coefficient to zero
  // or take the beta prior with alpha = 1 and beta = 3
  sel_cof_chn.col(0) = sel_cof;

  log_lik(0) = calculateLogLikelihood_arma(sel_cof_chn.col(0), rec_rat, pop_siz, smp_gen, smp_siz, ptl_cnt, ptn_num, pcl_num);

  double apt_cnt = 0;
  for (arma::uword i = 1; i < itn_num; i++) {
    cout << "iteration: " << i + 1 << endl;

    // draw the candidates of the selection coefficients from the random walk proposal
    U = arma::randn<arma::dcolvec>(2);
    sel_cof_chn.col(i) = sel_cof_chn.col(i - 1) + S * U;

    double alpha = 0;
    if (arma::any(sel_cof_chn.col(i) < -1)) {
      sel_cof_chn.col(i) = sel_cof_chn.col(i - 1);
      log_lik(1) = log_lik(0);
      // apt_cnt = apt_cnt + 0;
      cout << "acceptance: " << apt_cnt / i << endl;
    } else {
      // calculate the proposal
      // arma::drowvec log_psl = arma::zeros<arma::drowvec>(2);

      // calculate the likelihood
      log_lik(1) = calculateLogLikelihood_arma(sel_cof_chn.col(i), rec_rat, pop_siz, smp_gen, smp_siz, ptl_cnt, ptn_num, pcl_num);

      // calculate the acceptance ratio
      // double apt_rto = exp(log_lik(1) - log_lik(0));
      // double apt_rto = exp((log_pri(1) + log_lik(1) + log_psl(1)) - (log_pri(0) + log_lik(0) + log_psl(0)));

      alpha = exp(log_lik(1) - log_lik(0));
      alpha = (alpha > 1) ? 1 : alpha;
      if (arma::randu() > alpha) {
        sel_cof_chn.col(i) = sel_cof_chn.col(i - 1);
        log_lik(1) = log_lik(0);
        // apt_cnt = apt_cnt + 0;
        cout << "acceptance: " << apt_cnt / i << endl;
      } else {
        log_lik(0) = log_lik(1);
        apt_cnt = apt_cnt + 1;
        cout << "acceptance: " << apt_cnt / i << endl;
      }
    }

    U = arma::normalise(U);
    M = S * (I + stp_siz(i) * (alpha - apt_rto) * (U * U.t())) * S.t();
    S = arma::chol(M, "lower");
    cout << S << endl;
  }

  return sel_cof_chn;
}
/*************************/

