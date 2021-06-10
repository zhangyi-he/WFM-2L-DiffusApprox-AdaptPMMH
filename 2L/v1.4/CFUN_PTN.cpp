// Estimating selection coefficients and testing their changes from ancient DNA data
// Xiaoyang Dai, Mark Beaumont, Feng Yu, Zhangyi He

// version 1.4
// Phenotypes controlled by two genes (genetic linkage and epistatic interaction)
// Non-constant natural selection and non-constant demographic histories
// Prior knowledge from modern samples (gene polymorphism)
// Joint estimation of the underlying trajectory of haplotype frequencies

// Horse white coat patterns (KIT13 & KIT116)

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
  fts_mat(1, 0) = 1 + sel_cof(1);
  fts_mat(2, 0) = 1 + sel_cof(0);
  fts_mat(3, 0) = 1 + sel_cof(2);
  fts_mat(0, 1) = 1 + sel_cof(1);
  fts_mat(1, 1) = 1 + sel_cof(1);
  fts_mat(2, 1) = 1 + sel_cof(2);
  fts_mat(3, 1) = 1 + sel_cof(2);
  fts_mat(0, 2) = 1 + sel_cof(0);
  fts_mat(1, 2) = 1 + sel_cof(2);
  fts_mat(2, 2) = 1 + sel_cof(0);
  fts_mat(3, 2) = 1 + sel_cof(2);
  fts_mat(0, 3) = 1 + sel_cof(2);
  fts_mat(1, 3) = 1 + sel_cof(2);
  fts_mat(2, 3) = 1 + sel_cof(2);
  fts_mat(3, 3) = 1 + sel_cof(2);

  // return the fitness matrix for the Wright-Fisher model
  return fts_mat;
}

// Simulate the haplotype frequency trajectories according to the two-locus Wright-Fisher model with selection
// [[Rcpp::export]]
List simulateWFM_arma(const arma::dmat& fts_mat, const double& rec_rat, const arma::icolvec& pop_siz, const arma::dcolvec& int_frq, const int& int_gen, const int& lst_gen) {
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
    R::rmultinom(2 * pop_siz(k), prob.begin(), 4, hap_cnt.begin());
    hap_frq = as<arma::dcolvec>(hap_cnt) / 2 / pop_siz(k);
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
arma::dmat simulateWFD_arma(const arma::dcolvec& sel_cof, const double& rec_rat, const arma::icolvec& pop_siz, const int& ref_siz, const arma::dcolvec& int_frq, const int& int_gen, const int& lst_gen, const arma::uword& ptn_num) {
  // ensure RNG gets set/reset
  RNGScope scope;

  // rescale the selection coefficients and the recombination rate
  arma::dcolvec scl_sel_cof = 2 * ref_siz * sel_cof;
  double scl_rec_rat = 4 * ref_siz * rec_rat;

  // calculate the ratio of the population size to the reference population size
  arma::dcolvec siz_rto = arma::zeros<arma::dcolvec>(arma::uword(lst_gen - int_gen) * ptn_num);
  for (arma::uword k = 0; k < arma::uword(lst_gen - int_gen); k++) {
    arma::dcolvec siz_rto_tmp = arma::zeros<arma::dcolvec>(ptn_num);
    siz_rto_tmp.fill(pop_siz(k) / double(ref_siz));
    siz_rto.subvec(k * ptn_num, (k + 1) * ptn_num - 1) = siz_rto_tmp;
  }

  // calculate delta t
  double dt = 1.0 / (2 * ref_siz) / ptn_num;
  // generate delta W
  arma::dmat dW = pow(dt, 0.5) * arma::randn<arma::dmat>(6, arma::uword(lst_gen - int_gen) * ptn_num);

  // declare the haplotype frequency trajectories
  arma::dmat hap_frq_pth = arma::zeros<arma::dmat>(4, arma::uword(lst_gen - int_gen) * ptn_num + 1);

  // initialise the haplotype frequencies in generation 0
  hap_frq_pth.col(0) = int_frq;

  // simulate the haplotype frequency trajectories
  for (arma::uword t = 1; t < arma::uword(lst_gen - int_gen) * ptn_num + 1; t++) {
    // calculate the drift coefficient vector
    arma::dcolvec mu = arma::zeros<arma::dcolvec>(4);
    double LD = hap_frq_pth(0, t - 1) * hap_frq_pth(3, t - 1) - hap_frq_pth(1, t - 1) * hap_frq_pth(2, t - 1);
    mu(0) = -scl_sel_cof(0) * (hap_frq_pth(0, t - 1) * hap_frq_pth(2, t - 1) * (hap_frq_pth(0, t - 1) + hap_frq_pth(0, t - 1) + hap_frq_pth(2, t - 1)) - hap_frq_pth(2, t - 1) * hap_frq_pth(0, t - 1)) -
      scl_sel_cof(1) * (hap_frq_pth(0, t - 1) * hap_frq_pth(1, t - 1) * (hap_frq_pth(0, t - 1) + hap_frq_pth(0, t - 1) + hap_frq_pth(1, t - 1)) - hap_frq_pth(1, t - 1) * hap_frq_pth(0, t - 1)) -
      scl_sel_cof(2) * hap_frq_pth(0, t - 1) * (2 * hap_frq_pth(1, t - 1) * hap_frq_pth(2, t - 1) + hap_frq_pth(3, t - 1) - hap_frq_pth(3, t - 1) * hap_frq_pth(3, t - 1)) -
      0.5 * scl_rec_rat * LD;
    mu(1) = -scl_sel_cof(0) * hap_frq_pth(1, t - 1) * hap_frq_pth(2, t - 1) * (hap_frq_pth(0, t - 1) + hap_frq_pth(0 , t - 1) + hap_frq_pth(2, t - 1)) -
      scl_sel_cof(1) * (hap_frq_pth(1, t - 1) * hap_frq_pth(1, t - 1) * (hap_frq_pth(0, t - 1) + hap_frq_pth(0 , t - 1) + hap_frq_pth(1, t - 1)) - hap_frq_pth(1, t - 1) * (hap_frq_pth(0, t - 1) + hap_frq_pth(1, t - 1))) -
      scl_sel_cof(2) * (hap_frq_pth(1, t - 1) * (2 * hap_frq_pth(1, t - 1) * hap_frq_pth(2, t - 1) + hap_frq_pth(3, t - 1) - hap_frq_pth(3, t - 1) * hap_frq_pth(3, t - 1)) - hap_frq_pth(1, t - 1) * hap_frq_pth(2, t - 1)) -
      0.5 * scl_rec_rat * LD;
    mu(2) = -scl_sel_cof(0) * (hap_frq_pth(2, t - 1) * hap_frq_pth(2, t - 1) * (hap_frq_pth(0, t - 1) + hap_frq_pth(0, t - 1) + hap_frq_pth(2, t - 1)) - hap_frq_pth(2, t - 1) * (hap_frq_pth(0, t - 1) + hap_frq_pth(2, t - 1))) -
      scl_sel_cof(1) * hap_frq_pth(2, t - 1) * hap_frq_pth(1, t - 1) * (hap_frq_pth(0, t - 1) + hap_frq_pth(0, t - 1) + hap_frq_pth(1, t - 1)) -
      scl_sel_cof(2) * (hap_frq_pth(2, t - 1) * (2 * hap_frq_pth(1, t - 1) * hap_frq_pth(2, t - 1) + hap_frq_pth(3, t - 1) - hap_frq_pth(3, t - 1) * hap_frq_pth(3, t - 1)) - hap_frq_pth(1, t - 1) * hap_frq_pth(2, t - 1)) -
      0.5 * scl_rec_rat * LD;
    mu(3) = -scl_sel_cof(0) * hap_frq_pth(3, t - 1) * hap_frq_pth(2, t - 1) * (hap_frq_pth(0, t - 1) + hap_frq_pth(0, t - 1) + hap_frq_pth(2, t - 1)) -
      scl_sel_cof(1) * hap_frq_pth(3, t - 1) * hap_frq_pth(1, t - 1) * (hap_frq_pth(0, t - 1) + hap_frq_pth(0, t - 1) + hap_frq_pth(1, t - 1)) -
      scl_sel_cof(2) * (hap_frq_pth(3, t - 1) * (2 * hap_frq_pth(1, t - 1) * hap_frq_pth(2, t - 1) + hap_frq_pth(3, t - 1) - hap_frq_pth(3, t - 1) * hap_frq_pth(3, t - 1)) - hap_frq_pth(3, t - 1) + hap_frq_pth(3, t - 1) * hap_frq_pth(3, t - 1)) -
      0.5 * scl_rec_rat * LD;

    // calculate the diffusion coefficient matrix
    arma::dmat sigma = arma::zeros<arma::dmat>(4, 6);
    sigma(0, 0) = pow(hap_frq_pth(0, t - 1) * hap_frq_pth(1, t - 1), 0.5);
    sigma(0, 1) = pow(hap_frq_pth(0, t - 1) * hap_frq_pth(2, t - 1), 0.5);
    sigma(0, 2) = pow(hap_frq_pth(0, t - 1) * hap_frq_pth(3, t - 1), 0.5);
    // sigma(0, 3) = 0;
    // sigma(0, 4) = 0;
    // sigma(0, 5) = 0;
    sigma(1, 0) = -pow(hap_frq_pth(1, t - 1) * hap_frq_pth(0, t - 1), 0.5);
    // sigma(1, 1) = 0;
    // sigma(1, 2) = 0;
    sigma(1, 3) = pow(hap_frq_pth(1, t - 1) * hap_frq_pth(2, t - 1), 0.5);
    sigma(1, 4) = pow(hap_frq_pth(1, t - 1) * hap_frq_pth(3, t - 1), 0.5);
    // sigma(1, 5) = 0;
    // sigma(2, 0) = 0;
    sigma(2, 1) = -pow(hap_frq_pth(2, t - 1) * hap_frq_pth(0, t - 1), 0.5);
    // sigma(2, 2) = 0;
    sigma(2, 3) = -pow(hap_frq_pth(2, t - 1) * hap_frq_pth(1, t - 1), 0.5);
    sigma(2, 4) = 0;
    sigma(2, 5) = pow(hap_frq_pth(2, t - 1) * hap_frq_pth(3, t - 1), 0.5);
    // sigma(3, 0) = 0;
    // sigma(3, 1) = 0;
    sigma(3, 2) = -pow(hap_frq_pth(3, t - 1) * hap_frq_pth(0, t - 1), 0.5);
    // sigma(3, 3) = 0;
    sigma(3, 4) = -pow(hap_frq_pth(3, t - 1) * hap_frq_pth(1, t - 1), 0.5);
    sigma(3, 5) = -pow(hap_frq_pth(3, t - 1) * hap_frq_pth(2, t - 1), 0.5);
    sigma = pow(siz_rto(t - 1), -0.5) * sigma;

    // proceed the Euler-Maruyama scheme
    hap_frq_pth.col(t) = hap_frq_pth.col(t - 1) + mu * dt + sigma * dW.col(t - 1);

    // remove the noise from the numerical techniques
    for (arma::uword i = 0; i < 4; i++) {
      if (hap_frq_pth(i, t) < 0) {
        hap_frq_pth(i, t) = 0;
      }
      if (hap_frq_pth(i, t) > 1) {
        hap_frq_pth(i, t) = 1;
      }
    }
    hap_frq_pth.col(t) = hap_frq_pth.col(t) / sum(hap_frq_pth.col(t));
  }

  // return the haplotype frequency trajectories under the Wright-Fisher diffusion
  return hap_frq_pth;
}
/*************************/


/********** BPF **********/
// Calculate the possible genotype counts in the sample
// [[Rcpp::export]]
arma::imat calculateGenoCnt_arma(const arma::icolvec& smp_cnt) {
  // ensure RNG gets set/reset
  RNGScope scope;

  if (smp_cnt.n_elem == 9) {
    arma::imat gen_cnt = arma::zeros<arma::imat>(10, smp_cnt(5) + 1);

    for (int j = 0; j <= smp_cnt(5); j++) {
      gen_cnt(0, j) = smp_cnt(0);
      gen_cnt(1, j) = smp_cnt(1); // sabino
      gen_cnt(2, j) = smp_cnt(2); // sabino
      gen_cnt(3, j) = smp_cnt(3); // tobinao
      gen_cnt(4, j) = j;
      gen_cnt(5, j) = smp_cnt(5); // tobiano
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

  if (hap_frq(1) + hap_frq(3) == 0 || hap_frq(2) + hap_frq(3) == 0) {
    prob = 0; // the mutant alleles are observed in modern samples although they may not be found in ancient samples
  }

  return prob;
}

// Initialise the particles in the particle filter
// [[Rcpp::export]]
arma::dmat initialiseParticle_arma(const arma::uword& pcl_num){
  // ensure RNG gets set/reset
  RNGScope scope;

  arma::dmat part = arma::zeros<arma::dmat>(4, pcl_num);
  arma::dmat mut_frq = arma::randu<arma::dmat>(2, pcl_num);
  arma::drowvec ld = arma::randu<arma::drowvec>(pcl_num);
  for (arma::uword i = 0; i < pcl_num; i++) {
    double a = -mut_frq(0, i) * mut_frq(1, i);
    a = (a >= -(1 - mut_frq(0, i)) * (1 - mut_frq(1, i)))? a : -(1 - mut_frq(0, i)) * (1 - mut_frq(1, i));
    double b = mut_frq(0, i) * (1 - mut_frq(1, i));
    b = (b <= (1 - mut_frq(0, i)) * mut_frq(1, i))? b : (1 - mut_frq(0, i)) * mut_frq(1, i);
    ld(i) = a + (b - a) * ld(i);
  }
  part.row(0) = (1 - mut_frq.row(0)) % (1 - mut_frq.row(1)) + ld;
  part.row(1) = (1 - mut_frq.row(0)) % mut_frq.row(1) - ld;
  part.row(2) = mut_frq.row(0) % (1 - mut_frq.row(1)) - ld;
  part.row(3) = mut_frq.row(0) % mut_frq.row(1) + ld;

  return part;

  // NumericMatrix part(4, pcl_num);
  // for (int i = 0; i < 4; i++) {
  //   part(i, _) = rgamma(pcl_num, 1.0, 1.0);
  // }
  // for (int j = 0; j < pcl_num; j++) {
  //   part(_, j) = part(_, j) / sum(part(_, j));
  // }
  //
  // return as<arma::dmat>(part);
}

// Run the bootstrap particle filter
// [[Rcpp::export]]
List runBPF_arma(const arma::dmat& sel_cof, const double& rec_rat, const arma::icolvec& pop_siz, const int& ref_siz, const arma::irowvec& smp_gen, const arma::irowvec& smp_siz, const arma::imat& smp_cnt, const arma::uword& ptn_num, const arma::uword& pcl_num) {
  // ensure RNG gets set/reset
  RNGScope scope;

  double lik = 1;

  arma::uword evt_ind = arma::as_scalar(arma::find(smp_siz == 0));

  arma::dmat wght = arma::zeros<arma::dmat>(pcl_num, smp_gen.n_elem);
  arma::dcube hap_frq_pth = arma::zeros<arma::dcube>(4, arma::uword(arma::max(smp_gen) - arma::min(smp_gen)) * ptn_num + 1, pcl_num);
  arma::dcube hap_frq_pre = arma::zeros<arma::dcube>(4, pcl_num, smp_gen.n_elem);
  arma::dcube hap_frq_pst = arma::zeros<arma::dcube>(4, pcl_num, smp_gen.n_elem);
  arma::dcube gen_frq_pre = arma::zeros<arma::dcube>(10, pcl_num, smp_gen.n_elem);
  arma::dcube gen_frq_pst = arma::zeros<arma::dcube>(10, pcl_num, smp_gen.n_elem);

  arma::dcolvec wght_tmp = arma::zeros<arma::dcolvec>(pcl_num);
  arma::dmat hap_frq_tmp = arma::zeros<arma::dmat>(4, pcl_num);
  arma::dmat gen_frq_tmp = arma::zeros<arma::dmat>(10, pcl_num);

  // before the event of interest
  arma::dmat fts_mat = calculateFitnessMat_arma(sel_cof.col(0));

  // initialise the particles
  cout << "generation: " << smp_gen(0) << endl;
  hap_frq_tmp = initialiseParticle_arma(pcl_num);
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
                        Named("hap_frq_pth", hap_frq_pth),
                        Named("hap_frq_pre_resmp", hap_frq_pre),
                        Named("hap_frq_pst_resmp", hap_frq_pst),
                        Named("gen_frq_pre_resmp", gen_frq_pre),
                        Named("gen_frq_pst_resmp", gen_frq_pst));
  }

  // run the bootstrap particle filter
  for (arma::uword k = 1; k < evt_ind; k++) {
    cout << "generation: " << smp_gen(k) << endl;
    wght_tmp = arma::zeros<arma::dcolvec>(pcl_num);
    hap_frq_tmp = hap_frq_pst.slice(k - 1);
    arma::imat gen_cnt = calculateGenoCnt_arma(smp_cnt.col(k));
    for (arma::uword i = 0; i < pcl_num; i++) {
      arma::dmat path = simulateWFD_arma(sel_cof.col(0), rec_rat, pop_siz.subvec(smp_gen(k - 1), smp_gen(k)), ref_siz, hap_frq_tmp.col(i), smp_gen(k - 1), smp_gen(k), ptn_num);
      hap_frq_pth.subcube(0, (smp_gen(k - 1) - smp_gen(0)) * ptn_num, i, 3, (smp_gen(k) - smp_gen(0)) * ptn_num, i) = path;
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
      hap_frq_pth = hap_frq_pth.slices(indx);
      gen_frq_pre.slice(k) = gen_frq_tmp;
      gen_frq_pst.slice(k) = gen_frq_tmp.cols(indx);
    } else {
      lik = 0;
      wght.shed_cols(k, smp_gen.n_elem - 1);
      hap_frq_pre.shed_slices(k, smp_gen.n_elem - 1);
      hap_frq_pst.shed_slices(k, smp_gen.n_elem - 1);
      gen_frq_pre.shed_slices(k, smp_gen.n_elem - 1);
      gen_frq_pst.shed_slices(k, smp_gen.n_elem - 1);

      return List::create(Named("lik", lik),
                    Named("wght", wght),
                    Named("hap_frq_pth", hap_frq_pth),
                    Named("hap_frq_pre_resmp", hap_frq_pre),
                    Named("hap_frq_pst_resmp", hap_frq_pst),
                    Named("gen_frq_pre_resmp", gen_frq_pre),
                    Named("gen_frq_pst_resmp", gen_frq_pst));
    }
  }

  // after the event of interest
  fts_mat = calculateFitnessMat_arma(sel_cof.col(1));

  if (smp_gen(evt_ind) == smp_gen(evt_ind - 1)) {
    hap_frq_pre.slice(evt_ind) = hap_frq_pre.slice(evt_ind - 1);
    hap_frq_pst.slice(evt_ind) = hap_frq_pst.slice(evt_ind - 1);
    gen_frq_pre.slice(evt_ind) = gen_frq_pre.slice(evt_ind - 1);
    gen_frq_pst.slice(evt_ind) = gen_frq_pst.slice(evt_ind - 1);
  } else {
    cout << "generation: " << smp_gen(evt_ind) << endl;
    hap_frq_tmp = hap_frq_pst.slice(evt_ind - 1);
    for (arma::uword i = 0; i < pcl_num; i++) {
      arma::dmat path = simulateWFD_arma(sel_cof.col(1), rec_rat, pop_siz.subvec(smp_gen(evt_ind - 1), smp_gen(evt_ind)), ref_siz, hap_frq_tmp.col(i), smp_gen(evt_ind - 1), smp_gen(evt_ind), ptn_num);
      hap_frq_pth.subcube(0, (smp_gen(evt_ind - 1) - smp_gen(0)) * ptn_num, i, 3, (smp_gen(evt_ind) - smp_gen(0)) * ptn_num, i) = path;
      hap_frq_tmp.col(i) = arma::vectorise(path.tail_cols(1), 0);
      gen_frq_tmp.col(i) = calculateGenoFrq_arma(fts_mat, hap_frq_tmp.col(i));
    }
    hap_frq_pre.slice(evt_ind) = hap_frq_tmp;
    hap_frq_pst.slice(evt_ind) = hap_frq_tmp;
    gen_frq_pre.slice(evt_ind) = gen_frq_tmp;
    gen_frq_pst.slice(evt_ind) = gen_frq_tmp;
  }

  for (arma::uword k = evt_ind + 1; k < smp_gen.n_elem; k++) {
    cout << "generation: " << smp_gen(k) << endl;
    wght_tmp = arma::zeros<arma::dcolvec>(pcl_num);
    hap_frq_tmp = hap_frq_pst.slice(k - 1);
    arma::imat gen_cnt = calculateGenoCnt_arma(smp_cnt.col(k));
    for (arma::uword i = 0; i < pcl_num; i++) {
      arma::dmat path = simulateWFD_arma(sel_cof.col(1), rec_rat, pop_siz.subvec(smp_gen(k - 1), smp_gen(k)), ref_siz, hap_frq_tmp.col(i), smp_gen(k - 1), smp_gen(k), ptn_num);
      hap_frq_pth.subcube(0, (smp_gen(k - 1) - smp_gen(0)) * ptn_num, i, 3, (smp_gen(k) - smp_gen(0)) * ptn_num, i) = path;
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
      hap_frq_pth = hap_frq_pth.slices(indx);
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

  if (smp_gen(evt_ind) != smp_gen(evt_ind - 1)) {
    wght.shed_col(evt_ind);
    hap_frq_pre.shed_slice(evt_ind);
    hap_frq_pst.shed_slice(evt_ind);
    gen_frq_pre.shed_slice(evt_ind);
    gen_frq_pst.shed_slice(evt_ind);
  }

  return List::create(Named("lik", lik),
                      Named("wght", wght),
                      Named("hap_frq_pth", hap_frq_pth),
                      Named("hap_frq_pre_resmp", hap_frq_pre),
                      Named("hap_frq_pst_resmp", hap_frq_pst),
                      Named("gen_frq_pre_resmp", gen_frq_pre),
                      Named("gen_frq_pst_resmp", gen_frq_pst));
}
/*************************/


/********** PMMH **********/
// Calculate the log-likelihood using the bootstrap particle filter
// [[Rcpp::export]]
void calculateLogLikelihood_arma(double& log_lik, arma::dmat& frq_pth, const arma::dmat& sel_cof, const double& rec_rat, const arma::icolvec& pop_siz, const int& ref_siz, const arma::irowvec& smp_gen, const arma::irowvec& smp_siz, const arma::field<arma::imat>& ptl_cnt, const arma::uword& ptn_num, const arma::uword& pcl_num) {
  // ensure RNG gets set/reset
  RNGScope scope;

  log_lik = 0;
  frq_pth = arma::zeros<arma::dmat>(4, arma::uword(arma::max(smp_gen) - arma::min(smp_gen)) * ptn_num + 1);

  arma::uword evt_ind = arma::as_scalar(arma::find(smp_siz == 0));

  arma::dcolvec wght = arma::zeros<arma::dcolvec>(pcl_num);
  arma::dcube hap_frq_pth = arma::zeros<arma::dcube>(4, arma::uword(arma::max(smp_gen) - arma::min(smp_gen)) * ptn_num + 1, pcl_num);
  arma::dmat hap_frq_pre = arma::zeros<arma::dmat>(4, pcl_num);
  arma::dmat hap_frq_pst = arma::zeros<arma::dmat>(4, pcl_num);

  // before the event of interest
  arma::dmat fts_mat = calculateFitnessMat_arma(sel_cof.col(0));

  // initialise the particles
  hap_frq_pre = initialiseParticle_arma(pcl_num);
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

    return;
  }

  // run the bootstrap particle filter
  for (arma::uword k = 1; k < evt_ind; k++) {
    wght = arma::zeros<arma::dcolvec>(pcl_num);
    arma::imat gen_cnt = ptl_cnt(k);
    for (arma::uword i = 0; i < pcl_num; i++) {
      arma::dmat path = simulateWFD_arma(sel_cof.col(0), rec_rat, pop_siz.subvec(smp_gen(k - 1), smp_gen(k)), ref_siz, hap_frq_pst.col(i), smp_gen(k - 1), smp_gen(k), ptn_num);
      hap_frq_pth.subcube(0, (smp_gen(k - 1) - smp_gen(0)) * ptn_num, i, 3, (smp_gen(k) - smp_gen(0)) * ptn_num, i) = path;
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
      hap_frq_pth = hap_frq_pth.slices(indx);
    } else {
      log_lik = -(arma::datum::inf);

      return;
    }
  }

  // after the event of interest
  fts_mat = calculateFitnessMat_arma(sel_cof.col(1));

  if (smp_gen(evt_ind) != smp_gen(evt_ind - 1)) {
    for (arma::uword i = 0; i < pcl_num; i++) {
      arma::dmat path = simulateWFD_arma(sel_cof.col(1), rec_rat, pop_siz.subvec(smp_gen(evt_ind - 1), smp_gen(evt_ind)), ref_siz, hap_frq_pst.col(i), smp_gen(evt_ind - 1), smp_gen(evt_ind), ptn_num);
      hap_frq_pth.subcube(0, (smp_gen(evt_ind - 1) - smp_gen(0)) * ptn_num, i, 3, (smp_gen(evt_ind) - smp_gen(0)) * ptn_num, i) = path;
      hap_frq_pre.col(i) = arma::vectorise(path.tail_cols(1), 0);
    }
    hap_frq_pst = hap_frq_pre;
  }

  for (arma::uword k = evt_ind + 1; k < smp_gen.n_elem; k++) {
    wght = arma::zeros<arma::dcolvec>(pcl_num);
    arma::imat gen_cnt = ptl_cnt(k);
    for (arma::uword i = 0; i < pcl_num; i++) {
      arma::dmat path = simulateWFD_arma(sel_cof.col(1), rec_rat, pop_siz.subvec(smp_gen(k - 1), smp_gen(k)), ref_siz, hap_frq_pst.col(i), smp_gen(k - 1), smp_gen(k), ptn_num);
      hap_frq_pth.subcube(0, (smp_gen(k - 1) - smp_gen(0)) * ptn_num, i, 3, (smp_gen(k) - smp_gen(0)) * ptn_num, i) = path;
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
      hap_frq_pth = hap_frq_pth.slices(indx);
    } else {
      log_lik = -(arma::datum::inf);

      return;
    }
  }

  frq_pth = hap_frq_pth.slice(0);

  return;
}

// Calculate the optimal particle number in the particle marginal Metropolis-Hastings
// [[Rcpp::export]]
List calculateOptimalParticleNum_arma(const arma::dmat& sel_cof, const double& rec_rat, const arma::icolvec& pop_siz, const int& ref_siz, const arma::irowvec& smp_gen, const arma::irowvec& smp_siz, const arma::imat& smp_cnt, const arma::uword& ptn_num, const arma::uword& pcl_num, const arma::uword& gap_num) {
  // ensure RNG gets set/reset
  RNGScope scope;

  arma::field<arma::imat> ptl_cnt(smp_gen.n_elem);
  for (arma::uword k = 0; k < smp_gen.n_elem; k++) {
    ptl_cnt(k) = calculateGenoCnt_arma(smp_cnt.col(k));
  }

  arma::drowvec log_lik(300);
  arma::dmat frq_pth = arma::zeros<arma::dmat>(4, arma::uword(arma::max(smp_gen) - arma::min(smp_gen)) * ptn_num + 1);
  for (arma::uword i = 0; i < 300; i++) {
    calculateLogLikelihood_arma(log_lik(i), frq_pth, sel_cof, rec_rat, pop_siz, ref_siz, smp_gen, smp_siz, ptl_cnt, ptn_num, pcl_num);
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
        calculateLogLikelihood_arma(log_lik(i), frq_pth, sel_cof, rec_rat, pop_siz, ref_siz, smp_gen, smp_siz, ptl_cnt, ptn_num, opt_pcl_num(0));
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
        calculateLogLikelihood_arma(log_lik(i), frq_pth, sel_cof, rec_rat, pop_siz, ref_siz, smp_gen, smp_siz, ptl_cnt, ptn_num, opt_pcl_num(0));
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
        calculateLogLikelihood_arma(log_lik(i), frq_pth, sel_cof, rec_rat, pop_siz, ref_siz, smp_gen, smp_siz, ptl_cnt, ptn_num, opt_pcl_num(0));
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
        calculateLogLikelihood_arma(log_lik(i), frq_pth, sel_cof, rec_rat, pop_siz, ref_siz, smp_gen, smp_siz, ptl_cnt, ptn_num, opt_pcl_num(0));
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
List runPMMH_arma(const arma::dmat& sel_cof, const double& rec_rat, const arma::icolvec& pop_siz, const int& ref_siz, const arma::irowvec& smp_gen, const arma::irowvec& smp_siz, const arma::imat& smp_cnt, const arma::uword& ptn_num, const arma::uword& pcl_num, const arma::uword& itn_num) {
  // ensure RNG gets set/reset
  RNGScope scope;

  arma::field<arma::imat> ptl_cnt(smp_gen.n_elem);
  for (arma::uword k = 0; k < smp_gen.n_elem; k++) {
    ptl_cnt(k) = calculateGenoCnt_arma(smp_cnt.col(k));
  }

  arma::dcube sel_cof_chn = arma::zeros<arma::dcube>(3, 2, itn_num);
  arma::dcube frq_pth_chn = arma::zeros<arma::dcube>(4, arma::uword(max(smp_gen) - min(smp_gen)) * ptn_num + 1, itn_num);

  //arma::drowvec log_pri = arma::zeros<arma::drowvec>(2);
  arma::drowvec log_lik = arma::zeros<arma::drowvec>(2);
  arma::dmat frq_pth = arma::zeros<arma::dmat>(4, arma::uword(arma::max(smp_gen) - arma::min(smp_gen)) * ptn_num + 1);

  arma::dmat sel_cof_sd = {{5e-03, 5e-03},
                           {5e-03, 1e-02},
                           {1e-02, 1e-02}};

  // initialise the population genetic parameters
  cout << "iteration: " << 1 << endl;
  // take the uniform prior and fix the initial value of the selection coefficient to zero
  // or take the beta prior with alpha = 1 and beta = 3
  sel_cof_chn.slice(0) = sel_cof;

  calculateLogLikelihood_arma(log_lik(0), frq_pth, sel_cof_chn.slice(0), rec_rat, pop_siz, ref_siz, smp_gen, smp_siz, ptl_cnt, ptn_num, pcl_num);
  frq_pth_chn.slice(0) = frq_pth;

  double apt_cnt = 0;
  for (arma::uword i = 1; i < itn_num; i++) {
    cout << "iteration: " << i + 1 << endl;

    // draw the candidates of the selection coefficients from the random walk proposal
    sel_cof_chn.slice(i) = sel_cof_chn.slice(i - 1) + sel_cof_sd % arma::randn<arma::dmat>(3, 2);

    if (arma::any(arma::any(sel_cof_chn.slice(i) < -1, 1))) {
      sel_cof_chn.slice(i) = sel_cof_chn.slice(i - 1);
      log_lik(1) = log_lik(0);
      // apt_cnt = apt_cnt + 0;
      cout << "acceptance: " << apt_cnt / i << endl;
    } else {
      // calculate the proposal
      //arma::drowvec log_psl = arma::zeros<arma::drowvec>(2);

      // calculate the likelihood
      calculateLogLikelihood_arma(log_lik(1), frq_pth, sel_cof_chn.slice(i), rec_rat, pop_siz, ref_siz, smp_gen, smp_siz, ptl_cnt, ptn_num, pcl_num);
      frq_pth_chn.slice(i) = frq_pth;

      // calculate the acceptance ratio
      // double apt_rto = exp(log_lik(1) - log_lik(0));
      // double apt_rto = exp((log_pri(1) + log_lik(1) + log_psl(1)) - (log_pri(0) + log_lik(0) + log_psl(0)));

      if (arma::randu() > exp(log_lik(1) - log_lik(0))) {
        sel_cof_chn.slice(i) = sel_cof_chn.slice(i - 1);
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

  return List::create(Named("sel_cof_chn", sel_cof_chn),
                      Named("frq_pth_chn", frq_pth_chn));
}

// Run the adaptive particle marginal Metropolis-Hastings
//[[Rcpp::export]]
List runAdaptPMMH_arma(const arma::dmat& sel_cof, const double& rec_rat, const arma::icolvec& pop_siz, const int& ref_siz, const arma::irowvec& smp_gen, const arma::irowvec& smp_siz, const arma::imat& smp_cnt, const arma::uword& ptn_num, const arma::uword& pcl_num, const arma::uword& itn_num, const arma::drowvec& stp_siz, double& apt_rto) {
  // ensure RNG gets set/reset
  RNGScope scope;

  arma::field<arma::imat> ptl_cnt(smp_gen.n_elem);
  for (arma::uword k = 0; k < smp_gen.n_elem; k++) {
    ptl_cnt(k) = calculateGenoCnt_arma(smp_cnt.col(k));
  }

  arma::dcube sel_cof_chn = arma::zeros<arma::dcube>(3, 2, itn_num);
  arma::dcolvec sel_cof_tmp = arma::vectorise(sel_cof);
  arma::dcube frq_pth_chn = arma::zeros<arma::dcube>(4, arma::uword(max(smp_gen) - min(smp_gen)) * ptn_num + 1, itn_num);

  //arma::drowvec log_pri = arma::zeros<arma::drowvec>(2);
  arma::drowvec log_lik = arma::zeros<arma::drowvec>(2);
  arma::dmat frq_pth = arma::zeros<arma::dmat>(4, arma::uword(arma::max(smp_gen) - arma::min(smp_gen)) * ptn_num + 1);

  arma::dcolvec U = arma::zeros<arma::dcolvec>(6);
  arma::dmat S = {{5e-03, 0e-00, 0e-00, 0e-00, 0e-00, 0e-00}, 
                  {0e-00, 5e-03, 0e-00, 0e-00, 0e-00, 0e-00}, 
                  {0e-00, 0e-00, 1e-02, 0e-00, 0e-00, 0e-00}, 
                  {0e-00, 0e-00, 0e-00, 5e-03, 0e-00, 0e-00},
                  {0e-00, 0e-00, 0e-00, 0e-00, 1e-02, 0e-00},
                  {0e-00, 0e-00, 0e-00, 0e-00, 0e-00, 1e-02}};
  arma::dmat M = arma::zeros<arma::dmat>(6, 6);
  arma::dmat I = arma::eye<arma::dmat>(6, 6);

  // initialise the population genetic parameters
  cout << "iteration: " << 1 << endl;
  // take the uniform prior and fix the initial value of the selection coefficient to zero
  // or take the beta prior with alpha = 1 and beta = 3
  sel_cof_chn.slice(0) = sel_cof;
  // sel_cof_chn.slice(0) = arma::reshape(sel_cof_tmp, 2, 2);

  calculateLogLikelihood_arma(log_lik(0), frq_pth, sel_cof_chn.slice(0), rec_rat, pop_siz, ref_siz, smp_gen, smp_siz, ptl_cnt, ptn_num, pcl_num);
  frq_pth_chn.slice(0) = frq_pth;

  double apt_cnt = 0;
  double alpha = 0;
  for (arma::uword i = 1; i < itn_num; i++) {
    cout << "iteration: " << i + 1 << endl;

    // draw the candidates of the selection coefficients from the random walk proposal
    U = arma::randn<arma::dcolvec>(6);
    sel_cof_tmp = sel_cof_tmp + S * U;
    sel_cof_chn.slice(i) = arma::reshape(sel_cof_tmp, 3, 2);

    alpha = 0;
    if (arma::any(arma::any(sel_cof_chn.slice(i) < -1, 1))) {
      sel_cof_chn.slice(i) = sel_cof_chn.slice(i - 1);
      log_lik(1) = log_lik(0);
      // apt_cnt = apt_cnt + 0;
      cout << "acceptance: " << apt_cnt / i << endl;
    } else {
      // calculate the proposal
      // arma::drowvec log_psl = arma::zeros<arma::drowvec>(2);

      // calculate the likelihood
      calculateLogLikelihood_arma(log_lik(1), frq_pth, sel_cof_chn.slice(i), rec_rat, pop_siz, ref_siz, smp_gen, smp_siz, ptl_cnt, ptn_num, pcl_num);
      frq_pth_chn.slice(i) = frq_pth;

      // calculate the acceptance ratio
      // double apt_rto = exp(log_lik(1) - log_lik(0));
      // double apt_rto = exp((log_pri(1) + log_lik(1) + log_psl(1)) - (log_pri(0) + log_lik(0) + log_psl(0)));

      alpha = arma::find_finite(log_lik).is_empty() ? 0 : exp(log_lik(1) - log_lik(0));
      alpha = (alpha > 1) ? 1 : alpha;
      if (arma::randu() > alpha) {
        sel_cof_chn.slice(i) = sel_cof_chn.slice(i - 1);
        log_lik(1) = log_lik(0);
        // apt_cnt = apt_cnt + 0;
        cout << "acceptance: " << apt_cnt / i << endl;
      } else {
        log_lik(0) = log_lik(1);
        apt_cnt = apt_cnt + 1;
        cout << "acceptance: " << apt_cnt / i << endl;
      }
    }
    sel_cof_tmp = arma::vectorise(sel_cof_chn.slice(i));

    U = arma::normalise(U);
    M = S * (I + stp_siz(i) * (alpha - apt_rto) * (U * U.t())) * S.t();
    S = arma::chol(M, "lower");
    cout << M << endl;
  }

  return List::create(Named("sel_cof_chn", sel_cof_chn),
                      Named("frq_pth_chn", frq_pth_chn));
}
/*************************/
