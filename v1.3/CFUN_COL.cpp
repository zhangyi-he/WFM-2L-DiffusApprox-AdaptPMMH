// Inferring natural selection acting on horse coat colours and patterns during the process of domestication from ancient DNA data
// Xiaoyang Dai, Sile Hu, Mark Beaumont, Feng Yu, Zhangyi He

// version 1.3
// Horse coat colours (ASIP & MC1R) under non-constant natural selection and non-constant demographic histories (N/A is allowed)

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
    sigma = pow(siz_rto(t - 1), -0.5) * sigma;

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
// Calculate the possible genotype counts in the sample due to missingness
// [[Rcpp::export]]
arma::imat calculateGenoCnt_Mis_arma(const arma::icolvec& smp_cnt, const arma::icolvec& mis_cnt) {
  // ensure RNG gets set/reset
  RNGScope scope;

  arma::imat ptl_smp_cnt(9, 1);
  ptl_smp_cnt.col(0) = smp_cnt;

  arma::icolvec ptl_mis_cnt(9);

  // A1A1/?B1
  if (mis_cnt(0) > 0) {
    arma::uword cnt = 0;
    arma::imat ptl_smp_cnt_1101 = ptl_smp_cnt;
    for (arma::uword k = 0; k <= mis_cnt(0); k++) {
      ptl_mis_cnt.zeros();
      ptl_mis_cnt(0) = k;
      ptl_mis_cnt(1) = mis_cnt(0) - k;
      arma::imat ptl_smp_cnt_tmp = ptl_smp_cnt_1101;
      ptl_smp_cnt_tmp.each_col() += ptl_mis_cnt;
      if (cnt > 0) {
        ptl_smp_cnt.insert_cols(0, ptl_smp_cnt_tmp);
      } else {
        ptl_smp_cnt = ptl_smp_cnt_tmp;
      }
      cnt += 1;
    }
  }

  // A1A1/?B2
  if (mis_cnt(1) > 0) {
    arma::uword cnt = 0;
    arma::imat ptl_smp_cnt_1102 = ptl_smp_cnt;
    for (arma::uword k = 0; k <= mis_cnt(1); k++) {
      ptl_mis_cnt.zeros();
      ptl_mis_cnt(1) = k;
      ptl_mis_cnt(2) = mis_cnt(1) - k;
      arma::imat ptl_smp_cnt_tmp = ptl_smp_cnt_1102;
      ptl_smp_cnt_tmp.each_col() += ptl_mis_cnt;
      if (cnt > 0) {
        ptl_smp_cnt.insert_cols(0, ptl_smp_cnt_tmp);
      } else {
        ptl_smp_cnt = ptl_smp_cnt_tmp;
      }
      cnt += 1;
    }
  }

  // A1A2/?B1
  if (mis_cnt(2) > 0) {
    arma::uword cnt = 0;
    arma::imat ptl_smp_cnt_1201 = ptl_smp_cnt;
    for (arma::uword k = 0; k <= mis_cnt(2); k++) {
      ptl_mis_cnt.zeros();
      ptl_mis_cnt(3) = k;
      ptl_mis_cnt(4) = mis_cnt(2) - k;
      arma::imat ptl_smp_cnt_tmp = ptl_smp_cnt_1201;
      ptl_smp_cnt_tmp.each_col() += ptl_mis_cnt;
      if (cnt > 0) {
        ptl_smp_cnt.insert_cols(0, ptl_smp_cnt_tmp);
      } else {
        ptl_smp_cnt = ptl_smp_cnt_tmp;
      }
      cnt += 1;
    }
  }

  // A1A2/?B2
  if (mis_cnt(3) > 0) {
    arma::uword cnt = 0;
    arma::imat ptl_smp_cnt_1202 = ptl_smp_cnt;
    for (arma::uword k = 0; k <= mis_cnt(3); k++) {
      ptl_mis_cnt.zeros();
      ptl_mis_cnt(4) = k;
      ptl_mis_cnt(6) = mis_cnt(3) - k;
      arma::imat ptl_smp_cnt_tmp = ptl_smp_cnt_1202;
      ptl_smp_cnt_tmp.each_col() += ptl_mis_cnt;
      if (cnt > 0) {
        ptl_smp_cnt.insert_cols(0, ptl_smp_cnt_tmp);
      } else {
        ptl_smp_cnt = ptl_smp_cnt_tmp;
      }
      cnt += 1;
    }
  }

  // A2A2/?B1
  if (mis_cnt(4) > 0) {
    arma::uword cnt = 0;
    arma::imat ptl_smp_cnt_2201 = ptl_smp_cnt;
    for (arma::uword k = 0; k <= mis_cnt(4); k++) {
      ptl_mis_cnt.zeros();
      ptl_mis_cnt(5) = k;
      ptl_mis_cnt(7) = mis_cnt(4) - k;
      arma::imat ptl_smp_cnt_tmp = ptl_smp_cnt_2201;
      ptl_smp_cnt_tmp.each_col() += ptl_mis_cnt;
      if (cnt > 0) {
        ptl_smp_cnt.insert_cols(0, ptl_smp_cnt_tmp);
      } else {
        ptl_smp_cnt = ptl_smp_cnt_tmp;
      }
      cnt += 1;
    }
  }

  // A2A2/?B2
  if (mis_cnt(5) > 0) {
    arma::uword cnt = 0;
    arma::imat ptl_smp_cnt_2202 = ptl_smp_cnt;
    for (arma::uword k = 0; k <= mis_cnt(5); k++) {
      ptl_mis_cnt.zeros();
      ptl_mis_cnt(7) = k;
      ptl_mis_cnt(8) = mis_cnt(5) - k;
      arma::imat ptl_smp_cnt_tmp = ptl_smp_cnt_2202;
      ptl_smp_cnt_tmp.each_col() += ptl_mis_cnt;
      if (cnt > 0) {
        ptl_smp_cnt.insert_cols(0, ptl_smp_cnt_tmp);
      } else {
        ptl_smp_cnt = ptl_smp_cnt_tmp;
      }
      cnt += 1;
    }
  }

  // ?A1/B1B1
  if (mis_cnt(6) > 0) {
    arma::uword cnt = 0;
    arma::imat ptl_smp_cnt_0111 = ptl_smp_cnt;
    for (arma::uword k = 0; k <= mis_cnt(6); k++) {
      ptl_mis_cnt.zeros();
      ptl_mis_cnt(0) = k;
      ptl_mis_cnt(3) = mis_cnt(6) - k;
      arma::imat ptl_smp_cnt_tmp = ptl_smp_cnt_0111;
      ptl_smp_cnt_tmp.each_col() += ptl_mis_cnt;
      if (cnt > 0) {
        ptl_smp_cnt.insert_cols(0, ptl_smp_cnt_tmp);
      } else {
        ptl_smp_cnt = ptl_smp_cnt_tmp;
      }
      cnt += 1;
    }
  }

  // ?A2/B1B1
  if (mis_cnt(7) > 0) {
    arma::uword cnt = 0;
    arma::imat ptl_smp_cnt_0211 = ptl_smp_cnt;
    for (arma::uword k = 0; k <= mis_cnt(7); k++) {
      ptl_mis_cnt.zeros();
      ptl_mis_cnt(3) = k;
      ptl_mis_cnt(5) = mis_cnt(7) - k;
      arma::imat ptl_smp_cnt_tmp = ptl_smp_cnt_0211;
      ptl_smp_cnt_tmp.each_col() += ptl_mis_cnt;
      if (cnt > 0) {
        ptl_smp_cnt.insert_cols(0, ptl_smp_cnt_tmp);
      } else {
        ptl_smp_cnt = ptl_smp_cnt_tmp;
      }
      cnt += 1;
    }
  }

  // ?A1/B1B2
  if (mis_cnt(8) > 0) {
    arma::uword cnt = 0;
    arma::imat ptl_smp_cnt_0112 = ptl_smp_cnt;
    for (arma::uword k = 0; k <= mis_cnt(8); k++) {
      ptl_mis_cnt.zeros();
      ptl_mis_cnt(1) = k;
      ptl_mis_cnt(4) = mis_cnt(8) - k;
      arma::imat ptl_smp_cnt_tmp = ptl_smp_cnt_0112;
      ptl_smp_cnt_tmp.each_col() += ptl_mis_cnt;
      if (cnt > 0) {
        ptl_smp_cnt.insert_cols(0, ptl_smp_cnt_tmp);
      } else {
        ptl_smp_cnt = ptl_smp_cnt_tmp;
      }
      cnt += 1;
    }
  }

  // ?A2/B1B2
  if (mis_cnt(9) > 0) {
    arma::uword cnt = 0;
    arma::imat ptl_smp_cnt_0212 = ptl_smp_cnt;
    for (arma::uword k = 0; k <= mis_cnt(9); k++) {
      ptl_mis_cnt.zeros();
      ptl_mis_cnt(4) = k;
      ptl_mis_cnt(7) = mis_cnt(9) - k;
      arma::imat ptl_smp_cnt_tmp = ptl_smp_cnt_0212;
      ptl_smp_cnt_tmp.each_col() += ptl_mis_cnt;
      if (cnt > 0) {
        ptl_smp_cnt.insert_cols(0, ptl_smp_cnt_tmp);
      } else {
        ptl_smp_cnt = ptl_smp_cnt_tmp;
      }
      cnt += 1;
    }
  }

  // ?A1/B2B2
  if (mis_cnt(10) > 0) {
    arma::uword cnt = 0;
    arma::imat ptl_smp_cnt_0122 = ptl_smp_cnt;
    for (arma::uword k = 0; k <= mis_cnt(10); k++) {
      ptl_mis_cnt.zeros();
      ptl_mis_cnt(2) = k;
      ptl_mis_cnt(6) = mis_cnt(10) - k;
      arma::imat ptl_smp_cnt_tmp = ptl_smp_cnt_0122;
      ptl_smp_cnt_tmp.each_col() += ptl_mis_cnt;
      if (cnt > 0) {
        ptl_smp_cnt.insert_cols(0, ptl_smp_cnt_tmp);
      } else {
        ptl_smp_cnt = ptl_smp_cnt_tmp;
      }
      cnt += 1;
    }
  }

  // ?A2//B2B2
  if (mis_cnt(11) > 0) {
    arma::uword cnt = 0;
    arma::imat ptl_smp_cnt_0222 = ptl_smp_cnt;
    for (arma::uword k = 0; k <= mis_cnt(11); k++) {
      ptl_mis_cnt.zeros();
      ptl_mis_cnt(6) = k;
      ptl_mis_cnt(8) = mis_cnt(11) - k;
      arma::imat ptl_smp_cnt_tmp = ptl_smp_cnt_0222;
      ptl_smp_cnt_tmp.each_col() += ptl_mis_cnt;
      if (cnt > 0) {
        ptl_smp_cnt.insert_cols(0, ptl_smp_cnt_tmp);
      } else {
        ptl_smp_cnt = ptl_smp_cnt_tmp;
      }
      cnt += 1;
    }
  }

  // A1A1/??
  if (mis_cnt(12) > 0) {
    arma::uword cnt = 0;
    arma::imat ptl_smp_cnt_1100 = ptl_smp_cnt;
    for (arma::sword i = 0; i <= mis_cnt(12); i++) {
      for (arma::sword j = 0; j <= mis_cnt(12) - i; j++) {
        ptl_mis_cnt.zeros();
        ptl_mis_cnt(0) = i;
        ptl_mis_cnt(1) = j;
        ptl_mis_cnt(2) = mis_cnt(12) - i - j;
        arma::imat ptl_smp_cnt_tmp = ptl_smp_cnt_1100;
        ptl_smp_cnt_tmp.each_col() += ptl_mis_cnt;
        if (cnt > 0) {
          ptl_smp_cnt.insert_cols(0, ptl_smp_cnt_tmp);
        } else {
          ptl_smp_cnt = ptl_smp_cnt_tmp;
        }
        cnt += 1;
      }
    }
  }

  // A1A2/??
  if (mis_cnt(13) > 0) {
    arma::uword cnt = 0;
    arma::imat ptl_smp_cnt_1200 = ptl_smp_cnt;
    for (arma::sword i = 0; i <= mis_cnt(13); i++) {
      for (arma::sword j = 0; j <= mis_cnt(13) - i; j++) {
        ptl_mis_cnt.zeros();
        ptl_mis_cnt(3) = i;
        ptl_mis_cnt(4) = j;
        ptl_mis_cnt(6) = mis_cnt(13) - i - j;
        arma::imat ptl_smp_cnt_tmp = ptl_smp_cnt_1200;
        ptl_smp_cnt_tmp.each_col() += ptl_mis_cnt;
        if (cnt > 0) {
          ptl_smp_cnt.insert_cols(0, ptl_smp_cnt_tmp);
        } else {
          ptl_smp_cnt = ptl_smp_cnt_tmp;
        }
        cnt += 1;
      }
    }
  }

  // A2A2/??
  if (mis_cnt(14) > 0) {
    arma::uword cnt = 0;
    arma::imat ptl_smp_cnt_2200 = ptl_smp_cnt;
    for (arma::sword i = 0; i <= mis_cnt(14); i++) {
      for (arma::sword j = 0; j <= mis_cnt(14) - i; j++) {
        ptl_mis_cnt.zeros();
        ptl_mis_cnt(5) = i;
        ptl_mis_cnt(7) = j;
        ptl_mis_cnt(8) = mis_cnt(14) - i - j;
        arma::imat ptl_smp_cnt_tmp = ptl_smp_cnt_2200;
        ptl_smp_cnt_tmp.each_col() += ptl_mis_cnt;
        if (cnt > 0) {
          ptl_smp_cnt.insert_cols(0, ptl_smp_cnt_tmp);
        } else {
          ptl_smp_cnt = ptl_smp_cnt_tmp;
        }
        cnt += 1;
      }
    }
  }

  // ??/B1B1
  if (mis_cnt(15) > 0) {
    arma::uword cnt = 0;
    arma::imat ptl_smp_cnt_0011 = ptl_smp_cnt;
    for (arma::sword i = 0; i <= mis_cnt(15); i++) {
      for (arma::sword j = 0; j <= mis_cnt(15) - i; j++) {
        ptl_mis_cnt.zeros();
        ptl_mis_cnt(0) = i;
        ptl_mis_cnt(3) = j;
        ptl_mis_cnt(5) = mis_cnt(15) - i - j;
        arma::imat ptl_smp_cnt_tmp = ptl_smp_cnt_0011;
        ptl_smp_cnt_tmp.each_col() += ptl_mis_cnt;
        if (cnt > 0) {
          ptl_smp_cnt.insert_cols(0, ptl_smp_cnt_tmp);
        } else {
          ptl_smp_cnt = ptl_smp_cnt_tmp;
        }
        cnt += 1;
      }
    }
  }

  // ??/B1B2
  if (mis_cnt(16) > 0) {
    arma::uword cnt = 0;
    arma::imat ptl_smp_cnt_0012 = ptl_smp_cnt;
    for (arma::sword i = 0; i <= mis_cnt(16); i++) {
      for (arma::sword j = 0; j <= mis_cnt(16) - i; j++) {
        ptl_mis_cnt.zeros();
        ptl_mis_cnt(1) = i;
        ptl_mis_cnt(4) = j;
        ptl_mis_cnt(7) = mis_cnt(16) - i - j;
        arma::imat ptl_smp_cnt_tmp = ptl_smp_cnt_0012;
        ptl_smp_cnt_tmp.each_col() += ptl_mis_cnt;
        if (cnt > 0) {
          ptl_smp_cnt.insert_cols(0, ptl_smp_cnt_tmp);
        } else {
          ptl_smp_cnt = ptl_smp_cnt_tmp;
        }
        cnt += 1;
      }
    }
  }

  // ??/B2B2
  if (mis_cnt(17) > 0) {
    arma::uword cnt = 0;
    arma::imat ptl_smp_cnt_0022 = ptl_smp_cnt;
    for (arma::sword i = 0; i <= mis_cnt(17); i++) {
      for (arma::sword j = 0; j <= mis_cnt(17) - i; j++) {
        ptl_mis_cnt.zeros();
        ptl_mis_cnt(2) = i;
        ptl_mis_cnt(6) = j;
        ptl_mis_cnt(8) = mis_cnt(17) - i - j;
        arma::imat ptl_smp_cnt_tmp = ptl_smp_cnt_0022;
        ptl_smp_cnt_tmp.each_col() += ptl_mis_cnt;
        if (cnt > 0) {
          ptl_smp_cnt.insert_cols(0, ptl_smp_cnt_tmp);
        } else {
          ptl_smp_cnt = ptl_smp_cnt_tmp;
        }
        cnt += 1;
      }
    }
  }

  // ?A1/?B1
  if (mis_cnt(18) > 0) {
    arma::uword cnt = 0;
    arma::imat ptl_smp_cnt_0101 = ptl_smp_cnt;
    for (arma::sword i = 0; i <= mis_cnt(18); i++) {
      for (arma::sword j = 0; j <= mis_cnt(18) - i; j++) {
        for (arma::sword k = 0; k <= mis_cnt(18) - i - j; k++) {
          ptl_mis_cnt.zeros();
          ptl_mis_cnt(0) = i;
          ptl_mis_cnt(1) = j;
          ptl_mis_cnt(3) = k;
          ptl_mis_cnt(4) = mis_cnt(18) - i - j - k;
          arma::imat ptl_smp_cnt_tmp = ptl_smp_cnt_0101;
          ptl_smp_cnt_tmp.each_col() += ptl_mis_cnt;
          if (cnt > 0) {
            ptl_smp_cnt.insert_cols(0, ptl_smp_cnt_tmp);
          } else {
            ptl_smp_cnt = ptl_smp_cnt_tmp;
          }
          cnt += 1;
        }
      }
    }
  }

  // ?A1/?B2
  if (mis_cnt(19) > 0) {
    arma::uword cnt = 0;
    arma::imat ptl_smp_cnt_0102 = ptl_smp_cnt;
    for (arma::sword i = 0; i <= mis_cnt(19); i++) {
      for (arma::sword j = 0; j <= mis_cnt(19) - i; j++) {
        for (arma::sword k = 0; k <= mis_cnt(19) - i - j; k++) {
          ptl_mis_cnt.zeros();
          ptl_mis_cnt(1) = i;
          ptl_mis_cnt(2) = j;
          ptl_mis_cnt(4) = k;
          ptl_mis_cnt(6) = mis_cnt(19) - i - j - k;
          arma::imat ptl_smp_cnt_tmp = ptl_smp_cnt_0102;
          ptl_smp_cnt_tmp.each_col() += ptl_mis_cnt;
          if (cnt > 0) {
            ptl_smp_cnt.insert_cols(0, ptl_smp_cnt_tmp);
          } else {
            ptl_smp_cnt = ptl_smp_cnt_tmp;
          }
          cnt += 1;
        }
      }
    }
  }

  // ?A2/?B1
  if (mis_cnt(20) > 0) {
    arma::uword cnt = 0;
    arma::imat ptl_smp_cnt_0201 = ptl_smp_cnt;
    for (arma::sword i = 0; i <= mis_cnt(20); i++) {
      for (arma::sword j = 0; j <= mis_cnt(20) - i; j++) {
        for (arma::sword k = 0; k <= mis_cnt(20) - i - j; k++) {
          ptl_mis_cnt.zeros();
          ptl_mis_cnt(3) = i;
          ptl_mis_cnt(4) = j;
          ptl_mis_cnt(5) = k;
          ptl_mis_cnt(7) = mis_cnt(20) - i - j - k;
          arma::imat ptl_smp_cnt_tmp = ptl_smp_cnt_0201;
          ptl_smp_cnt_tmp.each_col() += ptl_mis_cnt;
          if (cnt > 0) {
            ptl_smp_cnt.insert_cols(0, ptl_smp_cnt_tmp);
          } else {
            ptl_smp_cnt = ptl_smp_cnt_tmp;
          }
          cnt += 1;
        }
      }
    }
  }

  // ?A2/?B2
  if (mis_cnt(21) > 0) {
    arma::uword cnt = 0;
    arma::imat ptl_smp_cnt_0202 = ptl_smp_cnt;
    for (arma::sword i = 0; i <= mis_cnt(21); i++) {
      for (arma::sword j = 0; j <= mis_cnt(21) - i; j++) {
        for (arma::sword k = 0; k <= mis_cnt(21) - i - j; k++) {
          ptl_mis_cnt.zeros();
          ptl_mis_cnt(4) = i;
          ptl_mis_cnt(6) = j;
          ptl_mis_cnt(7) = k;
          ptl_mis_cnt(8) = mis_cnt(21) - i - j - k;
          arma::imat ptl_smp_cnt_tmp = ptl_smp_cnt_0202;
          ptl_smp_cnt_tmp.each_col() += ptl_mis_cnt;
          if (cnt > 0) {
            ptl_smp_cnt.insert_cols(0, ptl_smp_cnt_tmp);
          } else {
            ptl_smp_cnt = ptl_smp_cnt_tmp;
          }
          cnt += 1;
        }
      }
    }
  }

  // ?A1/??
  if (mis_cnt(22) > 0) {
    arma::uword cnt = 0;
    arma::imat ptl_smp_cnt_0100 = ptl_smp_cnt;
    for (arma::sword i = 0; i <= mis_cnt(22); i++) {
      for (arma::sword j = 0; j <= mis_cnt(22) - i; j++) {
        for (arma::sword k = 0; k <= mis_cnt(22) - i - j; k++) {
          for (arma::sword l = 0; l <= mis_cnt(22) - i - j - k; l++) {
            for (arma::sword m = 0; m <= mis_cnt(22) - i - j - k - l; m++) {
              ptl_mis_cnt.zeros();
              ptl_mis_cnt(0) = i;
              ptl_mis_cnt(1) = j;
              ptl_mis_cnt(2) = k;
              ptl_mis_cnt(3) = l;
              ptl_mis_cnt(4) = m;
              ptl_mis_cnt(6) = mis_cnt(22) - i - j - k - l - m;
              arma::imat ptl_smp_cnt_tmp = ptl_smp_cnt_0100;
              ptl_smp_cnt_tmp.each_col() += ptl_mis_cnt;
              if (cnt > 0) {
                ptl_smp_cnt.insert_cols(0, ptl_smp_cnt_tmp);
              } else {
                ptl_smp_cnt = ptl_smp_cnt_tmp;
              }
              cnt += 1;
            }
          }
        }
      }
    }
  }

  // ?A2/??
  if (mis_cnt(23) > 0) {
    arma::uword cnt = 0;
    arma::imat ptl_smp_cnt_0200 = ptl_smp_cnt;
    for (arma::sword i = 0; i <= mis_cnt(23); i++) {
      for (arma::sword j = 0; j <= mis_cnt(23) - i; j++) {
        for (arma::sword k = 0; k <= mis_cnt(23) - i - j; k++) {
          for (arma::sword l = 0; l <= mis_cnt(23) - i - j - k; l++) {
            for (arma::sword m = 0; m <= mis_cnt(23) - i - j - k - l; m++) {
              ptl_mis_cnt.zeros();
              ptl_mis_cnt(3) = i;
              ptl_mis_cnt(4) = j;
              ptl_mis_cnt(5) = k;
              ptl_mis_cnt(6) = l;
              ptl_mis_cnt(7) = m;
              ptl_mis_cnt(8) = mis_cnt(23) - i - j - k - l - m;
              arma::imat ptl_smp_cnt_tmp = ptl_smp_cnt_0200;
              ptl_smp_cnt_tmp.each_col() += ptl_mis_cnt;
              if (cnt > 0) {
                ptl_smp_cnt.insert_cols(0, ptl_smp_cnt_tmp);
              } else {
                ptl_smp_cnt = ptl_smp_cnt_tmp;
              }
              cnt += 1;
            }
          }
        }
      }
    }
  }

  // ??/?B1
  if (mis_cnt(24) > 0) {
    arma::uword cnt = 0;
    arma::imat ptl_smp_cnt_0001 = ptl_smp_cnt;
    for (arma::sword i = 0; i <= mis_cnt(24); i++) {
      for (arma::sword j = 0; j <= mis_cnt(24) - i; j++) {
        for (arma::sword k = 0; k <= mis_cnt(24) - i - j; k++) {
          for (arma::sword l = 0; l <= mis_cnt(24) - i - j - k; l++) {
            for (arma::sword m = 0; m <= mis_cnt(24) - i - j - k - l; m++) {
              ptl_mis_cnt.zeros();
              ptl_mis_cnt(0) = i;
              ptl_mis_cnt(1) = j;
              ptl_mis_cnt(3) = k;
              ptl_mis_cnt(4) = l;
              ptl_mis_cnt(5)= m;
              ptl_mis_cnt(7) = mis_cnt(24) - i - j - k - l - m;
              arma::imat ptl_smp_cnt_tmp = ptl_smp_cnt_0001;
              ptl_smp_cnt_tmp.each_col() += ptl_mis_cnt;
              if (cnt > 0) {
                ptl_smp_cnt.insert_cols(0, ptl_smp_cnt_tmp);
              } else {
                ptl_smp_cnt = ptl_smp_cnt_tmp;
              }
              cnt += 1;
            }
          }
        }
      }
    }
  }

  // ??/?B2
  if (mis_cnt(25) > 0) {
    arma::uword cnt = 0;
    arma::imat ptl_smp_cnt_0002 = ptl_smp_cnt;
    for (arma::sword i = 0; i <= mis_cnt(25); i++) {
      for (arma::sword j = 0; j <= mis_cnt(25) - i; j++) {
        for (arma::sword k = 0; k <= mis_cnt(25) - i - j; k++) {
          for (arma::sword l = 0; l <= mis_cnt(25) - i - j - k; l++) {
            for (arma::sword m = 0; m <= mis_cnt(25) - i - j - k - l; m++) {
              ptl_mis_cnt.zeros();
              ptl_mis_cnt(1) = i;
              ptl_mis_cnt(2) = j;
              ptl_mis_cnt(4) = k;
              ptl_mis_cnt(6) = l;
              ptl_mis_cnt(7) = m;
              ptl_mis_cnt(8) = mis_cnt(25) - i - j - k - l - m;
              arma::imat ptl_smp_cnt_tmp = ptl_smp_cnt_0002;
              ptl_smp_cnt_tmp.each_col() += ptl_mis_cnt;
              if (cnt > 0) {
                ptl_smp_cnt.insert_cols(0, ptl_smp_cnt_tmp);
              } else {
                ptl_smp_cnt = ptl_smp_cnt_tmp;
              }
              cnt += 1;
            }
          }
        }
      }
    }
  }

  // ??/??
  if (mis_cnt(26) > 0) {
    arma::uword cnt = 0;
    arma::imat ptl_smp_cnt_0000 = ptl_smp_cnt;
    for (arma::sword i = 0; i <= mis_cnt(26); i++) {
      for (arma::sword j = 0; j <= mis_cnt(26) - i; j++) {
        for (arma::sword k = 0; k <= mis_cnt(26) - i - j; k++) {
          for (arma::sword l = 0; l <= mis_cnt(26) - i - j - k; l++) {
            for (arma::sword m = 0; m <= mis_cnt(26) - i - j - k - l; m++) {
              for (arma::sword n = 0; n <= mis_cnt(26) - i - j - k - l - m; n++) {
                for (arma::sword o = 0; o <= mis_cnt(26) - i - j - k - l - m - n; o++) {
                  for (arma::sword p = 0; p <= mis_cnt(26) - i - j - k - l - m - n - o; p++) {
                    ptl_mis_cnt.zeros();
                    ptl_mis_cnt(0) = i;
                    ptl_mis_cnt(1) = j;
                    ptl_mis_cnt(2) = k;
                    ptl_mis_cnt(3) = l;
                    ptl_mis_cnt(4) = m;
                    ptl_mis_cnt(5) = n;
                    ptl_mis_cnt(6) = o;
                    ptl_mis_cnt(7) = p;
                    ptl_mis_cnt(8) = mis_cnt(26) - i - j - k - l - m - n - o - p;
                    arma::imat ptl_smp_cnt_tmp = ptl_smp_cnt_0000;
                    ptl_smp_cnt_tmp.each_col() += ptl_mis_cnt;
                    if (cnt > 0) {
                      ptl_smp_cnt.insert_cols(0, ptl_smp_cnt_tmp);
                    } else {
                      ptl_smp_cnt = ptl_smp_cnt_tmp;
                    }
                    cnt += 1;
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  return ptl_smp_cnt;
}

// Calculate the possible genotype counts in the sample due to heterozygosity
// [[Rcpp::export]]
arma::imat calculateGenoCnt_Het_arma(const arma::icolvec& smp_cnt) {
  // ensure RNG gets set/reset
  RNGScope scope;

  if (smp_cnt.n_elem == 9) {
    arma::imat gen_cnt = arma::zeros<arma::imat>(10, smp_cnt(5) + 1);

    for (int j = 0; j <= smp_cnt(5); j++) {
      gen_cnt(0, j) = smp_cnt(0);
      gen_cnt(1, j) = smp_cnt(1);
      gen_cnt(2, j) = smp_cnt(2);
      gen_cnt(3, j) = smp_cnt(3);
      gen_cnt(4, j) = smp_cnt(4);
      gen_cnt(5, j) = j;
      gen_cnt(6, j) = smp_cnt(6);
      gen_cnt(7, j) = smp_cnt(5) - j;
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

// Calculate the possible genotype counts in the sample
// [[Rcpp::export]]
arma::imat calculateGenoCnt_arma(const arma::icolvec& smp_cnt, const arma::icolvec& mis_cnt) {
  // ensure RNG gets set/reset
  RNGScope scope;

  arma::imat gen_cnt = calculateGenoCnt_Mis_arma(smp_cnt, mis_cnt);

  arma::imat ptl_cnt = calculateGenoCnt_Het_arma(gen_cnt.col(0));
  for (arma::uword k = 1; k < gen_cnt.n_cols; k++) {
    arma::imat gen_cnt_tmp = calculateGenoCnt_Het_arma(gen_cnt.col(k));
    ptl_cnt.insert_cols(0, gen_cnt_tmp);
  }

  return ptl_cnt;
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
double calculateEmissionProb_arma(const arma::icolvec& smp_cnt, const arma::icolvec& mis_cnt, const int& smp_siz, const arma::dmat& fts_mat, const arma::dcolvec& hap_frq) {
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

  NumericMatrix part(pcl_num, 4);
  for (int j = 0; j < 4; j++) {
    part(_, j) = rgamma(pcl_num, 1.0, 1.0);
  }
  for (int i = 0; i < pcl_num; i++) {
    part(i, _) = part(i, _) / sum(part(i, _));
  }

  return as<arma::dmat>(transpose(part));
}

// Run the bootstrap particle filter
// [[Rcpp::export]]
List runBPF_arma(const arma::dmat& sel_cof, const double& rec_rat, const arma::icolvec& pop_siz, const int& ref_siz, const arma::irowvec& smp_gen, const arma::irowvec& smp_siz, const arma::imat& smp_cnt, const arma::imat& mis_cnt, const arma::uword& ptn_num, const arma::uword& pcl_num) {
  // ensure RNG gets set/reset
  RNGScope scope;

  double lik = 1;

  arma::uword evt_ind = arma::as_scalar(arma::find(smp_siz == 0));

  arma::dmat wght = arma::zeros<arma::dmat>(pcl_num, smp_gen.n_elem);
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
  hap_frq_tmp = initialiseParticle(pcl_num);
  arma::imat gen_cnt = calculateGenoCnt_arma(smp_cnt.col(0), mis_cnt.col(0));
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
  for (arma::uword k = 1; k < evt_ind; k++) {
    cout << "generation: " << smp_gen(k) << endl;
    wght_tmp = arma::zeros<arma::dcolvec>(pcl_num);
    hap_frq_tmp = hap_frq_pst.slice(k - 1);
    arma::imat gen_cnt = calculateGenoCnt_arma(smp_cnt.col(k), mis_cnt.col(k));
    for (arma::uword i = 0; i < pcl_num; i++) {
      arma::dmat path = simulateWFD_arma(sel_cof.col(0), rec_rat, pop_siz.subvec(smp_gen(k - 1), smp_gen(k)), ref_siz, hap_frq_tmp.col(i), smp_gen(k - 1), smp_gen(k), ptn_num);
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

      return List::create(Named("lik", lik),
                    Named("wght", wght),
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
    arma::imat gen_cnt = calculateGenoCnt_arma(smp_cnt.col(k), mis_cnt.col(k));
    for (arma::uword i = 0; i < pcl_num; i++) {
      arma::dmat path = simulateWFD_arma(sel_cof.col(1), rec_rat, pop_siz.subvec(smp_gen(k - 1), smp_gen(k)), ref_siz, hap_frq_tmp.col(i), smp_gen(k - 1), smp_gen(k), ptn_num);
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

  if (smp_gen(evt_ind) != smp_gen(evt_ind - 1)) {
    wght.shed_col(evt_ind);
    hap_frq_pre.shed_slice(evt_ind);
    hap_frq_pst.shed_slice(evt_ind);
    gen_frq_pre.shed_slice(evt_ind);
    gen_frq_pst.shed_slice(evt_ind);
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
double calculateLogLikelihood_arma(const arma::dmat& sel_cof, const double& rec_rat, const arma::icolvec& pop_siz, const int& ref_siz, const arma::irowvec& smp_gen, const arma::irowvec& smp_siz, const arma::field<arma::imat>& ptl_cnt, const arma::uword& ptn_num, const arma::uword& pcl_num) {
  // ensure RNG gets set/reset
  RNGScope scope;

  double log_lik = 0;

  arma::uword evt_ind = arma::as_scalar(arma::find(smp_siz == 0));

  arma::dcolvec wght = arma::zeros<arma::dcolvec>(pcl_num);
  arma::dmat hap_frq_pre = arma::zeros<arma::dmat>(4, pcl_num);
  arma::dmat hap_frq_pst = arma::zeros<arma::dmat>(4, pcl_num);

  // before the event of interest
  arma::dmat fts_mat = calculateFitnessMat_arma(sel_cof.col(0));

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
  for (arma::uword k = 1; k < evt_ind; k++) {
    wght = arma::zeros<arma::dcolvec>(pcl_num);
    arma::imat gen_cnt = ptl_cnt(k);
    for (arma::uword i = 0; i < pcl_num; i++) {
      arma::dmat path = simulateWFD_arma(sel_cof.col(0), rec_rat, pop_siz.subvec(smp_gen(k - 1), smp_gen(k)), ref_siz, hap_frq_pst.col(i), smp_gen(k - 1), smp_gen(k), ptn_num);
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

      return log_lik;
    }
  }

  // after the event of interest
  fts_mat = calculateFitnessMat_arma(sel_cof.col(1));

  if (smp_gen(evt_ind) != smp_gen(evt_ind - 1)) {
    for (arma::uword i = 0; i < pcl_num; i++) {
      arma::dmat path = simulateWFD_arma(sel_cof.col(1), rec_rat, pop_siz.subvec(smp_gen(evt_ind - 1), smp_gen(evt_ind)), ref_siz, hap_frq_pst.col(i), smp_gen(evt_ind - 1), smp_gen(evt_ind), ptn_num);
      hap_frq_pre.col(i) = arma::vectorise(path.tail_cols(1), 0);
    }
    hap_frq_pst = hap_frq_pre;
  }

  for (arma::uword k = evt_ind + 1; k < smp_gen.n_elem; k++) {
    wght = arma::zeros<arma::dcolvec>(pcl_num);
    arma::imat gen_cnt = ptl_cnt(k);
    for (arma::uword i = 0; i < pcl_num; i++) {
      arma::dmat path = simulateWFD_arma(sel_cof.col(1), rec_rat, pop_siz.subvec(smp_gen(k - 1), smp_gen(k)), ref_siz, hap_frq_pst.col(i), smp_gen(k - 1), smp_gen(k), ptn_num);
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

      return log_lik;
    }
  }

  return log_lik;
}

// Calculate the optimal particle number in the particle marginal Metropolis-Hastings
// [[Rcpp::export]]
List calculateOptimalParticleNum_arma(const arma::dmat& sel_cof, const double& rec_rat, const arma::icolvec& pop_siz, const int& ref_siz, const arma::irowvec& smp_gen, const arma::irowvec& smp_siz, const arma::imat& smp_cnt, const arma::imat& mis_cnt, const arma::uword& ptn_num, const arma::uword& pcl_num, const arma::uword& gap_num) {
  // ensure RNG gets set/reset
  RNGScope scope;

  arma::field<arma::imat> ptl_cnt(smp_gen.n_elem);
  for (arma::uword k = 0; k < smp_gen.n_elem; k++) {
    ptl_cnt(k) = calculateGenoCnt_arma(smp_cnt.col(k), mis_cnt.col(k));
  }

  arma::drowvec log_lik(300);
  for (arma::uword i = 0; i < 300; i++) {
    log_lik(i) = calculateLogLikelihood_arma(sel_cof, rec_rat, pop_siz, ref_siz, smp_gen, smp_siz, ptl_cnt, ptn_num, pcl_num);
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
        log_lik(i) = calculateLogLikelihood_arma(sel_cof, rec_rat, pop_siz, ref_siz, smp_gen, smp_siz, ptl_cnt, ptn_num, opt_pcl_num(0));
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
        log_lik(i) = calculateLogLikelihood_arma(sel_cof, rec_rat, pop_siz, ref_siz, smp_gen, smp_siz, ptl_cnt, ptn_num, opt_pcl_num(0));
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
        log_lik(i) = calculateLogLikelihood_arma(sel_cof, rec_rat, pop_siz, ref_siz, smp_gen, smp_siz, ptl_cnt, ptn_num, opt_pcl_num(0));
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
        log_lik(i) = calculateLogLikelihood_arma(sel_cof, rec_rat, pop_siz, ref_siz, smp_gen, smp_siz, ptl_cnt, ptn_num, opt_pcl_num(0));
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
arma::dcube runPMMH_arma(const arma::dmat& sel_cof, const double& rec_rat, const arma::icolvec& pop_siz, const int& ref_siz, const arma::irowvec& smp_gen, const arma::irowvec& smp_siz, const arma::imat& smp_cnt, const arma::imat& mis_cnt, const arma::uword& ptn_num, const arma::uword& pcl_num, const arma::uword& itn_num) {
  // ensure RNG gets set/reset
  RNGScope scope;

  arma::field<arma::imat> ptl_cnt(smp_gen.n_elem);
  for (arma::uword k = 0; k < smp_gen.n_elem; k++) {
    ptl_cnt(k) = calculateGenoCnt_arma(smp_cnt.col(k), mis_cnt.col(k));
  }

  arma::dcube sel_cof_chn = arma::zeros<arma::dcube>(2, 2, itn_num);

  //arma::drowvec log_pri_chn = arma::zeros<arma::drowvec>(itn_num);
  arma::drowvec log_lik_chn = arma::zeros<arma::drowvec>(itn_num);

  arma::dmat sel_cof_sd = {{5e-03, 5e-03},
                           {5e-03, 5e-03}};

  // initialise the population genetic parameters
  cout << "iteration: " << 1 << endl;
  // take the uniform prior and fix the initial value of the selection coefficient to zero
  // or take the beta prior with alpha = 1 and beta = 3
  sel_cof_chn.slice(0) = sel_cof;

  log_lik_chn(0) = calculateLogLikelihood_arma(sel_cof_chn.slice(0), rec_rat, pop_siz, ref_siz, smp_gen, smp_siz, ptl_cnt, ptn_num, pcl_num);

  double apt_rto = 0;
  for (arma::uword i = 1; i < itn_num; i++) {
    cout << "iteration: " << i + 1 << endl;

    // draw the candidates of the selection coefficients from the random walk proposal
    sel_cof_chn.slice(i) = sel_cof_chn.slice(i - 1) + sel_cof_sd % arma::randn<arma::dmat>(2, 2);

    if (arma::any(arma::any(sel_cof_chn.slice(i) < -1, 1))) {
      sel_cof_chn.slice(i) = sel_cof_chn.slice(i - 1);
      log_lik_chn(i) = log_lik_chn(i - 1);
    } else {
      // calculate the proposal
      //double log_psl_old_new
      //double log_psl_new_old

      // calculate the likelihood
      log_lik_chn(i) = calculateLogLikelihood_arma(sel_cof_chn.slice(i), rec_rat, pop_siz, ref_siz, smp_gen, smp_siz, ptl_cnt, ptn_num, pcl_num);

      // calculate the acceptance ratio
      apt_rto = exp(log_lik_chn(i) - log_lik_chn(i - 1));
      //apt_rto = exp((log_pri_chn(i) + log_lik_chn(i) + log_psl_old_new) - (log_pri_chn(i - 1) + log_lik_chn(i - 1) + log_psl_new_old));

      if (arma::randu() > apt_rto) {
        sel_cof_chn.slice(i) = sel_cof_chn.slice(i - 1);
        log_lik_chn(i) = log_lik_chn(i - 1);
      }
    }
  }

  return sel_cof_chn;
}
/*************************/
