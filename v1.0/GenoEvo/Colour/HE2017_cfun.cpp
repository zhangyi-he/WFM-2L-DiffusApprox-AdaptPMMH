// Inferring fluctuating selection acting on horse coat colours and patterns from ancient DNA data
// Zhangyi He, Xiaoyang Dai, Mark Beaumont and Feng Yu

// Genotype evolution (Wright-Fisher diffusion)

// Horse base coat colour: ASIP and MC1R (fluctuating selection caused by domestication)

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

/********** SIM **********/
// Simulate the phased genotype frequency trajectories according to the two-locus Wright-Fisher model with selection
// [[Rcpp::export]]
arma::dmat simulateTLWFMS_arma(const arma::dcolvec& sel_cof, const double& rec_rat, const int& pop_siz, const arma::dcolvec& int_frq, const int& int_gen, const int& lst_gen) {
  // ensure RNG gets set/reset
  RNGScope scope;
  
  // declare the fitness
  arma::dmat fts = arma::ones<arma::dmat>(4, 4);
  fts(1, 1) = 1 + sel_cof(1);
  fts(3, 1) = 1 + sel_cof(1);
  fts(2, 2) = 1 + sel_cof(0);
  fts(3, 2) = 1 + sel_cof(0);
  fts(1, 3) = 1 + sel_cof(1);
  fts(2, 3) = 1 + sel_cof(0);
  fts(3, 3) = 1 + sel_cof(1);
  
  // declare eta
  arma::dcolvec eta(4);
  eta(0) = -1;
  eta(1) = 1;
  eta(2) = 1;
  eta(3) = -1;
  
  // declare the genotype frequency trajectories
  arma::dmat frq_pth(10, arma::uword(lst_gen - int_gen) + 1);
  
  // initialise the genotype frequencies in generation 0
  frq_pth.col(0) = int_frq;
  
  for(arma::uword t = 1; t < arma::uword(lst_gen - int_gen) + 1; t++) {
    // meiosis (genetic recombination)
    arma::dcolvec hap_frq(4);
    hap_frq(0) = frq_pth(0, t - 1) + frq_pth(1, t - 1) / 2 + frq_pth(3, t - 1) / 2 + frq_pth(9, t - 1) / 2;
    hap_frq(1) = frq_pth(2, t - 1) + frq_pth(1, t - 1) / 2 + frq_pth(4, t - 1) / 2 + frq_pth(5, t - 1) / 2;
    hap_frq(2) = frq_pth(6, t - 1) + frq_pth(3, t - 1) / 2 + frq_pth(4, t - 1) / 2 + frq_pth(7, t - 1) / 2;
    hap_frq(3) = frq_pth(8, t - 1) + frq_pth(9, t - 1) / 2 + frq_pth(5, t - 1) / 2 + frq_pth(7, t - 1) / 2;
    hap_frq = hap_frq + eta * rec_rat * (hap_frq(0) * hap_frq(3) - hap_frq(1) * hap_frq(2));
    
    // reproduction (the Wright-Fisher sampling)
    IntegerVector hap_cnt(4);
    R::rmultinom(2 * pop_siz, hap_frq.begin(), 4, hap_cnt.begin());
    hap_frq = as<arma::dcolvec>(hap_cnt) / 2 / pop_siz;
    
    // random union of gametes
    arma::dmat gen_frq = hap_frq * hap_frq.t();
    
    // viability selection
    gen_frq = (fts % gen_frq) / arma::accu(fts % gen_frq);
    frq_pth(0, t) = gen_frq(0, 0);
    frq_pth(1, t) = gen_frq(0, 1) + gen_frq(1, 0);
    frq_pth(2, t) = gen_frq(1, 1);
    frq_pth(3, t) = gen_frq(0, 2) + gen_frq(2, 0);
    frq_pth(4, t) = gen_frq(1, 2) + gen_frq(2, 1);
    frq_pth(5, t) = gen_frq(1, 3) + gen_frq(3, 1);
    frq_pth(6, t) = gen_frq(2, 2);
    frq_pth(7, t) = gen_frq(2, 3) + gen_frq(3, 2);
    frq_pth(8, t) = gen_frq(3, 3);
    frq_pth(9, t) = gen_frq(0, 3) + gen_frq(3, 0);
  }
  
  return frq_pth;
}

// Simulate the phased genotype frequency trajectories according to the two-locus Wright-Fisher diffusion with selection using the Euler-Maruyama method
// [[Rcpp::export]]
arma::dmat simulateTLWFDS_arma(const arma::dcolvec& sel_cof, const double& rec_rat, const int& pop_siz, const arma::dcolvec& int_frq, const int& int_gen, const int& lst_gen, const arma::uword& ptn_num) {
  // ensure RNG gets set/reset
  RNGScope scope;
  
  // declare the fitness
  arma::dmat fts = arma::ones<arma::dmat>(4, 4);
  fts(1, 1) = 1 + sel_cof(1);
  fts(3, 1) = 1 + sel_cof(1);
  fts(2, 2) = 1 + sel_cof(0);
  fts(3, 2) = 1 + sel_cof(0);
  fts(1, 3) = 1 + sel_cof(1);
  fts(2, 3) = 1 + sel_cof(0);
  fts(3, 3) = 1 + sel_cof(1);
  
  // rescale the recombination rate
  double scl_rec_rat = 4 * pop_siz * rec_rat;
  
  // declare eta
  arma::dcolvec eta(4);
  eta(0) = -1;
  eta(1) = 1;
  eta(2) = 1;
  eta(3) = -1;
  
  // declare delta t
  double dt = 1.0 / (2 * pop_siz) / ptn_num;
  
  // declare the genotype frequency trajectories
  arma::dmat frq_pth(10, arma::uword(lst_gen - int_gen) + 1);
  
  // initialise the genotype frequencies at generation 0
  frq_pth.col(0) = int_frq;
  
  for(arma::uword t = 1; t < arma::uword(lst_gen - int_gen) + 1; t++) {
    // declare delta W
    arma::dmat dW = pow(dt, 0.5) * arma::randn<arma::dmat>(6, ptn_num);
    
    arma::dcolvec hap_frq = arma::zeros<arma::dcolvec>(4);
    hap_frq(0) = frq_pth(0, t - 1) + frq_pth(1, t - 1) / 2 + frq_pth(3, t - 1) / 2 + frq_pth(9, t - 1) / 2;
    hap_frq(1) = frq_pth(2, t - 1) + frq_pth(1, t - 1) / 2 + frq_pth(4, t - 1) / 2 + frq_pth(5, t - 1) / 2;
    hap_frq(2) = frq_pth(6, t - 1) + frq_pth(3, t - 1) / 2 + frq_pth(4, t - 1) / 2 + frq_pth(7, t - 1) / 2;
    hap_frq(3) = frq_pth(8, t - 1) + frq_pth(9, t - 1) / 2 + frq_pth(5, t - 1) / 2 + frq_pth(7, t - 1) / 2;
    
    for(arma::uword k = 0; k < ptn_num; k++) {
      // calculate the drift coefficient vector
      arma::dcolvec mu = eta * 0.5 * scl_rec_rat * (hap_frq(0) * hap_frq(3) - hap_frq(1) * hap_frq(2));
      
      // calculate the diffusion coefficient matrix
      arma::dmat sigma = arma::zeros<arma::dmat>(4, 6);
      sigma(0, 0) = pow(hap_frq(0) * hap_frq(1), 0.5);
      sigma(0, 1) = pow(hap_frq(0) * hap_frq(2), 0.5);
      sigma(0, 2) = pow(hap_frq(0) * hap_frq(3), 0.5);
      //sigma(0, 3) = 0;
      //sigma(0, 4) = 0;
      //sigma(0, 5) = 0;
      sigma(1, 0) = -pow(hap_frq(1) * hap_frq(0), 0.5);
      //sigma(1, 1) = 0;
      //sigma(1, 2) = 0;
      sigma(1, 3) = pow(hap_frq(1) * hap_frq(2), 0.5);
      sigma(1, 4) = pow(hap_frq(1) * hap_frq(3), 0.5);
      //sigma(1, 5) = 0;
      //sigma(2, 0) = 0;
      sigma(2, 1) = -pow(hap_frq(2) * hap_frq(0), 0.5);
      //sigma(2, 2) = 0;
      sigma(2, 3) = -pow(hap_frq(2) * hap_frq(1), 0.5);
      //sigma(2, 4) = 0;
      sigma(2, 5) = pow(hap_frq(2) * hap_frq(3), 0.5);
      //sigma(3, 0) = 0;
      //sigma(3, 1) = 0;
      sigma(3, 2) = -pow(hap_frq(3) * hap_frq(0), 0.5);
      //sigma(3, 3) = 0;
      sigma(3, 4) = -pow(hap_frq(3) * hap_frq(1), 0.5);
      sigma(3, 5) = -pow(hap_frq(3) * hap_frq(2), 0.5);
      
      // proceed the Euler-Maruyama scheme
      hap_frq = hap_frq + mu * dt + sigma * dW.col(k);
      
      // remove the noise from the numerical techniques
      for(arma::uword i = 0; i < 4; i++) {
        if(hap_frq(i) < 0) {
          hap_frq(i) = 0;
        } 
        if(hap_frq(i) > 1) {
          hap_frq(i) = 1;
        }
      }
      hap_frq = hap_frq / arma::accu(hap_frq);
    }
    
    // random union of gametes
    arma::dmat gen_frq = hap_frq * hap_frq.t();
    
    // viability selection
    gen_frq = (fts % gen_frq) / arma::accu(fts % gen_frq);
    frq_pth(0, t) = gen_frq(0, 0);
    frq_pth(1, t) = gen_frq(0, 1) + gen_frq(1, 0);
    frq_pth(2, t) = gen_frq(1, 1);
    frq_pth(3, t) = gen_frq(0, 2) + gen_frq(2, 0);
    frq_pth(4, t) = gen_frq(1, 2) + gen_frq(2, 1);
    frq_pth(5, t) = gen_frq(1, 3) + gen_frq(3, 1);
    frq_pth(6, t) = gen_frq(2, 2);
    frq_pth(7, t) = gen_frq(2, 3) + gen_frq(3, 2);
    frq_pth(8, t) = gen_frq(3, 3);
    frq_pth(9, t) = gen_frq(0, 3) + gen_frq(3, 0);
  }
  
  return frq_pth;
}
/*************************/


/********** CVT **********/
// Convert the phased genotype frequency trajectories to the unphased genotype frequency trajectories
// [[Rcpp::export]]
arma::dmat convertPhasedGenoPath2UnphasedGenoPath_arma(const arma::dmat& frq_pth) {
  // ensure RNG gets set/reset
  RNGScope scope;
  
  // declare the phased genotype frequency trajectories
  arma::dmat gen_frq_pth(9, frq_pth.n_cols);
  gen_frq_pth.row(0) = frq_pth.row(0);
  gen_frq_pth.row(1) = frq_pth.row(1);
  gen_frq_pth.row(2) = frq_pth.row(2);
  gen_frq_pth.row(3) = frq_pth.row(3);
  gen_frq_pth.row(4) = frq_pth.row(4) + frq_pth.row(9);
  gen_frq_pth.row(5) = frq_pth.row(5);
  gen_frq_pth.row(6) = frq_pth.row(6);
  gen_frq_pth.row(7) = frq_pth.row(7);
  gen_frq_pth.row(8) = frq_pth.row(8);
  
  return gen_frq_pth;
}

// Convert the phased genotype frequency trajectories to the haplotype frequency trajectories
// [[Rcpp::export]]
arma::dmat convertPhasedGenoPath2HaploPath_arma(const arma::dmat& frq_pth) {
  // ensure RNG gets set/reset
  RNGScope scope;
  
  // declare the haplotype frequency trajectories
  arma::dmat hap_frq_pth(4, frq_pth.n_cols);
  hap_frq_pth.row(0) = frq_pth.row(0) + frq_pth.row(1) / 2 + frq_pth.row(3) / 2 + frq_pth.row(9) / 2;
  hap_frq_pth.row(1) = frq_pth.row(2) + frq_pth.row(1) / 2 + frq_pth.row(4) / 2 + frq_pth.row(5) / 2;
  hap_frq_pth.row(2) = frq_pth.row(6) + frq_pth.row(3) / 2 + frq_pth.row(4) / 2 + frq_pth.row(7) / 2;
  hap_frq_pth.row(3) = frq_pth.row(8) + frq_pth.row(9) / 2 + frq_pth.row(5) / 2 + frq_pth.row(7) / 2;
  
  return hap_frq_pth;
}

// Convert the phased genotype frequency trajectories to the allele frequency trajectories
// [[Rcpp::export]]
arma::dmat convertPhasedGenoPath2AllelePath_arma(const arma::dmat& frq_pth) {
  // ensure RNG gets set/reset
  RNGScope scope;
  
  // declare the allele frequency trajectories
  arma::dmat ale_frq_pth(2, frq_pth.n_cols);
  ale_frq_pth.row(0) = frq_pth.row(6) + frq_pth.row(3) / 2 + frq_pth.row(4) / 2 + frq_pth.row(7) / 2 + frq_pth.row(8) + frq_pth.row(9) / 2 + frq_pth.row(5) / 2 + frq_pth.row(7) / 2;
  ale_frq_pth.row(1) = frq_pth.row(2) + frq_pth.row(1) / 2 + frq_pth.row(4) / 2 + frq_pth.row(5) / 2 + frq_pth.row(8) + frq_pth.row(9) / 2 + frq_pth.row(5) / 2 + frq_pth.row(7) / 2;
  
  return ale_frq_pth;
}
/*************************/


/********** SMP **********/
// Calculate the possible unobserved sample unphased genotype counts from the observed sample unphased genotype counts with uncertainties
// [[Rcpp::export]]
arma::imat calculateUnphasedGenoCount_arma(const arma::icolvec& smp_gen_cnt, const arma::icolvec& mis_gen_cnt) {
  // ensure RNG gets set/reset
  RNGScope scope;
  
  arma::imat ptl_smp_gen_cnt(9, 1);
  ptl_smp_gen_cnt.col(0) = smp_gen_cnt;
  
  arma::icolvec ptl_mis_gen_cnt(9);
  
  // A1A1/?B1
  if (mis_gen_cnt(0) > 0) {
    arma::uword cnt = 0;
    arma::imat ptl_smp_gen_cnt_1101 = ptl_smp_gen_cnt;
    for (arma::uword k = 0; k <= mis_gen_cnt(0); k++) {
      ptl_mis_gen_cnt.zeros();
      ptl_mis_gen_cnt(0) = k;
      ptl_mis_gen_cnt(1) = mis_gen_cnt(0) - k;
      arma::imat ptl_smp_gen_cnt_tmp = ptl_smp_gen_cnt_1101;
      ptl_smp_gen_cnt_tmp.each_col() += ptl_mis_gen_cnt;
      if (cnt > 0) {
        ptl_smp_gen_cnt.insert_cols(0, ptl_smp_gen_cnt_tmp);
      } else {
        ptl_smp_gen_cnt = ptl_smp_gen_cnt_tmp;
      }
      cnt += 1;
    }
  }
  
  // A1A1/?B2
  if (mis_gen_cnt(1) > 0) {
    arma::uword cnt = 0;
    arma::imat ptl_smp_gen_cnt_1102 = ptl_smp_gen_cnt;
    for (arma::uword k = 0; k <= mis_gen_cnt(1); k++) {
      ptl_mis_gen_cnt.zeros();
      ptl_mis_gen_cnt(1) = k;
      ptl_mis_gen_cnt(2) = mis_gen_cnt(1) - k;
      arma::imat ptl_smp_gen_cnt_tmp = ptl_smp_gen_cnt_1102;
      ptl_smp_gen_cnt_tmp.each_col() += ptl_mis_gen_cnt;
      if (cnt > 0) {
        ptl_smp_gen_cnt.insert_cols(0, ptl_smp_gen_cnt_tmp);
      } else {
        ptl_smp_gen_cnt = ptl_smp_gen_cnt_tmp;
      }
      cnt += 1;
    }
  }
  
  // A1A2/?B1
  if (mis_gen_cnt(2) > 0) {
    arma::uword cnt = 0;
    arma::imat ptl_smp_gen_cnt_1201 = ptl_smp_gen_cnt;
    for (arma::uword k = 0; k <= mis_gen_cnt(2); k++) {
      ptl_mis_gen_cnt.zeros();
      ptl_mis_gen_cnt(3) = k;
      ptl_mis_gen_cnt(4) = mis_gen_cnt(2) - k;
      arma::imat ptl_smp_gen_cnt_tmp = ptl_smp_gen_cnt_1201;
      ptl_smp_gen_cnt_tmp.each_col() += ptl_mis_gen_cnt;
      if (cnt > 0) {
        ptl_smp_gen_cnt.insert_cols(0, ptl_smp_gen_cnt_tmp);
      } else {
        ptl_smp_gen_cnt = ptl_smp_gen_cnt_tmp;
      }
      cnt += 1;
    }
  }
  
  // A1A2/?B2
  if (mis_gen_cnt(3) > 0) {
    arma::uword cnt = 0;
    arma::imat ptl_smp_gen_cnt_1202 = ptl_smp_gen_cnt;
    for (arma::uword k = 0; k <= mis_gen_cnt(3); k++) {
      ptl_mis_gen_cnt.zeros();
      ptl_mis_gen_cnt(4) = k;
      ptl_mis_gen_cnt(5) = mis_gen_cnt(3) - k;
      arma::imat ptl_smp_gen_cnt_tmp = ptl_smp_gen_cnt_1202;
      ptl_smp_gen_cnt_tmp.each_col() += ptl_mis_gen_cnt;
      if (cnt > 0) {
        ptl_smp_gen_cnt.insert_cols(0, ptl_smp_gen_cnt_tmp);
      } else {
        ptl_smp_gen_cnt = ptl_smp_gen_cnt_tmp;
      }
      cnt += 1;
    }
  }
  
  // A2A2/?B1
  if (mis_gen_cnt(4) > 0) {
    arma::uword cnt = 0;
    arma::imat ptl_smp_gen_cnt_2201 = ptl_smp_gen_cnt;
    for (arma::uword k = 0; k <= mis_gen_cnt(4); k++) {
      ptl_mis_gen_cnt.zeros();
      ptl_mis_gen_cnt(6) = k;
      ptl_mis_gen_cnt(7) = mis_gen_cnt(4) - k;
      arma::imat ptl_smp_gen_cnt_tmp = ptl_smp_gen_cnt_2201;
      ptl_smp_gen_cnt_tmp.each_col() += ptl_mis_gen_cnt;
      if (cnt > 0) {
        ptl_smp_gen_cnt.insert_cols(0, ptl_smp_gen_cnt_tmp);
      } else {
        ptl_smp_gen_cnt = ptl_smp_gen_cnt_tmp;
      }
      cnt += 1;
    }
  }
  
  // A2A2/?B2
  if (mis_gen_cnt(5) > 0) {
    arma::uword cnt = 0;
    arma::imat ptl_smp_gen_cnt_2202 = ptl_smp_gen_cnt;
    for (arma::uword k = 0; k <= mis_gen_cnt(5); k++) {
      ptl_mis_gen_cnt.zeros();
      ptl_mis_gen_cnt(7) = k;
      ptl_mis_gen_cnt(8) = mis_gen_cnt(5) - k;
      arma::imat ptl_smp_gen_cnt_tmp = ptl_smp_gen_cnt_2202;
      ptl_smp_gen_cnt_tmp.each_col() += ptl_mis_gen_cnt;
      if (cnt > 0) {
        ptl_smp_gen_cnt.insert_cols(0, ptl_smp_gen_cnt_tmp);
      } else {
        ptl_smp_gen_cnt = ptl_smp_gen_cnt_tmp;
      }
      cnt += 1;
    }
  }
  
  // ?A1/B1B1
  if (mis_gen_cnt(6) > 0) {
    arma::uword cnt = 0;
    arma::imat ptl_smp_gen_cnt_0111 = ptl_smp_gen_cnt;
    for (arma::uword k = 0; k <= mis_gen_cnt(6); k++) {
      ptl_mis_gen_cnt.zeros();
      ptl_mis_gen_cnt(0) = k;
      ptl_mis_gen_cnt(3) = mis_gen_cnt(6) - k;
      arma::imat ptl_smp_gen_cnt_tmp = ptl_smp_gen_cnt_0111;
      ptl_smp_gen_cnt_tmp.each_col() += ptl_mis_gen_cnt;
      if (cnt > 0) {
        ptl_smp_gen_cnt.insert_cols(0, ptl_smp_gen_cnt_tmp);
      } else {
        ptl_smp_gen_cnt = ptl_smp_gen_cnt_tmp;
      }
      cnt += 1;
    }
  }
  
  // ?A2/B1B1
  if (mis_gen_cnt(7) > 0) {
    arma::uword cnt = 0;
    arma::imat ptl_smp_gen_cnt_0211 = ptl_smp_gen_cnt;
    for (arma::uword k = 0; k <= mis_gen_cnt(7); k++) {
      ptl_mis_gen_cnt.zeros();
      ptl_mis_gen_cnt(3) = k;
      ptl_mis_gen_cnt(6) = mis_gen_cnt(7) - k;
      arma::imat ptl_smp_gen_cnt_tmp = ptl_smp_gen_cnt_0211;
      ptl_smp_gen_cnt_tmp.each_col() += ptl_mis_gen_cnt;
      if (cnt > 0) {
        ptl_smp_gen_cnt.insert_cols(0, ptl_smp_gen_cnt_tmp);
      } else {
        ptl_smp_gen_cnt = ptl_smp_gen_cnt_tmp;
      }
      cnt += 1;
    }
  }
  
  // ?A1/B1B2
  if (mis_gen_cnt(8) > 0) {
    arma::uword cnt = 0;
    arma::imat ptl_smp_gen_cnt_0112 = ptl_smp_gen_cnt;
    for (arma::uword k = 0; k <= mis_gen_cnt(8); k++) {
      ptl_mis_gen_cnt.zeros();
      ptl_mis_gen_cnt(1) = k;
      ptl_mis_gen_cnt(4) = mis_gen_cnt(8) - k;
      arma::imat ptl_smp_gen_cnt_tmp = ptl_smp_gen_cnt_0112;
      ptl_smp_gen_cnt_tmp.each_col() += ptl_mis_gen_cnt;
      if (cnt > 0) {
        ptl_smp_gen_cnt.insert_cols(0, ptl_smp_gen_cnt_tmp);
      } else {
        ptl_smp_gen_cnt = ptl_smp_gen_cnt_tmp;
      }
      cnt += 1;
    }
  }
  
  // ?A2/B1B2
  if (mis_gen_cnt(9) > 0) {
    arma::uword cnt = 0;
    arma::imat ptl_smp_gen_cnt_0212 = ptl_smp_gen_cnt;
    for (arma::uword k = 0; k <= mis_gen_cnt(9); k++) {
      ptl_mis_gen_cnt.zeros();
      ptl_mis_gen_cnt(4) = k;
      ptl_mis_gen_cnt(7) = mis_gen_cnt(9) - k;
      arma::imat ptl_smp_gen_cnt_tmp = ptl_smp_gen_cnt_0212;
      ptl_smp_gen_cnt_tmp.each_col() += ptl_mis_gen_cnt;
      if (cnt > 0) {
        ptl_smp_gen_cnt.insert_cols(0, ptl_smp_gen_cnt_tmp);
      } else {
        ptl_smp_gen_cnt = ptl_smp_gen_cnt_tmp;
      }
      cnt += 1;
    }
  }
  
  // ?A1/B2B2
  if (mis_gen_cnt(10) > 0) {
    arma::uword cnt = 0;
    arma::imat ptl_smp_gen_cnt_0122 = ptl_smp_gen_cnt;
    for (arma::uword k = 0; k <= mis_gen_cnt(10); k++) {
      ptl_mis_gen_cnt.zeros();
      ptl_mis_gen_cnt(2) = k;
      ptl_mis_gen_cnt(5) = mis_gen_cnt(10) - k;
      arma::imat ptl_smp_gen_cnt_tmp = ptl_smp_gen_cnt_0122;
      ptl_smp_gen_cnt_tmp.each_col() += ptl_mis_gen_cnt;
      if (cnt > 0) {
        ptl_smp_gen_cnt.insert_cols(0, ptl_smp_gen_cnt_tmp);
      } else {
        ptl_smp_gen_cnt = ptl_smp_gen_cnt_tmp;
      }
      cnt += 1;
    }
  }
  
  // ?A2//B2B2
  if (mis_gen_cnt(11) > 0) {
    arma::uword cnt = 0;
    arma::imat ptl_smp_gen_cnt_0222 = ptl_smp_gen_cnt;
    for (arma::uword k = 0; k <= mis_gen_cnt(11); k++) {
      ptl_mis_gen_cnt.zeros();
      ptl_mis_gen_cnt(5) = k;
      ptl_mis_gen_cnt(8) = mis_gen_cnt(11) - k;
      arma::imat ptl_smp_gen_cnt_tmp = ptl_smp_gen_cnt_0222;
      ptl_smp_gen_cnt_tmp.each_col() += ptl_mis_gen_cnt;
      if (cnt > 0) {
        ptl_smp_gen_cnt.insert_cols(0, ptl_smp_gen_cnt_tmp);
      } else {
        ptl_smp_gen_cnt = ptl_smp_gen_cnt_tmp;
      }
      cnt += 1;
    }
  }
  
  // A1A1/??
  if (mis_gen_cnt(12) > 0) {
    arma::uword cnt = 0;
    arma::imat ptl_smp_gen_cnt_1100 = ptl_smp_gen_cnt;
    for (arma::sword i = 0; i <= mis_gen_cnt(12); i++) {
      for (arma::sword j = 0; j <= mis_gen_cnt(12) - i; j++) {
        ptl_mis_gen_cnt.zeros();
        ptl_mis_gen_cnt(0) = i;
        ptl_mis_gen_cnt(1) = j;
        ptl_mis_gen_cnt(2) = mis_gen_cnt(12) - i - j;
        arma::imat ptl_smp_gen_cnt_tmp = ptl_smp_gen_cnt_1100;
        ptl_smp_gen_cnt_tmp.each_col() += ptl_mis_gen_cnt;
        if (cnt > 0) {
          ptl_smp_gen_cnt.insert_cols(0, ptl_smp_gen_cnt_tmp);
        } else {
          ptl_smp_gen_cnt = ptl_smp_gen_cnt_tmp;
        }
        cnt += 1;
      }
    }
  }
  
  // A1A2/??
  if (mis_gen_cnt(13) > 0) {
    arma::uword cnt = 0;
    arma::imat ptl_smp_gen_cnt_1200 = ptl_smp_gen_cnt;
    for (arma::sword i = 0; i <= mis_gen_cnt(13); i++) {
      for (arma::sword j = 0; j <= mis_gen_cnt(13) - i; j++) {
        ptl_mis_gen_cnt.zeros();
        ptl_mis_gen_cnt(3) = i;
        ptl_mis_gen_cnt(4) = j;
        ptl_mis_gen_cnt(5) = mis_gen_cnt(13) - i - j;
        arma::imat ptl_smp_gen_cnt_tmp = ptl_smp_gen_cnt_1200;
        ptl_smp_gen_cnt_tmp.each_col() += ptl_mis_gen_cnt;
        if (cnt > 0) {
          ptl_smp_gen_cnt.insert_cols(0, ptl_smp_gen_cnt_tmp);
        } else {
          ptl_smp_gen_cnt = ptl_smp_gen_cnt_tmp;
        }
        cnt += 1;
      }
    }
  }
  
  // A2A2/??
  if (mis_gen_cnt(14) > 0) {
    arma::uword cnt = 0;
    arma::imat ptl_smp_gen_cnt_2200 = ptl_smp_gen_cnt;
    for (arma::sword i = 0; i <= mis_gen_cnt(14); i++) {
      for (arma::sword j = 0; j <= mis_gen_cnt(14) - i; j++) {
        ptl_mis_gen_cnt.zeros();
        ptl_mis_gen_cnt(6) = i;
        ptl_mis_gen_cnt(7) = j;
        ptl_mis_gen_cnt(8) = mis_gen_cnt(14) - i - j;
        arma::imat ptl_smp_gen_cnt_tmp = ptl_smp_gen_cnt_2200;
        ptl_smp_gen_cnt_tmp.each_col() += ptl_mis_gen_cnt;
        if (cnt > 0) {
          ptl_smp_gen_cnt.insert_cols(0, ptl_smp_gen_cnt_tmp);
        } else {
          ptl_smp_gen_cnt = ptl_smp_gen_cnt_tmp;
        }
        cnt += 1;
      }
    }
  }
  
  // ??/B1B1
  if (mis_gen_cnt(15) > 0) {
    arma::uword cnt = 0;
    arma::imat ptl_smp_gen_cnt_0011 = ptl_smp_gen_cnt;
    for (arma::sword i = 0; i <= mis_gen_cnt(15); i++) {
      for (arma::sword j = 0; j <= mis_gen_cnt(15) - i; j++) {
        ptl_mis_gen_cnt.zeros();
        ptl_mis_gen_cnt(0) = i;
        ptl_mis_gen_cnt(3) = j;
        ptl_mis_gen_cnt(6) = mis_gen_cnt(15) - i - j;
        arma::imat ptl_smp_gen_cnt_tmp = ptl_smp_gen_cnt_0011;
        ptl_smp_gen_cnt_tmp.each_col() += ptl_mis_gen_cnt;
        if (cnt > 0) {
          ptl_smp_gen_cnt.insert_cols(0, ptl_smp_gen_cnt_tmp);
        } else {
          ptl_smp_gen_cnt = ptl_smp_gen_cnt_tmp;
        }
        cnt += 1;
      }
    }
  }
  
  // ??/B1B2
  if (mis_gen_cnt(16) > 0) {
    arma::uword cnt = 0;
    arma::imat ptl_smp_gen_cnt_0012 = ptl_smp_gen_cnt;
    for (arma::sword i = 0; i <= mis_gen_cnt(16); i++) {
      for (arma::sword j = 0; j <= mis_gen_cnt(16) - i; j++) {
        ptl_mis_gen_cnt.zeros();
        ptl_mis_gen_cnt(1) = i;
        ptl_mis_gen_cnt(4) = j;
        ptl_mis_gen_cnt(7) = mis_gen_cnt(16) - i - j;
        arma::imat ptl_smp_gen_cnt_tmp = ptl_smp_gen_cnt_0012;
        ptl_smp_gen_cnt_tmp.each_col() += ptl_mis_gen_cnt;
        if (cnt > 0) {
          ptl_smp_gen_cnt.insert_cols(0, ptl_smp_gen_cnt_tmp);
        } else {
          ptl_smp_gen_cnt = ptl_smp_gen_cnt_tmp;
        }
        cnt += 1;
      }
    }
  }
  
  // ??/B2B2
  if (mis_gen_cnt(17) > 0) {
    arma::uword cnt = 0;
    arma::imat ptl_smp_gen_cnt_0022 = ptl_smp_gen_cnt;
    for (arma::sword i = 0; i <= mis_gen_cnt(17); i++) {
      for (arma::sword j = 0; j <= mis_gen_cnt(17) - i; j++) {
        ptl_mis_gen_cnt.zeros();
        ptl_mis_gen_cnt(2) = i;
        ptl_mis_gen_cnt(5) = j;
        ptl_mis_gen_cnt(8) = mis_gen_cnt(17) - i - j;
        arma::imat ptl_smp_gen_cnt_tmp = ptl_smp_gen_cnt_0022;
        ptl_smp_gen_cnt_tmp.each_col() += ptl_mis_gen_cnt;
        if (cnt > 0) {
          ptl_smp_gen_cnt.insert_cols(0, ptl_smp_gen_cnt_tmp);
        } else {
          ptl_smp_gen_cnt = ptl_smp_gen_cnt_tmp;
        }
        cnt += 1;
      }
    }
  }
  
  // ?A1/?B1
  if (mis_gen_cnt(18) > 0) {
    arma::uword cnt = 0;
    arma::imat ptl_smp_gen_cnt_0101 = ptl_smp_gen_cnt;
    for (arma::sword i = 0; i <= mis_gen_cnt(18); i++) {
      for (arma::sword j = 0; j <= mis_gen_cnt(18) - i; j++) {
        for (arma::sword k = 0; k <= mis_gen_cnt(18) - i - j; k++) {
          ptl_mis_gen_cnt.zeros();
          ptl_mis_gen_cnt(0) = i;
          ptl_mis_gen_cnt(1) = j;
          ptl_mis_gen_cnt(3) = k;
          ptl_mis_gen_cnt(4) = mis_gen_cnt(18) - i - j - k;
          arma::imat ptl_smp_gen_cnt_tmp = ptl_smp_gen_cnt_0101;
          ptl_smp_gen_cnt_tmp.each_col() += ptl_mis_gen_cnt;
          if (cnt > 0) {
            ptl_smp_gen_cnt.insert_cols(0, ptl_smp_gen_cnt_tmp);
          } else {
            ptl_smp_gen_cnt = ptl_smp_gen_cnt_tmp;
          }
          cnt += 1;
        }
      }
    }
  }
  
  // ?A1/?B2
  if (mis_gen_cnt(19) > 0) {
    arma::uword cnt = 0;
    arma::imat ptl_smp_gen_cnt_0102 = ptl_smp_gen_cnt;
    for (arma::sword i = 0; i <= mis_gen_cnt(19); i++) {
      for (arma::sword j = 0; j <= mis_gen_cnt(19) - i; j++) {
        for (arma::sword k = 0; k <= mis_gen_cnt(19) - i - j; k++) {
          ptl_mis_gen_cnt.zeros();
          ptl_mis_gen_cnt(1) = i;
          ptl_mis_gen_cnt(2) = j;
          ptl_mis_gen_cnt(4) = k;
          ptl_mis_gen_cnt(5) = mis_gen_cnt(19) - i - j - k;
          arma::imat ptl_smp_gen_cnt_tmp = ptl_smp_gen_cnt_0102;
          ptl_smp_gen_cnt_tmp.each_col() += ptl_mis_gen_cnt;
          if (cnt > 0) {
            ptl_smp_gen_cnt.insert_cols(0, ptl_smp_gen_cnt_tmp);
          } else {
            ptl_smp_gen_cnt = ptl_smp_gen_cnt_tmp;
          }
          cnt += 1;
        }
      }
    }
  }
  
  // ?A2/?B1
  if (mis_gen_cnt(20) > 0) {
    arma::uword cnt = 0;
    arma::imat ptl_smp_gen_cnt_0201 = ptl_smp_gen_cnt;
    for (arma::sword i = 0; i <= mis_gen_cnt(20); i++) {
      for (arma::sword j = 0; j <= mis_gen_cnt(20) - i; j++) {
        for (arma::sword k = 0; k <= mis_gen_cnt(20) - i - j; k++) {
          ptl_mis_gen_cnt.zeros();
          ptl_mis_gen_cnt(3) = i;
          ptl_mis_gen_cnt(4) = j;
          ptl_mis_gen_cnt(6) = k;
          ptl_mis_gen_cnt(7) = mis_gen_cnt(20) - i - j - k;
          arma::imat ptl_smp_gen_cnt_tmp = ptl_smp_gen_cnt_0201;
          ptl_smp_gen_cnt_tmp.each_col() += ptl_mis_gen_cnt;
          if (cnt > 0) {
            ptl_smp_gen_cnt.insert_cols(0, ptl_smp_gen_cnt_tmp);
          } else {
            ptl_smp_gen_cnt = ptl_smp_gen_cnt_tmp;
          }
          cnt += 1;
        }
      }
    }
  }
  
  // ?A2/?B2
  if (mis_gen_cnt(21) > 0) {
    arma::uword cnt = 0;
    arma::imat ptl_smp_gen_cnt_0202 = ptl_smp_gen_cnt;
    for (arma::sword i = 0; i <= mis_gen_cnt(21); i++) {
      for (arma::sword j = 0; j <= mis_gen_cnt(21) - i; j++) {
        for (arma::sword k = 0; k <= mis_gen_cnt(21) - i - j; k++) {
          ptl_mis_gen_cnt.zeros();
          ptl_mis_gen_cnt(4) = i;
          ptl_mis_gen_cnt(5) = j;
          ptl_mis_gen_cnt(7) = k;
          ptl_mis_gen_cnt(8) = mis_gen_cnt(21) - i - j - k;
          arma::imat ptl_smp_gen_cnt_tmp = ptl_smp_gen_cnt_0202;
          ptl_smp_gen_cnt_tmp.each_col() += ptl_mis_gen_cnt;
          if (cnt > 0) {
            ptl_smp_gen_cnt.insert_cols(0, ptl_smp_gen_cnt_tmp);
          } else {
            ptl_smp_gen_cnt = ptl_smp_gen_cnt_tmp;
          }
          cnt += 1;
        }
      }
    }
  }
  
  // ?A1/??
  if (mis_gen_cnt(22) > 0) {
    arma::uword cnt = 0;
    arma::imat ptl_smp_gen_cnt_0100 = ptl_smp_gen_cnt;
    for (arma::sword i = 0; i <= mis_gen_cnt(22); i++) {
      for (arma::sword j = 0; j <= mis_gen_cnt(22) - i; j++) {
        for (arma::sword k = 0; k <= mis_gen_cnt(22) - i - j; k++) {
          for (arma::sword l = 0; l <= mis_gen_cnt(22) - i - j - k; l++) {
            for (arma::sword m = 0; m <= mis_gen_cnt(22) - i - j - k - l; m++) {
              ptl_mis_gen_cnt.zeros();
              ptl_mis_gen_cnt(0) = i;
              ptl_mis_gen_cnt(1) = j;
              ptl_mis_gen_cnt(2) = k;
              ptl_mis_gen_cnt(3) = l;
              ptl_mis_gen_cnt(4) = m;
              ptl_mis_gen_cnt(5) = mis_gen_cnt(22) - i - j - k - l - m;
              arma::imat ptl_smp_gen_cnt_tmp = ptl_smp_gen_cnt_0100;
              ptl_smp_gen_cnt_tmp.each_col() += ptl_mis_gen_cnt;
              if (cnt > 0) {
                ptl_smp_gen_cnt.insert_cols(0, ptl_smp_gen_cnt_tmp);
              } else {
                ptl_smp_gen_cnt = ptl_smp_gen_cnt_tmp;
              }
              cnt += 1;
            }
          }
        }
      }
    }
  }
  
  // ?A2/??
  if (mis_gen_cnt(23) > 0) {
    arma::uword cnt = 0;
    arma::imat ptl_smp_gen_cnt_0200 = ptl_smp_gen_cnt;
    for (arma::sword i = 0; i <= mis_gen_cnt(23); i++) {
      for (arma::sword j = 0; j <= mis_gen_cnt(23) - i; j++) {
        for (arma::sword k = 0; k <= mis_gen_cnt(23) - i - j; k++) {
          for (arma::sword l = 0; l <= mis_gen_cnt(23) - i - j - k; l++) {
            for (arma::sword m = 0; m <= mis_gen_cnt(23) - i - j - k - l; m++) {
              ptl_mis_gen_cnt.zeros();
              ptl_mis_gen_cnt(3) = i;
              ptl_mis_gen_cnt(4) = j;
              ptl_mis_gen_cnt(5) = k;
              ptl_mis_gen_cnt(6) = l;
              ptl_mis_gen_cnt(7) = m;
              ptl_mis_gen_cnt(8) = mis_gen_cnt(23) - i - j - k - l - m;
              arma::imat ptl_smp_gen_cnt_tmp = ptl_smp_gen_cnt_0200;
              ptl_smp_gen_cnt_tmp.each_col() += ptl_mis_gen_cnt;
              if (cnt > 0) {
                ptl_smp_gen_cnt.insert_cols(0, ptl_smp_gen_cnt_tmp);
              } else {
                ptl_smp_gen_cnt = ptl_smp_gen_cnt_tmp;
              }
              cnt += 1;
            }
          }
        }
      }
    }
  }
  
  // ??/?B1
  if (mis_gen_cnt(24) > 0) {
    arma::uword cnt = 0;
    arma::imat ptl_smp_gen_cnt_0001 = ptl_smp_gen_cnt;
    for (arma::sword i = 0; i <= mis_gen_cnt(24); i++) {
      for (arma::sword j = 0; j <= mis_gen_cnt(24) - i; j++) {
        for (arma::sword k = 0; k <= mis_gen_cnt(24) - i - j; k++) {
          for (arma::sword l = 0; l <= mis_gen_cnt(24) - i - j - k; l++) {
            for (arma::sword m = 0; m <= mis_gen_cnt(24) - i - j - k - l; m++) {
              ptl_mis_gen_cnt.zeros();
              ptl_mis_gen_cnt(0) = i;
              ptl_mis_gen_cnt(1) = j;
              ptl_mis_gen_cnt(3) = k;
              ptl_mis_gen_cnt(4) = l;
              ptl_mis_gen_cnt(6)= m;
              ptl_mis_gen_cnt(7) = mis_gen_cnt(24) - i - j - k - l - m;
              arma::imat ptl_smp_gen_cnt_tmp = ptl_smp_gen_cnt_0001;
              ptl_smp_gen_cnt_tmp.each_col() += ptl_mis_gen_cnt;
              if (cnt > 0) {
                ptl_smp_gen_cnt.insert_cols(0, ptl_smp_gen_cnt_tmp);
              } else {
                ptl_smp_gen_cnt = ptl_smp_gen_cnt_tmp;
              }
              cnt += 1;
            }
          }
        }
      }
    }
  }
  
  // ??/?B2
  if (mis_gen_cnt(25) > 0) {
    arma::uword cnt = 0;
    arma::imat ptl_smp_gen_cnt_0002 = ptl_smp_gen_cnt;
    for (arma::sword i = 0; i <= mis_gen_cnt(25); i++) {
      for (arma::sword j = 0; j <= mis_gen_cnt(25) - i; j++) {
        for (arma::sword k = 0; k <= mis_gen_cnt(25) - i - j; k++) {
          for (arma::sword l = 0; l <= mis_gen_cnt(25) - i - j - k; l++) {
            for (arma::sword m = 0; m <= mis_gen_cnt(25) - i - j - k - l; m++) {
              ptl_mis_gen_cnt.zeros();
              ptl_mis_gen_cnt(1) = i;
              ptl_mis_gen_cnt(2) = j;
              ptl_mis_gen_cnt(4) = k;
              ptl_mis_gen_cnt(5) = l;
              ptl_mis_gen_cnt(7) = m;
              ptl_mis_gen_cnt(8) = mis_gen_cnt(25) - i - j - k - l - m;
              arma::imat ptl_smp_gen_cnt_tmp = ptl_smp_gen_cnt_0002;
              ptl_smp_gen_cnt_tmp.each_col() += ptl_mis_gen_cnt;
              if (cnt > 0) {
                ptl_smp_gen_cnt.insert_cols(0, ptl_smp_gen_cnt_tmp);
              } else {
                ptl_smp_gen_cnt = ptl_smp_gen_cnt_tmp;
              }
              cnt += 1;
            }
          }
        }
      }
    }
  }
  
  // ??/??
  if (mis_gen_cnt(26) > 0) {
    arma::uword cnt = 0;
    arma::imat ptl_smp_gen_cnt_0000 = ptl_smp_gen_cnt;
    for (arma::sword i = 0; i <= mis_gen_cnt(26); i++) {
      for (arma::sword j = 0; j <= mis_gen_cnt(26) - i; j++) {
        for (arma::sword k = 0; k <= mis_gen_cnt(26) - i - j; k++) {
          for (arma::sword l = 0; l <= mis_gen_cnt(26) - i - j - k; l++) {
            for (arma::sword m = 0; m <= mis_gen_cnt(26) - i - j - k - l; m++) {
              for (arma::sword n = 0; n <= mis_gen_cnt(26) - i - j - k - l - m; n++) {
                for (arma::sword o = 0; o <= mis_gen_cnt(26) - i - j - k - l - m - n; o++) {
                  for (arma::sword p = 0; p <= mis_gen_cnt(26) - i - j - k - l - m - n - o; p++) {
                    ptl_mis_gen_cnt.zeros();
                    ptl_mis_gen_cnt(0) = i;
                    ptl_mis_gen_cnt(1) = j;
                    ptl_mis_gen_cnt(2) = k;
                    ptl_mis_gen_cnt(3) = l;
                    ptl_mis_gen_cnt(4) = m;
                    ptl_mis_gen_cnt(5) = n;
                    ptl_mis_gen_cnt(6) = o;
                    ptl_mis_gen_cnt(7) = p;
                    ptl_mis_gen_cnt(8) = mis_gen_cnt(26) - i - j - k - l - m - n - o - p;
                    arma::imat ptl_smp_gen_cnt_tmp = ptl_smp_gen_cnt_0000;
                    ptl_smp_gen_cnt_tmp.each_col() += ptl_mis_gen_cnt;
                    if (cnt > 0) {
                      ptl_smp_gen_cnt.insert_cols(0, ptl_smp_gen_cnt_tmp);
                    } else {
                      ptl_smp_gen_cnt = ptl_smp_gen_cnt_tmp;
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
  
  return ptl_smp_gen_cnt;
}

// Convert the unphased genotype counts to the phased genotype counts
// [[Rcpp::export]]
arma::imat convertUnphasedGenoCount2PhasedGenoCount_arma(const arma::icolvec& smp_gen_cnt) {
  // ensure RNG gets set/reset
  RNGScope scope;
  
  // declare the phased genotype counts
  arma::imat ptl_smp_cnt(10, smp_gen_cnt(4) + 1);
  
  for (arma::uword k = 0; k <= smp_gen_cnt(4); k++) {
    ptl_smp_cnt(0, k) = smp_gen_cnt(0);
    ptl_smp_cnt(1, k) = smp_gen_cnt(1);
    ptl_smp_cnt(2, k) = smp_gen_cnt(2);
    ptl_smp_cnt(3, k) = smp_gen_cnt(3);
    ptl_smp_cnt(4, k) = k;
    ptl_smp_cnt(5, k) = smp_gen_cnt(5);
    ptl_smp_cnt(6, k) = smp_gen_cnt(6);
    ptl_smp_cnt(7, k) = smp_gen_cnt(7);
    ptl_smp_cnt(8, k) = smp_gen_cnt(8);
    ptl_smp_cnt(9, k) = smp_gen_cnt(4) - k;
  }
  
  return ptl_smp_cnt;
}

// Calculate the possible unobserved sample phased genotype counts from the observed sample unphased genotype counts with uncertainties
// [[Rcpp::export]]
arma::imat calculatePhasedGenoCount_arma(const arma::icolvec& smp_gen_cnt, const arma::icolvec& mis_gen_cnt) {
  // ensure RNG gets set/reset
  RNGScope scope;
  
  arma::imat ptl_smp_gen_cnt = calculateUnphasedGenoCount_arma(smp_gen_cnt, mis_gen_cnt);
  
  arma::imat ptl_smp_cnt = convertUnphasedGenoCount2PhasedGenoCount_arma(ptl_smp_gen_cnt.col(0));
  for (arma::uword k = 1; k < ptl_smp_gen_cnt.n_cols; k++) {
    arma::imat tmp_smp_cnt = convertUnphasedGenoCount2PhasedGenoCount_arma(ptl_smp_gen_cnt.col(k));
    ptl_smp_cnt.insert_cols(0, tmp_smp_cnt);
  }
  
  return ptl_smp_cnt;
}
/*************************/


/********** EMN **********/
// Calculate the multinomial probabilities
// [[Rcpp::export]]
double calculateMultinomialProb_arma(const arma::icolvec& smp_cnt, const int& smp_siz, const arma::dcolvec& pop_frq) {
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
/*************************/


/********** BPF **********/
// Run the particle filter
// [[Rcpp::export]]
List runBPF_arma(const arma::dcolvec& sel_cof, const double& rec_rat, const int& pop_siz, const arma::irowvec& smp_gen, const arma::irowvec& smp_siz, const arma::imat& smp_cnt, const arma::imat& mis_cnt, const arma::uword& ptn_num, const arma::uword& pcl_num) {
  // ensure RNG gets set/reset
  RNGScope scope;
  
  double lik = 1;
  
  arma::dmat wght = arma::zeros<arma::dmat>(pcl_num, smp_gen.n_elem);
  arma::dcube part_pre = arma::zeros<arma::dcube>(10, pcl_num, smp_gen.n_elem);
  arma::dcube part_pst = arma::zeros<arma::dcube>(10, pcl_num, smp_gen.n_elem);
  
  arma::dcolvec wght_tmp = arma::zeros<arma::dcolvec>(pcl_num);
  arma::dmat part_tmp = arma::zeros<arma::dmat>(10, pcl_num);
  
  arma::uword dom_ind = arma::as_scalar(arma::find(smp_siz == 0));
  
  arma::dcolvec sel_cof_bd = sel_cof.rows(0, 1);
  arma::dcolvec sel_cof_ad = sel_cof.rows(2, 3);
  
  // initialise the particles
  cout << "generation: " << smp_gen(0) << endl;
  part_tmp = arma::normalise(arma::randu<arma::dmat>(10, pcl_num), 1, 0);
  arma::imat ptl_smp_cnt = calculatePhasedGenoCount_arma(smp_cnt.col(0), mis_cnt.col(0));
  for (arma::uword i = 0; i < pcl_num; i++) {
    for (arma::uword j = 0; j < ptl_smp_cnt.n_cols; j++) {
      wght_tmp(i) = wght_tmp(i) + calculateMultinomialProb_arma(ptl_smp_cnt.col(j), smp_siz(0), part_tmp.col(i)) / smp_cnt.n_cols;
    }
  }
  
  if (arma::any(wght_tmp > 0)) {
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
  for (arma::uword k = 1; k < dom_ind; k++) {
    cout << "generation: " << smp_gen(k) << endl;
    wght_tmp = arma::zeros<arma::dcolvec>(pcl_num);
    part_tmp = part_pst.slice(k - 1);
    arma::imat ptl_smp_cnt = calculatePhasedGenoCount_arma(smp_cnt.col(k), mis_cnt.col(k));
    for (arma::uword i = 0; i < pcl_num; i++) {
      arma::dmat path = simulateTLWFDS_arma(sel_cof_bd, rec_rat, pop_siz, part_tmp.col(i), smp_gen(k - 1), smp_gen(k), ptn_num);
      part_tmp.col(i) = arma::vectorise(path.tail_cols(1), 0);
      for (arma::uword j = 0; j < ptl_smp_cnt.n_cols; j++) {
        wght_tmp(i) = wght_tmp(i) + calculateMultinomialProb_arma(ptl_smp_cnt.col(j), smp_siz(k), part_tmp.col(i)) / smp_cnt.n_cols;
      }
    }
    
    if (arma::any(wght_tmp > 0)) {
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
      
      return List::create(Named("lik", lik), 
                          Named("wght", wght), 
                          Named("part_pre_resmp", part_pre), 
                          Named("part_pst_resmp", part_pst));
    }
  }
  
  if (smp_gen(dom_ind) == smp_gen(dom_ind - 1)) {
    part_pre.slice(dom_ind) = part_pre.slice(dom_ind - 1);
    part_pst.slice(dom_ind) = part_pst.slice(dom_ind - 1);
  } else {
    cout << "generation: " << smp_gen(dom_ind) << endl;
    part_tmp = part_pst.slice(dom_ind - 1);
    for (arma::uword i = 0; i < pcl_num; i++) {
      arma::dmat path = simulateTLWFDS_arma(sel_cof_bd, rec_rat, pop_siz, part_tmp.col(i), smp_gen(dom_ind - 1), smp_gen(dom_ind), ptn_num);
      part_tmp.col(i) = arma::vectorise(path.tail_cols(1), 0);
    }
    part_pre.slice(dom_ind) = part_tmp;
    part_pst.slice(dom_ind) = part_tmp;
  }
  
  for (arma::uword k = dom_ind + 1; k < smp_gen.n_elem; k++) {
    cout << "generation: " << smp_gen(k) << endl;
    wght_tmp = arma::zeros<arma::dcolvec>(pcl_num);
    part_tmp = part_pst.slice(k - 1);
    arma::imat ptl_smp_cnt = calculatePhasedGenoCount_arma(smp_cnt.col(k), mis_cnt.col(k));
    for (arma::uword i = 0; i < pcl_num; i++) {
      arma::dmat path = simulateTLWFDS_arma(sel_cof_ad, rec_rat, pop_siz, part_tmp.col(i), smp_gen(k - 1), smp_gen(k), ptn_num);
      part_tmp.col(i) = arma::vectorise(path.tail_cols(1), 0);
      for (arma::uword j = 0; j < ptl_smp_cnt.n_cols; j++) {
        wght_tmp(i) = wght_tmp(i) + calculateMultinomialProb_arma(ptl_smp_cnt.col(j), smp_siz(k), part_tmp.col(i)) / smp_cnt.n_cols;
      }
    }
    
    if (arma::any(wght_tmp > 0)) {
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
  
  wght.shed_col(dom_ind);
  part_pre.shed_slice(dom_ind);
  part_pst.shed_slice(dom_ind);
  
  return List::create(Named("lik", lik),
                      Named("wght", wght), 
                      Named("part_pre_resmp", part_pre), 
                      Named("part_pst_resmp", part_pst));
}
/*************************/


/********** PMMH **********/
// Simulate the phased genotype frequency trajectories in a pre-determined generation according to the two-locus Wright-Fisher diffusion with selection using the Euler-Maruyama method
// [[Rcpp::export]]
arma::dcolvec simulateWFD_arma(const arma::dmat& fts, const double& scl_rec_rat, const arma::dcolvec& int_frq, const double& dt, const int& int_gen, const int& lst_gen, const arma::uword& ptn_num) {
  // ensure RNG gets set/reset
  RNGScope scope;
  
  // declare eta
  arma::dcolvec eta(4);
  eta(0) = -1;
  eta(1) = 1;
  eta(2) = 1;
  eta(3) = -1;
  
  // declare the genotype frequency trajectories
  arma::dcolvec frq_pth = int_frq;
  
  for(arma::uword t = 1; t < arma::uword(lst_gen - int_gen) + 1; t++) {
    // declare delta W
    arma::dmat dW = pow(dt, 0.5) * arma::randn<arma::dmat>(6, ptn_num);
    
    arma::dcolvec hap_frq = arma::zeros<arma::dcolvec>(4);
    hap_frq(0) = frq_pth(0) + frq_pth(1) / 2 + frq_pth(3) / 2 + frq_pth(9) / 2;
    hap_frq(1) = frq_pth(2) + frq_pth(1) / 2 + frq_pth(4) / 2 + frq_pth(5) / 2;
    hap_frq(2) = frq_pth(6) + frq_pth(3) / 2 + frq_pth(4) / 2 + frq_pth(7) / 2;
    hap_frq(3) = frq_pth(8) + frq_pth(9) / 2 + frq_pth(5) / 2 + frq_pth(7) / 2;
    
    for(arma::uword k = 0; k < ptn_num; k++) {
      // calculate the drift coefficient vector
      arma::dcolvec mu = eta * 0.5 * scl_rec_rat * (hap_frq(0) * hap_frq(3) - hap_frq(1) * hap_frq(2));
      
      // calculate the diffusion coefficient matrix
      arma::dmat sigma = arma::zeros<arma::dmat>(4, 6);
      sigma(0, 0) = pow(hap_frq(0) * hap_frq(1), 0.5);
      sigma(0, 1) = pow(hap_frq(0) * hap_frq(2), 0.5);
      sigma(0, 2) = pow(hap_frq(0) * hap_frq(3), 0.5);
      //sigma(0, 3) = 0;
      //sigma(0, 4) = 0;
      //sigma(0, 5) = 0;
      sigma(1, 0) = -pow(hap_frq(1) * hap_frq(0), 0.5);
      //sigma(1, 1) = 0;
      //sigma(1, 2) = 0;
      sigma(1, 3) = pow(hap_frq(1) * hap_frq(2), 0.5);
      sigma(1, 4) = pow(hap_frq(1) * hap_frq(3), 0.5);
      //sigma(1, 5) = 0;
      //sigma(2, 0) = 0;
      sigma(2, 1) = -pow(hap_frq(2) * hap_frq(0), 0.5);
      //sigma(2, 2) = 0;
      sigma(2, 3) = -pow(hap_frq(2) * hap_frq(1), 0.5);
      //sigma(2, 4) = 0;
      sigma(2, 5) = pow(hap_frq(2) * hap_frq(3), 0.5);
      //sigma(3, 0) = 0;
      //sigma(3, 1) = 0;
      sigma(3, 2) = -pow(hap_frq(3) * hap_frq(0), 0.5);
      //sigma(3, 3) = 0;
      sigma(3, 4) = -pow(hap_frq(3) * hap_frq(1), 0.5);
      sigma(3, 5) = -pow(hap_frq(3) * hap_frq(2), 0.5);
      
      // proceed the Euler-Maruyama scheme
      hap_frq = hap_frq + mu * dt + sigma * dW.col(k);
      
      // remove the noise from the numerical techniques
      for(arma::uword i = 0; i < 4; i++) {
        if(hap_frq(i) < 0) {
          hap_frq(i) = 0;
        } 
        if(hap_frq(i) > 1) {
          hap_frq(i) = 1;
        }
      }
      hap_frq = hap_frq / arma::accu(hap_frq);
    }
    
    // random union of gametes
    arma::dmat gen_frq = hap_frq * hap_frq.t();
    
    // viability selection
    gen_frq = (fts % gen_frq) / arma::accu(fts % gen_frq);
    frq_pth(0) = gen_frq(0, 0);
    frq_pth(1) = gen_frq(0, 1) + gen_frq(1, 0);
    frq_pth(2) = gen_frq(1, 1);
    frq_pth(3) = gen_frq(0, 2) + gen_frq(2, 0);
    frq_pth(4) = gen_frq(1, 2) + gen_frq(2, 1);
    frq_pth(5) = gen_frq(1, 3) + gen_frq(3, 1);
    frq_pth(6) = gen_frq(2, 2);
    frq_pth(7) = gen_frq(2, 3) + gen_frq(3, 2);
    frq_pth(8) = gen_frq(3, 3);
    frq_pth(9) = gen_frq(0, 3) + gen_frq(3, 0);
  }
  
  return frq_pth;
}

// Calculate the log-likelihood using the bootstrap particle filter
// [[Rcpp::export]]
void calculateLogLikelihood_arma(double& log_lik, const arma::dcolvec& sel_cof, const double& scl_rec_rat, const arma::irowvec& smp_gen, const arma::irowvec& smp_siz, const arma::field<arma::imat>& ptl_smp_cnt, const double& dt, const arma::uword& ptn_num, const arma::uword& pcl_num) {
  // ensure RNG gets set/reset
  RNGScope scope;
  
  arma::dcolvec wght = arma::zeros<arma::dcolvec>(pcl_num);
  arma::dmat part_pre = arma::zeros<arma::dmat>(10, pcl_num);
  arma::dmat part_pst = arma::zeros<arma::dmat>(10, pcl_num);
  
  arma::uword dom_ind = arma::as_scalar(arma::find(smp_siz == 0));
 
  arma::dmat fts_bd = arma::ones<arma::dmat>(4, 4);
  fts_bd(1, 1) = 1 + sel_cof(1);
  fts_bd(3, 1) = 1 + sel_cof(1);
  fts_bd(2, 2) = 1 + sel_cof(0);
  fts_bd(3, 2) = 1 + sel_cof(0);
  fts_bd(1, 3) = 1 + sel_cof(1);
  fts_bd(2, 3) = 1 + sel_cof(0);
  fts_bd(3, 3) = 1 + sel_cof(1);
  
  arma::dmat fts_ad = arma::ones<arma::dmat>(4, 4);
  fts_ad(1, 1) = 1 + sel_cof(3);
  fts_ad(3, 1) = 1 + sel_cof(3);
  fts_ad(2, 2) = 1 + sel_cof(2);
  fts_ad(3, 2) = 1 + sel_cof(2);
  fts_ad(1, 3) = 1 + sel_cof(3);
  fts_ad(2, 3) = 1 + sel_cof(2);
  fts_ad(3, 3) = 1 + sel_cof(3);
  
  // initialise the particles
  part_pre = arma::normalise(arma::randu<arma::dmat>(10, pcl_num), 1, 0);
  arma::imat smp_cnt = ptl_smp_cnt(0);
  for (arma::uword i = 0; i < pcl_num; i++) {
    for (arma::uword j = 0; j < smp_cnt.n_cols; j++) {
      wght(i) = wght(i) + calculateMultinomialProb_arma(smp_cnt.col(j), smp_siz(0), part_pre.col(i));
    }
  }
  
  if (arma::any(wght > 0)) {
    log_lik = log(arma::mean(wght) / smp_cnt.n_cols);
    arma::dcolvec prob = arma::normalise(wght, 1);
    arma::ucolvec elem = arma::linspace<arma::ucolvec>(0, pcl_num - 1, pcl_num);
    arma::ucolvec indx = RcppArmadillo::sample(elem, pcl_num, TRUE, prob);
    part_pst = part_pre.cols(indx);
  } else {
    log_lik = -(arma::datum::inf);
    return;
  }
  
  // run the bootstrap particle filter
  for (arma::uword k = 1; k < dom_ind; k++) {
    wght = arma::zeros<arma::dcolvec>(pcl_num);
    arma::imat smp_cnt = ptl_smp_cnt(k);
    for (arma::uword i = 0; i < pcl_num; i++) {
      part_pre.col(i) = simulateWFD_arma(fts_bd, scl_rec_rat, part_pst.col(i), dt, smp_gen(k - 1), smp_gen(k), ptn_num);
      for (arma::uword j = 0; j < smp_cnt.n_cols; j++) {
        wght(i) = wght(i) + calculateMultinomialProb_arma(smp_cnt.col(j), smp_siz(k), part_pre.col(i));
      }
    }
    
    if (arma::any(wght > 0)) {
      log_lik = log_lik + log(arma::mean(wght) / smp_cnt.n_cols);
      arma::dcolvec prob = arma::normalise(wght, 1);
      arma::ucolvec elem = arma::linspace<arma::ucolvec>(0, pcl_num - 1, pcl_num);
      arma::ucolvec indx = RcppArmadillo::sample(elem, pcl_num, TRUE, prob);
      part_pst = part_pre.cols(indx);
    } else {
      log_lik = -(arma::datum::inf);
      return;
    }
  }
  
  if (smp_gen(dom_ind) != smp_gen(dom_ind - 1)) {
    for (arma::uword i = 0; i < pcl_num; i++) {
      part_pre.col(i) = simulateWFD_arma(fts_bd, scl_rec_rat, part_pst.col(i), dt, smp_gen(dom_ind - 1), smp_gen(dom_ind), ptn_num);
    }
    part_pst = part_pre;
  }
  
  for (arma::uword k = dom_ind + 1; k < smp_gen.n_elem; k++) {
    wght = arma::zeros<arma::dcolvec>(pcl_num);
    arma::imat smp_cnt = ptl_smp_cnt(k);
    for (arma::uword i = 0; i < pcl_num; i++) {
      part_pre.col(i) = simulateWFD_arma(fts_ad, scl_rec_rat, part_pst.col(i), dt, smp_gen(k - 1), smp_gen(k), ptn_num);
      for (arma::uword j = 0; j < smp_cnt.n_cols; j++) {
        wght(i) = wght(i) + calculateMultinomialProb_arma(smp_cnt.col(j), smp_siz(k), part_pre.col(i));
      }
    }
    
    if (arma::any(wght > 0)) {
      log_lik = log_lik + log(arma::mean(wght) / smp_cnt.n_cols);
      arma::dcolvec prob = arma::normalise(wght, 1);
      arma::ucolvec elem = arma::linspace<arma::ucolvec>(0, pcl_num - 1, pcl_num);
      arma::ucolvec indx = RcppArmadillo::sample(elem, pcl_num, TRUE, prob);
      part_pst = part_pre.cols(indx);
    } else {
      log_lik = -(arma::datum::inf);
      return;
    }
  }
}

// Calculate the optimal particle number in the particle marginal Metropolis-Hastings
// [[Rcpp::export]]
List calculateOptimalParticleNum_arma(const arma::dcolvec& sel_cof, const double& rec_rat, const int& pop_siz, const arma::irowvec& smp_gen, const arma::irowvec& smp_siz, const arma::imat& smp_cnt, const arma::imat& mis_cnt, const arma::uword& ptn_num, const arma::uword& pcl_num, const arma::uword& gap_num) {
  // ensure RNG gets set/reset
  RNGScope scope;
  
  double scl_rec_rat = 4 * pop_siz * rec_rat;
  
  // declare delta t
  double dt = 1.0 / (2 * pop_siz) / ptn_num;
  
  arma::field<arma::imat> ptl_smp_cnt(smp_gen.n_elem);
  for (arma::uword k = 0; k < smp_gen.n_elem; k++) {
    ptl_smp_cnt(k) = calculatePhasedGenoCount_arma(smp_cnt.col(k), mis_cnt.col(k));
  }
  
  arma::drowvec log_lik(300);
  for (arma::uword i = 0; i < 300; i++) {
    calculateLogLikelihood_arma(log_lik(i), sel_cof, scl_rec_rat, smp_gen, smp_siz, ptl_smp_cnt, dt, ptn_num, pcl_num);
  }
  
  arma::drowvec log_lik_sdv(1);
  log_lik_sdv(0) = arma::stddev(log_lik);
  log_lik_sdv.print();
  arma::urowvec opt_pcl_num(1);
  opt_pcl_num(0) = pcl_num;
  
  if (log_lik_sdv(0) > 1.7) {
    while (log_lik_sdv(0) > 1.0 && opt_pcl_num(0) > gap_num) {
      opt_pcl_num.insert_cols(0, 1);
      opt_pcl_num(0) = opt_pcl_num(1) + gap_num;
      for (arma::uword i = 0; i < 300; i++) {
        calculateLogLikelihood_arma(log_lik(i), sel_cof, scl_rec_rat, smp_gen, smp_siz, ptl_smp_cnt, dt, ptn_num, opt_pcl_num(0));
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
        calculateLogLikelihood_arma(log_lik(i), sel_cof, scl_rec_rat, smp_gen, smp_siz, ptl_smp_cnt, dt, ptn_num, opt_pcl_num(0));
      }
      log_lik_sdv.insert_cols(0, 1);
      log_lik_sdv(0) = arma::stddev(log_lik);
      log_lik_sdv.print();
    }
  } else {
    while (log_lik_sdv(0) > 1.0 && opt_pcl_num(0) > gap_num) {
      opt_pcl_num.insert_cols(0, 1);
      opt_pcl_num(0) = opt_pcl_num(1) + gap_num;
      for (arma::uword i = 0; i < 300; i++) {
        calculateLogLikelihood_arma(log_lik(i), sel_cof, scl_rec_rat, smp_gen, smp_siz, ptl_smp_cnt, dt, ptn_num, opt_pcl_num(0));
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
        calculateLogLikelihood_arma(log_lik(i), sel_cof, scl_rec_rat, smp_gen, smp_siz, ptl_smp_cnt, dt, ptn_num, opt_pcl_num(0));
      }
      log_lik_sdv.insert_cols(0, 1);
      log_lik_sdv(0) = arma::stddev(log_lik);
      log_lik_sdv.print();
    }
  }
  
  return List::create(Named("opt_pcl_num", opt_pcl_num), 
                      Named("log_lik_sdv", log_lik_sdv));
}

// Run the Metropolis-Hastings
// [[Rcpp::export]]
void runMH_arma(bool& apt, double& log_lik, const double& log_lik_pre, const arma::dcolvec& sel_cof, const double& scl_rec_rat, const arma::irowvec& smp_gen, const arma::irowvec& smp_siz, const arma::field<arma::imat>& ptl_smp_cnt, const double& dt, const arma::uword& ptn_num, const arma::uword& pcl_num) {
  // ensure RNG gets set/reset
  RNGScope scope;
  
  double log_apt_rto = log(arma::as_scalar(arma::randu()));
  
  arma::dcolvec wght = arma::zeros<arma::dcolvec>(pcl_num);
  arma::dmat part_pre = arma::zeros<arma::dmat>(10, pcl_num);
  arma::dmat part_pst = arma::zeros<arma::dmat>(10, pcl_num);
    
  arma::uword dom_ind = arma::as_scalar(arma::find(smp_siz == 0));
  
  arma::dmat fts_bd = arma::ones<arma::dmat>(4, 4);
  fts_bd(1, 1) = 1 + sel_cof(1);
  fts_bd(3, 1) = 1 + sel_cof(1);
  fts_bd(2, 2) = 1 + sel_cof(0);
  fts_bd(3, 2) = 1 + sel_cof(0);
  fts_bd(1, 3) = 1 + sel_cof(1);
  fts_bd(2, 3) = 1 + sel_cof(0);
  fts_bd(3, 3) = 1 + sel_cof(1);
  
  arma::dmat fts_ad = arma::ones<arma::dmat>(4, 4);
  fts_ad(1, 1) = 1 + sel_cof(3);
  fts_ad(3, 1) = 1 + sel_cof(3);
  fts_ad(2, 2) = 1 + sel_cof(2);
  fts_ad(3, 2) = 1 + sel_cof(2);
  fts_ad(1, 3) = 1 + sel_cof(3);
  fts_ad(2, 3) = 1 + sel_cof(2);
  fts_ad(3, 3) = 1 + sel_cof(3);
  
  // initialise the particles
  part_pre = arma::normalise(arma::randu<arma::dmat>(10, pcl_num), 1, 0);
  arma::imat smp_cnt = ptl_smp_cnt(0);
  for (arma::uword i = 0; i < pcl_num; i++) {
    for (arma::uword j = 0; j < smp_cnt.n_cols; j++) {
      wght(i) = wght(i) + calculateMultinomialProb_arma(smp_cnt.col(j), smp_siz(0), part_pre.col(i));
    }
  }
  
  if (arma::any(wght > 0)) {
    log_lik = log(arma::mean(wght) / smp_cnt.n_cols);
    arma::dcolvec prob = arma::normalise(wght, 1);
    arma::ucolvec elem = arma::linspace<arma::ucolvec>(0, pcl_num - 1, pcl_num);
    arma::ucolvec indx = RcppArmadillo::sample(elem, pcl_num, TRUE, prob);
    part_pst = part_pre.cols(indx);
  } else {
    log_lik = log_lik_pre;
    apt = false;
    return;
  }
  
  // run the bootstrap particle filter
  for (arma::uword k = 1; k < dom_ind; k++) {
    wght = arma::zeros<arma::dcolvec>(pcl_num);
    arma::imat smp_cnt = ptl_smp_cnt(k);
    for (arma::uword i = 0; i < pcl_num; i++) {
      part_pre.col(i) = simulateWFD_arma(fts_bd, scl_rec_rat, part_pst.col(i), dt, smp_gen(k - 1), smp_gen(k), ptn_num);
      for (arma::uword j = 0; j < smp_cnt.n_cols; j++) {
        wght(i) = wght(i) + calculateMultinomialProb_arma(smp_cnt.col(j), smp_siz(k), part_pre.col(i));
      }
    }
    
    if (arma::any(wght > 0)) {
      log_lik = log_lik + log(arma::mean(wght) / smp_cnt.n_cols);
      arma::dcolvec prob = arma::normalise(wght, 1);
      arma::ucolvec elem = arma::linspace<arma::ucolvec>(0, pcl_num - 1, pcl_num);
      arma::ucolvec indx = RcppArmadillo::sample(elem, pcl_num, TRUE, prob);
      part_pst = part_pre.cols(indx);
    } else {
      log_lik = log_lik_pre;
      apt = false;
      return;
    }
  }
  
  if (smp_gen(dom_ind) != smp_gen(dom_ind - 1)) {
    for (arma::uword i = 0; i < pcl_num; i++) {
      part_pre.col(i) = simulateWFD_arma(fts_bd, scl_rec_rat, part_pst.col(i), dt, smp_gen(dom_ind - 1), smp_gen(dom_ind), ptn_num);
    }
    part_pst = part_pre;
  }
  
  for (arma::uword k = dom_ind + 1; k < smp_gen.n_elem; k++) {
    wght = arma::zeros<arma::dcolvec>(pcl_num);
    arma::imat smp_cnt = ptl_smp_cnt(k);
    for (arma::uword i = 0; i < pcl_num; i++) {
      part_pre.col(i) = simulateWFD_arma(fts_ad, scl_rec_rat, part_pst.col(i), dt, smp_gen(k - 1), smp_gen(k), ptn_num);
      for (arma::uword j = 0; j < smp_cnt.n_cols; j++) {
        wght(i) = wght(i) + calculateMultinomialProb_arma(smp_cnt.col(j), smp_siz(k), part_pre.col(i));
      }
    }
    
    if (arma::any(wght > 0)) {
      log_lik = log_lik + log(arma::mean(wght) / smp_cnt.n_cols);
      arma::dcolvec prob = arma::normalise(wght, 1);
      arma::ucolvec elem = arma::linspace<arma::ucolvec>(0, pcl_num - 1, pcl_num);
      arma::ucolvec indx = RcppArmadillo::sample(elem, pcl_num, TRUE, prob);
      part_pst = part_pre.cols(indx);
    } else {
      log_lik = log_lik_pre;
      apt = false;
      return;
    }
  }
  
  if (log_apt_rto > log_lik - log_lik_pre) {
    log_lik = log_lik_pre;
    apt = false;
  } else {
    apt = true;
  }
}

// Run the particle marginal Metropolis-Hastings
//[[Rcpp::export]]
arma::dmat runPMMH_arma(const arma::dcolvec& sel_cof, const double& rec_rat, const int& pop_siz, const arma::irowvec& smp_gen, const arma::irowvec& smp_siz, const arma::imat& smp_cnt, const arma::imat& mis_cnt, const arma::uword& ptn_num, const arma::uword& pcl_num, const arma::uword& itn_num) {
  // ensure RNG gets set/reset
  RNGScope scope;
  
  arma::field<arma::imat> ptl_smp_cnt(smp_gen.n_elem);
  for (arma::uword k = 0; k < smp_gen.n_elem; k++) {
    ptl_smp_cnt(k) = calculatePhasedGenoCount_arma(smp_cnt.col(k), mis_cnt.col(k));
  }
  
  arma::dmat sel_cof_chn = arma::zeros<arma::dmat>(4, itn_num);
  
  //arma::drowvec log_pri_chn = arma::zeros<arma::drowvec>(itn_num);
  arma::drowvec log_lik_chn = arma::zeros<arma::drowvec>(itn_num);
  
  // initialise the population genetic parameters
  cout << "iteration: " << 1 << endl;
  // take the uniform prior and fix the initial value of the selection coefficient to zero
  // or take the beta prior with alpha = 1 and beta = 3
  sel_cof_chn.col(0) = sel_cof;
  
  // rescale the recombination rate
  double scl_rec_rat = 4 * pop_siz * rec_rat;
  
  // declare delta t
  double dt = 1.0 / (2 * pop_siz) / ptn_num;
  
  calculateLogLikelihood_arma(log_lik_chn(0), sel_cof_chn.col(0), scl_rec_rat, smp_gen, smp_siz, ptl_smp_cnt, dt, ptn_num, pcl_num);
  
  double sel_cof_sd = 2e-03;
  arma::dmat sel_cof_cov = arma::eye<arma::dmat>(4, 4) * pow(sel_cof_sd, 2);
  
  bool apt = true;
  for (arma::uword i = 1; i < itn_num; i++) {
    cout << "iteration: " << i + 1 << endl;
    apt = true;
    
    // draw the candidates of the selection coefficients from the proposals
    sel_cof_chn.col(i) = arma::mvnrnd(sel_cof_chn.col(i - 1), sel_cof_cov);
    if (arma::any(sel_cof_chn.col(i) < -1)) {
      sel_cof_chn.col(i) = sel_cof_chn.col(i - 1);
      log_lik_chn(i) = log_lik_chn(i - 1);
    } else {
      //double log_psl_old_new = log(arma::normpdf(sel_cof_chn.col(i - 1), sel_cof_chn.col(i), sel_cof_cov));
      //double log_psl_new_old = log(arma::normpdf(sel_cof_chn.col(i), sel_cof_chn.col(i - 1), sel_cof_cov));
      
      // calculate the acceptance ratio
      //log_apt_rto = arma::as_scalar(log(arma::randu())) - log_pri_chn(i) + log_pri_chn(i - 1);
      //log_apt_rto = arma::as_scalar(log(arma::randu())) - log_pri_chn(i) - log_psl_old_new + log_pri_chn(i - 1) + log_psl_new_old;
      runMH_arma(apt, log_lik_chn(i), log_lik_chn(i - 1), sel_cof_chn.col(i), scl_rec_rat, smp_gen, smp_siz, ptl_smp_cnt, dt, ptn_num, pcl_num);
      if (apt == false) {
        sel_cof_chn.col(i) = sel_cof_chn.col(i - 1);
      }
    }
  }
  
  return sel_cof_chn;
}

// Run the parallel adaptive particle marginal Metropolis-Hastings
//[[Rcpp::export]]
arma::dmat runPAPMMH_arma(const arma::dcolvec& sel_cof, const double& rec_rat, const int& pop_siz, const arma::irowvec& smp_gen, const arma::irowvec& smp_siz, const arma::imat& smp_cnt, const arma::imat& mis_cnt, const arma::uword& ptn_num, const arma::uword& pcl_num, const arma::uword& itn_num, const arma::uword nap_num) {
  // ensure RNG gets set/reset
  RNGScope scope;
  
  arma::field<arma::imat> ptl_smp_cnt(smp_gen.n_elem);
  for (arma::uword k = 0; k < smp_gen.n_elem; k++) {
    ptl_smp_cnt(k) = calculatePhasedGenoCount_arma(smp_cnt.col(k), mis_cnt.col(k));
  }
  
  arma::dmat sel_cof_chn = arma::zeros<arma::dmat>(4, itn_num);
  
  //arma::drowvec log_pri_chn = arma::zeros<arma::drowvec>(itn_num);
  arma::drowvec log_lik_chn = arma::zeros<arma::drowvec>(itn_num);
  
  // initialise the population genetic parameters
  cout << "iteration: " << 1 << endl;
  // take the uniform prior and fix the initial value of the selection coefficient to zero
  // or take the beta prior with alpha = 1 and beta = 3
  sel_cof_chn.col(0) = sel_cof;
  
  // rescale the recombination rate
  double scl_rec_rat = 4 * pop_siz * rec_rat;
  
  // declare delta t
  double dt = 1.0 / (2 * pop_siz) / ptn_num;
  
  calculateLogLikelihood_arma(log_lik_chn(0), sel_cof_chn.col(0), scl_rec_rat, smp_gen, smp_siz, ptl_smp_cnt, dt, ptn_num, pcl_num);
  
  double sel_cof_sd = 2e-03;
  arma::dmat sel_cof_cov = arma::eye<arma::dmat>(4, 4) * pow(sel_cof_sd, 2);
  
  bool apt = true;
  // non-adaptation period
  for (arma::uword i = 1; i < nap_num; i++) {
    cout << "iteration: " << i + 1 << endl;
    apt = true;
    
    // draw the candidates of the selection coefficients from the proposals
    sel_cof_chn.col(i) = arma::mvnrnd(sel_cof_chn.col(i - 1), sel_cof_cov);
    if (arma::any(sel_cof_chn.col(i) < -1)) {
      sel_cof_chn.col(i) = sel_cof_chn.col(i - 1);
      log_lik_chn(i) = log_lik_chn(i - 1);
    } else {
      //double log_psl_old_new = log(arma::normpdf(sel_cof_chn.col(i - 1), sel_cof_chn.col(i), sel_cof_cov));
      //double log_psl_new_old = log(arma::normpdf(sel_cof_chn.col(i), sel_cof_chn.col(i - 1), sel_cof_cov));
      
      // calculate the acceptance ratio
      //log_apt_rto = arma::as_scalar(log(arma::randu())) - log_pri_chn(i) + log_pri_chn(i - 1);
      //log_apt_rto = arma::as_scalar(log(arma::randu())) - log_pri_chn(i) - log_psl_old_new + log_pri_chn(i - 1) + log_psl_new_old;
      runMH_arma(apt, log_lik_chn(i), log_lik_chn(i - 1), sel_cof_chn.col(i), scl_rec_rat, smp_gen, smp_siz, ptl_smp_cnt, dt, ptn_num, pcl_num);
      if (apt == false) {
        sel_cof_chn.col(i) = sel_cof_chn.col(i - 1);
      }
    }
  }
  
  // adaptation period
  arma::dmat sel_cof_mat_pre = arma::mean(sel_cof_chn.cols(0, nap_num - 2), 1) * arma::trans(arma::mean(sel_cof_chn.cols(0, nap_num - 2), 1));
  arma::dmat sel_cof_mat_pst = arma::mean(sel_cof_chn.cols(0, nap_num - 1), 1) * arma::trans(arma::mean(sel_cof_chn.cols(0, nap_num - 1), 1));
  for (arma::uword i = nap_num; i < itn_num; i++) {
    cout << "iteration: " << i + 1 << endl;
    apt = true;
    
    // draw the candidates of the selection coefficients from the proposals
    sel_cof_cov = (1 - 1.0 / (i - 1)) * sel_cof_cov + (2.88 / (i - 1)) * ((i - 2) * sel_cof_mat_pre - (i - 1) * sel_cof_mat_pst + sel_cof_chn.col(i - 1) * arma::trans(sel_cof_chn.col(i - 1)) + 1e-06 * arma::eye<arma::dmat>(4, 4));
    sel_cof_chn.col(i) = arma::mvnrnd(sel_cof_chn.col(i - 1), sel_cof_cov);
    if (arma::any(sel_cof_chn.col(i) < -1)) {
      sel_cof_chn.col(i) = sel_cof_chn.col(i - 1);
      log_lik_chn(i) = log_lik_chn(i - 1);
    } else {
      //double log_psl_old_new = log(arma::normpdf(sel_cof_chn.col(i - 1), sel_cof_chn.col(i), sel_cof_cov));
      //double log_psl_new_old = log(arma::normpdf(sel_cof_chn.col(i), sel_cof_chn.col(i - 1), sel_cof_cov));
      
      // calculate the acceptance ratio
      //log_apt_rto = arma::as_scalar(log(arma::randu())) - log_pri_chn(i) + log_pri_chn(i - 1);
      //log_apt_rto = arma::as_scalar(log(arma::randu())) - log_pri_chn(i) - log_psl_old_new + log_pri_chn(i - 1) + log_psl_new_old;
      runMH_arma(apt, log_lik_chn(i), log_lik_chn(i - 1), sel_cof_chn.col(i), scl_rec_rat, smp_gen, smp_siz, ptl_smp_cnt, dt, ptn_num, pcl_num);
      if (apt == false) {
        sel_cof_chn.col(i) = sel_cof_chn.col(i - 1);
      }
    }
    sel_cof_mat_pre = sel_cof_mat_pst;
    sel_cof_mat_pst = arma::mean(sel_cof_chn.cols(0, i), 1) * arma::trans(arma::mean(sel_cof_chn.cols(0, i), 1));
  }
  
  return sel_cof_chn;
}

// Run the early rejection Metropolis-Hastings
// [[Rcpp::export]]
void runERMH_arma(bool& apt, double& log_lik, const double& log_lik_pre, const arma::dcolvec& sel_cof, const double& scl_rec_rat, const arma::irowvec& smp_gen, const arma::irowvec& smp_siz, const arma::field<arma::imat>& ptl_smp_cnt, const double& dt, const arma::uword& ptn_num, const arma::uword& pcl_num) {
  // ensure RNG gets set/reset
  RNGScope scope;
  
  double log_apt_rto = log(arma::as_scalar(arma::randu()));
  
  arma::dcolvec wght = arma::zeros<arma::dcolvec>(pcl_num);
  arma::dmat part_pre = arma::zeros<arma::dmat>(10, pcl_num);
  arma::dmat part_pst = arma::zeros<arma::dmat>(10, pcl_num);
  
  arma::uword dom_ind = arma::as_scalar(arma::find(smp_siz == 0));
  
  arma::dmat fts_bd = arma::ones<arma::dmat>(4, 4);
  fts_bd(1, 1) = 1 + sel_cof(1);
  fts_bd(3, 1) = 1 + sel_cof(1);
  fts_bd(2, 2) = 1 + sel_cof(0);
  fts_bd(3, 2) = 1 + sel_cof(0);
  fts_bd(1, 3) = 1 + sel_cof(1);
  fts_bd(2, 3) = 1 + sel_cof(0);
  fts_bd(3, 3) = 1 + sel_cof(1);
  
  arma::dmat fts_ad = arma::ones<arma::dmat>(4, 4);
  fts_ad(1, 1) = 1 + sel_cof(3);
  fts_ad(3, 1) = 1 + sel_cof(3);
  fts_ad(2, 2) = 1 + sel_cof(2);
  fts_ad(3, 2) = 1 + sel_cof(2);
  fts_ad(1, 3) = 1 + sel_cof(3);
  fts_ad(2, 3) = 1 + sel_cof(2);
  fts_ad(3, 3) = 1 + sel_cof(3);
  
  // initialise the particles
  part_pre = arma::normalise(arma::randu<arma::dmat>(10, pcl_num), 1, 0);
  arma::imat smp_cnt = ptl_smp_cnt(0);
  for (arma::uword i = 0; i < pcl_num; i++) {
    for (arma::uword j = 0; j < smp_cnt.n_cols; j++) {
      wght(i) = wght(i) + calculateMultinomialProb_arma(smp_cnt.col(j), smp_siz(0), part_pre.col(i));
    }
  }
  
  if (arma::any(wght > 0)) {
    log_lik = log(arma::mean(wght) / smp_cnt.n_cols);
    if (log_apt_rto > log_lik - log_lik_pre) {
      log_lik = log_lik_pre;
      apt = false;
      //cout << "early rejection: " << 1 << endl;
      return;
    } else {
      arma::dcolvec prob = arma::normalise(wght, 1);
      arma::ucolvec elem = arma::linspace<arma::ucolvec>(0, pcl_num - 1, pcl_num);
      arma::ucolvec indx = RcppArmadillo::sample(elem, pcl_num, TRUE, prob);
      part_pst = part_pre.cols(indx);
    }
  } else {
    log_lik = log_lik_pre;
    apt = false;
    //cout << "early rejection: " << 1 << endl;
    return;
  }
  
  // run the bootstrap particle filter
  for (arma::uword k = 1; k < dom_ind; k++) {
    wght = arma::zeros<arma::dcolvec>(pcl_num);
    arma::imat smp_cnt = ptl_smp_cnt(k);
    for (arma::uword i = 0; i < pcl_num; i++) {
      part_pre.col(i) = simulateWFD_arma(fts_bd, scl_rec_rat, part_pst.col(i), dt, smp_gen(k - 1), smp_gen(k), ptn_num);
      for (arma::uword j = 0; j < smp_cnt.n_cols; j++) {
        wght(i) = wght(i) + calculateMultinomialProb_arma(smp_cnt.col(j), smp_siz(k), part_pre.col(i));
      }
    }
    
    if (arma::any(wght > 0)) {
      log_lik = log_lik + log(arma::mean(wght) / smp_cnt.n_cols);
      if (log_apt_rto > log_lik - log_lik_pre) {
        log_lik = log_lik_pre;
        apt = false;
        //cout << "early rejection: " << k + 1 << endl;
        return;
      } else {
        arma::dcolvec prob = arma::normalise(wght, 1);
        arma::ucolvec elem = arma::linspace<arma::ucolvec>(0, pcl_num - 1, pcl_num);
        arma::ucolvec indx = RcppArmadillo::sample(elem, pcl_num, TRUE, prob);
        part_pst = part_pre.cols(indx);
      }
    } else {
      log_lik = log_lik_pre;
      apt = false;
      //cout << "early rejection: " << k + 1 << endl;
      return;
    }
  }
  
  if (smp_gen(dom_ind) != smp_gen(dom_ind - 1)) {
    for (arma::uword i = 0; i < pcl_num; i++) {
      part_pre.col(i) = simulateWFD_arma(fts_bd, scl_rec_rat, part_pst.col(i), dt, smp_gen(dom_ind - 1), smp_gen(dom_ind), ptn_num);
    }
    part_pst = part_pre;
  }
  
  for (arma::uword k = dom_ind + 1; k < smp_gen.n_elem; k++) {
    wght = arma::zeros<arma::dcolvec>(pcl_num);
    arma::imat smp_cnt = ptl_smp_cnt(k);
    for (arma::uword i = 0; i < pcl_num; i++) {
      part_pre.col(i) = simulateWFD_arma(fts_ad, scl_rec_rat, part_pst.col(i), dt, smp_gen(k - 1), smp_gen(k), ptn_num);
      for (arma::uword j = 0; j < smp_cnt.n_cols; j++) {
        wght(i) = wght(i) + calculateMultinomialProb_arma(smp_cnt.col(j), smp_siz(k), part_pre.col(i));
      }
    }
    
    if (arma::any(wght > 0)) {
      log_lik = log_lik + log(arma::mean(wght) / smp_cnt.n_cols);
      if (log_apt_rto > log_lik - log_lik_pre) {
        log_lik = log_lik_pre;
        apt = false;
        //cout << "early rejection: " << k << endl;
        return;
      } else {
        arma::dcolvec prob = arma::normalise(wght, 1);
        arma::ucolvec elem = arma::linspace<arma::ucolvec>(0, pcl_num - 1, pcl_num);
        arma::ucolvec indx = RcppArmadillo::sample(elem, pcl_num, TRUE, prob);
        part_pst = part_pre.cols(indx);
      }
    } else {
      log_lik = log_lik_pre;
      apt = false;
      //cout << "early rejection: " << k << endl;
      return;
    }
  }
  
  if (log_apt_rto > log_lik - log_lik_pre) {
    log_lik = log_lik_pre;
    apt = false;
  } else {
    apt = true;
  }
}

// Run the early rejection particle marginal Metropolis-Hastings
//[[Rcpp::export]]
arma::dmat runERPMMH_arma(const arma::dcolvec& sel_cof, const double& rec_rat, const int& pop_siz, const arma::irowvec& smp_gen, const arma::irowvec& smp_siz, const arma::imat& smp_cnt, const arma::imat& mis_cnt, const arma::uword& ptn_num, const arma::uword& pcl_num, const arma::uword& itn_num) {
  // ensure RNG gets set/reset
  RNGScope scope;
  
  arma::field<arma::imat> ptl_smp_cnt(smp_gen.n_elem);
  for (arma::uword k = 0; k < smp_gen.n_elem; k++) {
    ptl_smp_cnt(k) = calculatePhasedGenoCount_arma(smp_cnt.col(k), mis_cnt.col(k));
  }
  
  arma::dmat sel_cof_chn = arma::zeros<arma::dmat>(4, itn_num);
  
  //arma::drowvec log_pri_chn = arma::zeros<arma::drowvec>(itn_num);
  arma::drowvec log_lik_chn = arma::zeros<arma::drowvec>(itn_num);
  
  // initialise the population genetic parameters
  cout << "iteration: " << 1 << endl;
  // take the uniform prior and fix the initial value of the selection coefficient to zero
  // or take the beta prior with alpha = 1 and beta = 3
  sel_cof_chn.col(0) = sel_cof;
  
  // rescale the recombination rate
  double scl_rec_rat = 4 * pop_siz * rec_rat;
  
  // declare delta t
  double dt = 1.0 / (2 * pop_siz) / ptn_num;
  
  calculateLogLikelihood_arma(log_lik_chn(0), sel_cof_chn.col(0), scl_rec_rat, smp_gen, smp_siz, ptl_smp_cnt, dt, ptn_num, pcl_num);
  
  double sel_cof_sd = 2e-03;
  arma::dmat sel_cof_cov = arma::eye<arma::dmat>(4, 4) * pow(sel_cof_sd, 2);
  
  bool apt = true;
  for (arma::uword i = 1; i < itn_num; i++) {
    cout << "iteration: " << i + 1 << endl;
    apt = true;
    
    // draw the candidates of the selection coefficients from the proposals
    sel_cof_chn.col(i) = arma::mvnrnd(sel_cof_chn.col(i - 1), sel_cof_cov);
    if (arma::any(sel_cof_chn.col(i) < -1)) {
      sel_cof_chn.col(i) = sel_cof_chn.col(i - 1);
      log_lik_chn(i) = log_lik_chn(i - 1);
    } else {
      //double log_psl_old_new = log(arma::normpdf(sel_cof_chn.col(i - 1), sel_cof_chn.col(i), sel_cof_cov));
      //double log_psl_new_old = log(arma::normpdf(sel_cof_chn.col(i), sel_cof_chn.col(i - 1), sel_cof_cov));
      
      // calculate the acceptance ratio
      //log_apt_rto = arma::as_scalar(log(arma::randu())) - log_pri_chn(i) + log_pri_chn(i - 1);
      //log_apt_rto = arma::as_scalar(log(arma::randu())) - log_pri_chn(i) - log_psl_old_new + log_pri_chn(i - 1) + log_psl_new_old;
      runERMH_arma(apt, log_lik_chn(i), log_lik_chn(i - 1), sel_cof_chn.col(i), scl_rec_rat, smp_gen, smp_siz, ptl_smp_cnt, dt, ptn_num, pcl_num);
      if (apt == false) {
        sel_cof_chn.col(i) = sel_cof_chn.col(i - 1);
      }
    }
  }
  
  return sel_cof_chn;
}

// Run the early rejection parallel adaptive particle marginal Metropolis-Hastings
//[[Rcpp::export]]
arma::dmat runERPAPMMH_arma(const arma::dcolvec& sel_cof, const double& rec_rat, const int& pop_siz, const arma::irowvec& smp_gen, const arma::irowvec& smp_siz, const arma::imat& smp_cnt, const arma::imat& mis_cnt, const arma::uword& ptn_num, const arma::uword& pcl_num, const arma::uword& itn_num, const arma::uword nap_num) {
  // ensure RNG gets set/reset
  RNGScope scope;
  
  arma::field<arma::imat> ptl_smp_cnt(smp_gen.n_elem);
  for (arma::uword k = 0; k < smp_gen.n_elem; k++) {
    ptl_smp_cnt(k) = calculatePhasedGenoCount_arma(smp_cnt.col(k), mis_cnt.col(k));
  }
  
  arma::dmat sel_cof_chn = arma::zeros<arma::dmat>(4, itn_num);
  
  //arma::drowvec log_pri_chn = arma::zeros<arma::drowvec>(itn_num);
  arma::drowvec log_lik_chn = arma::zeros<arma::drowvec>(itn_num);
  
  // initialise the population genetic parameters
  cout << "iteration: " << 1 << endl;
  // take the uniform prior and fix the initial value of the selection coefficient to zero
  // or take the beta prior with alpha = 1 and beta = 3
  sel_cof_chn.col(0) = sel_cof;
  
  // rescale the recombination rate
  double scl_rec_rat = 4 * pop_siz * rec_rat;
  
  // declare delta t
  double dt = 1.0 / (2 * pop_siz) / ptn_num;
  
  calculateLogLikelihood_arma(log_lik_chn(0), sel_cof_chn.col(0), scl_rec_rat, smp_gen, smp_siz, ptl_smp_cnt, dt, ptn_num, pcl_num);
  
  double sel_cof_sd = 2e-03;
  arma::dmat sel_cof_cov = arma::eye<arma::dmat>(4, 4) * pow(sel_cof_sd, 2);
  
  bool apt = true;
  // non-adaptation period
  for (arma::uword i = 1; i < nap_num; i++) {
    cout << "iteration: " << i + 1 << endl;
    apt = true;
    
    // draw the candidates of the selection coefficients from the proposals
    sel_cof_chn.col(i) = arma::mvnrnd(sel_cof_chn.col(i - 1), sel_cof_cov);
    if (arma::any(sel_cof_chn.col(i) < -1)) {
      sel_cof_chn.col(i) = sel_cof_chn.col(i - 1);
      log_lik_chn(i) = log_lik_chn(i - 1);
    } else {
      //double log_psl_old_new = log(arma::normpdf(sel_cof_chn.col(i - 1), sel_cof_chn.col(i), sel_cof_cov));
      //double log_psl_new_old = log(arma::normpdf(sel_cof_chn.col(i), sel_cof_chn.col(i - 1), sel_cof_cov));
      
      // calculate the acceptance ratio
      //log_apt_rto = arma::as_scalar(log(arma::randu())) - log_pri_chn(i) + log_pri_chn(i - 1);
      //log_apt_rto = arma::as_scalar(log(arma::randu())) - log_pri_chn(i) - log_psl_old_new + log_pri_chn(i - 1) + log_psl_new_old;
      runERMH_arma(apt, log_lik_chn(i), log_lik_chn(i - 1), sel_cof_chn.col(i), scl_rec_rat, smp_gen, smp_siz, ptl_smp_cnt, dt, ptn_num, pcl_num);
      if (apt == false) {
        sel_cof_chn.col(i) = sel_cof_chn.col(i - 1);
      }
    }
  }
  
  // adaptation period
  arma::dmat sel_cof_mat_pre = arma::mean(sel_cof_chn.cols(0, nap_num - 2), 1) * arma::trans(arma::mean(sel_cof_chn.cols(0, nap_num - 2), 1));
  arma::dmat sel_cof_mat_pst = arma::mean(sel_cof_chn.cols(0, nap_num - 1), 1) * arma::trans(arma::mean(sel_cof_chn.cols(0, nap_num - 1), 1));
  for (arma::uword i = nap_num; i < itn_num; i++) {
    cout << "iteration: " << i + 1 << endl;
    apt = true;
    
    // draw the candidates of the selection coefficients from the proposals
    sel_cof_cov = (1 - 1.0 / (i - 1)) * sel_cof_cov + (2.88 / (i - 1)) * ((i - 2) * sel_cof_mat_pre - (i - 1) * sel_cof_mat_pst + sel_cof_chn.col(i - 1) * arma::trans(sel_cof_chn.col(i - 1)) + 1e-06 * arma::eye<arma::dmat>(4, 4));
    sel_cof_chn.col(i) = arma::mvnrnd(sel_cof_chn.col(i - 1), sel_cof_cov);
    if (arma::any(sel_cof_chn.col(i) < -1)) {
      sel_cof_chn.col(i) = sel_cof_chn.col(i - 1);
      log_lik_chn(i) = log_lik_chn(i - 1);
    } else {
      //double log_psl_old_new = log(arma::normpdf(sel_cof_chn.col(i - 1), sel_cof_chn.col(i), sel_cof_cov));
      //double log_psl_new_old = log(arma::normpdf(sel_cof_chn.col(i), sel_cof_chn.col(i - 1), sel_cof_cov));
      
      // calculate the acceptance ratio
      //log_apt_rto = arma::as_scalar(log(arma::randu())) - log_pri_chn(i) + log_pri_chn(i - 1);
      //log_apt_rto = arma::as_scalar(log(arma::randu())) - log_pri_chn(i) - log_psl_old_new + log_pri_chn(i - 1) + log_psl_new_old;
      runERMH_arma(apt, log_lik_chn(i), log_lik_chn(i - 1), sel_cof_chn.col(i), scl_rec_rat, smp_gen, smp_siz, ptl_smp_cnt, dt, ptn_num, pcl_num);
      if (apt == false) {
        sel_cof_chn.col(i) = sel_cof_chn.col(i - 1);
      }
    }
    sel_cof_mat_pre = sel_cof_mat_pst;
    sel_cof_mat_pst = arma::mean(sel_cof_chn.cols(0, i), 1) * arma::trans(arma::mean(sel_cof_chn.cols(0, i), 1));
  }
  
  return sel_cof_chn;
}
/*************************/
