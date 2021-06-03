// Estimating selection coefficients and testing their changes from ancient DNA data
// Xiaoyang Dai, Mark Beaumont, Feng Yu, Zhangyi He

// version 1.1
// Single-gene phenotypes under constant natural selection and non-constant demographic histories
// Time series data of allele frequencies

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
// Calculate the fitness matrix for the single-locus Wright-Fisher model with selection
// [[Rcpp::export]]
arma::dmat calculateFitnessMat_arma(const double& sel_cof, const double& dom_par) {
  // ensure RNG gets set/reset
  RNGScope scope;

  // calculate the fitness matrix
  arma::dmat fts_mat = arma::ones<arma::dmat>(2, 2);
  // fts_mat(0, 0) = 1;
  fts_mat(1, 0) = 1 + sel_cof * dom_par;
  fts_mat(0, 1) = 1 + sel_cof * dom_par;
  fts_mat(1, 1) = 1 + sel_cof;

  // return the fitness matrix for the Wright-Fisher model
  return fts_mat;
}

// Simulate the mutant allele frequency trajectory according to the single-locus Wright-Fisher model with selection
// [[Rcpp::export]]
arma::drowvec simulateWFM_arma(const arma::dmat& fts_mat, const arma::icolvec& pop_siz, const double& int_frq, const int& int_gen, const int& lst_gen) {
  // ensure RNG gets set/reset
  RNGScope scope;

  // declare the mutant allele frequency trajectory
  arma::drowvec ale_frq_pth = arma::zeros<arma::drowvec>(arma::uword(lst_gen - int_gen) + 1);

  // initialise the mutant allele frequency in generation 0
  arma::dcolvec ale_frq = {1 - int_frq, int_frq};
  ale_frq_pth(0) = ale_frq(1);
  arma::dmat gen_frq = fts_mat % (ale_frq * ale_frq.t()) / arma::as_scalar(ale_frq.t() * fts_mat * ale_frq);
  gen_frq = arma::diagmat(gen_frq) + 2 * arma::trimatu(gen_frq, 1);

  // simulate the mutant allele frequency trajectory
  for (arma::uword k = 1; k < arma::uword(lst_gen - int_gen) + 1; k++) {
    // calculate the sampling probabilities
    arma::dcolvec prob = ale_frq;
    prob = ale_frq % (fts_mat * ale_frq) / arma::as_scalar(ale_frq.t() * fts_mat * ale_frq);

    // proceed the Wright-Fisher sampling
    ale_frq(0) = R::rbinom(2 * pop_siz(k), prob(0)) / 2 / pop_siz(k);
    ale_frq(1) = 1 - ale_frq(0);
    ale_frq_pth(k) = ale_frq(1);

    gen_frq = fts_mat % (ale_frq * ale_frq.t()) / arma::as_scalar(ale_frq.t() * fts_mat * ale_frq);
    gen_frq = arma::diagmat(gen_frq) + 2 * arma::trimatu(gen_frq, 1);
  }

  // return the mutant allele frequency trajectory under the Wright-Fisher model
  return ale_frq_pth;
}
/*************************/


/********** WFD **********/
// Simulate the mutant allele frequency trajectory according to the single-locus Wright-Fisher diffusion with selection using the Euler-Maruyama method
// [[Rcpp::export]]
arma::drowvec simulateWFD_arma(const double& sel_cof, const double& dom_par, const arma::icolvec& pop_siz, const int& ref_siz, const double& int_frq, const int& int_gen, const int& lst_gen, const arma::uword& ptn_num) {
  // ensure RNG gets set/reset
  RNGScope scope;

  // rescale the selection coefficient
  double scl_sel_cof = 2 * ref_siz * sel_cof;

  // calculate the ratio of the population size to the reference population size
  arma::dcolvec siz_rto = arma::zeros<arma::dcolvec>(arma::uword(lst_gen - int_gen) * ptn_num);
  for (arma::uword k = 0; k < arma::uword(lst_gen - int_gen); k++) {
    arma::dcolvec siz_rto_tmp = arma::zeros<arma::dcolvec>(ptn_num);
    siz_rto_tmp.fill(pop_siz(k) / double(ref_siz));
    siz_rto.subvec(k * ptn_num, (k + 1) * ptn_num - 1) = siz_rto_tmp;
  }

  // declare delta t
  double dt = 1.0 / (2 * ref_siz) / ptn_num;
  // declare delta W
  arma::drowvec dW = pow(dt, 0.5) * arma::randn<arma::drowvec>(arma::uword(lst_gen - int_gen) * ptn_num);

  // declare the mutant allele frequency trajectory
  arma::drowvec frq_pth = arma::zeros<arma::drowvec>(arma::uword(lst_gen - int_gen) * ptn_num + 1);

  // initialise the mutant allele frequency in generation 0
  frq_pth(0) = int_frq;

  // simulate the mutant allele frequency trajectory
  for (arma::uword t = 1; t < arma::uword(lst_gen - int_gen) * ptn_num + 1; t++) {
    // calculate the drift coefficient
    double mu = scl_sel_cof * frq_pth(t - 1) * (1 - frq_pth(t - 1)) * (dom_par + (1 - 2 * dom_par) * frq_pth(t - 1));

    // calculate the diffusion coefficient
    double sigma = pow(frq_pth(t - 1) * (1 - frq_pth(t - 1)) / siz_rto(t - 1), 0.5);

    // proceed the Euler-Maruyama scheme
    frq_pth(t) = frq_pth(t - 1) + mu * dt + sigma * dW(t - 1);

    // remove the noise from the numerical techniques
    if (frq_pth(t) < 0) {
      frq_pth(t) = 0;
    }
    if (frq_pth(t) > 1) {
      frq_pth(t) = 1;
    }
  }

  // return the mutant allele frequency trajectory under the Wright-Fisher diffusion
  return frq_pth;
}
/*************************/


/********** BPF **********/
// Run the bootstrap particle filter
// [[Rcpp::export]]
List runBPF_arma(const double& sel_cof, const double& dom_par, const arma::icolvec& pop_siz, const int& ref_siz, const arma::irowvec& smp_gen, const arma::irowvec& smp_siz, const arma::irowvec& smp_cnt, const arma::uword& ptn_num, const arma::uword& pcl_num) {
  // ensure RNG gets set/reset
  RNGScope scope;

  double lik = 1;

  arma::dmat wght = arma::zeros<arma::dmat>(pcl_num, smp_gen.n_elem);
  arma::dmat ale_frq_pre = arma::zeros<arma::dmat>(pcl_num, smp_gen.n_elem);
  arma::dmat ale_frq_pst = arma::zeros<arma::dmat>(pcl_num, smp_gen.n_elem);

  arma::dcolvec wght_tmp = arma::zeros<arma::dcolvec>(pcl_num);
  arma::dcolvec ale_frq_tmp = arma::zeros<arma::dcolvec>(pcl_num);

  // initialise the particles
  cout << "generation: " << smp_gen(0) << endl;
  ale_frq_tmp = arma::randu<arma::dcolvec>(pcl_num);
  for (arma::uword i = 0; i < pcl_num; i++) {
    wght_tmp(i) = R::dbinom(smp_cnt(0), smp_siz(0), ale_frq_tmp(i), false);
  }

  if (arma::sum(wght_tmp) > 0) {
    arma::dcolvec prob = arma::normalise(wght_tmp, 1);
    arma::ucolvec elem = arma::linspace<arma::ucolvec>(0, pcl_num - 1, pcl_num);
    arma::ucolvec indx = RcppArmadillo::sample(elem, pcl_num, true, prob);

    lik = lik * arma::mean(wght_tmp);
    wght.col(0) = wght_tmp;
    ale_frq_pre.col(0) = ale_frq_tmp;
    ale_frq_pst.col(0) = ale_frq_tmp.elem(indx);
  } else {
    lik = 0;
    wght.shed_cols(0, smp_gen.n_elem - 1);
    ale_frq_pre.shed_cols(0, smp_gen.n_elem - 1);
    ale_frq_pst.shed_cols(0, smp_gen.n_elem - 1);

    return List::create(Named("lik", lik),
                        Named("wght", wght),
                        Named("ale_frq_pre_resmp", ale_frq_pre),
                        Named("ale_frq_pst_resmp", ale_frq_pst));
  }

  // run the bootstrap particle filter
  for (arma::uword k = 1; k < smp_gen.n_elem; k++) {
    cout << "generation: " << smp_gen(k) << endl;
    wght_tmp = arma::zeros<arma::dcolvec>(pcl_num);
    ale_frq_tmp = ale_frq_pst.col(k - 1);
    for (arma::uword i = 0; i < pcl_num; i++) {
      arma::drowvec path = simulateWFD_arma(sel_cof, dom_par, pop_siz.subvec(smp_gen(k - 1), smp_gen(k)), ref_siz, ale_frq_tmp(i), smp_gen(k - 1), smp_gen(k), ptn_num);
      ale_frq_tmp(i) = arma::as_scalar(path.tail(1));
      wght_tmp(i) = R::dbinom(smp_cnt(k), smp_siz(k), ale_frq_tmp(i), false);
    }

    if (arma::sum(wght_tmp) > 0) {
      arma::dcolvec prob = arma::normalise(wght_tmp, 1);
      arma::ucolvec elem = arma::linspace<arma::ucolvec>(0, pcl_num - 1, pcl_num);
      arma::ucolvec indx = RcppArmadillo::sample(elem, pcl_num, true, prob);

      lik = lik * arma::mean(wght_tmp);
      wght.col(k) = wght_tmp;
      ale_frq_pre.col(k) = ale_frq_tmp;
      ale_frq_pst.col(k) = ale_frq_tmp.elem(indx);
    } else {
      lik = 0;
      wght.shed_cols(k, smp_gen.n_elem - 1);
      ale_frq_pre.shed_cols(k, smp_gen.n_elem - 1);
      ale_frq_pst.shed_cols(k, smp_gen.n_elem - 1);

      break;
    }
  }

  return List::create(Named("lik", lik),
                      Named("wght", wght),
                      Named("ale_frq_pre_resmp", ale_frq_pre),
                      Named("ale_frq_pst_resmp", ale_frq_pst));
}
/*************************/


/********** PMMH **********/
// Calculate the log-likelihood using the bootstrap particle filter
// [[Rcpp::export]]
double calculateLogLikelihood_arma(const double& sel_cof, const double& dom_par, const arma::icolvec& pop_siz, const int& ref_siz, const arma::irowvec& smp_gen, const arma::irowvec& smp_siz, const arma::irowvec& smp_cnt, const arma::uword& ptn_num, const arma::uword& pcl_num) {
  // ensure RNG gets set/reset
  RNGScope scope;

  double log_lik = 0;

  arma::dcolvec wght = arma::zeros<arma::dcolvec>(pcl_num);
  arma::dcolvec ale_frq_pre = arma::zeros<arma::dcolvec>(pcl_num);
  arma::dcolvec ale_frq_pst = arma::zeros<arma::dcolvec>(pcl_num);

  // initialise the particles
  ale_frq_pre = arma::randu<arma::dcolvec>(pcl_num);;
  for (arma::uword i = 0; i < pcl_num; i++) {
    wght(i) = R::dbinom(smp_cnt(0), smp_siz(0), ale_frq_pre(i), false);
  }

  if (arma::mean(wght) > 0) {
    log_lik = log_lik + log(arma::mean(wght));
    arma::dcolvec prob = arma::normalise(wght, 1);
    arma::ucolvec elem = arma::linspace<arma::ucolvec>(0, pcl_num - 1, pcl_num);
    arma::ucolvec indx = RcppArmadillo::sample(elem, pcl_num, true, prob);
    ale_frq_pst = ale_frq_pre.elem(indx);
  } else {
    log_lik = -(arma::datum::inf);

    return log_lik;
  }

  // run the bootstrap particle filter
  for (arma::uword k = 1; k < smp_gen.n_elem; k++) {
    wght = arma::zeros<arma::dcolvec>(pcl_num);
    for (arma::uword i = 0; i < pcl_num; i++) {
      arma::drowvec path = simulateWFD_arma(sel_cof, dom_par, pop_siz.subvec(smp_gen(k - 1), smp_gen(k)), ref_siz, ale_frq_pst(i), smp_gen(k - 1), smp_gen(k), ptn_num);
      ale_frq_pre(i) = arma::as_scalar(path.tail(1));
      wght(i) = R::dbinom(smp_cnt(k), smp_siz(k), ale_frq_pre(i), false);
    }

    if (arma::mean(wght) > 0) {
      log_lik = log_lik + log(arma::mean(wght));
      arma::dcolvec prob = arma::normalise(wght, 1);
      arma::ucolvec elem = arma::linspace<arma::ucolvec>(0, pcl_num - 1, pcl_num);
      arma::ucolvec indx = RcppArmadillo::sample(elem, pcl_num, true, prob);
      ale_frq_pst = ale_frq_pre.elem(indx);
    } else {
      log_lik = -(arma::datum::inf);

      return log_lik;
    }
  }

  return log_lik;
}

// Calculate the optimal particle number in the particle marginal Metropolis-Hastings
// [[Rcpp::export]]
List calculateOptimalParticleNum_arma(const double& sel_cof, const double& dom_par, const arma::icolvec& pop_siz, const int& ref_siz, const arma::irowvec& smp_gen, const arma::irowvec& smp_siz, const arma::irowvec& smp_cnt, const arma::uword& ptn_num, const arma::uword& pcl_num, const arma::uword& gap_num) {
  // ensure RNG gets set/reset
  RNGScope scope;

  arma::drowvec log_lik(300);
  for (arma::uword i = 0; i < 300; i++) {
    log_lik(i) = calculateLogLikelihood_arma(sel_cof, dom_par, pop_siz, ref_siz, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num);
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
        log_lik(i) = calculateLogLikelihood_arma(sel_cof, dom_par, pop_siz, ref_siz, smp_gen, smp_siz, smp_cnt, ptn_num, opt_pcl_num(0));
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
        log_lik(i) = calculateLogLikelihood_arma(sel_cof, dom_par, pop_siz, ref_siz, smp_gen, smp_siz, smp_cnt, ptn_num, opt_pcl_num(0));
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
        log_lik(i) = calculateLogLikelihood_arma(sel_cof, dom_par, pop_siz, ref_siz, smp_gen, smp_siz, smp_cnt, ptn_num, opt_pcl_num(0));
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
        log_lik(i) = calculateLogLikelihood_arma(sel_cof, dom_par, pop_siz, ref_siz, smp_gen, smp_siz, smp_cnt, ptn_num, opt_pcl_num(0));
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
arma::drowvec runPMMH_arma(const double& sel_cof, const double& dom_par, const arma::icolvec& pop_siz, const int& ref_siz, const arma::irowvec& smp_gen, const arma::irowvec& smp_siz, const arma::irowvec& smp_cnt, const arma::uword& ptn_num, const arma::uword& pcl_num, const arma::uword& itn_num) {
  // ensure RNG gets set/reset
  RNGScope scope;

  arma::drowvec sel_cof_chn = arma::zeros<arma::drowvec>(itn_num);

  //arma::drowvec log_pri_chn = arma::zeros<arma::drowvec>(2);
  arma::drowvec log_lik = arma::zeros<arma::drowvec>(2);

  double sel_cof_sd = 5e-03;

  // initialise the population genetic parameters
  cout << "iteration: " << 1 << endl;
  // take the uniform prior and fix the initial value of the selection coefficient to zero
  // or take the beta prior with alpha = 1 and beta = 3
  sel_cof_chn(0) = sel_cof;

  log_lik(0) = calculateLogLikelihood_arma(sel_cof_chn(0), dom_par, pop_siz, ref_siz, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num);

  double apt_cnt = 0;
  for (arma::uword i = 1; i < itn_num; i++) {
    cout << "iteration: " << i + 1 << endl;

    // draw the candidate of the selection coefficient from the random walk proposal
    sel_cof_chn(i) = sel_cof_chn(i - 1) + sel_cof_sd * arma::randn();

    if (sel_cof_chn(i) < -1) {
      sel_cof_chn(i) = sel_cof_chn(i - 1);
      log_lik(1) = log_lik(0);
      // apt_cnt = apt_cnt + 0;
      cout << "acceptance: " << apt_cnt / i << endl;
    } else {
      // calculate the proposal
      // arma::drowvec log_psl = arma::zeros<arma::drowvec>(2);

      // calculate the likelihood
      log_lik(1) = calculateLogLikelihood_arma(sel_cof_chn(i), dom_par, pop_siz, ref_siz, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num);

      // calculate the acceptance ratio
      // double apt_rto = exp(log_lik(1) - log_lik(0));
      // double apt_rto = exp((log_pri(1) + log_lik(1) + log_psl(1)) - (log_pri(0) + log_lik(0) + log_psl(0)));

      // calculate the acceptance ratio
      // double apt_rto = exp(log_lik(1) - log_lik(0));
      // double apt_rto = exp((log_pri(1) + log_lik(1) + log_psl(1)) - (log_pri(0) + log_lik(0) + log_psl(0)));

      if (arma::randu() > exp(log_lik(1) - log_lik(0))) {
        sel_cof_chn(i) = sel_cof_chn(i - 1);
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
/*************************/