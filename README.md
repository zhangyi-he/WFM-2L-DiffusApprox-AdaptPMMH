# Estimating temporally variable selection intensity from ancient DNA data with the flexibility of modelling linkage and epistasis
The source code is implemented for an MCMC-based method for inferring temporally variable selection from ancient DNA data (in the format of genotype likelihoods) with the flexibility of modelling linkage and epistasis, and the preprint has been submitted to Molecular Ecology Resources, available at https://doi.org/10.1101/2022.08.02.502360.

[Code v1.0](https://github.com/zhangyi-he/WFM-2L-DiffusApprox-AdaptPMMH/tree/master/Code%20v1.0) / [Code v1.0.1](https://github.com/zhangyi-he/WFM-2L-DiffusApprox-AdaptPMMH/tree/master/Code%20v1.0.1) includes the source code implemented for horse base coat colours (ASIP and MC1R with epistasis) and horse pinto coat patterns (KIT13 and KIT16 with linkage), where the linkage disequilibrium is fixed to be 0 for horse base coat colours and the selection coefficient of mixed against solid is fixed to be -1 for horse pinto coat patterns. Uniform priors for population mutant allele frequencies and linkage disequilibrium are adopted. Note that [Code v1.0.1](https://github.com/zhangyi-he/WFM-2L-DiffusApprox-AdaptPMMH/tree/master/Code%20v1.0.1) provides a version that no loss or fixation event occurred in the underlying population.

[Code v1.1](https://github.com/zhangyi-he/WFM-2L-DiffusApprox-AdaptPMMH/tree/master/Code%20v1.1) includes the source code implemented for horse base coat colours (ASIP and MC1R with epistasis) and horse pinto coat patterns (KIT13 and KIT16 with linkage). Uniform priors for population mutant allele frequencies and linkage disequilibrium are adopted.

[Code v1.2](https://github.com/zhangyi-he/WFM-2L-DiffusApprox-AdaptPMMH/tree/master/Code%20v1.2) / [Code v1.2.1](https://github.com/zhangyi-he/WFM-2L-DiffusApprox-AdaptPMMH/tree/master/Code%20v1.2.1) includes the source code implemented for horse base coat colours (ASIP and MC1R with epistasis) and horse pinto coat patterns (KIT13 and KIT16 with linkage). A Dirichlet flat prior for population gamete frequencies is adopted. Note that [Code v1.2.1](https://github.com/zhangyi-he/WFM-2L-DiffusApprox-AdaptPMMH/tree/master/Code%20v1.2.1) provides a version that the selection coefficient of mixed against solid is fixed to be -1 for horse pinto coat patterns.

[Data](https://github.com/zhangyi-he/WFM-2L-DiffusApprox-AdaptPMMH/tree/master/Data) includes the ancient horse samples genotyped at the loci for coat colouration.

Note that source code currently available on GitHub is implemented for horse base coat colours and pinto coat patterns. We will release a version that allows for a user-specified genotype-phenotype map soon.
