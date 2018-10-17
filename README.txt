This repository contains code and data for the paper, "Discovering de novo peptide substrates for enzymes using machine learning," by Lorillee Tallorin, JiaLei Wang, Woojoo E. Kim, Swagat Sahu, Nicolas M. Kosa, Pu Yang, Matthew Thompsona, Michael K. Gilson, Peter I. Frazier, Nathan C. Gianneschia, and Michael D. Burkart.  

The paper proposes an optimal learning approach for optimizing peptides, called Peptide Optimization with Optimal Learning or POOL.  POOL's goal is to produce sets of peptides to test ("recommendations") to find short peptides that have a desired chemical activity together with each of a collection of enzymes.  We call these peptides "short hits".  POOL proceeds iteratively, consuming a dataset of the results from past experiments and producing a recommendation of which set of peptides to produce next.

This repository is organized as follows:

gen_recommendation.sh
Generate recommendations of peptides to test using POOL, based on each round of training data, and using the discovery goals used in the experiments conducted for the paper.  In the experiments conducted for the paper, recommendations were obtained from several variants of POOL (and mutation in the early rounds) as described in the Supplemental Information.  This script gives an example of how to generate peptide using a single variant of POOL, which samples theta from its posterior distribution and uses three classifiers combined together (one for each enzyme) rather than a single classifier.  This script uses the following additional script:
gen_recommendation.R

gen_figures.sh
Reproduce figures in the paper.  This script uses the following additional scripts:
  gen_roc_data.R
  gen_benchmark.R
  gen_simulation_coord.py
  gen_figures.py


data/
Contains data produced by experiments done in the process of writing this paper.  It also contains data taken from the literature used to initalize POOL in the first round.  In this folder, "TS" indicates "training set" (i.e., rounds of recommendations), with larger numbers indicating later training sets.

recommendation/
Contains recommendations produced by the POOL algorithm during the course of the optimization process.

src/
Contains the software implementation of the POOL algorithm in src/core. Also contains custom software written to normalize images of membranes created during experiments (src/normalization).

Requirements to run the code
- R
- MCMCpack (an R package)
- Python 3
