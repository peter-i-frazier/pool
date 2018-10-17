This repository contains code and data for the paper, "Discovering de novo peptide substrates for enzymes using machine learning," by Lorillee Tallorin, JiaLei Wang, Woojoo E. Kim, Swagat Sahu, Nicolas M. Kosa, Pu Yang, Matthew Thompsona, Michael K. Gilson, Peter I. Frazier, Nathan C. Gianneschia, and Michael D. Burkart.  

The paper proposes an optimal learning approach for optimizing peptides, called Peptide Optimization with Optimal Learning or POOL.  POOL's goal is to produce sets of peptides to test ("recommendations") to find short peptides that have a desired chemical activity together with each of a collection of enzymes.  We call these peptides "short hits".  POOL proceeds iteratively, consuming a dataset of the results from past experiments and producing a recommendation of which set of peptides to produce next.

This repository is organized as follows:

generate_all_peptides.sh
Generate all 5 rounds of the peptides using POOL as described in the paper. Note that due to stochasticity of the algorithm, the generated peptides will be different from time to time.  This script uses the following additional script:
gen_recommendation.R

generate_figures.sh
Reproduce figures in the paper.  This script uses the following additional scripts:
  gen_roc_data.R
  gen_benchmark.R
  gen_simulation_coord.py
  make_figures.py


data/
Contains data produced by experiments done in the process of writing this paper, in both raw and processed form.  It also contains data taken from the literature used to initalize POOL in the first round

recommendation/
Contains recommendations produced by the POOL algorithm during the course of the optimization process.

src/
Contains the software implementation of the POOL algorithm in src/core. Also contains custom software written to normalize images of membranes created during experiments (src/normalization).

Requirements to run the code
- R
- MCMCpack (an R package)
- Python 3
