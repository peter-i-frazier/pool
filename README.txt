This repository contains code and data for the paper, "Discovering de novo peptide substrates for enzymes using machine learning," by Lorillee Tallorin, JiaLei Wang, Woojoo E. Kim, Swagat Sahu, Nicolas M. Kosa, Matthew Thompsona, Michael K. Gilson, Peter I. Frazier, Nathan C. Gianneschia, and Michael D. Burkart.  

The paper proposes an optimal learning approach for optimizing peptides, called Peptide Optimization with Optimal Learning or POOL.  POOL's goal is to produce sets of peptides to test ("recommendations") to find short peptides that have a desired chemical activity together with each of a collection of enzymes.  We call these peptides "short hits".  POOL proceeds iteratively, consuming a dataset of the results from past experiments and producing a recommendation of which set of peptides to produce next.

This repository is organized as follows:

data/
Contains data produced by experiments done in the process of writing this paper, in both raw and processed form.  It also contains data taken from the literature used to initalize POOL in the first round

recommendation/
Contains recommendations produced by the POOL algorithm during the course of the optimizatin process.

src/
Contains the software implementation of the POOL algorithm in src/core and an earlier implementation in src/archive_code. Also contains custom software written to normalize images of membranes created during experiments (src/normalization). Also contains code for making plots, both in the paper (src/make_plots_for_paper) and in presentations (src/make_plots_for_presentation).
