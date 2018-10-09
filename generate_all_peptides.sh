#!/usr/bin/env bash

# Generate peptides for Round 1 sfp with unlabeling
Rscript gen_recommendation.R 0 267 sfp null PfAcpH round1.csv

# Generate peptides for Round 2 sfp with unlabeling
Rscript gen_recommendation.R 1 580 sfp null PfAcpH round2.csv

# Generate peptides for Round 3 sfp specific with unlabeling
Rscript gen_recommendation.R 2 285 sfp AcpS PfAcpH round3_sfp.csv

# Generate peptides for Round 3 AcpS specific with unlabeling
Rscript gen_recommendation.R 2 285 AcpS sfp PfAcpH round3_AcpS.csv

# Generate peptides for Round 4 sfp specific with unlabeling
Rscript gen_recommendation.R 3 300 sfp AcpS PfAcpH round4_sfp.csv

# Generate peptides for Round 4 AcpS specific with unlabeling
Rscript gen_recommendation.R 3 300 AcpS sfp PfAcpH round4_AcpS.csv

# Generate peptides for Round 5 sfp specific with unlabeling
Rscript gen_recommendation.R 4 300 sfp AcpS PfAcpH round5_sfp.csv

# Generate peptides for Round 5 AcpS specific with unlabeling
Rscript gen_recommendation.R 4 300 AcpS sfp PfAcpH round5_AcpS.csv
