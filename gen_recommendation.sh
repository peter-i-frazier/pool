#!/usr/bin/env bash

# Generate peptides based on training set 0 for sfp with unlabeling
Rscript gen_recommendation.R 0 267 sfp null PfAcpH round1.csv

# Generate peptides based on training sets 0-1 for sfp with unlabeling
Rscript gen_recommendation.R 1 580 sfp null PfAcpH round2.csv

# Generate peptides based on training sets 0-2 for sfp specific with unlabeling
Rscript gen_recommendation.R 2 285 sfp AcpS PfAcpH round3_sfp.csv

# Generate peptides based on training sets 0-2 for AcpS specific with unlabeling
Rscript gen_recommendation.R 2 285 AcpS sfp PfAcpH round3_AcpS.csv

# Generate peptides based on training sets 0-3 for sfp specific with unlabeling
Rscript gen_recommendation.R 3 300 sfp AcpS PfAcpH round4_sfp.csv

# Generate peptides based on training sets 0-3 or AcpS specific with unlabeling
Rscript gen_recommendation.R 3 300 AcpS sfp PfAcpH round4_AcpS.csv

# Generate peptides based on training sets 0-4 for sfp specific with unlabeling
Rscript gen_recommendation.R 4 300 sfp AcpS PfAcpH round5_sfp.csv

# Generate peptides based on training sets 0-4 for AcpS specific with unlabeling
Rscript gen_recommendation.R 4 300 AcpS sfp PfAcpH round5_AcpS.csv
