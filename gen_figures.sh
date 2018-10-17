#!/usr/bin/env bash
Rscript gen_roc_data.R

Rscript gen_benchmark.R 100 sfp AcpS sfp_specific
Rscript gen_benchmark.R 100 AcpS sfp acps_specific

python gen_simulation_coord.py sfp
python gen_simulation_coord.py acps

python gen_figures.py
