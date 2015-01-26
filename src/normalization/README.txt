Jan 16, 2015

The raw data we are going to normalize and analyze is from root_path/data/all_TS_raw_reading.csv

script.R showed the latest method to do normalization and make histogram plots

May 28, 2014.

This directory contains files that Jialei and Peter put together to address a
question about the right normalization method for experimental data gathered by
Lori.  She had used one particular normalization scheme to prepare the slide
presentation presentations/2014-05-15_BAC Update.pptx (gzipped in the
repository).  This normalization scheme did not properly shift everything to
the same scale, and so Jialei proposed a new normalization scheme. 
The results are summarized in presentations/orthogonal_peptide_result_analysis.pptx

To run Jialei's analysis, run the script README.R.

This script creates the following files:

sfp_diff.csv
sfp_category.csv
AcpS_diff.csv
AcpS_category.csv
Lori_result_sfp_analysis.csv
Lori_result_AcpS_analysis.csv

The first four files contain peptides specifically labeled by sfp and AcpS,
according to two different normalization schemes -- "diff" and "category".  The
one that we recommend is "diff".

The last two files, Lori_result_sfp_analysis.csv and
Lori_result_AcpS_analysis.csv, contain information about each of the peptides
that Lori listed as being labeled by sfp.  The column "unlabeling" is whether
or not Lori listed it as being unlabeled.  Then it contains normalized values
for sfp labeling, sfp unlabeling, AcpS labeling, and AcpS unlabeling.

presentation_preparation/
contains files that Peter used to prepare histograms of normalized light
intensity using Lori's normalization scheme, included in
orthogonal_peptide_result_analysis.pptx


analysis.R contains functions used by README.R
