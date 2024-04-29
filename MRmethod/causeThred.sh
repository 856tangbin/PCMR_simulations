# ----------------------
# Copyright 2023 PMG Lab
# Author: bintang
# Licence: MIT
# Version: 20240429
# ----------------------

export MKL_NUM_THREADS=128;
export NUMEXPR_NUM_THREADS=128;
export OMP_NUM_THREADS=128;

source /app/conda/etc/profile.d/conda.sh;


Rscript ./MRmethod/causeThred.R