#!/bin/sh
#Torque directives
#PBS -W group_list=yetistats
#PBS -N make_rng_data 
#PBS -l nodes=1,walltime=01:00:00,mem=2gb
#PBS -V
#PBS -M smm2253@columbia.edu
#PBS -m a
#PBS -j oe
#PBS -t 1
#PBS -o localhost:/vega/stats/users/smm2253/cluster_sampling/output/simulation/outfiles
#PBS -e localhost:/vega/stats/users/smm2253/cluster_sampling/output/simulation/outfiles
export CCACHE_DISABLE=1
OUTDIR="/vega/stats/users/smm2253/cluster_sampling/output/simulation/outfiles"
Rscript --no-save --vanilla ${PBS_O_WORKDIR}/make_rng_data.R > ${OUTDIR}/make_rng_data.routput 2>&1
