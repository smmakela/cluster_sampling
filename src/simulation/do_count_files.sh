#!/bin/sh
#Torque directives
#PBS -W group_list=yetistats
#PBS -l nodes=1,walltime=00:30:00,mem=1gb
#PBS -V
#PBS -N count_files
#PBS -M smm2253@columbia.edu
#PBS -m a
#PBS -j oe
#PBS -o localhost:/vega/stats/users/smm2253/cluster_sampling/output/simulation/outfiles
#PBS -e localhost:/vega/stats/users/smm2253/cluster_sampling/output/simulation/outfiles
export CCACHE_DISABLE=1
OUTDIR="/vega/stats/users/smm2253/cluster_sampling/output/simulation/outfiles"
R CMD BATCH --no-save --vanilla ${PBS_O_WORKDIR}/count_files.R ${OUTDIR}/count_files.routput 2>&1
