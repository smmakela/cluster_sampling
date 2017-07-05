#!/bin/sh
#Torque directives
#PBS -W group_list=yetistats
#PBS -N compile_svy_results
#PBS -l nodes=1,walltime=00:05:00,mem=1gb
#PBS -V
#PBS -M smm2253@columbia.edu
#PBS -m a
#PBS -j oe
#PBS -t 1-128
#PBS -o localhost:/vega/stats/users/smm2253/cluster_sampling/output/simulation/outfiles
#PBS -e localhost:/vega/stats/users/smm2253/cluster_sampling/output/simulation/outfiles
export CCACHE_DISABLE=1
OUTDIR="/vega/stats/users/smm2253/cluster_sampling/output/simulation/outfiles"
R CMD BATCH --no-save --vanilla -${PBS_ARRAYID} ${PBS_O_WORKDIR}/compile_svy_results.R ${OUTDIR}/${PBS_JOBNAME}.routput
