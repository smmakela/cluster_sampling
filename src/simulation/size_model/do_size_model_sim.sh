#!/bin/sh
#Torque directives
#PBS -W group_list=yetistats
#PBS -N size_model_sim
#PBS -l nodes=1,walltime=02:00:00,mem=3gb
#PBS -V
#PBS -M smm2253@columbia.edu
#PBS -m a
#PBS -j oe
#PBS -o localhost:/vega/stats/users/smm2253/cluster_sampling/output/simulation/size_model
#PBS -e localhost:/vega/stats/users/smm2253/cluster_sampling/output/simulation/size_model
export CCACHE_DISABLE=1
OUTDIR="/vega/stats/users/smm2253/cluster_sampling/output/simulation/size_model"
Rscript --no-save --vanilla ${PBS_O_WORKDIR}/master.R\
 > ${OUTDIR}/master.routput 2>&1
