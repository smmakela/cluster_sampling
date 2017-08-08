#!/bin/sh
#Torque directives
#PBS -W group_list=yetistats
#PBS -l nodes=1,walltime=00:30:00,mem=2gb
#PBS -V
#PBS -M smm2253@columbia.edu
#PBS -m a
#PBS -j oe
#PBS -o localhost:/vega/stats/users/smm2253/cluster_sampling/output/simulation/outfiles
#PBS -e localhost:/vega/stats/users/smm2253/cluster_sampling/output/simulation/outfiles
export CCACHE_DISABLE=1
OUTDIR="/vega/stats/users/smm2253/cluster_sampling/output/simulation/outfiles"
JOB_NAME="negbin_usz_1_continuous_ff"
Rscript --no-save --vanilla ${PBS_O_WORKDIR}/sim_master_stan.R\
 --simno=1 --use_sizes=1\
 --numclusters=100 --outcome_type=continuous --size_model=ff_strat --model_name=negbin > tmp.routput 2>&1
