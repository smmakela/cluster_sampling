#!/bin/sh
#Torque directives
#PBS -W group_list=yetistats
#PBS -l nodes=1,walltime=02:00:00,mem=12gb
#PBS -V
#PBS -M smm2253@columbia.edu
#PBS -m a
#PBS -j oe
#PBS -o localhost:/vega/stats/users/smm2253/cluster_sampling/output/simulation/outfiles
#PBS -e localhost:/vega/stats/users/smm2253/cluster_sampling/output/simulation/outfiles
export CCACHE_DISABLE=1
OUTDIR="/vega/stats/users/smm2253/cluster_sampling/output/simulation/outfiles"
JOB_NAME="db_usz_${US_VAL}_${OT_VAL}_${SM_VAL}"
Rscript --no-save --vanilla ${PBS_O_WORKDIR}/sim_master_db.R\
 --simno=${PBS_ARRAYID} --use_sizes=${US_VAL}\
 --outcome_type=${OT_VAL} --size_model=${SM_VAL} --numclusters=${J} > ${OUTDIR}/${JOB_NAME}-${PBS_ARRAYID}.routput 2>&1
