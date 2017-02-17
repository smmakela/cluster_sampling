#!/bin/sh
#Torque directives
#PBS -W group_list=yetistats
#PBS -l nodes=1,walltime=12:00:00,mem=24gb
#PBS -V
#PBS -M smm2253@columbia.edu
#PBS -m a
#PBS -j oe
#PBS -o localhost:/vega/stats/users/smm2253/cluster_sampling/output/simulation/outfiles
#PBS -e localhost:/vega/stats/users/smm2253/cluster_sampling/output/simulation/outfiles
export CCACHE_DISABLE=1
OUTDIR="/vega/stats/users/smm2253/cluster_sampling/output/simulation/outfiles"
JOB_NAME="ccc_${US_VAL}_${OT_VAL}_${MOD_NAME}"
Rscript --no-save --vanilla ${PBS_O_WORKDIR}/ccc.R\
 --use_sizes=${US_VAL} --outcome_type=${OT_VAL} --mod_name=${MOD_NAME} >\
 ${OUTDIR}/${JOB_NAME}.routput 2>&1
