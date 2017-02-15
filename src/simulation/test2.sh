#!/bin/sh
#Torque directives
#PBS -N testing
#PBS -W group_list=yetistats
#PBS -l nodes=1,walltime=12:00:00,mem=12gb
#PBS -V
#PBS -M smm2253@columbia.edu
#PBS -m a
#PBS -j oe
#PBS -o localhost:/vega/stats/users/smm2253/cluster_sampling/output/simulation/outfiles
#PBS -e localhost:/vega/stats/users/smm2253/cluster_sampling/output/simulation/outfiles
echo $PBS_JOBNAME
export CCACHE_DISABLE=1
OUTDIR="/vega/stats/users/smm2253/cluster_sampling/output/simulation/outfiles"
Rscript --no-save --vanilla $PBS_O_WORKDIR/sim_master.R --simno=${PBS_ARRAYID} --use_sizes=${US} --outcome_type=${OT} --numclusters=$J > ${PBS_O_WORKDIR}/testing.routput 2>&1
#Rscript --no-save --vanilla ${PBS_O_WORKDIR}/sim_master.R --simno ${PBS_ARRAYID}\
# --use_sizes ${US_VAL} --outcome_type ${OT_VAL} --numclusters ${J} \ 
# >> ${OUTFILE_NAME} 2>&1
