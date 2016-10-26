#!/bin/sh
#Torque directives
#PBS -N master_sim
#PBS -W group_list=yetistats
#PBS -l nodes=1,walltime=12:00:00,mem=12gb
#PBS -V
#PBS -M smm2253@columbia.edu
#PBS -m a
#PBS -j oe
#PBS -o localhost:/vega/stats/users/smm2253/cluster_sampling/output/simulation/outfiles
#PBS -e localhost:/vega/stats/users/smm2253/cluster_sampling/output/simulation/outfiles
outdir="/vega/stats/users/smm2253/cluster_sampling/output/simulation/outfiles"
export CCACHE_DISABLE=1
Rscript --no-save --vanilla $PBS_O_WORKDIR/sim_master.R --simno ${PBS_ARRAYID} --numclusters $J > $outdir/$PBS_JOBNAME.routput 2>&1
