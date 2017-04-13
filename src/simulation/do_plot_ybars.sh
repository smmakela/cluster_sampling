#!/bin/sh
#Torque directives
#PBS -W group_list=yetistats
#PBS -N plot_ybars
#PBS -l nodes=1,walltime=04:00:00,mem=4gb
#PBS -V
#PBS -M smm2253@columbia.edu
#PBS -m a
#PBS -j oe
#PBS -o localhost:/vega/stats/users/smm2253/cluster_sampling/output/simulation/outfiles
#PBS -e localhost:/vega/stats/users/smm2253/cluster_sampling/output/simulation/outfiles
export CCACHE_DISABLE=1
OUTDIR="/vega/stats/users/smm2253/cluster_sampling/output/simulation/outfiles"
Rscript --no-save --vanilla ${PBS_O_WORKDIR}/plot_ybars.R > ${OUTDIR}/plot_ybars.routput 2>&1
