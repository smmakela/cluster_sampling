#!/bin/sh
#Torque directives
#PBS -N makepopdata
#PBS -W group_list=yetistats
#PBS -l nodes=1,walltime=04:00:00,mem=4gb
#PBS -V
#PBS -M smm2253@columbia.edu
#PBS -m a
#PBS -j oe
#PBS -o localhost:/vega/stats/users/smm2253/Projects/Cluster_Sampling/Code/Simplify/vary_K/outfiles
#PBS -e localhost:/vega/stats/users/smm2253/Projects/Cluster_Sampling/Code/Simplify/vary_K/outfiles
outdir="/vega/stats/users/smm2253/Projects/Cluster_Sampling/Code/Simplify/vary_K/outfiles"
export CCACHE_DISABLE=1
R CMD BATCH --no-save --vanilla '--args 1' $PBS_O_WORKDIR/makepopdata.r $outdir/$PBS_JOBNAME.routput

