#!/bin/sh
#Torque directives
#PBS -N process_results_varyK
#PBS -W group_list=yetistats
#PBS -l nodes=1,walltime=04:00:00,mem=8gb
#PBS -V
#PBS -j oe
#PBS -o localhost:/vega/stats/users/smm2253/Projects/Cluster_Sampling/Code/Simplify/vary_K/outfiles
#PBS -e localhost:/vega/stats/users/smm2253/Projects/Cluster_Sampling/Code/Simplify/vary_K/outfiles
outdir="/vega/stats/users/smm2253/Projects/Cluster_Sampling/Code/Simplify/vary_K/outfiles"
export CCACHE_DISABLE=1
R CMD BATCH --no-save --vanilla $PBS_O_WORKDIR/process_results_test.r $outdir/$PBS_JOBNAME.routput

