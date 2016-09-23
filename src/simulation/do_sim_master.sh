#!/bin/sh
#Torque directives
#PBS -N master_varyK
#PBS -W group_list=yetistats
#PBS -l nodes=1,walltime=12:00:00,mem=12gb
#PBS -V
#PBS -M smm2253@columbia.edu
#PBS -m a
#PBS -t 1-100
#PBS -j oe
#PBS -o localhost:/vega/stats/users/smm2253/Projects/Cluster_Sampling/Code/Simplify/vary_K/outfiles
#PBS -e localhost:/vega/stats/users/smm2253/Projects/Cluster_Sampling/Code/Simplify/vary_K/outfiles
outdir="/vega/stats/users/smm2253/Projects/Cluster_Sampling/Code/Simplify/vary_K/outfiles"
export CCACHE_DISABLE=1
R CMD BATCH --no-save --vanilla -${PBS_ARRAYID} $PBS_O_WORKDIR/do_sim_master.r $outdir/$PBS_JOBNAME.routput

