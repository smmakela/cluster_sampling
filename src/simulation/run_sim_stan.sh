#!/bin/bash
echo "arguments:"
echo "J = $1"
echo "US_VAL = $2"
echo "OT_VAL = $3"
echo "MOD_NAM = $4"
echo "NUM_SIMS = $5"
echo "JOBNAME = $6"
echo "FIRST QSUB COMMAND:"
FIRST=$(qsub -v J=$1,US_VAL=$2,OT_VAL=$3,MOD_NAM=$4\
            -t 1-$5 -N "$6" do_sim_master_stan.sh)
echo "qsub -v J=$1,US_VAL=$2,OT_VAL=$3,MOD_NAM=$4 -t 1-$5 -N "$6" do_sim_master_stan.sh"
echo $FIRST
