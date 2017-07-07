#!/bin/bash
echo "arguments:"
echo "J = $1"
echo "US_VAL = $2"
echo "OT_VAL = $3"
echo "SM_VAL = $4"
echo "MOD_NAM = $5"
echo "NUM_SIMS = $6"
echo "JOBNAME = $7"
echo "FIRST QSUB COMMAND:"
FIRST=$(qsub -v J=$1,US_VAL=$2,OT_VAL=$3,SM_VAL=$4,MOD_NAM=$5\
            -t 1-$6 -N "$7" do_sim_master_stan.sh)
echo "qsub -v J=$1,US_VAL=$2,OT_VAL=$3,SM_VAL=$4,MOD_NAM=$5 -t 1-$6 -N "$7" do_sim_master_stan.sh"
echo $FIRST
