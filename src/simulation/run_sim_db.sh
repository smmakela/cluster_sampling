#!/bin/bash
echo "arguments:"
echo "arguments:"
echo "J = $1"
echo "US_VAL = $2"
echo "OT_VAL = $3"
echo "SM_VAL = $4"
echo "NUM_SIMS = $5"
echo "JOBNAME = $6"
echo $1
echo $2
echo $3
echo $4
echo $5
echo $6
echo "FIRST QSUB COMMAND:"
FIRST=$(qsub -v J=$1,US_VAL=$2,OT_VAL=$3,SM_VAL=$4\
            -t 1-$5 -N "$6" do_sim_master_db.sh)
echo "qsub -v J=$1,US_VAL=$2,OT_VAL=$3,SM_VAL=$4 -t 1-$5 -N "$6" do_sim_master_db.sh"
echo $FIRST
