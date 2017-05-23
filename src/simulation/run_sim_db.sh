#!/bin/bash
echo "arguments:"
echo "arguments:"
echo "J = $1"
echo "US_VAL = $2"
echo "OT_VAL = $3"
echo "NUM_SIMS = $4"
echo "JOBNAME = $5"
echo $1
echo $2
echo $3
echo $4
echo $5
echo "FIRST QSUB COMMAND:"
FIRST=$(qsub -v J=$1,US_VAL=$2,OT_VAL=$3\
            -t 1-$4 -N "$5" do_sim_master_db.sh)
echo "qsub -v J=$1,US_VAL=$2,OT_VAL=$3 -t 1-$4 -N "$5" do_sim_master_db.sh"
echo $FIRST
