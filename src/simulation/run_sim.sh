#!/bin/bash
echo "arguments:"
echo $1
echo $2
echo $3
echo $4
echo $5
echo $6
FIRST=$(qsub -v J=$1,US_VAL=$4,OT_VAL=$5\
            -t 1-$6 -N "$2" do_sim_master.sh)
echo "FIRST QSUB COMMAND:"
echo $FIRST
#SECOND=$(qsub -W depend=afterokarray:$FIRST\
#             -v US_VAL=$4,OT_VAL=$5,NUM_SIMS=$6\
#             -N "$3" do_count_results_files.sh)
