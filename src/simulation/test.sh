#!/bin/sh
US_VAL=0
OT_VAL=continuous
J=100
PBS_JOBNAME=TEST
OUTDIR="/vega/stats/users/smm2253/cluster_sampling/output/simulation/outfiles"
OUTFILE_NAME="${OUTDIR}/${PBS_JOBNAME}_${US_VAL}_${OT_VAL}.routput"
echo "Rscript --no-save --vanilla sim_master.R --simno 1\
 --use_sizes ${US_VAL} --outcome_type ${OT_VAL} --numclusters ${J} \ 
 >> ${OUTFILE_NAME} 2>&1" > ${OUTFILE_NAME}
