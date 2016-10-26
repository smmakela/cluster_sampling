# This is the makefile for the simulation


# This little snippet will print out whatever variable is put in the % position
# useful for debugging: just call make print-VARNAME to display the value of
# VARNAME
print-%  : ; @echo $* = $($*)

# *****************************************************************************
# Define variables that will be used in this makefile
ROOT_DIR       = /vega/stats/users/smm2253/cluster_sampling
OUT_DIR       := $(ROOT_DIR)/output
SIM_OUT_DIR   := $(OUT_DIR)/simulation
FIG_OUT_DIR   := $(OUT_DIR)/figures
PAPER_OUT_DIR := $(OUT_DIR)/manuscript
SRC_DIR       := $(ROOT_DIR)/src
SIM_SRC_DIR   := $(SRC_DIR)/simulation
FIG_SRC_DIR   := $(SIM_OUT_DIR)/figures
MODEL_SRC_DIR := $(SRC_DIR)/analysis
R_OPTS         = --no-save --no-restore --no-init-file --no-site-file

# Hard-code highest-level simulation parameters
NUM_POP_CLUSTERS = 100
NUM_SIMS         = 1

# Variables that store filenames of population data
# e.g. <file path>/popdata_usesizes_<0,1>_<continuous, binary>.<rds, Rout>
USE_SIZE_VALS  = 0 1
# Here we create the stems of the file names:
# popdata_usesizes_<0,1>_<continuous,binary>
# First we add the prefix "popdata_usesizes_" to the values of USE_SIZE_VALS
# There's probably a more elegant way to do this in one step with vectors,
# but I couldn't figure out how to do it because addprefix/addsuffix will not
# take vectors for the prefix/suffix parameter
POP_FILE_STEMS := $(addprefix popdata_usesizes_, $(USE_SIZE_VALS))
POP_FILE_STEMS := $(addsuffix _continuous, $(POP_FILE_STEMS)) \
                  $(addsuffix _binary, $(POP_FILE_STEMS))
# Now add the file path in $(SIM_OUT_DIR) to the filename stems
POP_FILE_PATHS := $(addprefix $(SIM_OUT_DIR)/, $(POP_FILE_STEMS))
# Then add the .rds and .Rout file endings to the file paths
POP_DATA_FILES := $(addsuffix .rds,  $(POP_FILE_PATHS))
POP_ROUT_FILES := $(addsuffix .Rout, $(POP_FILE_PATHS))
# Concatenate all pop files; used in clean:
POP_FILES_ALL  := $(POP_DATA_FILES) $(POP_ROUT_FILES)

# Code files that are used as either prereqs or in recipes
# e.g. sampling code, stan code, other model code files, figure code
SAMPLE_SRC_FILE := $(SIM_SRC_DIR)/sampledata.R
STAN_SRC_FILES  := $(wildcard $(MODEL_SRC_DIR)/*.stan)
MODEL_SRC_FILES := $(STAN_SRC_FILES) $(MODEL_SRC_DIR)/lmer_compare.R \
                   $(MODEL_SRC_DIR)/svy_ests.R
FIG_SRC_FILES   := $(wildcard $(FIG_CODE_DIR)/figure_*.R)

# Compiled simulation results
SIM_RESULTS := $(OUT_DIR)/compiled_simulation_results.csv

# Prerequisites for $(SIM_OUT_DIR)/sim_res_check.txt
# - population files
# - model code (stan files plus lmer/svy code)
# If any of the above are newer than $(SIM_OUT_DIR)/sim_res_check.txt, rerun
# the simulation
SIM_RES_CHECK_PREREQS := $(POP_ROUT_FILES) $(MODEL_SRC_FILES)

# Prerequisites for $(SIM_OUT_DIR)/compiled_simulation_results.csv:
# - check file to make sure that all simulations finished successfully
# - code to compile simulation results
SIM_RES_PREREQS := $(sim_out_dir)/sim_res_check.txt \
                   $(SIM_SRC_DIR)/compile_simulation_results.R

# Figure pdf files
FIG_PDF_FILES := $(wildcard $(FIG_OUT_DIR)/figures_*.pdf)

# Paper file
PAPER_TEX_FILE := $(PAPER_OUT_DIR)/paper.tex
PAPER_BIB_FILE := $(PAPER_OUT_DIR)/clusterrefs.bib
PAPER_PDF_FILE := $(PAPER_OUT_DIR)/paper.pdf
# *****************************************************************************

.PHONY: all clean 
#all: $(POP_ROUT_FILES)
all: $(SIM_OUT_DIR)/sim_res_check.txt

# Define substring extracting functions
get_firstword  = $(firstword $(subst _, ,$1))
get_secondword = $(word 2, $(subst _, ,$1))

# GENERATE POPULATION DATA
# We only need to check whether the output files (the ones
# ending in .Rout) are older than makepopdata.r, since running makepopdata.r
# will *create* the .rds files
# The % is used to extract the <use_sizes>_<outcome_type> part of the filenames,
# e.g. 0_continuous for the case that use_sizes is 0 and outcome_type is
# "continuous". The $* in the recipe uses the stem, which for this example
# would be 0_continuous, and the arg1 and arg2 functions split 0_continuous
# into "0" and "continous"
$(SIM_OUT_DIR)/popdata_usesizes_%.Rout: $(SIM_SRC_DIR)/makepopdata.R
	Rscript $(R_OPTS) makepopdata.R\
	 --numclusters $(NUM_POP_CLUSTERS)\
	 --use_sizes $(call get_firstword, $*)\
	 --outcome_type $(call get_secondword, $*)\
	 > $@ 2>&1 

# RUN SIMULATION
# If the check file for the simulation is older than its prereqs (pop data,
# model code), then rerun entire simulation
$(SIM_OUT_DIR)/sim_res_check.txt: $(SIM_RES_CHECK_PREREQS)
	qsub -v J=$(NUM_POP_CLUSTERS) -t 1-$(NUM_SIMS)\
	 $(SIM_SRC_DIR)/do_sim_master.sh 

# COMPILE RESULTS
$(SIM_RESULTS): $(SIM_RES_PREREQS)
	Rscript $(R_OPTS) $<

# CREATE FIGURES
$(FIG_OUT_DIR)/figure_%.pdf: $(FIG_SRC_DIR)/figure_%.R
	Rscript $(R_OPTS) $<

# GENERATE PAPER
$(PAPER_PDF_FILE): $(FIG_PDF_FILES) $(PAPER_TEX_FILE) $(PAPER_BIB_FILE)
	pdflatex paper
	bibtex paper
	pdflatex paper
	pdflatex paper

# only remove simulation results (not pop data)
simclean:
	rm -f $(SIM_RESULTS) $(FIG_PDF_FILES)

clean:
	rm -f $(POP_FILES_ALL) $(SIM_RESULTS) $(FIG_PDF_FILES) $(PAPER_PDF_FILE)
