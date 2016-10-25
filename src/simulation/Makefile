# This is the makefile for the simulation


# This little snippet will print out whatever variable is put in the % position
# useful for debugging: just call make print-VARNAME to display the value of
# VARNAME
print-%  : ; @echo $* = $($*)

# *****************************************************************************
# Define variables that will be used in this makefile
ROOT_DIR        = /vega/stats/users/smm2253/cluster_sampling
SIM_OUT_DIR    := $(ROOT_DIR)/output/simulation
FIG_OUT_DIR    := $(SIM_OUT_DIR)/figures
SRC_DIR        := $(ROOT_DIR)/src
SIM_CODE_DIR   := $(SRC_DIR)/simulation
FIG_CODE_DIR   := $(SIM_CODE_DIR)/figures
FIG_CODE_FIlES := $(wildcard $(FIG_CODE_DIR)/figure_*.R)
STAN_DIR       := $(SRC_DIR)/analysis
R_OPTS          = --no-save --no-restore --no-init-file --no-site-file

# Variables needed for making population data
# e.g. <file path>/popdata_usesizes_<0,1>_<continuous, binary>.<rds, Rout>
USE_SIZE_VALS   = 0 1
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

# Variables needed for sampling and running stan
STAN_FILES := $(wildcard $(STAN_DIR)/*.stan)
SAMPLED_DATA_FILES := $(wildcard $(SIM_OUT_DIR)/sampledata_*.rds)

# Compiled simulation results
SIM_RESULTS := $(SIM_OUT_DIR)/compiled_results.csv

# Figure pdf files
FIG_PDF_FILES := $(wildcard $(FIG_OUT_DIR)/figures_*.pdf)

# Paper file
PAPER_TEX_FILE := $(ROOT_DIR)/output/paper.tex
PAPER_BIB_FILE := $(ROOT_DIR)/output/clusterrefs.bib
PAPER_PDF_FILE := $(ROOT_DIR)/output/paper.pdf
# *****************************************************************************


.PHONY: all clean 
all: $(POP_ROUT_FILES)

# Define substring extracting functions
get_firstword = $(firstword $(subst _, ,$1))
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
$(SIM_OUT_DIR)/popdata_usesizes_%.Rout: $(SIM_CODE_DIR)/makepopdata.r
	Rscript $(R_OPTS) makepopdata.r\
	 --use_sizes $(call get_firstword, $*)\
	 --outcome_type $(call get_secondword, $*)\
	 > $@ 2>&1 

# GENERATE SAMPLE DATA
$(SAMPLED_DATA_FILES): $(POP_FILES_ALL)
	Rscript $(R_OPTS) 

# COMPILE RESULTS
$(SIM_RESULTS): $(SIM_CODE_DIR)/compile_results.R
	Rscript $(R_OPTS) $<

# CREATE FIGURES
$(FIG_OUT_DIR)/figure_%.pdf: $(FIG_CODE_DIR)/figure_%.R
	Rscript $(R_OPTS) $<

# GENERATE PAPER
$(PAPER_PDF_FILE): $(FIG_PDF_FILES) $(PAPER_TEX_FILE) $(PAPER_BIB_FILE)
	pdflatex paper
	bibtex paper
	pdflatex paper
	pdflatex paper

clean:
	rm -f $(POP_FILES_ALL) $(FIG_PDF_FILES) $(PDF_FILE)
