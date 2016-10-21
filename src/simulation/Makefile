# This is the makefile for the simulation

# Define variables that will be used in this makefile
ROOT_DIR        = /vega/stats/users/smm2253/cluster_sampling
SIM_OUT_DIR    := $(ROOT_DIR)/output/simulation
SRC_DIR        := $(ROOT_DIR)/src
SIM_CODE_DIR   := $(SRC_DIR)/simulation
R_OPTS          = --no-save --no-restore --no-init-file --no-site-file

USE_SIZE_VALS   = 0 1
POP_FILE_STEMS := $(addprefix popdata_usesizes_, $(USE_SIZE_VALS))
POP_FILE_STEMS := $(addsuffix _continuous, $(POP_FILE_STEMS)) \
                  $(addsuffix _binary, $(POP_FILE_STEMS))
POP_FILE_PATHS := $(addprefix $(SIM_OUT_DIR)/, $(POP_FILE_STEMS))
POP_DATA_FILES := $(addsuffix .rds,  $(POP_FILE_PATHS))
POP_ROUT_FILES := $(addsuffix .Rout, $(POP_FILE_PATHS))
POP_FILES_ALL  := $(POP_FILES) $(POP_ROUT_FILES)

# This little snippet will print out whatever variable is put in the % position
# useful for debugging: just call make print-VARNAME to display the value of
# VARNAME
print-%  : ; @echo $* = $($*)

.PHONY: all clean
all: $(POP_ROUT_FILES)

# Define substring extracting functions
get_firstword = $(firstword $(subst _, ,$1))
get_secondword = $(word 2, $(subst _, ,$1))

# Generate pop data -- we only need to check whether the output files (the ones
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
 

clean:
	rm -f $(POP_TARGETS)
