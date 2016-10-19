# This is the makefile for the simulation

# Constants
ROOT_DIR     = /vega/stats/users/smm2253/cluster_sampling
R_OPTS       = --no-save --no-restore --no-init-file --no-site-file

SRC_DIR      := $(ROOT_DIR)/src
SIM_CODE_DIR := $(SRC_DIR)/simulation
SIM_OUT_DIR  := $(ROOT_DIR)/output/simulation
#POP_FILES     := $(SIM_OUT_DIR)/$(wildcard popdata_usesizes_*.rds)
POP_ROUT_FILES := $(SIM_OUT_DIR)/popdata_usesizes_%.Rout

#@echo $(POPFILES)
print-%  : ; @echo $* = $($*)
all: $(POP_ROUT_FILES)
.PHONY: all clean

# Define substring extracting functions
define get_firstword = 
$(firstword $(subst _, ,$1))
endef

define get_secondword =
$(word 2, $(subst _, ,$1))
endef

# Generate pop data
$(SIM_OUT_DIR)/popdata_usesizes_%.Rout: $(SIM_CODE_DIR)/makepopdata.r
	Rscript $(R_OPTS) makepopdata.r --use_sizes $(call get_firstword, $*)\
	 --outcome_type $(call get_secondword, $*) > $@< 2>&1 
 
# The % is used to extract the <use_sizes>_<outcome_type> part of the filenames,
# e.g. 0_continuous for the case that use_sizes is 0 and outcome_type is
# "continuous". The $* in the recipe uses the stem, which for this example
# would be 0_continuous, and the arg1 and arg2 functions split 0_continuous
# into "0" and "continous"

#$(SIM_OUT_DIR)/popdata_usesizes_%.rds: $(SIM_OUT_DIR)/popdata_usesizes_%.Rout
#	arg1 = $(firstword $(subst _, ,$1)) \
	arg2 = $(word 2, $(subst _, ,$1)) \
	Rscript $(R_OPTS) makepopdata.r --use_sizes $(call arg1, $*)\
	 --outcome_type $(call arg2, $*) > $< 2>&1 

#clean:
#	rm -rf $(POP_FILES) $(POP_ROUT_FILES)
