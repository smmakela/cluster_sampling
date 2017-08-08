# This is the makefile for the simulation


# This little snippet will print out whatever variable is put in the % position
# useful for debugging: just call make print-VARNAME to display the value of
# VARNAME
print-%  : ; @echo $* = $($*)

# *****************************************************************************
# Define variables that will be used in this makefile
# *****************************************************************************
# File paths
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
NUM_SIMS         = 100
USE_SIZE_VALS    = 1
OUTCOME_TYPES    = continuous
SIZE_MODEL_NAMES = multinomial ff ffstrat
STAN_MODEL_NAMES = knowsizes negbin lognormal

# Variables that store filenames of population data
# e.g. <file path>/popdata_usesizes_<0,1>_<continuous, binary>.<rds, Rout>
# Here we create the stems of the file names:
# popdata_usesizes_<0,1>_<continuous,binary>

# First we add the prefix "popdata_usesizes_" to the values of USE_SIZE_VALS
# There's probably a more elegant way to do this in one step with vectors,
# but I couldn't figure out how to do it because addprefix/addsuffix will not
# take vectors for the prefix/suffix parameter
POP_FILE_STEMS := $(addprefix popdata_usesizes_, $(USE_SIZE_VALS))
POP_FILE_PATHS := $(foreach v1, $(SIZE_MODEL_NAMES), $(foreach v2, $(OUTCOME_TYPES), $(foreach v3, $(POP_FILE_STEMS), $(v3)_$(v2)_$(v1))))
# Now add the file path in $(SIM_OUT_DIR) to the filename stems
POP_FILE_PATHS := $(addprefix $(SIM_OUT_DIR)/, $(POP_FILE_PATHS))
# Then add the .rds and .Rout file endings to the file paths
POP_DATA_FILES := $(addsuffix .rds,  $(POP_FILE_PATHS))
POP_ROUT_FILES := $(addsuffix .Rout, $(POP_FILE_PATHS))
# Concatenate all pop files; used in clean:
POP_FILES_ALL  := $(POP_DATA_FILES) $(POP_ROUT_FILES) 

# Code files that are used as either prereqs or in recipes
# e.g. sampling code, stan code, other model code files, figure code
SAMPLE_SRC_FILE := $(SIM_SRC_DIR)/sampledata.R
STAN_SRC_FILES  := $(wildcard $(MODEL_SRC_DIR)/*.stan)
STAN_MOD_FILES  := $(patsubst %.stan, %.rds, $(STAN_SRC_FILES))
MODEL_SRC_FILES := $(STAN_SRC_FILES) $(MODEL_SRC_DIR)/lmer_compare.R \
                   $(MODEL_SRC_DIR)/svy_ests.R
DB_SRC_FILES    := $(MODEL_SRC_DIR)/lmer_compare.R $(MODEL_SRC_DIR)/svy_ests.R

# Files used to check whether all simulations completed 
STAN_CHECK_STEMS := $(addprefix $(SIM_OUT_DIR)/total_stan_check_usesizes_, $(USE_SIZE_VALS))
STAN_CHECK_FILES := $(foreach v1, $(STAN_MODEL_NAMES), $(foreach v2, $(SIZE_MODEL_NAMES), $(foreach v3, $(OUTCOME_TYPES), $(foreach v4, $(STAN_CHECK_STEMS), $(v4)_$(v3)_$(v2)_$(v1).txt))))

DB_CHECK_STEMS  := $(addprefix $(SIM_OUT_DIR)/total_db_check_usesizes_, $(USE_SIZE_VALS))
DB_CHECK_FILES := $(foreach v1, $(SIZE_MODEL_NAMES), $(foreach v2, $(OUTCOME_TYPES), $(foreach v3, $(DB_CHECK_STEMS), $(v3)_$(v2)_$(v1).txt)))

STAN_PREREQS := $(STAN_MOD_FILES) \
                $(SIM_SRC_DIR)/sampledata.R \
                $(SIM_SRC_DIR)/sim_master_stan.R \
                $(SIM_SRC_DIR)/runstan.R
DB_PREREQS   := $(DB_SRC_FILES) \
                $(SIM_SRC_DIR)/sampledata.R \
                $(SIM_SRC_DIR)/sim_master_db.R
# Prerequisites for $(SIM_OUT_DIR)/compiled_simulation_results.rds:
# - check file to make sure that all simulations finished successfully
# - code to compile simulation results
COMPILED_PREREQS := $(STAN_CHECK_FILES) \
                    $(DB_CHECK_FILES) \
                    $(SIM_SRC_DIR)/compile_simulation_results.R 

# Define substring extracting functions -- these are used to pull out various
# parts of strings like 0_continuous_multinomial_bb
get_firstword  = $(firstword $(subst _, ,$1))      # pulls out 0
get_secondword = $(word 2, $(subst _, ,$1))        # pulls out continuous
get_thirdword = $(word 3, $(subst _, ,$1))         # pulls out multinomial
empty :=
space := $(empty) $(empty)
get_fourthword  = $(subst $(space),_,$(wordlist 4, 6, $(subst _, ,$1)))  # pulls out bb
# get_fourthword uses wordlist so we get all of cluster_inds_only

# *****************************************************************************
# Commands to build targets
# *****************************************************************************

.PHONY: all pop stan db stanmodels ccc
.SECONDARY:
#all: $(POP_ROUT_FILES)
#all: $(SIM_OUT_DIR)/sim_res_check.txt
all: $(STAN_CHECK_FILES) $(DB_CHECK_FILES)
pop: $(POP_ROUT_FILES)
stan: $(STAN_CHECK_FILES)
db: $(DB_CHECK_FILES)
stanmodels: $(STAN_MOD_FILES)
ccc: $(CCC_FILES)
#all: $(SIM_RES_CHECK_FILES) $(SIM_OUT_DIR)/compiled_simulation_results.rds

# GENERATE POPULATION DATA
# We only need to check whether the output files (the ones
# ending in .Rout) are older than makepopdata.r, since running makepopdata.r
# will *create* the .rds files
# The % is used to extract the <use_sizes>_<outcome_type> part of the filenames,
# e.g. 0_continuous for the case that use_sizes is 0 and outcome_type is
# "continuous". The $* in the recipe uses the stem, which for this example
# would be 0_continuous, and the functions get_firstword and get_secondword 
# split 0_continuous into "0" and "continous". The $@ part tells R to use the
# name of the target as the R output file filename (this is how the target
# gets created).
$(SIM_OUT_DIR)/popdata_usesizes_%.Rout: $(SIM_SRC_DIR)/makepopdata.R
	@echo $@
	@echo ----------------------- making pop data --------------------------
	Rscript $(R_OPTS) $<\
	 --numclusters $(NUM_POP_CLUSTERS)\
	 --use_sizes $(call get_firstword, $*)\
	 --outcome_type $(call get_secondword, $*)\
	 --size_model $(call get_thirdword, $*)\
	 > $@ 2>&1 

# COMPILE STAN MODELS
$(MODEL_SRC_DIR)/%.rds: $(MODEL_SRC_DIR)/%.stan
	@echo -------------------- compiling stan models -----------------------
	Rscript $(R_OPTS) compile_stan_models.R --model_name $* > $*.Rout 2>&1

# RUN SIMULATION
# If a check file for the simulation is older than its prereqs (pop data,
# model code), then rerun entire simulation for that value of ($US_VAL, $OT_VAL)
# the % will be of the from x_y, where x is either 0 or 1 (one of the possible
# values of $US_VAL) and y is either "continuous" or "discrete", the two
# possible values of $OT_VAL
# $(SIM_OUT_DIR)/total_stan_check_usesizes_%.txt: $(SIM_OUT_DIR)/popdata_usesizes_*.Rout $(STAN_PREREQS)
#$(SIM_OUT_DIR)/total_stan_check_usesizes_%.txt: $$(subst $$(space),, $$(SIM_OUT_DIR)/popdata_usesizes_$$(subst $$(to_remove),, $$*).Rout) $(STAN_PREREQS)
# FIRST RUN STAN
to_remove := $(addprefix _, $(STAN_MODEL_NAMES))
#TEST := $(subst $(space),, $(SIM_OUT_DIR)/popdata_usesizes_$(subst $(filter 1_continuous_ffstrat_negbin, $(to_remove)), $(to_remove),, 1_continuous_ffstrat_negbin).Rout)
TT := $(subst $(space),, $(SIM_OUT_DIR)/popdata_usesizes_$(subst $(strip $(foreach tt, $(to_remove), $(findstring $(tt), 1_continuous_ffstrat_negbin))),, 1_continuous_ffstrat_negbin).Rout)
#TT := $(foreach tt, $(to_remove), $(tt))
FF := $(foreach tt, $(to_remove), $(findstring $(tt), 1_continuous_ffstrat_negbin))
FF2 := $(strip $(foreach tt, $(to_remove), $(findstring $(tt), 1_continuous_ffstrat_negbin)))
TEST := $(subst $(space),, $(SIM_OUT_DIR)/popdata_usesizes_$(subst $(to_remove),, 1_continuous_ffstrat_negbin).Rout)
.SECONDEXPANSION:
$(SIM_OUT_DIR)/total_stan_check_usesizes_%.txt: $$(subst $$(space),, $$(SIM_OUT_DIR)/popdata_usesizes_$$(subst $$(strip $$(foreach tt, $$(to_remove), $$(findstring $$(tt), $$*))),, $$*).Rout) $(STAN_PREREQS)
	@echo ---------------------- doing run_sim_stan.sh -------------------------
	echo $@
	USV=$(call get_firstword, $*);\
	OTV=$(call get_secondword, $*);\
	SM=$(call get_thirdword, $*);\
	MN=$(call get_fourthword, $*);\
	JN="$${MN}_usz_$${USV}_$${OTV}_$${SM}";\
	$(SIM_SRC_DIR)/run_sim_stan.sh $(NUM_POP_CLUSTERS) $${USV} $${OTV} $${SM} $${MN} $(NUM_SIMS) $${JN};\
	#JN="$${MN}_usz_"; JN+=$${USV}; JN+="_"; JN+=$${OTV};\
	#Rscript $(R_OPTS) count_stan_check_files.R --use_sizes ${USV} \
	#--outcome_type ${OTV} --num_sims ${NUM_SIMS} --model_name ${MN};

# THEN RUN SVY, LMER
$(SIM_OUT_DIR)/total_db_check_usesizes_%.txt: $(SIM_OUT_DIR)/popdata_usesizes_*.Rout $(DB_PREREQS)
	@echo ---------------------- doing run_sim_db.sh -------------------------
	USV=$(call get_firstword, $*);\
	OTV=$(call get_secondword, $*);\
	SM=$(call get_thirdword, $*);\
	JN="db_usz_$${USV}_$${OTV}_$${SM}";\
	$(SIM_SRC_DIR)/run_sim_db.sh $(NUM_POP_CLUSTERS) $${USV} $${OTV} $${SM} $(NUM_SIMS) $${JN};\
	#JN="db_usz_"; JN+=$${USV}; JN+="_"; JN+=$${OTV};\
	#Rscript $(R_OPTS) count_db_check_files.R --use_sizes ${USV} \
	#--outcome_type ${OTV} --num_sims ${NUM_SIMS};

# only remove simulation results (not pop data)
#simclean:
#	rm -f $(SIM_RESULTS) $(FIG_PDF_FILES)

#clean:
#	rm -f $(POP_FILES_ALL) $(SIM_RESULTS) $(FIG_PDF_FILES) $(PAPER_PDF_FILE)
