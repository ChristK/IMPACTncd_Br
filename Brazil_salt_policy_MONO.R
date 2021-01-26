## Brazil Salt Policy model: a decision support tool for primary prevention of NCDs
## Copyright (C) 2019 Chris Kypridemos

## Brazil Sodium Policy model is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have a copy of the GNU General Public License
## along with this program; if not, see <http://www.gnu.org/licenses/>
## or write to the Free Software Foundation, Inc., 51 Franklin Street,
## Fifth Floor, Boston, MA 02110-1301  USA.

# ********************************************************************************
# IMPACT NCD Brazil ONE CORE
# ********************************************************************************
design <- list()
design$iteration_n    <- 10L # Overrides design saved in file. If changed in design.R all outputs will be deleted
design$clusternumber  <- 1L # Change to your number of CPU cores (explicit parallelisation)
design$n_cpus         <- 1L  # Change to your number of CPU cores (implicit parallelisation/Open MP)
design$use_fixed_pop  <- TRUE
design$logs           <- FALSE
design$process_output <- TRUE

source(file = "./Scenarios/design.R")
# *************************************************************************************
cat("Initialising Brazil Salt Policy model...\n\n")

# Main --------------------------------------------------------------------
source(file = "./initialisation.R")
setDTthreads(20L)
options(datatable.auto.index = TRUE)

source(file = "./lifetable_engine.R")
source(file = "./diseases_epidemiology.R")

mc_iter = 4L
setDTthreads(1L)
my_env <- environment() # get environment of this branch

# Actual simulation
sys.source(file = "./2dmc.R", my_env)

# Load synthetic population
sys.source(file = "./synth_pop_gen.R", my_env)

# Scenarios
lapply(scenarios_list, sys.source, envir = my_env)

# sys.source(file = "./partial_output.R", my_env)


setDTthreads(20L)
if (design$process_output == TRUE) {
		source(file = "./output.R")
}

