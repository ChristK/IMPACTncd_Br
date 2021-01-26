#!/usr/bin/Rscript
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
# IMPACT NCD Brazil
# ********************************************************************************
design <- list()
design$iteration_n    <- 100L # Overrides design saved in file. If changed in design.R all outputs will be deleted
design$clusternumber  <- 10L # Change to your number of CPU cores (explicit parallelisation)
design$n_cpus         <- 1L  # Change to your number of CPU cores (implicit parallelisation/Open MP)
design$use_fixed_pop  <- TRUE
design$logs           <- TRUE
design$process_output <- TRUE

source(file = "./Scenarios/design.R")
source(file = "./initialisation.R")
source(file = "./lifetable_engine.R")
source(file = "./diseases_epidemiology.R")

if (Sys.info()[1] == "Windows") {
  cl <- makeCluster(design$clusternumber) # used for clustering. Windows compatible
  registerDoParallel(cl)
} else {
  registerDoParallel(design$clusternumber)  # used for forking. Only Linux/OSX compatible
}

gc()
time_mark("start parallelisation")
foreach(
	mc_iter = (1L + mc_iter_from_previous_run):(design$iteration_n + mc_iter_from_previous_run),
	.inorder = FALSE,
	.verbose = TRUE,
	.packages = c("data.table", "BrazilSaltModelmisc", "fst"),
	.export = ls(),
	.noexport = c("scenarios_list", "time_mark")
) %dorng%
{
  my_env <- environment()  # get environment of this branch
  sys.source(file = "./simulation.R", my_env)
  gc()
  0L
}

if (exists("cl")) stopCluster(cl)
time_mark("End of parallelisation")
gc()

setDTthreads(20L)
file.rename(output_dir("simulation parameters temp.txt"),
            output_dir("simulation parameters.txt"))
while (sink.number() > 0L) sink()

if (design$process_output) source(file = "./output.R")

end_sim()

