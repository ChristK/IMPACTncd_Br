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

dir.create(output_dir("logs"), showWarnings = FALSE)


if (design$logs) {
sink(file = normalizePath(output_dir(paste0("logs/log", mc_iter, ".txt")), mustWork = FALSE),
     append = TRUE,
     type = "output",
     split = FALSE)
}

setDTthreads(1L)
gc()
time_mark("start simulation")
sys.source("./2dmc.R", my_env)
sys.source("./synth_pop_gen.R", my_env)
lapply(scenarios_list, sys.source, envir = my_env)
# sys.source("./partial_output.R", my_env)
rm(my_env)  # BE CAREFULL. DTs altered with in my_env, change universaly.
# You need copies of DTs to by handled within my_env

cat("\nSimulation ended successfully.\n")
if (design$logs) sink()
gc()
