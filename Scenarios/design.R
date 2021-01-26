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

# Simulation design
if (!exists("design")) {
  design <- list()
  design$iteration_n    <- 20L # override saved file. If changed in file all outputs will be deleted
  design$clusternumber  <- 20L # Change to your number of CPU cores (explicit parallelisation)
  design$n_cpus         <- 1L  # Change to your number of CPU cores (implicit parallelisation)
  design$use_fixed_pop  <- TRUE
  design$logs           <- FALSE
  design$process_output <- TRUE
}

design$strata_for_outputs <- c("year", "age", "sex")
design$exposures <- c("age", "sex", "sbp", "salt")
design$init_year   <- 2013L # min = 2010. We subtract 2000 in initialisation.R
design$n           <- 4e5L # Define the sample size
design$sim_horizon <- 28L  # NEED TO force >=1 and up to 50
design$ageL        <- 30L  # Define lower age limit to diseases-model simulation (min = 30)
design$ageH        <- 79L  # Define lower age limit to diseases-model simulation (max = 84)
design$cvd_lag     <- 5L
design$chd_incidence_ratio    <- 0.5
design$stroke_incidence_ratio <- 0.5
design$stochastic <- TRUE
design$kismet <- TRUE
design$jumpiness <- 0.03 # from 0 to 1
design$diseases <- c("chd", "stroke")
design$friendly_disease_names <- c("CHD", "Stroke")
design$fatal_diseases <- c("other", "chd", "stroke")
design$friendly_fatal_diseases_names <- c("Other COD", "CHD", "Stroke")
design$scenarios <- gsub(".R", "", sort(list.files("./Scenarios", pattern = "^sc.*\\.R$")))
design$max_prvl_for_outputs          <- 2L # Need to be > 0. If 1 incidence is
# not recorded, only prevalence. If 2, duration 1 is incidence and duration
# 2 is prevalence
