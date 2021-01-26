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


cat("Initialising Brazil Salt Policy model...\n\n")
if (!require(BrazilSaltModelmisc)) {
  if (!require(remotes)) install.packages("remotes")
remotes::install_local("./Rpackage/Brazil_salt_model_pkg/")
library(BrazilSaltModelmisc)
}

options(warn = 1)

design$init_year <- design$init_year - 2000L
output_dir <- function(x = character(0)) paste0("./Output/", x)

synthpop_dir <- "./SynthPop"


options(rgl.useNULL = TRUE)  # suppress error by demography in rstudio server
# Then try/install packages...

dependencies(
  c(
    # "gamlss", # only necesary when fitting the models
    # "gamlss.tr", # only necesary when fitting the models
    # "mc2d", # only necessary for generating fixed_mc
    "doParallel",
    "doRNG",
    "foreach",
    "qs",
    "fst",
    "data.table"
  ), TRUE, FALSE, FALSE, FALSE
)

options(datatable.verbose = FALSE)
options(datatable.showProgress = FALSE)

# Define RNG for parallel use with doRNG
RNGkind("L'Ecuyer-CMRG")
dice <- function(n = .N) my_runif(n, 0, 1)


# Define end() function to beep end print a message
end_sim <- function(...) {
	cat("All done! \a\n")
	sink(
		file = output_dir("simulation parameters.txt"),
		append = TRUE,
		type = "output",
		split = FALSE
	)
	cat(paste0("Simulation ended successfully at: ", Sys.time(), "\n"))
	sink()
	if (Sys.info()[1] == "Windows") {
		system("rundll32 user32.dll,MessageBeep -1")
		Sys.sleep(.5)
	}
}

# Function for timing log
time_mark <- function(x) {
	sink(
		file = output_dir("times.txt"),
		append = TRUE,
		type = "output",
		split = FALSE
	)
	cat(paste0(x, " at: ", Sys.time(), "\n"))
	sink()
}

# Logic to continue iterations ----
# from previous runs if the scenarios are completely the same. If you need a completely new
# run, please delete contents in the output_dir. If scenarios have changed, then
# contents in output_dir are deleted.
# TODO include fixed_mc.qs in the files to monitor

mc_iter_from_previous_run <- 0L
if (all(file.exists(output_dir(".scenarios_snapshot.qs")),
        file.exists(output_dir("disease_output.csv")),
        file.exists(output_dir("xps_output.csv")))) {
  # if scenarios snapshot and output files exist
  snapshot <-
    changedFiles(qread(output_dir(".scenarios_snapshot.qs")))

  if (!any(nzchar(snapshot$added) + nzchar(snapshot$deleted) + nzchar(snapshot$changed))) {
    # If no changes, get previous max mc_iter and update current run
    mc_iter_from_previous_run <-
      rbind(fread(output_dir("disease_output.csv"), select = "mc"),
            fread(output_dir("xps_output.csv"), select = "mc"))[, max(mc)]
  } else {
    # if scenarios have changed DELETE contents of output_dir and start from scratch
    cat("Model parameters have changed since last run. Previous model outputs have been deleted...\n")
    delete_output_files()
  }
  rm(snapshot)
} else {
  # if scenarios snapshot or some output files don't exist
  # DELETE contents of output_dir and start again
  cat("Model parameters have changed since last run. Previous model outputs have been deleted...\n")
  delete_output_files()

  qsave(
    fileSnapshot(
      "./Scenarios",
      timestamp = tempfile("timestamp"),
      md5sum = TRUE
    ),
    output_dir(".scenarios_snapshot.qs")
  )
}




# Find and load scenarios -------------------------------------------------
scenarios_list <- sort(list.files(
	path = "./Scenarios",
	pattern = glob2rx("sc*.R"),
	full.names = TRUE,
	recursive = FALSE
))

n_scenarios <- length(scenarios_list)
# Specify output.txt file for simulation parameters -----------------------
fileOut <- file(output_dir("simulation parameters temp.txt"))
writeLines(
  c(
    "Brazil Salt Policy model\nA dynamic microsimulation, by Dr Chris Kypridemos",
    "\n",
    paste0("Simulation started at: ", Sys.time(), "\n"),
    "Simulation parameters:\n",
    paste0(
      "Continue from previous run = ",
      ifelse(mc_iter_from_previous_run > 0L, "Yes", "No")
    ),
    paste0("Kismet = "                      , design$kismet),
    paste0("Jumpiness = "                   , design$jumpiness),
    paste0("First year of the simulation = ", design$init_year),
    paste0("Years to project = "            , design$sim_horizon),
    paste0("design$ageL = "                  , design$ageL),
    paste0("design$ageH = "                  , design$ageH),
    paste0("design$cvd_lag = "               , design$cvd_lag),
    paste0("diseases = "                    , design$diseases),
    paste0(
      "CHD incidence contribution to mortality decline = ",
      design$chd_incidence_ratio
    ),
    paste0(
      "Stroke incidence contribution to mortality decline = ",
      design$stroke_incidence_ratio
    ),
    paste0("Sample size = "                 , format(design$n, scientific = FALSE)),
    paste0(
      "Number of mc_iter = "        ,
      design$iteration_n + mc_iter_from_previous_run
    ),
    paste0("Stochastic = "                  , design$stochastic)
  ),
  fileOut
)
close(fileOut)

# Sample for parameter distributions --------------------------------------
# precalculate realisations of random variables

# Load file if exist, else create it from start
if ((design$iteration_n + mc_iter_from_previous_run) > 1e4L)
  stop(
    "Please generate a larger fixed_mc file.\nThe current supports up to 10K iterations.\nPlease delete file 'fixed_mc.qs' and increase max_iter to more than 10K for this.\n"
  )
if (file.exists("./Scenarios/fixed_mc.qs")) {
	fixed_mc <- qread("./Scenarios/fixed_mc.qs") # TODO will not work when more than 4k iterations are needed
} else {
  cat("Rebuilding fixed_mc... \n")
	max_iter <- 1e4L
	fixed_mc <- list()
	fixed_mc$cvd_lag_l <-
		#(design$cvd_lag - 1) / 9 # calculates p of binom for mean = user input design$cvd_lag
		1L + rbinom(max_iter, 9, (design$cvd_lag - 1) / 9) # max lag = 10 years

	fixed_mc$sodium_optim_l <-
		mc2d::rpert(max_iter, 614, 1500, 2391, 4) # optimal level for sodium around 1500 mg/day. Under which no risk from mozaffarian NEJM

	fixed_mc$salt_optim_l <-
	  mc2d::rpert(max_iter, 614*2.5/1000, 1500*2.5/1000, 2391*2.5/1000, 4) # optimal level for salt around 3.75 g/day. Under which no risk from mozaffarian NEJM

	fixed_mc$chd_incidence_ratio_l <-
		do.call(rbeta, c(
			list(n = max_iter),
			estim_beta_params(design$chd_incidence_ratio, 0.005),
			list(ncp = 0)
		))

	fixed_mc$stroke_incidence_ratio_l <-
		do.call(rbeta, c(
			list(n = max_iter),
			estim_beta_params(design$stroke_incidence_ratio, 0.005),
			list(ncp = 0)
		))

	#Inflate all_cause (proxy to non_cvd) mortality for hypertensives
	# RR (= 1.3) roughly based on stringhini_socioeconomic_2017 figure 4
	fixed_mc$death_sbp_rr_l <-  runif(max_iter, 1.2, 1.4)

	# lifetable stochasticity
	# need to apply after year 16 since for years 13 to 16 we have observed data
	fixed_mc$breakpoint_year_noncvd_l <-
		sample(16:(design$init_year + design$sim_horizon), max_iter, TRUE)
	fixed_mc$breakpoint_year_cvd_l <-
		sample(16:(design$init_year + design$sim_horizon), max_iter, TRUE)
	fixed_mc$breakpoint_quant_noncvd_l <- runif(max_iter)
	fixed_mc$breakpoint_quant_cvd_l <- runif(max_iter)

	# percentile of distribution for incidence, prevalence
	fixed_mc$chd_burden_l <- runif(max_iter)

	fixed_mc$chd_sbp_rr_l <- setkey(
		fread(
			"./CVD_statistics/sbp.rrchd.csv",
			stringsAsFactors = F,
			colClasses = c("factor", "factor",
										 "numeric", "numeric")),
		agegroup,	sex)
	fixed_mc$chd_sbp_rr_l[, `:=` (mean.rr = exp(mean.rr),
																ci.rr = exp(ci.rr))] # original table has coef of rr ie log(rr)

	fixed_mc$chd_sbp_rr_l <- clone_dt(fixed_mc$chd_sbp_rr_l, max_iter)
	fixed_mc$chd_sbp_rr_l[, rr := stochRRtabl(mean.rr, ci.rr, design$stochastic),
												by = .id]
	fixed_mc$chd_sbp_rr_l[is.na(rr) | rr < 1, rr := 1]
	# fixed_mc$chd_sbp_rr_l[, .(unique(mean.rr), median(rr)), by = .(agegroup, sex)]

	# theoretical minimum distribution.
	# level of sbp below no risk exist. From Singh et al
	fixed_mc$chd_sbp_tmred_l <-	rnorm(
		max_iter,
		mean = runif(max_iter, 110, 115),
		sd = runif(max_iter, 4, 6))

	# percentile of distribution for incidence, prevalence
	fixed_mc$stroke_burden_l <- runif(max_iter)

	fixed_mc$stroke_sbp_rr_l <- setkey(
		fread(
			"./CVD_statistics/sbp.rrstroke.csv",
			stringsAsFactors = F,
			colClasses = c("factor", "factor",
										 "numeric", "numeric")),
		agegroup,	sex)

	fixed_mc$stroke_sbp_rr_l[, `:=` (mean.rr = exp(mean.rr),
																	 ci.rr = exp(ci.rr))] # original table has coef of rr ie log(rr)
	fixed_mc$stroke_sbp_rr_l <- clone_dt(fixed_mc$stroke_sbp_rr_l, max_iter)
	fixed_mc$stroke_sbp_rr_l[, rr := stochRRtabl(mean.rr, ci.rr), by = .id]
	fixed_mc$stroke_sbp_rr_l[is.na(rr) | rr < 1, rr := 1]

	# theoretical minimum distribution.
	# level of sbp below no risk exist. From Singh et al
	fixed_mc$stroke_sbp_tmred_l <-
		rnorm(
			max_iter,
			mean = runif(max_iter, 110, 115),
			sd = runif(max_iter, 4, 6))

	tt <- CJ(age = design$ageL:design$ageH, sex = 1:2)
	to_agegrp(tt, agegroup_colname = "agegroup", max_age = design$ageH)
	tt[, `:=` (sex = as.character(sex))]
	fixed_mc$chd_sbp_rr_l <-
	  na.omit(tt[fixed_mc$chd_sbp_rr_l, on = c("agegroup", "sex"), allow.cartesian = TRUE])
	fixed_mc$chd_sbp_rr_l[, rr := predict(loess(rr ~ age, span = 0.75)), by = .(.id, sex)]
	fixed_mc$chd_sbp_rr_l[rr < 1, rr := 1]
	fixed_mc$chd_sbp_rr_l[, sex := factor(sex, c("1", "2"), c("men", "women"))]
	#chd_sbp_rr_l[sex == "1", plot(age, rr)]
	#chd_sbp_rr_l[sex == "1" & .id == 1, lines(age, rr.sm)]


	fixed_mc$stroke_sbp_rr_l <-
	  na.omit(tt[fixed_mc$stroke_sbp_rr_l, on = c("agegroup", "sex"), allow.cartesian = TRUE])
	fixed_mc$stroke_sbp_rr_l[, rr := predict(loess(rr ~ age, span = 0.75)), by = .(.id, sex)]
	fixed_mc$stroke_sbp_rr_l[rr < 1, rr := 1]
	fixed_mc$stroke_sbp_rr_l[, sex := factor(sex, c("1", "2"), c("men", "women"))]

	#stroke_sbp_rr_l[sex == "1", plot(age, rr)]
	#stroke_sbp_rr_l[sex == "1" & .id == 1, lines(age, rr.sm)]

	qsave(fixed_mc, file = "./Scenarios/fixed_mc.qs")
}

if (!design$stochastic) {
	fixed_mc$cvd_lag_l <- rep(design$cvd_lag, design$iteration_n)
	fixed_mc$sodium_optim_l <- rep(1500, design$iteration_n)
	fixed_mc$design$chd_incidence_ratio_l <- rep(design$chd_incidence_ratio, design$iteration_n)
	fixed_mc$design$stroke_incidence_ratio_l <- rep(design$stroke_incidence_ratio, design$iteration_n)
	fixed_mc$death_sbp_rr_l <-  rep(1.3, design$iteration_n)
	fixed_mc$breakpoint_year_noncvd_l <- rep(design$sim_horizon, design$iteration_n)
	fixed_mc$breakpoint_year_cvd_l <- rep(design$sim_horizon, design$iteration_n)
	fixed_mc$breakpoint_quant_noncvd_l <- rep(0.5, design$iteration_n)
	fixed_mc$breakpoint_quant_cvd_l <- rep(0.5, design$iteration_n)
	fixed_mc$chd_burden_l <- rep(0.5, design$iteration_n)
	fixed_mc$chd_sbp_rr_l[, rr := mean.rr]
	fixed_mc$chd_sbp_tmred_l <- rep(110, design$iteration_n)
	fixed_mc$stroke_burden_l <- rep(0.5, design$iteration_n)
	fixed_mc$stroke_sbp_rr_l[, rr := mean.rr]
	fixed_mc$stroke_sbp_tmred_l <- rep(110, design$iteration_n)
	fixed_mc$stroke_sbp_tmred_l <- rep(1.0, design$iteration_n)
}


# Generate fixed synthpop if file doesn't exists
if (design$use_fixed_pop &&
         !file.exists("./SynthPop/brazil_pop_1.fst")) {
  cat("Generating synthetic population...\n")
  ptm <- proc.time()

  POP <- copy(population_actual)
  POP <-
    POP[rep.int(seq_len(.N), pct)][, c("population", "pct") := NULL]
  gen_rn_cols(POP,
              c("rn_sbp", "rn_salt_added", "rn_salt_other", "rn_black_race"))

  POP[, `:=` (pid = .I,
              # create the correlarion (-0.132) observed between added and other salt
              rn_salt_other = ifelse(rbinom(.N, 1L, 0.132), (1 - rn_salt_added), rn_salt_other))]
  # Note: if I perform above after scrambling it introduces big jumbs of ranks
  # cor(POP$rn_salt_added, POP$rn_salt_other)

  # synthesise Black race
  # Black race gamlss model was only fitted to ages 20-79. So for ages <20
  # I have to assume age == 20.
  # TODO better approach. Fit gamlss model for younger ages
  tt <-
    fread("./Lifecourse_models/black_race_table.csv")[between(age, min(POP$age), max(POP$age))]
  tt1 <- CJ(age = min(POP$age):20L, sex = unique(tt$sex))
  tt1[tt[age == min(age),], on = "sex", mu := i.mu]
  tt <- rbind(tt, tt1)
  rm(tt1)

  POP[tt, on = intersect(names(POP), names(tt)),
      xps_black_race := qbinom(rn_black_race, 1L, mu)]
  # POP[, prop.table(table(xps_black_race))]

  # Clone POP
  # tt <- as.integer(do.call(max, design[grep("_lag$", names(design))]))
  tt <- as.integer(max(fixed_mc$cvd_lag_l))

  POP <- clone_dt(POP, design$sim_horizon + tt + 1L)
  POP[, .id := .id - tt - 1L]
  POP[, `:=` (age  = age + .id,
              year = design$init_year + .id)]

  POP <-
    POP[between(age, design$ageL - tt, design$ageH)] # delete unnecessary ages
  POP[, .id := NULL]
  setkey(POP, pid, year)
  POP[, rn_sbp := scramble_trajectories(rn_sbp, pid, design$jumpiness)]
  POP[, rn_salt_other := scramble_trajectories(rn_salt_other, pid, design$jumpiness)]
  POP[, rn_salt_added := scramble_trajectories(rn_salt_added, pid, design$jumpiness)]

  tt <-
    fread("./Lifecourse_models/sbp_table.csv")[between(age, POP[, min(age)], POP[, max(age)])]
  # POP[tt, on = intersect(names(POP), names(tt)), c("mu", "sigma", "nu", "tau") := .(i.mu, i.sigma, i.nu, i.tau)]
  # POP[, xps_sbp := qBCPEo_sbp(rn_sbp, TRUE, FALSE, mu, sigma, nu, tau), by = age] # by is for performance only (50% faster)
  POP[tt, on = intersect(names(POP), names(tt)), xps_sbp := my_qBCPEo_trunc(rn_sbp, mu, sigma, nu, tau, 75, 240, getDTthreads())] # OpenMP optimised

  tt <-
    fread("./Lifecourse_models/salt_added_table.csv")[between(age, POP[, min(age)], POP[, max(age)])]
  POP[tt, on = intersect(names(POP), names(tt)), xps_salt_added :=
        my_qBCT_trunc(rn_salt_added, mu, sigma, nu, tau, 0, 20, getDTthreads())] # OpenMP optimised

  tt <-
    fread("./Lifecourse_models/salt_othersources_table.csv")[between(age, POP[, min(age)], POP[, max(age)])]
  POP[tt, on = intersect(names(POP), names(tt)), xps_salt_other :=
        my_qBCPEo_trunc(rn_salt_other, mu, sigma, nu, tau, 0, 20, getDTthreads())] # OpenMP optimised

  POP[, grep("^rn_", names(POP), value = TRUE) := NULL] # delete rn cols

  setkey(POP, pid, year)
  POP[, new_simulant := mk_new_simulant_markers(pid)]
  write_fst(POP, "./SynthPop/brazil_pop_1.fst", 100L)
  print(proc.time() - ptm)
}
