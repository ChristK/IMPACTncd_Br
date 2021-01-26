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

if (design$use_fixed_pop &&
    file.exists("./SynthPop/brazil_pop_1.fst")) {
  cat("Loading precalculated synthetic population...\n")
  ptm <- proc.time()
  POP <-
    read_fst("./SynthPop/brazil_pop_1.fst", as.data.table = TRUE)[year >= (design$init_year - cvd_lag_mc)]
  gen_rn_cols(
    POP,
    c(
      "rn_chd_inc",
      "rn_chd_death",
      "rn_chd_death30",
      "rn_stroke_inc",
      "rn_stroke_death",
      "rn_stroke_death30",
      "rn_other_death",
      "rn_resolve_death"
    )
  )
  print(proc.time() - ptm)

} else if (design$use_fixed_pop &&
           !file.exists("./SynthPop/brazil_pop_1.fst")) {
  # This branch only entered in MONO if file missing. In parallel the file is created during initialisation
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
  POP[tt, on = intersect(names(POP), names(tt)), xps_sbp := my_qBCPEo_trunc(rn_sbp, mu, sigma, nu, tau, 75, 240, design$n_cpus)] # OpenMP optimised

  tt <-
    fread("./Lifecourse_models/salt_added_table.csv")[between(age, POP[, min(age)], POP[, max(age)])]
  POP[tt, on = intersect(names(POP), names(tt)), xps_salt_added :=
        my_qBCT_trunc(rn_salt_added, mu, sigma, nu, tau, 0, 20, design$n_cpus)] # OpenMP optimised

  tt <-
    fread("./Lifecourse_models/salt_othersources_table.csv")[between(age, POP[, min(age)], POP[, max(age)])]
  POP[tt, on = intersect(names(POP), names(tt)), xps_salt_other :=
        my_qBCPEo_trunc(rn_salt_other, mu, sigma, nu, tau, 0, 20, design$n_cpus)] # OpenMP optimised

  POP[, grep("^rn_", names(POP), value = TRUE) := NULL] # delete rn cols

  setkey(POP, pid, year)
  POP[, new_simulant := mk_new_simulant_markers(pid)]
  write_fst(POP, "./SynthPop/brazil_pop_1.fst", 100L)

  # aux
  gen_rn_cols(
    POP,
    c(
      "rn_chd_inc",
      "rn_chd_death",
      "rn_chd_death30",
      "rn_stroke_inc",
      "rn_stroke_death",
      "rn_stroke_death30",
      "rn_other_death",
      "rn_resolve_death"
    )
  )

  print(proc.time() - ptm)

} else {
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
  tt <- as.integer(cvd_lag_mc)
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

  # synthesise SBP
  # qBCPEo_sbp <- trun.q(
  #   par = c(75, 240),
  #   family = "BCPEo",
  #   type = "both"
  # )


  tt <-
    fread("./Lifecourse_models/sbp_table.csv")[between(age, POP[, min(age)], POP[, max(age)])]
  # POP[tt, on = intersect(names(POP), names(tt)), c("mu", "sigma", "nu", "tau") := .(i.mu, i.sigma, i.nu, i.tau)]
  # POP[, xps_sbp := qBCPEo_sbp(rn_sbp, TRUE, FALSE, mu, sigma, nu, tau), by = age] # by is for performance only (50% faster)
  POP[tt, on = intersect(names(POP), names(tt)), xps_sbp := my_qBCPEo_trunc(rn_sbp, mu, sigma, nu, tau, 75, 240, design$n_cpus)] # OpenMP optimised

  tt <-
    fread("./Lifecourse_models/salt_added_table.csv")[between(age, POP[, min(age)], POP[, max(age)])]
  POP[tt, on = intersect(names(POP), names(tt)), xps_salt_added :=
        my_qBCT_trunc(rn_salt_added, mu, sigma, nu, tau, 0, 20, design$n_cpus)] # OpenMP optimised

  tt <-
    fread("./Lifecourse_models/salt_othersources_table.csv")[between(age, POP[, min(age)], POP[, max(age)])]
  POP[tt, on = intersect(names(POP), names(tt)), xps_salt_other :=
        my_qBCPEo_trunc(rn_salt_other, mu, sigma, nu, tau, 0, 20, design$n_cpus)] # OpenMP optimised

  POP[, grep("^rn_", names(POP), value = TRUE) := NULL] # delete rn cols


  # aux
  setkey(POP, pid, year)
  POP[, new_simulant := mk_new_simulant_markers(pid)]
  gen_rn_cols(
    POP,
    c(
      "rn_chd_inc",
      "rn_chd_death",
      "rn_chd_death30",
      "rn_stroke_inc",
      "rn_stroke_death",
      "rn_stroke_death30",
      "rn_other_death",
      "rn_resolve_death"
    )
  )

  print(proc.time() - ptm)
}

setindexv(POP, c("age", "year", "sex", "xps_black_race"))
