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

cat("Loading stroke (I60-I69) model...\n")
ptm <- proc.time()
set(POP, NULL, "stroke_prvl_sc00",  0L)
set(POP, NULL, "stroke_mrtl_sc00", 0L)

# RR for SBP from Optimal SBP level at 115mmHg and RR(HR) of dying from stroke was
# taken from "Singh GM, Danaei G, Farzadfar F, Stevens GA, Woodward M, Wormser D,
# et al. The age-specific quantitative effects of metabolic risk factors on
# cardiovascular diseases and diabetes: a pooled analysis. PLOS ONE.
# 2013 Jul 30;8(7):e65174.

# cat("sbp RR\n")
# set(POP, NULL, "stroke_sbp_rr",  1)
# POP[sbp_rr_stroke_mc, stroke_sbp_rr := bound(rr ^ ((xps_sbp_cvdlag - stroke_sbp_tmred_mc) / 10), 1, 20),
#     on = c("age", "sex")]

set(POP, NULL, "rr",  1)
POP[sbp_rr_stroke_mc, rr := i.rr, on = c("age", "sex")]
POP[, stroke_sbp_rr := rr ^ ((xps_sbp_cvdlag - stroke_sbp_tmred_mc) / 10)]
# POP[, stroke_sbp_rr := bound(stroke_sbp_rr, 1, 20)]
POP[stroke_sbp_rr < 1, stroke_sbp_rr := 1]
POP[, rr := NULL]

if (!"stroke_prvl_init" %in% names(POP)) {
  # Estimate prevalence -------------------------------------------------
  set(POP, NULL, "stroke_prvl_init",  0L)

  #cat(paste0("Estimating stroke prevalence in ", init.year, " ...\n\n"))
  age_structure <- population_actual[, .(age, sex, pct)][between(age, design$ageL - cvd_lag_mc, design$ageH)]
  # POP[between(age, design$ageL - cvd_lag_mc, design$ageH) & year == design$init_year, .N,
  #            keyby = .(age, sex)]
  setkey(age_structure, age, sex)
  if (design$stochastic) {
    age_structure[stroke_epi_mc$prevalence[between(age, design$ageL - cvd_lag_mc, design$ageH), ],
                  Nprev := rbinom(.N, ifelse(is.na(pct), 0L, pct), prevalence)]
    age_structure[stroke_epi_mc$incidence[between(age, design$ageL - cvd_lag_mc, design$ageH)],
                  Nprev := ifelse(is.na(Nprev), 0L, Nprev) -
                    rbinom(.N,
                           ifelse(is.na(pct), 0L, pct) - ifelse(is.na(Nprev),
                                                                0L, Nprev),
                           incidence)]
  } else {
    age_structure[stroke_epi_mc$prevalence[between(age, design$ageL - cvd_lag_mc, design$ageH)],
                  Nprev := round(pct * prevalence)]
    age_structure[stroke_epi_mc$incidence[between(age, design$ageL - cvd_lag_mc, design$ageH)],
                  Nprev := Nprev - round((pct - Nprev) * incidence)]
  }

  age_structure[Nprev < 0, Nprev := 0]
  #age_structure[strokesurv  , Nprev := round(Nprev *(1-fatality * 1.03))]
  setnames(age_structure, "pct", "population")

  pid_stroke <-
    POP[between(age, design$ageL - cvd_lag_mc, design$ageH) &
          year == design$init_year,
        .(pid = resample(pid,
                         age_structure[age == .BY[[1]] &
                                         sex == .BY[[2]], Nprev],
                         replace = F)),
        by = .(age, sex, year)]

  pid_stroke[,
             duration :=
               as.integer(rpois(.N, bound(
                 rnorm(
                   .N,
                   stroke_epi_l$duration[1]$stroke_duration,
                   stroke_epi_l$duration[1]$se
                 ),
                 1,
                 100
               )))] # estimate duration
  pid_stroke[(age - duration) < 0, duration := age]
  POP[pid_stroke, stroke_prvl_init := duration, on = c("pid", "year")]

  rm(pid_stroke)

  # Estimate PAF --------------------------------------------------------
  #cat("Estimating stroke PAF...\n")
  strokepaf <-
    POP[between(age, design$ageL, design$ageH) &
          stroke_prvl_init == 0L &
          year == design$init_year,
        .(paf = 1 - 1 / (sum(stroke_sbp_rr) / .N)),
        keyby = .(age, sex)]

  # strokepaf[, plot(age, paf, ylim = c(0,1), main=paste0(.BY[[1]])), keyby = .(sex)]
  strokepaf[, paf := predict(loess(paf ~ age, span = 0.5)), by = .(sex)]
  setkey(stroke_epi_mc$incidence, age, sex)
  stroke_epi_mc$incidence[strokepaf, p0 := incidence * (1 - paf)]
  stroke_epi_mc$incidence[is.na(p0), p0 := incidence]


  # Estimate non SBP related incidence trends -----------------------------------------------
  stroke_epi_mc$incidence <-
    stroke_epi_mc$incidence[lifetable_mc[between(age, design$ageL, design$ageH),
                                         .(age, sex, year, stroke_prop_change)],
                            on = c("age", "sex")]

  stroke_epi_mc$incidence[, p0 := p0 * (1 + (stroke_prop_change - 1) * stroke_incidence_ratio_mc)]

  set(POP, NULL, "p0_stroke",  0)
  POP[stroke_epi_mc$incidence,
      p0_stroke := i.p0,
      on = c("age", "sex", "year")] # NA age > 84
  #POP[between(age, design$ageL, design$ageH) & year >= 0, summary(p0*100000)]
}

# Estimate stroke incidence -------------------------------
POP[, prb_stroke_inc := p0_stroke * stroke_sbp_rr]
if (anyNA(POP$prb_stroke_inc)) setnafill(POP, "const", 0, cols = "prb_stroke_inc")


if (!design$kismet)
  gen_rn_cols(POP,
              c("rn_stroke_inc", "rn_stroke_death", "rn_stroke_death30"))

setkey(POP, pid, year)
POP[, stroke_prvl_sc00 := incidence_type_2(
  year,
  new_simulant,
  prb_stroke_inc,
  rn_stroke_inc,
  stroke_prvl_init,
  design$init_year
)]


# Estimate stroke mortality ----
# Calibrated to match mortality projections. Only run in sc00
# qx * pop = prev * fatality rate
if (!"qx_stroke" %in% names(POP)) {
  tt <-
    POP[between(age, design$ageL, design$ageH) &
          year >= design$init_year &
          other_mrtl_sc00 == 0L & chd_mrtl_sc00 == 0L,
        .(N = as.numeric(.N), prev = as.numeric(sum(stroke_prvl_sc00 > 0L))),
        keyby = .(year, age, sex)]
  tt[lifetable_mc,
     qx := qx_stroke_mc,
     on = c("age", "sex", "year")]
  # tt[stroke_epi_l$day30_fatality,
  #    day30_fatality := i.day30_fatality, # day30_fatality is the probability of dying within 30 days from the event
  #    on = "age"]

  # account for stroke expected deaths (more accurate calculation of future prevalence
  for (ii in (design$init_year + head((seq_len(
    design$sim_horizon
  ) - 1L), -1L))) {
    ttt <- tt[year == ii, .(
      year = year + 1L,
      age = age + 1L,
      sex,
      "expected_deaths" = qx * N
    )]
    # subtract expected deaths from next year
    tt[ttt, `:=` (prev = prev - i.expected_deaths,
                  N    = N    - i.expected_deaths),
       on = c("age", "sex", "year")]
  }

  tt[prev <= 0, prev := 1L] # avoid div by 0
  # tt[, day30_deaths  := day30_fatality * incid]
  tt[, fatality      := qx * N / prev] # fatality is the probability of dying from the event

  POP[tt, `:=` (qx_stroke = i.fatality),
                # qx_stroke30 = i.day30_fatality),
                on = c("age", "sex", "year")]

  # Calibration. Needed because of rounding errors of the prevalence (0 and integers)
  POP[, qx_stroke := qx_stroke * 1.05]
}

setkey(POP, pid, year)
POP[, stroke_mrtl_sc00 := mortality_type_1(
  year,
  new_simulant,
  qx_stroke,
  rn_stroke_death,
  stroke_prvl_sc00,
  design$init_year,
  3L
)]

if (exists("tt")) rm(tt)
print(proc.time() - ptm)
