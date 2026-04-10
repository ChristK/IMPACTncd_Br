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

# This is the baseline scenario. All other scenarios are compared with this one
# Because we inform the model with only one time-point, calendar trends of sbp
# and salt have not been extracted from the data. Therefore, these trends are
# user assumptions and can be specified separately for men and women. All trends
# are centred to 2013.

## User inputs -----
sc_inputs <- list()
sc_inputs$names <- "sc00" # REMEMBER TO UPDATE for new scenarios. Baseline needs always be sc00
sc_inputs$men_sbp_annual_trend          <- 0.0 # Mean sbp for men will be reducing by i.e. 0.1 mmHg annualy
sc_inputs$women_sbp_annual_trend        <- 0.0 # Mean sbp for women will be reducing by i.e. 0.1 mmHg annualy

sc_inputs$men_salt_added_baseline_calibration   <- 0.0 # Mean added salt for men one-off change for init year by i.e. 2 g/d
sc_inputs$women_salt_added_baseline_calibration <- 0.0 # Mean added salt for women one-off change for init year by i.e. 2 g/d
sc_inputs$men_salt_other_baseline_calibration   <- 0.0 # Mean other salt for men one-off change for init year by i.e. 2 g/d
sc_inputs$women_salt_other_baseline_calibration <- 0.0 # Mean other salt for women one-off change for init year by i.e. 2 g/d

sc_inputs$men_salt_added_annual_trend   <- -0.04248  # Mean added salt for men will be reducing by i.e. 0.1 g/d annualy
sc_inputs$women_salt_added_annual_trend <- -0.04248  # Mean added salt for women will be reducing by i.e. 0.1 g/d annualy
sc_inputs$men_salt_other_annual_trend   <- +0.04248  # Mean other salt for men will be reducing by i.e. 0.1 g/d annualy
sc_inputs$women_salt_other_annual_trend <- +0.04248  # Mean other salt for women will be reducing by i.e. 0.1 g/d annualy


## Model logic. DO NOT ALTER ----
cat(paste0("Loading ", sc_inputs$names, "...\n"))
sc_inputs$names <- paste0("_", sc_inputs$names)

# added salt
POP[, xps_salt_added_sc00 := xps_salt_added]

POP[sex == "men", xps_salt_added_sc00 := ((
  mean(xps_salt_added) + sc_inputs$men_salt_added_baseline_calibration
) /
  mean(xps_salt_added)) * xps_salt_added]
POP[sex == "women", xps_salt_added_sc00 := ((
  mean(xps_salt_added) + sc_inputs$women_salt_added_baseline_calibration
) /
  mean(xps_salt_added)) * xps_salt_added]

POP[sex == "men", xps_salt_added_sc00 := ((
  mean(xps_salt_added) + sc_inputs$men_salt_added_annual_trend * (year - design$init_year)
) /
  mean(xps_salt_added)) * xps_salt_added]
POP[sex == "women", xps_salt_added_sc00 := ((
  mean(xps_salt_added) + sc_inputs$women_salt_added_annual_trend * (year - design$init_year)
) /
  mean(xps_salt_added)) * xps_salt_added]
POP[, xps_salt_added := NULL]

# other salt
POP[, xps_salt_other_sc00 := xps_salt_other]

POP[sex == "men", xps_salt_other_sc00 := ((
  mean(xps_salt_other) + sc_inputs$men_salt_other_baseline_calibration
) /
  mean(xps_salt_other)) * xps_salt_other]
POP[sex == "women", xps_salt_other_sc00 := ((
  mean(xps_salt_other) + sc_inputs$women_salt_other_baseline_calibration
) /
  mean(xps_salt_other)) * xps_salt_other]


POP[sex == "men", xps_salt_other_sc00 := ((
  mean(xps_salt_other) + sc_inputs$men_salt_other_annual_trend * (year - design$init_year)
) /
  mean(xps_salt_other)) * xps_salt_other]
POP[sex == "women", xps_salt_other_sc00 := ((
  mean(xps_salt_other) + sc_inputs$women_salt_other_annual_trend * (year - design$init_year)
) /
  mean(xps_salt_other)) * xps_salt_other]

POP[, xps_salt_sc00 := xps_salt_added_sc00 + xps_salt_other_sc00]
POP[, xps_salt_other := NULL]

# sbp
POP[, xps_sbp_sc00 := xps_sbp]
POP[sex == "men", xps_sbp_sc00 := ((mean(xps_sbp) + sc_inputs$men_sbp_annual_trend * (year - design$init_year)) /
                                                             mean(xps_sbp)) * xps_sbp]
POP[sex == "women", xps_sbp_sc00 := ((mean(xps_sbp) + sc_inputs$women_sbp_annual_trend * (year - design$init_year)) /
                                                               mean(xps_sbp)) * xps_sbp]
POP[, xps_sbp := NULL]

# create lagged variables
POP[, xps_sbp_cvdlag := shift_byid(xps_sbp_sc00,
                                   cvd_lag_mc, 115, pid)] # lag exposure


# calculate morbidity/mortalty
source(exprs = str2expression(gsub("_sc00",  sc_inputs$names,
                                 readLines(file.path("other_model.R")))), local = TRUE
)
source(exprs = str2expression(gsub("_sc00",  sc_inputs$names,
                                 readLines(file.path("chd_model.R")))), local = TRUE
)
source(exprs = str2expression(gsub("_sc00",  sc_inputs$names,
                                 readLines(file.path("stroke_model.R")))), local = TRUE
)

# resolve deaths from multiple causes at the same year
POP[, all_cause_mrtl_sc00 :=
      mortality_resolve(
        year,
        new_simulant,
        other_mrtl_sc00,
        chd_mrtl_sc00,
        stroke_mrtl_sc00,
        rn_resolve_death,
        design$init_year
      )]

# write to disk
cat(paste0("Writing to disk...\n"))
ptm <- proc.time()
sc_inputs$names <- gsub("_", "", sc_inputs$names)
prunned_POP <- prune_pop(POP, sc_inputs$names, design)

out_disease <- summarise_hlp(prunned_POP, sc_inputs$names, design)
out_disease[, `:=` (scenario = sc_inputs$names,
                    year = year + 2000L,
                    mc = mc_iter)]
setkeyv(out_disease, c(design$strata_for_outputs, "duration"))

to_agegrp(out_disease, grp_width = 20L, max_age = design$ageH)
out_disease[, age := NULL]
out_disease <- out_disease[, lapply(.SD, sum), by = c(
  setdiff(design$strata_for_outputs, "age"),
  "agegroup", "duration",
  "mc",
  "scenario"
)]

fwrite_safe(out_disease, output_dir("disease_output.csv"))
rm(out_disease)

prunned_POP[, lifecourse :=
              morbidity_resolve(
                year,
                # htn_prvl,
                # t2dm_prvl,
                chd_prvl,
                stroke_prvl,
                # lung_ca_prvl,
                all_cause_mrtl,
                design$friendly_disease_names,
                design$friendly_fatal_diseases_names,
                design$init_year
              )]
del_col_nam <-
  setdiff(grep("_mrtl$", names(prunned_POP), value = TRUE), "all_cause_mrtl")
prunned_POP[, (del_col_nam) := NULL]
prunned_POP[, `:=` (scenario = sc_inputs$names,
                    mc = mc_iter)]
setkeyv(prunned_POP, c("pid", "year"))

replace_from_table(
  prunned_POP,
  "age",
  design$ageL:design$ageH,
  agegrp_name(design$ageL, design$ageH, 20L, FALSE, TRUE, design$ageH),
  "agegroup"
)

out_xps <- groupingsets(
  prunned_POP,
  j = lapply(.SD, mean),
  by = c("year", "sex", "agegroup", "scenario", "mc"),
  .SDcols = grep("^xps_", names(prunned_POP), value = TRUE),
  sets = list(c("scenario", "mc", "year"), c("scenario", "mc", "year", "sex"), c("scenario", "mc", "year", "sex", "agegroup"))
)[, year := year + 2000L]
for (j in seq_len(ncol(out_xps)))
  set(out_xps, which(is.na(out_xps[[j]])), j, "All")

setkey(out_xps, year, scenario)
fwrite_safe(out_xps, output_dir("xps_output.csv"))
rm(out_xps, prunned_POP)
print(proc.time() - ptm)
