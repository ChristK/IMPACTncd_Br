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

# added salt
ptm <- proc.time()
POP[, xps_salt_added____changeme := xps_salt_added_sc00]
POP[sex == "men", xps_salt_added____changeme := (
  (mean(xps_salt_added_sc00) + linear_diffusion(
    sc_inputs$men_salt_added_change,
    year - sc_inputs$year_change_starts,
    sc_inputs$year_change_completed - sc_inputs$year_change_starts
  )) /
    mean(xps_salt_added_sc00)) * xps_salt_added_sc00]
POP[sex == "women", xps_salt_added____changeme := (
  (mean(xps_salt_added_sc00) + linear_diffusion(
    sc_inputs$women_salt_added_change,
    year - sc_inputs$year_change_starts,
    sc_inputs$year_change_completed - sc_inputs$year_change_starts
  )) /
    mean(xps_salt_added_sc00)) * xps_salt_added_sc00]

# other salt
POP[, xps_salt_other____changeme := xps_salt_other_sc00]
POP[sex == "men", xps_salt_other____changeme := (
  (mean(xps_salt_other_sc00) + linear_diffusion(
    sc_inputs$men_salt_other_change,
    year - sc_inputs$year_change_starts,
    sc_inputs$year_change_completed - sc_inputs$year_change_starts
  )) /
    mean(xps_salt_other_sc00)) * xps_salt_other_sc00]
POP[sex == "women", xps_salt_other____changeme := (
  (mean(xps_salt_other_sc00) + linear_diffusion(
    sc_inputs$women_salt_other_change,
    year - sc_inputs$year_change_starts,
    sc_inputs$year_change_completed - sc_inputs$year_change_starts
  )) /
    mean(xps_salt_other_sc00)) * xps_salt_other_sc00]

POP[, xps_salt____changeme := xps_salt_added____changeme + xps_salt_other____changeme]

# Calculate salt difference. Thisi is only relevant when salt is above the lowest
# limit after which no health effect is observed (salt_optim_mc)
set(POP, NULL, "xps_salt_diff", 0) # Covers case when bothe before and after salt < salt_optim_mc
POP[xps_salt____changeme > salt_optim_mc &
      xps_salt_sc00 > salt_optim_mc,
    xps_salt_diff := xps_salt____changeme - xps_salt_sc00]

POP[xps_salt____changeme < salt_optim_mc &
      xps_salt_sc00 > salt_optim_mc,
    xps_salt_diff := salt_optim_mc - xps_salt_sc00]

POP[xps_salt____changeme > salt_optim_mc &
      xps_salt_sc00 < salt_optim_mc,
    xps_salt_diff := xps_salt____changeme - salt_optim_mc]

# sbp
POP[, xps_sbp____changeme := xps_sbp_sc00 + salt_sbp_effect(xps_salt_diff, age, xps_sbp_sc00, xps_black_race, .N, design$stochastic)]

# create lagged variables
POP[, xps_sbp_cvdlag := shift_byid(xps_sbp____changeme,
                                   cvd_lag_mc, 115, pid)] # lag exposure

print(proc.time() - ptm)

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
POP[, all_cause_mrtl____changeme :=
      mortality_resolve(
        year,
        new_simulant,
        other_mrtl____changeme,
        chd_mrtl____changeme,
        stroke_mrtl____changeme,
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


POP[, (grep("____changeme$", names(POP), value = TRUE)) := NULL]


print(proc.time() - ptm)
