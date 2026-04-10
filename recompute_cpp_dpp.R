# =============================================================================
# recompute_cpp_dpp.R
# -----------------------------------------------------------------------------
# Statistically correct re-implementation of every calculation that the
# `Results abstract new scenarios.xlsx` workbook performs on the IMPACTncd_Br
# raw MC-level output. The goal is to reproduce the workbook's Summary,
# Summary UI, Total2, BR voluntary and CVD-aggregate figures from
# `Output/disease_output.csv` directly, fixing the methodological shortcuts
# that the workbook takes:
#
#   SHORTCUT 1  Sum-of-medians across years
#       The workbook's BR voluntary sheet computes 20-year cumulative CPPs as
#       SUM(D2:D27) where each row is the median over MC iterations of the
#       annual CPP. The median of a sum is NOT the sum of medians, so the
#       cumulative point estimate is biased and -- more importantly -- there
#       is no valid uncertainty interval around it. Fix: cumulate within each
#       MC iteration, then take quantiles across MCs at the cumulative level.
#
#   SHORTCUT 2  Persons = Men + Women on median values
#       The workbook's Summary tab reports "Persons" as Men + Women, where
#       both columns are MC medians. Median(M + W) is NOT median(M) + median(W).
#       Fix: per MC, sum men + women within each (year, agegroup), then take
#       the quantile across MCs. This also yields a correct UI for the
#       persons total, which the workbook does not provide.
#
#   SHORTCUT 3  CVD = CHD + Stroke on median values
#       The workbook's Summary tab reports CVD CPP as median(CHD) + median(Stroke).
#       Same problem. Fix: per MC, compute CHD + Stroke first, then quantile.
#
#   SHORTCUT 4  Mixed model runs pasted into Total2
#       The Summary cells are hardcoded paste-values from a model run that no
#       longer matches the embedded cpp_dpp_out tab. Fix: this script
#       regenerates everything from the on-disk disease_output.csv so the
#       provenance of every figure is the same single run.
#
#   SHORTCUT 5  Cumulative CPPs only at year 2038
#       The workbook focuses on the 2038 cumulative figure and does not always
#       carry the per-year cumulative trajectory consistently. Fix: full
#       cumulative time series 2019..2038 by (scenario, sex, agegroup) with
#       median, 2.5%, 97.5%, 10%, 90% quantiles.
#
# Inputs (relative to repo root, override with command-line args if needed):
#   Output/disease_output.csv   raw MC-level events and deaths by scenario
#   (optional) cost parameters in COSTS section below for the Costs tab
#
# Outputs (written to Output/Summaries/):
#   cpp_dpp_recomputed.csv      tidy long format, drop-in replacement for
#                               cpp_dpp_out.csv
#   summary_recomputed.csv      wide Summary-style table (one row per
#                               scenario/outcome with men, women, persons,
#                               and 95% UIs for each)
#   costs_recomputed.csv        per-MC cost aggregates with quantiles
#                               (only if cost parameters provided)
#
# Author: Chris (epidemiology), drafted 2026-04-09
# British English throughout. All quantiles are computed across MCs as the
# very last step; everything upstream is per-MC arithmetic.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# -----------------------------------------------------------------------------
# 0. CONFIGURATION
# -----------------------------------------------------------------------------
disease_out_path <- "Output/disease_output.csv"
sim_params_path  <- "Output/simulation parameters.txt"
out_dir          <- "Output/Summaries"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Counterfactual scenario name (must match disease_output.csv exactly)
baseline_scenario <- "Baseline"

# Years over which to cumulate
year_first <- 2019
year_last  <- 2038

# Quantile probabilities (median plus 95% and 80% UIs)
quant_probs <- c(0.025, 0.10, 0.50, 0.90, 0.975)

# Age groupings produced in the workbook
age_groups <- c("30-49", "50-69", "70+")

# Outcomes to track. The script aggregates these and then derives:
#   cumul_cvd_cpp    = cumul_chd_cpp + cumul_stroke_cpp           (per MC)
#   cumul_all_dpp    = cumul_chd_mrtl + cumul_stroke_mrtl + cumul_other_mrtl
event_cols    <- c("chd_prvl", "stroke_prvl")               # incident events -> CPP
death_cols    <- c("chd_mrtl", "stroke_mrtl", "other_mrtl") # deaths -> DPP

# Optional cost parameters (set to NULL to skip costs entirely).
# Per-event costs in 2019 Reais; replace with the audited figures from the
# Costs tab of the workbook. If supplied, costs are computed per MC and then
# quantiled, again avoiding the median-of-products shortcut.
cost_params <- list(
  chd_hosp_men   = NA_real_,
  chd_hosp_women = NA_real_,
  chd_phc_men    = NA_real_,
  chd_phc_women  = NA_real_,
  stroke_hosp_men   = NA_real_,
  stroke_hosp_women = NA_real_,
  stroke_phc_men    = NA_real_,
  stroke_phc_women  = NA_real_
)
do_costs <- !any(is.na(unlist(cost_params)))

# -----------------------------------------------------------------------------
# 1. LOAD AND VALIDATE
# -----------------------------------------------------------------------------
dis <- fread(disease_out_path)
setnames(dis, tolower(names(dis)))

# Upscale event and death counts from simulation-sample scale to full
# Brazilian-population scale, matching what post_simulation_functions.R's
# process_disease_out() does at the melt stage. Without this step the
# recomputed numbers are ~1/246 of the figures in cpp_dpp_out.csv and the
# workbook, because the simulation cohort is a ~400k sample of the ~100M
# Brazilian adults aged 30-79. pop_fraction is written into simulation
# parameters.txt by lifetable_engine.R; we parse it with the same regex
# used in post_simulation_functions.R so the two stay in sync.
sim_params <- readLines(sim_params_path)
pop_fraction <- as.numeric(
  sub("^Population fraction [^=]*= *", "",
      grep("^Population fraction ", sim_params, value = TRUE))
)
if (length(pop_fraction) != 1L || is.na(pop_fraction) || pop_fraction <= 0)
  stop("Could not parse pop_fraction from ", sim_params_path)
message("Upscaling counts by 1/pop_fraction = ", round(1 / pop_fraction, 2),
        "x (pop_fraction = ", signif(pop_fraction, 6), ")")
for (col in c(event_cols, death_cols)) {
  set(dis, j = col, value = as.integer(round(dis[[col]] / pop_fraction)))
}

# The raw disease_output.csv uses scenario CODES (sc00..sc12). Rename them to
# the friendly scenario labels used by the workbook/Total2 convention so that
# cpp_dpp_recomputed.csv is a true drop-in replacement for cpp_dpp_out.csv,
# and so the baseline_scenario config above ("Baseline") matches.
#
# design$scenarios is auto-derived from list.files("./Scenarios/sc*.R") in
# Scenarios/design.R, so new scenarios added later are picked up automatically.
# The inline friendly_scenario_names vector below MUST be kept in sync with
# the one in output.R:21-27 (currently 13 entries, last added: "Potassium Salt
# 10% (Na+K)" for the post-hoc potassium CRA via Scenarios/sc12.R).
#
# We ALWAYS re-source design.R and ALWAYS overwrite friendly_scenario_names,
# rather than guarding on exists(). Long-lived R sessions (e.g. re-running
# this script interactively after adding sc12.R to disk) otherwise carry a
# stale design$scenarios from before the new file existed, producing a
# length mismatch with the hardcoded 13-entry name list below. Re-sourcing
# design.R is safe because design.R's own `if (!exists("design"))` guard
# protects iteration_n / clusternumber / etc., while its unconditional
# `design$scenarios <- list.files(...)` assignment always refreshes.
source("./Scenarios/design.R")
friendly_scenario_names <- c(
  "Baseline", "Br Voluntary Targets", "Br Regulatory Targets",
  "PAHO Regulatory Targets", "Lowest World Regulatory Targets",
  "Potassium Salt 10%", "FOPL Warning 100%", "FOPL Warning 70%",
  "FOPL Traffic Light 100%", "FOPL Traffic light 70%",
  "Media Campaigns", "5% Reduction in 10 Years",
  "Potassium Salt 10% (Na+K)"
)
if (length(design$scenarios) != length(friendly_scenario_names))
  stop("Scenario name mismatch: design$scenarios has ",
       length(design$scenarios), " entries (",
       paste(design$scenarios, collapse = ", "),
       ") but friendly_scenario_names has ",
       length(friendly_scenario_names),
       ". Add the new scenario's friendly label to recompute_cpp_dpp.R and ",
       "output.R:21-27, or remove the dangling Scenarios/sc*.R file.")
scenario_rename <- setNames(friendly_scenario_names, design$scenarios)
dis[, scenario := scenario_rename[scenario]]
if (anyNA(dis$scenario))
  stop("Unmapped scenario code(s) in disease_output.csv. Check that ",
       "design$scenarios and friendly_scenario_names are aligned with the ",
       "current file contents.")

needed <- c("year", "sex", "agegroup", "mc", "scenario",
            event_cols, death_cols)
missing <- setdiff(needed, names(dis))
if (length(missing) > 0)
  stop("disease_output.csv is missing required columns: ",
       paste(missing, collapse = ", "))

if (!(baseline_scenario %in% dis$scenario))
  stop("Baseline scenario '", baseline_scenario, "' not found in disease_output.csv")

dis <- dis[year >= year_first & year <= year_last]
n_mc <- dis[, uniqueN(mc)]
message("Loaded ", nrow(dis), " rows; ", n_mc, " MC iterations; ",
        dis[, uniqueN(scenario)], " scenarios.")

# -----------------------------------------------------------------------------
# 2. EVENTS PREVENTED (PER MC) VERSUS BASELINE
# -----------------------------------------------------------------------------
# For each (year, sex, agegroup, mc) we compute:
#   events_prevented[scenario] = events[Baseline] - events[scenario]
# Positive values mean fewer events under the policy than baseline.
# IMPORTANT: this difference is taken WITHIN the same MC index so that
# correlated parameter draws cancel correctly -- this is the standard
# probabilistic-sensitivity-analysis convention.
# duration MUST be part of the join key. The raw disease_output.csv stores
# mortality (at duration=0), incidence (at duration=1), and prevalence (at
# duration>=2) as separate rows sharing the same (year, sex, agegroup, mc)
# tuple. Omitting duration from the merge produces a Cartesian explosion
# (each policy row joins every baseline duration row and vice versa).
# After the subtract we can safely sum across duration in step 3 because
# the zero cells for "wrong" duration levels fall out automatically
# (e.g. mrtl columns are 0 at duration>0; prvl columns are 0 at duration=0).
key_cols <- c("year", "sex", "agegroup", "mc", "duration")

base <- dis[scenario == baseline_scenario,
            c(key_cols, event_cols, death_cols), with = FALSE]
setnames(base, c(event_cols, death_cols),
         paste0(c(event_cols, death_cols), "_base"))

policy <- dis[scenario != baseline_scenario,
              c(key_cols, "scenario", event_cols, death_cols), with = FALSE]

prev <- merge(policy, base, by = key_cols, all.x = TRUE)
for (col in c(event_cols, death_cols)) {
  prev[, (paste0(col, "_prev")) := get(paste0(col, "_base")) - get(col)]
}
prev[, paste0(c(event_cols, death_cols), "_base") := NULL]
prev[, (c(event_cols, death_cols)) := NULL]

# Tidy long format -------------------------------------------------------------
prev_long <- melt(
  prev,
  id.vars = c("scenario", key_cols),
  variable.name = "outcome",
  value.name = "prevented"
)
prev_long[, outcome := sub("_prev$", "", as.character(outcome))]

# -----------------------------------------------------------------------------
# 3. AGGREGATE TO THE STRATA THE WORKBOOK REPORTS, PER MC
# -----------------------------------------------------------------------------
# We need every combination of:
#   sex      in {men, women, All}
#   agegroup in {30-49, 50-69, 70+, All}
# Each combination is built by SUMMING the prevented events within each MC,
# never by averaging or medianing. Quantiles are deferred to step 5.

# Normalise sex labels (some runs use Men/Women)
prev_long[, sex := tolower(sex)]
prev_long <- prev_long[sex %in% c("men", "women")]
prev_long[, agegroup := as.character(agegroup)]

agg_one <- function(d, by_sex, by_age) {
  out <- copy(d)
  if (!by_sex) out[, sex := "All"]
  if (!by_age) out[, agegroup := "All"]
  out[, .(prevented = sum(prevented)),
      by = .(scenario, year, sex, agegroup, mc, outcome)]
}

per_mc <- rbindlist(list(
  agg_one(prev_long, by_sex = TRUE,  by_age = TRUE ),
  agg_one(prev_long, by_sex = TRUE,  by_age = FALSE),
  agg_one(prev_long, by_sex = FALSE, by_age = TRUE ),
  agg_one(prev_long, by_sex = FALSE, by_age = FALSE)
), use.names = TRUE)
per_mc <- unique(per_mc)

# -----------------------------------------------------------------------------
# 4. PER-MC CUMULATIVE SUMS OVER YEARS (FIX FOR SHORTCUT 1)
# -----------------------------------------------------------------------------
setorder(per_mc, scenario, sex, agegroup, mc, outcome, year)
per_mc[, cumulative := cumsum(prevented),
       by = .(scenario, sex, agegroup, mc, outcome)]

# -----------------------------------------------------------------------------
# 4b. DERIVED PER-MC OUTCOMES (FIX FOR SHORTCUT 3)
# -----------------------------------------------------------------------------
# CVD CPP = CHD CPP + Stroke CPP, computed per MC. All-cause DPP = sum of all
# three death-cause DPPs per MC.
make_derived <- function(dt, components, new_name, var) {
  d <- dt[outcome %in% components,
          .(value = sum(get(var))),
          by = .(scenario, year, sex, agegroup, mc)]
  d[, outcome := new_name]
  setnames(d, "value", var)
  d
}

cvd_cumul <- make_derived(per_mc, c("chd_prvl", "stroke_prvl"),
                          "cvd_prvl", "cumulative")
cvd_annual <- make_derived(per_mc, c("chd_prvl", "stroke_prvl"),
                           "cvd_prvl", "prevented")
all_death_cumul <- make_derived(per_mc, c("chd_mrtl", "stroke_mrtl", "other_mrtl"),
                                "all_mrtl", "cumulative")
all_death_annual <- make_derived(per_mc, c("chd_mrtl", "stroke_mrtl", "other_mrtl"),
                                 "all_mrtl", "prevented")

derived <- merge(
  rbind(cvd_annual, all_death_annual),
  rbind(cvd_cumul,  all_death_cumul),
  by = c("scenario", "year", "sex", "agegroup", "mc", "outcome")
)
per_mc <- rbind(per_mc, derived, use.names = TRUE)

# CVD DPP = all-cause DPP minus non-CVD (other-cause) DPP, per MC. Arithmetically
# equivalent to chd_mrtl + stroke_mrtl in the current set-up (because all_mrtl is
# itself derived as the sum of chd + stroke + other), but computed as
# "total - non-CVD" to match the standard epidemiological convention and to
# remain correct if additional CVD sub-causes (e.g., heart failure) are ever
# added to the disease model without updating every call site.
all_mc   <- per_mc[outcome == "all_mrtl"]
other_mc <- per_mc[outcome == "other_mrtl",
                   .(scenario, year, sex, agegroup, mc,
                     other_prev = prevented, other_cum = cumulative)]
cvd_dpp  <- merge(all_mc, other_mc,
                  by = c("scenario", "year", "sex", "agegroup", "mc"))
cvd_dpp[, `:=`(prevented  = prevented  - other_prev,
               cumulative = cumulative - other_cum,
               outcome    = "cvd_mrtl")]
cvd_dpp[, c("other_prev", "other_cum") := NULL]
per_mc <- rbind(per_mc, cvd_dpp, use.names = TRUE)

# -----------------------------------------------------------------------------
# 5. QUANTILES ACROSS MCS  --  THE ONLY PLACE QUANTILES ARE TAKEN
# -----------------------------------------------------------------------------
qfun <- function(x) as.list(quantile(x, probs = quant_probs,
                                     names = FALSE, na.rm = TRUE))
# Format each percentile individually so integer values (10, 50, 90) stay
# un-padded. Vectorised format(quant_probs * 100, trim = TRUE) would align
# to 1 decimal place because the vector contains 2.5 and 97.5, producing
# names like p50_0 / p10_0 / p90_0 that would mismatch the hardcoded
# references in section 6 (wide summary) below.
qnames <- paste0("p",
                 vapply(quant_probs * 100,
                        function(x) sub("\\.", "_", format(x, trim = TRUE)),
                        character(1)))

summarise_mc <- function(value_col) {
  s <- per_mc[, qfun(get(value_col)),
              by = .(scenario, year, sex, agegroup, outcome)]
  setnames(s, paste0("V", seq_along(quant_probs)), qnames)
  s[, statistic := value_col]
  s
}

annual_summary <- summarise_mc("prevented")
cumul_summary  <- summarise_mc("cumulative")

# Re-label outcomes to the workbook's naming convention so the output is a
# drop-in replacement for cpp_dpp_out.csv.
outcome_relabel <- c(
  chd_prvl    = "cumul_chd_cpp",
  stroke_prvl = "cumul_stroke_cpp",
  cvd_prvl    = "cumul_cvd_cpp",
  chd_mrtl    = "cumul_chd_dpp",
  stroke_mrtl = "cumul_stroke_dpp",
  cvd_mrtl    = "cumul_cvd_dpp",
  other_mrtl  = "cumul_other_dpp",
  all_mrtl    = "cumul_all_dpp"
)
cumul_summary[, variable := outcome_relabel[outcome]]
cumul_summary[, outcome  := NULL]

annual_relabel <- sub("^cumul_", "annual_", outcome_relabel)
annual_summary[, variable := annual_relabel[outcome]]
annual_summary[, outcome  := NULL]

cpp_dpp_recomputed <- rbind(cumul_summary, annual_summary, use.names = TRUE)
setcolorder(cpp_dpp_recomputed,
            c("scenario", "year", "variable", qnames, "sex", "agegroup",
              "statistic"))
fwrite(cpp_dpp_recomputed,
       file.path(out_dir, "cpp_dpp_recomputed.csv"))

# -----------------------------------------------------------------------------
# 6. WIDE SUMMARY TABLE (FIX FOR SHORTCUT 2)
# -----------------------------------------------------------------------------
# One row per (scenario, outcome) at year_last, all/all stratification, with
# Men, Women, and Persons columns each carrying median + 95% UI. Persons is
# computed PER MC as men + women, then quantiled -- not as men_median +
# women_median.
final_year <- year_last
key_outcomes <- c("cumul_chd_cpp", "cumul_stroke_cpp", "cumul_cvd_cpp",
                  "cumul_chd_dpp", "cumul_stroke_dpp", "cumul_cvd_dpp",
                  "cumul_other_dpp", "cumul_all_dpp")

wide_block <- cpp_dpp_recomputed[
  year == final_year &
  agegroup == "All" &
  variable %in% key_outcomes
]

# Reshape so each row is (scenario, outcome) with men/women/Persons quantiles.
wide <- dcast(
  wide_block,
  scenario + variable ~ sex,
  value.var = c("p2_5", "p50", "p97_5")
)
setnames(wide,
         c("p2_5_men",  "p50_men",  "p97_5_men",
           "p2_5_women","p50_women","p97_5_women",
           "p2_5_All",  "p50_All",  "p97_5_All"),
         c("men_lo",    "men_med",  "men_hi",
           "women_lo",  "women_med","women_hi",
           "persons_lo","persons_med","persons_hi"))
setorder(wide, scenario, variable)
fwrite(wide, file.path(out_dir, "summary_recomputed.csv"))

# -----------------------------------------------------------------------------
# 7. COSTS (OPTIONAL, FIX FOR SHORTCUT 3 APPLIED TO MONETARY AGGREGATES)
# -----------------------------------------------------------------------------
if (do_costs) {
  # Per-MC cost = events_prevented * unit_cost, summed over years and ages,
  # then quantiled across MCs. Costs are kept sex-specific because the unit
  # costs differ by sex in the workbook.
  cost_long <- per_mc[
    outcome %in% c("chd_prvl", "stroke_prvl") &
    year == final_year & agegroup == "All" & sex %in% c("men", "women")
  ]
  cp <- cost_params
  cost_long[, unit_hosp := fcase(
    outcome == "chd_prvl"    & sex == "men",   cp$chd_hosp_men,
    outcome == "chd_prvl"    & sex == "women", cp$chd_hosp_women,
    outcome == "stroke_prvl" & sex == "men",   cp$stroke_hosp_men,
    outcome == "stroke_prvl" & sex == "women", cp$stroke_hosp_women
  )]
  cost_long[, unit_phc := fcase(
    outcome == "chd_prvl"    & sex == "men",   cp$chd_phc_men,
    outcome == "chd_prvl"    & sex == "women", cp$chd_phc_women,
    outcome == "stroke_prvl" & sex == "men",   cp$stroke_phc_men,
    outcome == "stroke_prvl" & sex == "women", cp$stroke_phc_women
  )]
  cost_long[, cost_per_mc := cumulative * (unit_hosp + unit_phc)]
  cost_mc <- cost_long[, .(cost_per_mc = sum(cost_per_mc)),
                       by = .(scenario, mc)]
  cost_summary <- cost_mc[, qfun(cost_per_mc), by = scenario]
  setnames(cost_summary, paste0("V", seq_along(quant_probs)), qnames)
  fwrite(cost_summary, file.path(out_dir, "costs_recomputed.csv"))
}

# -----------------------------------------------------------------------------
# 8. DONE
# -----------------------------------------------------------------------------
message("Wrote: ", file.path(out_dir, "cpp_dpp_recomputed.csv"))
message("Wrote: ", file.path(out_dir, "summary_recomputed.csv"))
if (do_costs)
  message("Wrote: ", file.path(out_dir, "costs_recomputed.csv"))
message("All quantiles computed across MC iterations as the FINAL step. ",
        "Persons = sum(men + women) per MC. ",
        "CVD = sum(CHD + Stroke) per MC. ",
        "Cumulative = cumsum within MC then quantiled. ",
        "No medians-of-medians, no sums-of-medians.")
# =============================================================================
