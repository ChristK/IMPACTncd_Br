# =============================================================================
# paper_tables.R
# -----------------------------------------------------------------------------
# Build the policy-impact tables that are reported in the paper:
#   * "Br Voluntary Targets" alone
#   * "Lowest World Regulatory Targets" alone
#   * "Br Voluntary + FOPL Warning 100%"
#   * "Br Voluntary + Potassium Salt 10% (Na+K)"
#
# Combination logic
# -----------------
# IMPACTncd_Br computes each policy scenario relative to the same Baseline
# (sc00, no policy). Two of the reported arms are standalone policies (Br
# Voluntary Targets, Lowest World Regulatory Targets); the other two stack
# FOPL Warning 100% and Potassium Salt 10% (Na+K) on top of the Br Voluntary
# Targets that are already in place. Because every scenario shares Baseline
# as its counterfactual, the additive total benefit of "Br Voluntary + X"
# is, per MC iteration:
#
#     prevented[BrVol + X]_mc  =  prevented[BrVol]_mc  +  prevented[X]_mc
#
# This assumes the two policies act independently (no interaction in the
# dose-response or in saturation of the prevented-event pool). For the salt
# reduction policies modelled here this is the standard convention used in
# the IMPACTncd literature; the supplementary text should state it explicitly.
# All arithmetic is done WITHIN each MC index so that correlated parameter
# draws cancel correctly when quantiles are taken at the very end.
#
# Methodology
# -----------
# Same correctness rules as recompute_cpp_dpp.R:
#   * cumulate within each MC, then quantile across MCs (no sum-of-medians)
#   * Persons = sum(men + women) per MC, then quantile (no median+median)
#   * CVD    = sum(CHD + Stroke) per MC, then quantile (no median+median)
#   * All quantiles are the FINAL step
#
# Output
# ------
# A single CSV `paper_tables.csv` and a printed table for each combined
# scenario, with median (95% UI) for CHD CPP, Stroke CPP, CVD CPP, CHD DPP,
# Stroke DPP, other-cause DPP, all-cause DPP, by sex (men, women, persons),
# at the end-of-horizon year (default 2038), all ages combined. Numeric
# values are rounded to 2 significant figures using signif().
#
# Author: Chris (epidemiology), drafted 2026-04-09. British English.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# -----------------------------------------------------------------------------
# 0. CONFIGURATION
# -----------------------------------------------------------------------------
disease_out_path <- "Output/disease_output.csv"
sim_params_path  <- "Output/simulation parameters.txt"
out_path         <- "Output/Tables/paper_tables.csv"
dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)

baseline_scenario <- "Baseline"
year_first        <- 2019
year_last         <- 2038
sig_digits        <- 2          # rounding for the paper tables

# --- Cost tariffs ------------------------------------------------------------
# 20-year per-case costs in Brazilian Reais (R$, 2019), taken directly from
# the Costs tab of `Results abstract new scenarios.xlsx` (rows 13-14,
# "Individual Costs (R$) -- 20y costs"). The components sum to the per-case
# total used in the workbook's cost columns.
#
# These unit costs are NOT sex-specific in the workbook -- the sex breakdown
# in the Costs tab arises only because the prevented event counts are
# sex-specific. We multiply per-MC prevented events by these tariffs and
# then quantile across MCs (no median-of-products).
#
# Conversion to 2019 PPP international dollars uses the workbook's stated
# factor 2.218 R$ per US$ (Costs tab cell J6).
unit_cost_chd_R  <- 444.38902708406613 +   # PHC
                    283.2336334857545  +   # Outpatient
                    1622.1562645093213 +   # Medications
                    2059.275123461133  +   # Informal
                    4306.06                # Inpatient (hospitalisation)
unit_cost_stk_R  <- 107.74848090124978 +
                    200.59431790177786 +
                    127.85522267206478 +
                    854.5365252596373  +
                    2015.76
ppp_R_per_USD    <- 2.218
unit_cost_chd_USD <- unit_cost_chd_R / ppp_R_per_USD
unit_cost_stk_USD <- unit_cost_stk_R / ppp_R_per_USD

# Combinations to report. Each list element is an arm of the paper's results:
# the name (list key) is how the arm will appear in paper_tables.csv, and the
# character vector (list value) is the set of scenarios whose prevented events
# are summed per MC iteration to build that arm. Single-element vectors are
# standalone policies (reported as-is); multi-element vectors stack policies
# on top of Br Voluntary Targets under the independence assumption documented
# in the "Combination logic" section at the top of this file.
combinations <- list(
  "Br Voluntary Targets"                   = c("Br Voluntary Targets"),
  "Lowest World Regulatory Targets"        = c("Lowest World Regulatory Targets"),
  "Br Voluntary + FOPL Warning 100%"       = c("Br Voluntary Targets",
                                               "FOPL Warning 100%"),
  "Br Voluntary + Potassium Salt 10% (Na+K)" = c("Br Voluntary Targets",
                                                 "Potassium Salt 10% (Na+K)")
)

event_cols <- c("chd_prvl", "stroke_prvl")
death_cols <- c("chd_mrtl", "stroke_mrtl", "other_mrtl")

# -----------------------------------------------------------------------------
# 1. LOAD AND VALIDATE
# -----------------------------------------------------------------------------
dis <- fread(disease_out_path)
setnames(dis, tolower(names(dis)))
needed <- c("year", "sex", "agegroup", "mc", "scenario",
            event_cols, death_cols)
missing <- setdiff(needed, names(dis))
if (length(missing) > 0)
  stop("disease_output.csv is missing required columns: ",
       paste(missing, collapse = ", "))

# Upscale event and death counts from simulation-sample scale to full
# Brazilian-population scale, matching what post_simulation_functions.R's
# process_disease_out() does at the melt stage. Without this step the
# recomputed numbers (and the cost aggregates built on top of them) are
# ~1/246 of the figures in cpp_dpp_out.csv and the workbook, because the
# simulation cohort is a ~400k sample of the ~100M Brazilian adults aged
# 30-79. pop_fraction is written into simulation parameters.txt by
# lifetable_engine.R; parse with the same regex used in
# post_simulation_functions.R so the two stay in sync.
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

# The raw disease_output.csv uses scenario CODES (sc00..sc12). Rename them
# to the friendly scenario labels used by the workbook/paper convention so
# that baseline_scenario ("Baseline") and the `combinations` map above
# (which reference friendly names like "Br Voluntary Targets" and
# "Potassium Salt 10% (Na+K)") all resolve correctly. design$scenarios is
# auto-derived from list.files("./Scenarios/sc*.R") in Scenarios/design.R.
#
# We ALWAYS re-source design.R and ALWAYS overwrite friendly_scenario_names,
# rather than guarding on exists(). Long-lived R sessions (e.g. re-running
# this script interactively after adding sc12.R to disk) otherwise carry a
# stale design$scenarios from before the new file existed, producing a
# length mismatch with the hardcoded 13-entry name list below. Re-sourcing
# design.R is safe because design.R's own `if (!exists("design"))` guard
# protects iteration_n / clusternumber / etc., while its unconditional
# `design$scenarios <- list.files(...)` assignment always refreshes.
# The friendly names MUST stay in sync with output.R:21-27.
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
       ". Add the new scenario's friendly label to paper_tables.R and ",
       "output.R:21-27, or remove the dangling Scenarios/sc*.R file.")
scenario_rename <- setNames(friendly_scenario_names, design$scenarios)
dis[, scenario := scenario_rename[scenario]]
if (anyNA(dis$scenario))
  stop("Unmapped scenario code(s) in disease_output.csv. Check that ",
       "design$scenarios and friendly_scenario_names are aligned with the ",
       "current file contents.")

dis[, sex := tolower(sex)]
dis <- dis[year >= year_first & year <= year_last & sex %in% c("men", "women")]

scen_have <- unique(dis$scenario)
needed_scen <- unique(unlist(combinations))
miss_scen <- setdiff(c(baseline_scenario, needed_scen), scen_have)
if (length(miss_scen) > 0)
  stop("disease_output.csv is missing scenarios: ",
       paste(miss_scen, collapse = ", "),
       "\nDid you run potassium_cra.R first to add 'Potassium Salt 10% (Na+K)'?")

# -----------------------------------------------------------------------------
# 2. PREVENTED EVENTS PER MC VERSUS BASELINE
# -----------------------------------------------------------------------------
# duration MUST be part of the join key. The raw disease_output.csv stores
# mortality (at duration=0), incidence (at duration=1), and prevalence (at
# duration>=2) as separate rows sharing the same (year, sex, agegroup, mc)
# tuple. Omitting duration from the merge produces a Cartesian explosion
# (each policy row joins every baseline duration row and vice versa).
# After the subtract, the step 4 all_ages aggregation collapses across
# both agegroup AND duration, which is correct because wrong-duration
# cells are zero in the raw CSV (mrtl=0 at duration>0; prvl=0 at duration=0).
key_cols <- c("year", "sex", "agegroup", "mc", "duration")
base <- dis[scenario == baseline_scenario,
            c(key_cols, event_cols, death_cols), with = FALSE]
setnames(base, c(event_cols, death_cols),
         paste0(c(event_cols, death_cols), "_base"))

policy <- dis[scenario != baseline_scenario,
              c(key_cols, "scenario", event_cols, death_cols), with = FALSE]
prev <- merge(policy, base, by = key_cols, all.x = TRUE)
for (col in c(event_cols, death_cols))
  prev[, (col) := get(paste0(col, "_base")) - get(col)]
prev[, paste0(c(event_cols, death_cols), "_base") := NULL]

# -----------------------------------------------------------------------------
# 3. BUILD COMBINED-SCENARIO PER-MC PREVENTED EVENTS
# -----------------------------------------------------------------------------
combine_one <- function(combo_label, component_scens) {
  d <- prev[scenario %in% component_scens]
  d <- d[, lapply(.SD, sum),
         by = key_cols,
         .SDcols = c(event_cols, death_cols)]
  d[, scenario := combo_label]
  d
}
combined <- rbindlist(Map(combine_one,
                          names(combinations), combinations),
                      use.names = TRUE)

# -----------------------------------------------------------------------------
# 4. AGGREGATE TO ALL AGES (KEEPING SEX), THEN ADD A PERSONS ROW PER MC
# -----------------------------------------------------------------------------
all_ages <- combined[, lapply(.SD, sum),
                     by = .(scenario, year, sex, mc),
                     .SDcols = c(event_cols, death_cols)]

# Persons per MC: sum men + women within the same MC, then treat as another sex
persons <- all_ages[, lapply(.SD, sum),
                    by = .(scenario, year, mc),
                    .SDcols = c(event_cols, death_cols)]
persons[, sex := "Persons"]

mc_table <- rbindlist(list(all_ages, persons), use.names = TRUE)

# -----------------------------------------------------------------------------
# 5. PER-MC CUMULATIVE OVER YEARS, AND DERIVED OUTCOMES
# -----------------------------------------------------------------------------
setorder(mc_table, scenario, sex, mc, year)
cumul <- mc_table[, lapply(.SD, cumsum),
                  by = .(scenario, sex, mc),
                  .SDcols = c(event_cols, death_cols)]
cumul[, year := mc_table$year]

# CVD CPP and all-cause DPP per MC (sums, not medians of sums)
cumul[, cvd_cpp := chd_prvl + stroke_prvl]
cumul[, all_dpp := chd_mrtl + stroke_mrtl + other_mrtl]

# CVD DPP = all-cause DPP minus non-CVD (other-cause) DPP, per MC.
# Arithmetically equivalent to chd_mrtl + stroke_mrtl in the current set-up
# (all_dpp is itself defined as the sum above), but computed as
# "total - non-CVD" to match the standard epidemiological convention and to
# remain correct if additional CVD sub-causes (e.g., heart failure) are ever
# added to the disease model without updating every call site.
cumul[, cvd_dpp := all_dpp - other_mrtl]

# -----------------------------------------------------------------------------
# 5b. PER-MC COST SAVINGS (Reais and 2019 PPP US$)
# -----------------------------------------------------------------------------
# Cost savings = prevented events x 20-year per-case cost. Computed per MC so
# the 95% UI on cost savings inherits the same correlation structure as the
# event UIs (avoids the workbook's "median(events) x unit cost" shortcut).
cumul[, cost_chd_R     := chd_prvl    * unit_cost_chd_R]
cumul[, cost_stroke_R  := stroke_prvl * unit_cost_stk_R]
cumul[, cost_cvd_R     := cost_chd_R  + cost_stroke_R]
cumul[, cost_chd_USD   := chd_prvl    * unit_cost_chd_USD]
cumul[, cost_stroke_USD:= stroke_prvl * unit_cost_stk_USD]
cumul[, cost_cvd_USD   := cost_chd_USD + cost_stroke_USD]

# -----------------------------------------------------------------------------
# 6. QUANTILES ACROSS MCS  --  THE ONLY PLACE QUANTILES ARE TAKEN
# -----------------------------------------------------------------------------
end <- cumul[year == year_last]
quant_cols <- c(chd_cpp        = "chd_prvl",
                stroke_cpp     = "stroke_prvl",
                cvd_cpp        = "cvd_cpp",
                chd_dpp        = "chd_mrtl",
                stroke_dpp     = "stroke_mrtl",
                cvd_dpp        = "cvd_dpp",
                other_dpp      = "other_mrtl",
                all_dpp        = "all_dpp",
                cost_chd_R     = "cost_chd_R",
                cost_stroke_R  = "cost_stroke_R",
                cost_cvd_R     = "cost_cvd_R",
                cost_chd_USD   = "cost_chd_USD",
                cost_stroke_USD= "cost_stroke_USD",
                cost_cvd_USD   = "cost_cvd_USD")

q <- function(x) {
  qs <- quantile(x, probs = c(0.025, 0.5, 0.975), names = FALSE, na.rm = TRUE)
  list(median = qs[2], lo = qs[1], hi = qs[3])
}

out_rows <- list()
for (lbl in names(quant_cols)) {
  src <- quant_cols[[lbl]]
  agg <- end[, q(get(src)), by = .(scenario, sex)]
  agg[, outcome := lbl]
  out_rows[[lbl]] <- agg
}
paper_long <- rbindlist(out_rows, use.names = TRUE)

# -----------------------------------------------------------------------------
# 7. ROUND TO 2 SIGNIFICANT FIGURES AND FORMAT AS "median (lo-hi)"
# -----------------------------------------------------------------------------
fmt <- function(m, l, h) {
  m <- signif(m, sig_digits)
  l <- signif(l, sig_digits)
  h <- signif(h, sig_digits)
  sprintf("%s (%s to %s)",
          format(m, big.mark = ",", scientific = FALSE, trim = TRUE),
          format(l, big.mark = ",", scientific = FALSE, trim = TRUE),
          format(h, big.mark = ",", scientific = FALSE, trim = TRUE))
}
paper_long[, value := fmt(median, lo, hi)]

paper_wide <- dcast(
  paper_long,
  scenario + outcome ~ sex,
  value.var = "value"
)

# Order outcomes the way the paper presents them
outcome_order <- c("chd_cpp", "stroke_cpp", "cvd_cpp",
                   "chd_dpp", "stroke_dpp", "cvd_dpp", "other_dpp", "all_dpp",
                   "cost_chd_R", "cost_stroke_R", "cost_cvd_R",
                   "cost_chd_USD", "cost_stroke_USD", "cost_cvd_USD")
paper_wide[, outcome := factor(outcome, levels = outcome_order)]
scenario_order <- names(combinations)
paper_wide[, scenario := factor(scenario, levels = scenario_order)]
setorder(paper_wide, scenario, outcome)
setcolorder(paper_wide, c("scenario", "outcome", "men", "women", "Persons"))
setnames(paper_wide,
         c("scenario", "outcome", "men", "women", "Persons"),
         c("Scenario", "Outcome", "Men",  "Women", "Persons"))

fwrite(paper_wide, out_path)

# -----------------------------------------------------------------------------
# 8. PRINT TO CONSOLE
# -----------------------------------------------------------------------------
cat("\n========================================================================\n")
cat("Paper tables  --  cumulative ", year_first, "-", year_last,
    ", all ages, median (95% UI), rounded to ",
    sig_digits, " significant figures\n", sep = "")
cat("========================================================================\n\n")
for (sc in scenario_order) {
  cat("### ", sc, "\n", sep = "")
  print(paper_wide[Scenario == sc,
                   .(Outcome, Men, Women, Persons)],
        row.names = FALSE)
  cat("\n")
}
message("Wrote: ", out_path)
# =============================================================================
