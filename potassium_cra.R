# =============================================================================
# potassium_cra.R
# -----------------------------------------------------------------------------
# Adds the direct potassium-on-CVD effect to the "Potassium Salt 10%" scenario
# of IMPACTncd_Br, using a comparative risk assessment (CRA) layered on top of
# the existing sodium-mediated SBP pathway already produced by the simulation.
#
# This implements the supplementary description (Supplementary file SR, para.
# 79): "The impact of increasing potassium intake was estimated using a
# comparative risk assessment approach in parallel to the sodium model,
# considering the meta-analysis by D'Elia et al (29) and similar lag times to
# sodium for the impact on cases and deaths prevented or postponed."
#
# Run AFTER the main simulation, before post_simulation_functions.R aggregates.
# It rewrites the rows of disease_output.csv corresponding to the Potassium
# Salt 10% scenario, leaving every other scenario untouched.
#
# Author: Chris (epidemiology), drafted 2026-04-09
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(qs)
})

# -----------------------------------------------------------------------------
# 1. PARAMETERS  --  edit here, do NOT hard-code further down
# -----------------------------------------------------------------------------

# --- File locations ----------------------------------------------------------
disease_out_path <- "Output/disease_output.csv"   # raw MC-level model output
xps_out_path     <- "Output/xps_output.csv"        # exposure outputs (sodium)
fixed_mc_path    <- "Scenarios/fixed_mc.qs"        # MC parameter draws
out_backup_path  <- "Output/disease_output_pre_kCRA.csv"

# --- Scenario label ----------------------------------------------------------
# The RAW disease_output.csv uses scenario CODES (sc00..sc11), not friendly
# names. The friendly-name mapping (sc05 -> "Potassium Salt 10%", etc.) is
# applied later by post_simulation_functions.R's process_disease_out() via
# replace_from_table(dt, "scenario", design$scenarios, friendly_scenario_names).
# This script operates BEFORE that rename, so we match on the scenario code.
# The new Na+K arm is assigned a fresh code "sc12" which is registered by a
# stub file Scenarios/sc12.R and by a 13th friendly-name entry in output.R.
k_scenario_name      <- "sc05"   # existing sodium-only arm (kept intact; becomes "Potassium Salt 10%")
k_scenario_name_new  <- "sc12"   # new combined arm added by this script     (becomes "Potassium Salt 10% (Na+K)")
sodium_baseline      <- "sc00"   # counterfactual scenario                   (becomes "Baseline")

# --- Potassium exposure shift ------------------------------------------------
# Sex-specific NaCl REDUCTIONS under the 10% KCl substitute arm, taken from
# Scenarios/sc05.R (men_salt_added_change = -0.706, women_salt_added_change =
# -0.509 g/d). IMPORTANT: these values are already the 10% fraction of the
# original discretionary salt -- i.e. the 10% w/w substitution factor is
# ALREADY baked in. Under a 1:1 mass substitution of NaCl by KCl in the
# substitute, these are therefore ALSO the masses of KCl added per day:
#   KCl added (men)   = 0.706 g/d
#   KCl added (women) = 0.509 g/d
# Do NOT re-apply the 10% factor below -- only convert KCl mass to elemental K.
discr_salt_g_men   <- 0.706   # NaCl removed = KCl added (men),   g/d
discr_salt_g_women <- 0.509   # NaCl removed = KCl added (women), g/d

# K elemental mass fraction of KCl = 39.10 / 74.55 = 0.5244
k_fraction_of_substitute <- 39.10 / 74.55   # ~0.5244

# Resulting added K intake (g/day) after full diffusion
k_added_g_men   <- discr_salt_g_men   * k_fraction_of_substitute
k_added_g_women <- discr_salt_g_women * k_fraction_of_substitute

# Diffusion ramp -- matches linear_diffusion(2019, 2025) used by the sodium
# template_policy_scenario.R for all reformulation scenarios.
diffusion_start_year <- 2019
diffusion_end_year   <- 2025

# --- D'Elia et al. 2011 meta-analytic effect sizes ---------------------------
# Stroke: RR 0.79 (95% CI 0.68-0.90) per 1.64 g/day higher K
# CHD:    RR ~0.93 (95% CI 0.87-0.99) per 1.64 g/day higher K
# Convert to log-RR per gram K and a normal SE on the log scale.
delia_unit_g            <- 1.64
delia_logrr_stroke_mean <- log(0.79) / delia_unit_g
delia_logrr_stroke_se   <- ((log(0.90) - log(0.68)) / (2 * 1.96)) / delia_unit_g
delia_logrr_chd_mean    <- log(0.93) / delia_unit_g
delia_logrr_chd_se      <- ((log(0.99) - log(0.87)) / (2 * 1.96)) / delia_unit_g
include_chd_pathway     <- TRUE   # set FALSE for stroke-only conservative run

# --- Lag distribution (matches the sodium arm) -------------------------------
# Supplementary file SR, paragraph 94: "1 + Binomial(9, (5-1)/9) to vary the
# lag time between 1 and 10 years (median 5 years)". We hardcode the resulting
# probability mass function so the potassium-CVD lag exactly mirrors the
# SBP-CVD lag used in the sodium arm. Index k of lag_kernel corresponds to a
# lag of k years (k = 1..10); index 0 is held at zero so apply_lag() can be
# indexed from 1.
.lag_p <- (5 - 1) / 9                             # = 4/9
.lag_pmf <- dbinom(0:9, size = 9, prob = .lag_p)  # for lag values 1..10
lag_kernel <- c(0, .lag_pmf)                      # length 11; lag_kernel[k+1] = P(lag = k)
stopifnot(abs(sum(lag_kernel) - 1) < 1e-6)

set.seed(20260409)

# -----------------------------------------------------------------------------
# 2. LOAD DATA
# -----------------------------------------------------------------------------
dis <- fread(disease_out_path)
setnames(dis, tolower(names(dis)))

# Backup once, never overwrite an existing backup
if (!file.exists(out_backup_path)) fwrite(dis, out_backup_path)

# Idempotence: drop any pre-existing rows for the post-hoc Na+K arm before
# we append fresh ones. This matters when this script is called from the
# standard pipeline (Brazil_salt_policy.R / Brazil_salt_policy_MONO.R) and
# the user re-runs the simulation. initialisation.R's continue-from-previous
# mechanism preserves disease_output.csv between runs that share the same
# scenarios snapshot, which means the previous run's sc12 rows are still in
# the file when this script is called the second time. Without this filter,
# the plain rbind at the end of section 8 would duplicate sc12 rows for the
# overlap MC range and produce double-counted aggregates downstream.
if (k_scenario_name_new %in% dis$scenario) {
  n_before <- nrow(dis)
  dis <- dis[scenario != k_scenario_name_new]
  message("Idempotence: dropped ", n_before - nrow(dis),
          " pre-existing '", k_scenario_name_new, "' rows before re-running CRA.")
}

# MC iteration values in this pipeline are NOT guaranteed to be 1..n_mc --
# Brazil_salt_policy.R's foreach writes whatever mc_iter integer the loop
# counter is at, and a partial / single-iter run (e.g. Brazil_salt_policy_MONO.R
# which uses mc_iter = 4L) can produce a CSV containing a single non-1 value.
# So we take the actual unique mc values from the data and use them as the
# join key, rather than assuming they are contiguous starting from 1.
mc_vals   <- dis[, sort(unique(mc))]
n_mc      <- length(mc_vals)
years_all <- sort(unique(dis$year))

# -----------------------------------------------------------------------------
# 3. PER-MC POTASSIUM LOG-RR DRAWS
# -----------------------------------------------------------------------------
# One draw per MC iteration so the K parameter uncertainty propagates through
# CPP/DPP UIs in the same way the sodium parameters do.
mc_draws <- data.table(
  mc           = mc_vals,
  beta_stroke  = rnorm(n_mc, delia_logrr_stroke_mean, delia_logrr_stroke_se),
  beta_chd     = if (include_chd_pathway)
                   rnorm(n_mc, delia_logrr_chd_mean, delia_logrr_chd_se)
                 else 0
)

# -----------------------------------------------------------------------------
# 4. POTASSIUM EXPOSURE SHIFT BY YEAR / SEX
# -----------------------------------------------------------------------------
# Linear diffusion ramp from 0 in (start-1) to 1 in end_year
diffusion <- function(y) {
  pmin(1, pmax(0, (y - (diffusion_start_year - 1)) /
                    (diffusion_end_year - diffusion_start_year + 1)))
}

k_shift <- CJ(year = years_all, sex = c("men", "women"))
k_shift[, dk_g := fifelse(sex == "men",
                          k_added_g_men   * diffusion(year),
                          k_added_g_women * diffusion(year))]

# -----------------------------------------------------------------------------
# 5. APPLY DISTRIBUTED LAG TO THE EXPOSURE SHIFT
# -----------------------------------------------------------------------------
# Convolve the K shift time series with the lag kernel, separately by sex.
apply_lag <- function(x, kernel) {
  out <- numeric(length(x))
  for (k in seq_along(kernel)) {
    if (kernel[k] == 0) next
    shifted <- c(rep(0, k - 1), x[seq_len(length(x) - (k - 1))])
    out <- out + kernel[k] * shifted
  }
  out
}

setorder(k_shift, sex, year)
k_shift[, dk_lagged := apply_lag(dk_g, lag_kernel), by = sex]

# -----------------------------------------------------------------------------
# 6. POTASSIUM PIF PER (year, sex, mc) FOR STROKE AND CHD
# -----------------------------------------------------------------------------
# PIF = 1 - exp(beta * delta_K)  (collapses to a population shift form because
# the K shift is the same for everyone in a stratum). For an individual-level
# CRA, swap this for sum(p_i * exp(beta * k_i_new)) / sum(p_i * exp(...)).
pif <- CJ(year = years_all, sex = c("men", "women"), mc = mc_vals)
pif <- merge(pif, k_shift[, .(year, sex, dk_lagged)], by = c("year", "sex"))
pif <- merge(pif, mc_draws, by = "mc")

pif[, pif_stroke := 1 - exp(beta_stroke * dk_lagged)]
pif[, pif_chd    := 1 - exp(beta_chd    * dk_lagged)]
pif[, c("beta_stroke", "beta_chd", "dk_lagged") := NULL]

# -----------------------------------------------------------------------------
# 7. APPLY THE PIF TO A COPY OF THE POTASSIUM SALT 10% ROWS
# -----------------------------------------------------------------------------
# IMPORTANT: we do NOT overwrite the existing Potassium Salt 10% rows. Instead
# we copy them, apply the PIF, and append them under a new scenario label
# (k_scenario_name_new) so the original sodium-only outputs remain untouched
# and downstream code (cpp_dpp aggregation, Total2 paste, supplementary tables)
# can still report the sodium-only Potassium Salt 10% values alongside the new
# combined Na+K version.
k_rows <- copy(dis[scenario == k_scenario_name])
stopifnot(nrow(k_rows) > 0)

k_rows <- merge(k_rows, pif, by = c("year", "sex", "mc"), all.x = TRUE)
k_rows[is.na(pif_stroke), pif_stroke := 0]
k_rows[is.na(pif_chd),    pif_chd    := 0]

# The raw disease_output.csv stores outcomes disjointly by `duration`. The
# authoritative source is summarise_hlp() in
# Rpackage/Brazil_salt_model_pkg/R/functions.R:1506-1533, which assembles the
# file: every *_mrtl column is hardcoded to duration = 0 ("summarise mortality
# (always recorded as duration 0)"), while *_prvl rows inherit their `duration`
# field from the POP column value itself -- and POP$chd_prvl / POP$stroke_prvl
# store YEARS SINCE DIAGNOSIS (see chd_model.R:80-86), so chd_prvl == 1 is a
# case diagnosed THIS year and chd_prvl == k is a case diagnosed k years ago.
# post_simulation_functions.R:67-79 confirms this layout downstream by dropping
# *_prvl at duration == 0 and dropping *_mrtl / *_size at duration > 0.
#
#     duration == 0  -> *_mrtl populated (deaths that year); *_prvl == 0
#     duration == 1  -> *_prvl populated (NEW / incident cases);  *_mrtl == 0
#     duration >= 2  -> *_prvl populated (CARRY-OVER prevalent cases from
#                       prior years, still alive);                *_mrtl == 0
#
# D'Elia's RR is a HAZARD ratio on incident CVD events, so the PIF must be
# applied ONLY to:
#   (a) incidence        -- *_prvl at duration == 1
#   (b) cause-specific mortality -- *_mrtl at duration == 0
#
# We must NOT multiply *_prvl at duration >= 2 by (1 - PIF): those rows are
# people who already had the disease in an earlier year, and raising K cannot
# retroactively un-cause their past stroke/CHD. The previous version of this
# script applied the PIF to ALL prvl rows, which silently deflated the
# carry-over stock and then got summed by recompute_cpp_dpp.R into the
# cumulative CPP -- double-counting the incidence reduction across every
# subsequent duration stratum. That was the dominant error in the prior run.
#
# The explicit `duration == k` gates below are self-documenting: the
# "wrong-duration" cells are zero in the raw CSV anyway (so an ungated multiply
# would arithmetically collapse to the same thing for duration 0 vs 1), but
# gating encodes the epidemiological intent and protects against future
# changes in the upstream CSV shape.
#
# Rounding back to integer preserves the storage class of the original column
# (see paper_tables.R / post_simulation_functions.R, which assume integer
# event counts); banker's rounding avoids systematic upward bias.

# (a) Incidence: *_prvl at duration == 1
k_rows[duration == 1L,
       stroke_prvl := as.integer(round(stroke_prvl * (1 - pif_stroke)))]
if (include_chd_pathway) {
  k_rows[duration == 1L,
         chd_prvl := as.integer(round(chd_prvl * (1 - pif_chd)))]
}

# (b) Cause-specific mortality: *_mrtl at duration == 0
k_rows[duration == 0L,
       stroke_mrtl := as.integer(round(stroke_mrtl * (1 - pif_stroke)))]
if (include_chd_pathway) {
  k_rows[duration == 0L,
         chd_mrtl := as.integer(round(chd_mrtl * (1 - pif_chd)))]
}

# Carry-over prevalence (duration >= 2) and non-CVD mortality are left
# untouched. D'Elia is CVD-specific, and K does not reduce the risk held by
# people already diagnosed.
k_rows[, c("pif_stroke", "pif_chd") := NULL]

# Re-label as the new combined-effect scenario
k_rows[, scenario := k_scenario_name_new]

# -----------------------------------------------------------------------------
# 8. WRITE BACK -- original rows preserved, combined arm appended
# -----------------------------------------------------------------------------
dis_out <- rbind(dis, k_rows)
setorder(dis_out, scenario, year, sex, agegroup, duration, mc)
fwrite(dis_out, disease_out_path)

message("Potassium CRA applied to '", k_scenario_name,
        "' (n_mc=", n_mc,
        ", CHD pathway=", include_chd_pathway,
        "). Backup at: ", out_backup_path)

# -----------------------------------------------------------------------------
# 9. NEXT STEPS (run separately, not in this script)
# -----------------------------------------------------------------------------
# After this script writes the new disease_output.csv, re-run the existing
# post-simulation pipeline to regenerate cpp_dpp_out.csv:
#
#   source("post_simulation_functions.R")
#   process_disease_out(...)
#   calculate_cpp_dpp(...)
#
# Then paste the updated Potassium Salt 10% rows into the Total2 sheet of
# Results abstract new scenarios.xlsx -- the Summary and Summary UI tabs do
# not need formula changes because they reference Total2 by row.
#
# VALIDATION CHECKS to run after re-aggregation:
#   * Stroke CPPs for Potassium Salt 10% should be larger than the sodium-only
#     run by ~10-25% (driven by D'Elia stroke beta and ~0.22-0.30 g/day K shift)
#   * CHD CPPs should rise by less (or not at all if include_chd_pathway=FALSE)
#   * 95% UIs should remain plausibly skewed; 2.5% bound > 0
#   * Non-CVD DPPs should be UNCHANGED versus the pre-CRA backup
# =============================================================================
