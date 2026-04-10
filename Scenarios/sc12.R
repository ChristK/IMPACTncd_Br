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

# =============================================================================
# sc12.R -- STUB for "Potassium Salt 10% (Na+K)"
# =============================================================================
#
# This file exists ONLY to register the scenario code "sc12" with the pipeline
# so that:
#   (a) initialisation.R's scenarios_list glob "sc*.R" picks it up and
#       lapply(scenarios_list, sys.source, ...) sources it during simulation,
#   (b) Scenarios/design.R's design$scenarios (which is built from
#       list.files("./Scenarios", pattern = "^sc.*\\.R$")) automatically
#       includes "sc12" so that process_disease_out()'s positional rename
#       `replace_from_table(dt, "scenario", design$scenarios, friendly_scenario_names)`
#       has a slot to map "sc12" -> "Potassium Salt 10% (Na+K)".
#
# It DOES NOT run any simulation logic. The rows for scenario "sc12" are not
# produced by the main simulation at all -- they are generated POST HOC by
# potassium_cra.R, which reads the "sc05" (Potassium Salt 10%) rows from the
# raw disease_output.csv, layers a comparative risk assessment for direct
# potassium-on-CVD effects (D'Elia et al 2011 meta-analysis) on top of the
# sodium-mediated SBP pathway, re-labels the modified rows as "sc12", and
# appends them to disease_output.csv.
#
# So the intended flow is:
#
#   1. Run Brazil_salt_policy.R (or Brazil_salt_policy_MONO.R). This sources
#      every file matching Scenarios/sc*.R including THIS FILE, which is a
#      no-op. disease_output.csv is written with rows for sc00..sc11 only.
#
#   2. Run potassium_cra.R. It reads the sc05 rows from disease_output.csv,
#      applies the potassium CRA, writes the modified rows back with
#      scenario = "sc12", and leaves every other scenario untouched.
#
#   3. Run output.R / post_simulation_functions.R as usual. The positional
#      code->friendly-name mapping handles "sc12" because this file exists
#      and friendly_scenario_names in output.R has a matching 13th entry.
#
# If you ever want "sc12" to be a real simulation scenario (i.e. have the
# main simulation populate it directly instead of the post-hoc CRA route),
# replace this stub with the real scenario logic -- follow sc05.R as a
# template. You will also then need to remove or rewrite potassium_cra.R.
#
# Loading this file must be cheap and side-effect-free. Do NOT add any POP
# manipulation, disease-model sourcing, or fwrite calls here.
# =============================================================================

cat("Loading sc12 (stub - populated post-hoc by potassium_cra.R)...\n")
