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

cat("Summarising output...\n")

friendly_scenario_names <-
  c("Baseline", "Br Voluntary Targets", "Br Regulatory Targets",
    "PAHO Regulatory Targets", "Lowest World Regulatory Targets",
    "Potassium Salt 10%", "FOPL Warning 100%", "FOPL Warning 70%",
    "FOPL Traffic Light 100%", "FOPL Traffic light 70%",
    "Media Campaigns", "5% Reduction in 10 Years",
    "Potassium Salt 10% (Na+K)") # in the order of sc00, sc01, sc02, ..., sc12

source(file = "./post_simulation_functions.R")
options(future.fork.enable = TRUE) # enable fork in Rstudio
plan(multicore) # was plan(multiprocess); multiprocess removed in future 1.32.0
# plan(sequential)

# Muffle benign graphics-device hygiene warnings from future workers.
# cowplot::ggsave2 necessarily opens a graphics device to render each plot;
# future's worker-hygiene check flags this as "opened the default graphics
# device" even though the files are written correctly. We match only the
# two specific known-benign messages, so real warnings still surface.
muffle_device_warnings <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    if (grepl("opened the default graphics device|must close any opened devices",
              conditionMessage(w))) {
      invokeRestart("muffleWarning")
    }
  })
}

disease_out <- fread(output_dir("disease_output.csv"))
disease_out[, all_mrtl := other_mrtl + stroke_mrtl + chd_mrtl]
disease_out <- process_disease_out(disease_out, design,
                                   pop_fraction, friendly_scenario_names)

# Calculate cases & deaths --------------------------------------------
disease_out_cmpct <-
  summarise_distr(disease_out, design, pr, TRUE, TRUE)
fwrite(disease_out_cmpct,
            output_dir("Summaries/disease_out_cmpct.csv"))


# Calculate CPPs, DPPs -------------------------------------------------
cpp_dpp_out <-
  calculate_cpp_dpp(disease_out, design, pr,
                    friendly_scenario_names, TRUE)
fwrite(cpp_dpp_out, output_dir("Summaries/cpp_dpp_out.csv"))


# Extract exposures trends -------------------------------------
xps_out <- fread(output_dir("xps_output.csv"))
xps_out <-
  melt(xps_out, c(
    setdiff(design$strata_for_outputs, "age"),
    "agegroup",
    "scenario",
    "mc"
  ))
replace_from_table(xps_out,
                   "scenario",
                   design$scenarios,
                   friendly_scenario_names)

xps_out <- summarise_distr(xps_out, design, pr, FALSE, FALSE)
fwrite(xps_out, output_dir("Exposures/rf.csv"))


# Graphs -------------------------------------------------------------------
# Render std_graph across the three (sex, agegroup) slices.
# Replaces 9 near-identical future_lapply blocks from the pre-modernisation version.
graph_slices <- function(dt, subdir) {
  slices <- list(
    overall     = dt[sex == "All"  & agegroup == "All"],
    by_sex      = dt[sex != "All"  & agegroup == "All"],
    by_sex_age  = dt[sex != "All"  & agegroup != "All"]
  )
  for (sl in slices) {
    if (nrow(sl) == 0L) next
    muffle_device_warnings(
      future_lapply(unique(sl$variable), std_graph, sl, subdir,
                    future.packages = c("ggplot2", "cowplot"))
    )
  }
}

graph_slices(cpp_dpp_out,       "Summaries") # CPP & DPP
graph_slices(disease_out_cmpct, "Graphs")    # Cases & Deaths for inspection
graph_slices(xps_out,           "Exposures") # Exposures for inspection


# Validation ------------------------------------------------
# Load observed mortality once, in both absolute and per-100k forms.
# Replaces two copy-pasted blocks that each re-read the same 3 CSVs.
load_observed_deaths <- function() {
  load_one <- function(path) {
    fread(path, header = TRUE)[
      agegroup %in% agegrp_name(design$ageL, design$ageH, 5L, FALSE, FALSE, design$ageH)
    ]
  }
  all_cause <- load_one("./Population_statistics/mortality_all_cause_estimates_2000_2016_age_sex_race_educ.csv")
  chd       <- load_one("./Population_statistics/mortality_chd_estimates_2000_2016_age_sex_race_educ.csv")
  stroke    <- load_one("./Population_statistics/mortality_stroke_estimates_2000_2016_age_sex_race_educ.csv")

  # absolute counts
  abs_dt <- all_cause[, .(deaths = sum(deaths), pop_size = sum(pop)), by = .(year, sex)]
  abs_dt[chd[,    .(deaths = sum(deaths)), by = .(year, sex)], chd_mrtl    := i.deaths, on = c("year", "sex")]
  abs_dt[stroke[, .(deaths = sum(deaths)), by = .(year, sex)], stroke_mrtl := i.deaths, on = c("year", "sex")]
  abs_dt[, other_mrtl := deaths - chd_mrtl - stroke_mrtl]

  # rates per 100k
  rate_dt <- all_cause[, .(deaths = 1e5 * sum(deaths) / sum(pop)), by = .(year, sex)]
  rate_dt[chd[,    .(deaths = 1e5 * sum(deaths) / sum(pop)), by = .(year, sex)], chd_mrtl_rate    := i.deaths, on = c("year", "sex")]
  rate_dt[stroke[, .(deaths = 1e5 * sum(deaths) / sum(pop)), by = .(year, sex)], stroke_mrtl_rate := i.deaths, on = c("year", "sex")]
  rate_dt[, other_mrtl_rate := deaths - chd_mrtl_rate - stroke_mrtl_rate]

  list(abs = abs_dt, rate = rate_dt)
}

# Given a wide-form observed table and a modelled subset, melt + join + zero-NAs
# into the shape make_validation_plot expects.
finalise_validation_frame <- function(observed_wide, modelled, deaths_rename) {
  obs <- melt(observed_wide, 1:2, value.name = "50.0%")
  obs[, type := "Observed"]
  replace_from_table(obs, "sex", 1:2, c("men", "women"))
  obs[variable == "deaths", variable := deaths_rename]

  mod <- copy(modelled)[, type := "Modelled"]
  out <- rbind(mod, obs, use.names = TRUE, fill = TRUE)
  for (j in seq_len(ncol(out)))
    if (is.numeric(out[[j]])) set(out, which(is.na(out[[j]])), j, 0)
  out
}

make_validation_plot <- function(x, dt, filename_prefix, y_label_suffix = "") {
  # See comment in std_graph: close any leaked devices before the worker
  # returns so future's device-hygiene check stays quiet.
  on.exit(while (!is.null(dev.list())) dev.off(), add = TRUE)

  gg <- ggplot(dt[variable == x],
               aes(
                 x = year,
                 y =    `50.0%`,
                 ymin = `2.5%`,
                 ymax = `97.5%`,
                 col = type,
                 fill = type
               )) +
    geom_point(size = 1, alpha = 1, show.legend = FALSE) +
    geom_line(linewidth = 1, alpha = 1) +
    geom_ribbon(alpha = 1 / 5, linetype = 0, show.legend = FALSE) +
    scale_x_continuous(name = "Year") +
    scale_y_continuous(name = paste0(x, y_label_suffix)) +
    facet_grid(sex ~ .) +
    ggtitle(x)

  cowplot::ggsave2(
    filename = paste0(filename_prefix, x, "_validation.png"),
    gg,
    height = 9,
    width = 16,
    units = "in",
    path = output_dir("Validation")
  )
}

observed <- load_observed_deaths()

base_modelled <- disease_out_cmpct[sex != "All" &
                                     agegroup == "All" &
                                     scenario == friendly_scenario_names[[1]]]

# absolute counts validation
abs_dt <- finalise_validation_frame(observed$abs, base_modelled, "all_mrtl")
muffle_device_warnings(
  future_lapply(
    grep("_mrtl$|_size$", unique(abs_dt$variable), value = TRUE),
    make_validation_plot, abs_dt, "absolute_",
    future.packages = c("ggplot2", "cowplot")
  )
)

# rates validation
# TODO rates need to calculate per each iteration
rate_dt <- finalise_validation_frame(
  observed$rate,
  base_modelled[like(variable, "_mrtl_rate$")],
  "all_mrtl_rate"
)
muffle_device_warnings(
  future_lapply(
    unique(rate_dt$variable), make_validation_plot, rate_dt, "rate_", " per 100,000",
    future.packages = c("ggplot2", "cowplot")
  )
)
