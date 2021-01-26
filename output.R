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
    "Media Campaigns", "5% Reduction in 10 Years") # in the order of sc00, sc01, sc02....

source(file = "./post_simulation_functions.R")
options(future.fork.enable = TRUE) # enable fork in Rstudio
plan(multiprocess)
# plan(sequential)

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
# CPP & DPP
# sex == "All" & agegroup == "All"
tt <- cpp_dpp_out[sex == "All" & agegroup == "All"]
future_lapply(unique(tt$variable), std_graph, tt, "Summaries", future.packages = c("ggplot2", "cowplot"))

# sex != "All" & agegroup == "All"
tt <- cpp_dpp_out[sex != "All" & agegroup == "All"]
future_lapply(unique(tt$variable), std_graph, tt, "Summaries", future.packages = c("ggplot2", "cowplot"))

# sex != "All" & agegroup != "All"
tt <- cpp_dpp_out[sex != "All" & agegroup != "All"]
future_lapply(unique(tt$variable), std_graph, tt, "Summaries", future.packages = c("ggplot2", "cowplot"))

# Cases & Deaths graphs for inspection
# sex == "All" & agegroup == "All"
tt <- disease_out_cmpct[sex == "All" & agegroup == "All"]
future_lapply(unique(tt$variable), std_graph, tt, "Graphs", future.packages = c("ggplot2", "cowplot"))

# sex != "All" & agegroup == "All"
tt <- disease_out_cmpct[sex != "All" & agegroup == "All"]
future_lapply(unique(tt$variable), std_graph, tt, "Graphs", future.packages = c("ggplot2", "cowplot"))

# sex != "All" & agegroup != "All"
tt <- disease_out_cmpct[sex != "All" & agegroup != "All"]
future_lapply(unique(tt$variable), std_graph, tt, "Graphs", future.packages = c("ggplot2", "cowplot"))

# Exposure graphs for inspection
# sex == "All" & agegroup == "All"
tt <- xps_out[sex == "All" & agegroup == "All"]
future_lapply(unique(tt$variable), std_graph, tt, "Exposures", future.packages = c("ggplot2", "cowplot"))

# sex != "All" & agegroup == "All"
tt <- xps_out[sex != "All" & agegroup == "All"]
future_lapply(unique(tt$variable), std_graph, tt, "Exposures", future.packages = c("ggplot2", "cowplot"))

# sex != "All" & agegroup != "All"
tt <- xps_out[sex != "All" & agegroup != "All"]
future_lapply(unique(tt$variable), std_graph, tt, "Exposures", future.packages = c("ggplot2", "cowplot"))


# Validation ------------------------------------------------
# absolute
deaths <-
  fread(
    "./Population_statistics/mortality_all_cause_estimates_2000_2016_age_sex_race_educ.csv",
    header = TRUE,
    skip = 0
  )[agegroup %in% agegrp_name(design$ageL, design$ageH, 5L, FALSE, FALSE, design$ageH)]
deaths <-
  deaths[, .(deaths = sum(deaths), pop_size = sum(pop)), by = .(year, sex)]
deaths_chd <-
  fread(
    "./Population_statistics/mortality_chd_estimates_2000_2016_age_sex_race_educ.csv",
    header = TRUE,
    skip = 0
  )[agegroup %in% agegrp_name(design$ageL, design$ageH, 5L, FALSE, FALSE, design$ageH)]
deaths_chd <-
  deaths_chd[, .(deaths = sum(deaths)), by = .(year, sex)]
deaths_stroke <-
  fread(
    "./Population_statistics/mortality_stroke_estimates_2000_2016_age_sex_race_educ.csv",
    header = TRUE,
    skip = 0
  )[agegroup %in% agegrp_name(design$ageL, design$ageH, 5L, FALSE, FALSE, design$ageH)]
deaths_stroke <-
  deaths_stroke[, .(deaths = sum(deaths)), by = .(year, sex)]

strata <- c("year", "sex")
deaths[deaths_chd, chd_mrtl := i.deaths, on = strata]
deaths[deaths_stroke, stroke_mrtl := i.deaths, on = strata]
deaths[, other_mrtl := deaths - chd_mrtl - stroke_mrtl]
deaths <- melt(deaths, 1:2, value.name = "50.0%")
deaths[, type := "Observed"]
replace_from_table(deaths, "sex", 1:2, c("men", "women"))
deaths[variable == "deaths", variable := "all_mrtl"]
tt <-
  disease_out_cmpct[sex != "All" &
                      agegroup == "All" &
                      scenario == friendly_scenario_names[[1]]]
tt[, type := "Modelled"]


tt <- rbind(tt, deaths, use.names = TRUE, fill = TRUE)
for (j in seq_len(ncol(tt))) {
  if (is.numeric(tt[[j]]))
    set(tt, which(is.na(tt[[j]])), j, 0)
}

future_lapply(unique(tt$variable), function(x) {
  if (grepl("_mrtl$|_size$", x)) {
    gg <- ggplot(tt[variable == x],
                 aes(
                   x = year,
                   y =    `50.0%`,
                   ymin = `2.5%`,
                   ymax = `97.5%`,
                   col = type,
                   fill = type
                 )) +
      geom_point(size = 1,
                 alpha = 5 / 5,
                 show.legend = F) +
      geom_line(size = 1, alpha = 5 / 5) +
      geom_ribbon(alpha = 1 / 5,
                  linetype = 0,
                  show.legend = F) +
      scale_x_continuous(name = "Year") +
      scale_y_continuous(name = x) +
      facet_grid(sex ~ .)+
      ggtitle(x)

    cowplot::ggsave2(
      filename = paste0("absolute_", x, "_validation.png"),
      gg,
      height = 9,
      width = 16,
      units = "in",
      path = output_dir("Validation")
    )
  }
}, future.packages = c("ggplot2", "cowplot"))

# rates
deaths <-
  fread(
    "./Population_statistics/mortality_all_cause_estimates_2000_2016_age_sex_race_educ.csv",
    header = TRUE,
    skip = 0
  )[agegroup %in% agegrp_name(design$ageL, design$ageH, 5L, FALSE, FALSE, design$ageH)]
deaths <-
  deaths[, .(deaths = 1e5 * sum(deaths) / sum(pop)), by = .(year, sex)]
deaths_chd <-
  fread(
    "./Population_statistics/mortality_chd_estimates_2000_2016_age_sex_race_educ.csv",
    header = TRUE,
    skip = 0
  )[agegroup %in% agegrp_name(design$ageL, design$ageH, 5L, FALSE, FALSE, design$ageH)]
deaths_chd <-
  deaths_chd[, .(deaths = 1e5 * sum(deaths) / sum(pop)), by = .(year, sex)]
deaths_stroke <-
  fread(
    "./Population_statistics/mortality_stroke_estimates_2000_2016_age_sex_race_educ.csv",
    header = TRUE,
    skip = 0
  )[agegroup %in% agegrp_name(design$ageL, design$ageH, 5L, FALSE, FALSE, design$ageH)]
deaths_stroke <-
  deaths_stroke[, .(deaths = 1e5 * sum(deaths) / sum(pop)), by = .(year, sex)]

strata <- c("year", "sex")
deaths[deaths_chd, chd_mrtl_rate := i.deaths, on = strata]
deaths[deaths_stroke, stroke_mrtl_rate := i.deaths, on = strata]
deaths[, other_mrtl_rate := deaths - chd_mrtl_rate - stroke_mrtl_rate]
# deaths[, deaths := NULL]
deaths <- melt(deaths, 1:2, value.name = "50.0%")
deaths[, type := "Observed"]
replace_from_table(deaths, "sex", 1:2, c("men", "women"))
deaths[variable=="deaths", variable := "all_mrtl_rate"]
tt <-
  disease_out_cmpct[sex != "All" &
                      agegroup == "All" &
                      scenario == friendly_scenario_names[[1]] &
                      like(variable, "_mrtl_rate$")]
tt[, type := "Modelled"]

# TODO rates need to calculate per each iteration
tt <- rbind(tt, deaths, use.names = TRUE, fill = TRUE)
for (j in seq_len(ncol(tt))) {
  if (is.numeric(tt[[j]]))
    set(tt, which(is.na(tt[[j]])), j, 0)
}

future_lapply(unique(tt$variable), function(x) {
  gg <- ggplot(tt[variable == x],
               aes(
                 x = year,
                 y =    `50.0%`,
                 ymin = `2.5%`,
                 ymax = `97.5%`,
                 col = type,
                 fill = type
               )) +
    geom_point(size = 1,
               alpha = 5 / 5,
               show.legend = FALSE) +
    geom_line(size = 1, alpha = 5 / 5) +
    geom_ribbon(alpha = 1 / 5,
                linetype = 0,
                show.legend = FALSE) +
    scale_x_continuous(name = "Year") +
    scale_y_continuous(name = paste0(x, " per 100,000")) +
    facet_grid(sex ~ .) +
    ggtitle(x)

  cowplot::ggsave2(
    filename = paste0("rate_", x, "_validation.png"),
    gg,
    height = 9,
    width = 16,
    units = "in",
    path = output_dir("Validation")
  )
}, future.packages = c("ggplot2", "cowplot"))
# Export tables --------------------------------------------------

# if (all(c(value(f1), value(f2), value(f3), value(f4)) == rep("finish", 4))) print("Finish exporting summaries")
