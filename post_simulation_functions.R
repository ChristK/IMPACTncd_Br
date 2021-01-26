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

pr <-
  c(0.5, 0.025, 0.975, 0.1, 0.9, 0.2, 0.8) # percentiles for summaries

if (!exists("design"))
  source("./Scenarios/design.R")
if (exists("out_xps"))
  rm(out_xps)
if (exists("out_lifecourse"))
  rm(out_lifecourse)
if (exists("POP"))
  rm(POP)

cat("Loading post-simulation functions...\n\n")
if (!require(BrazilSaltModelmisc)) {
  if (!require(remotes))
    install.packages("remotes")
  remotes::install_local("./Rpackage/Brazil_salt_model_pkg/")
  library(BrazilSaltModelmisc)
}

output_dir <- function(x = character(0))
  paste0("./Output/", x)

tt <- readLines(output_dir("simulation parameters.txt"))
pop_fraction <-
  as.numeric(substring(tt[[grep("^Population fraction ", tt)]], 40))

dependencies(c("future",
               "future.apply",
               "scales",
               "ggplot2",
               "cowplot",
               "data.table"),
             TRUE,
             FALSE,
             FALSE,
             FALSE)
theme_set(theme_cowplot())

process_disease_out <-
  function(dt,
           design,
           pop_fraction,
           friendly_scenario_names,
           per = 1e5) {
    # Some columns pop_size and _mrtl cols are only contain information for
    # duration =0. _prvl only contain info for duration > 0.
    # For efficiency, I am splitting on duration, delete unnecessary cols conditional
    # on duration, before melting.
    dt <- split(dt, by = "duration")

    # delete unnecessary cols
    dt <- lapply(dt, function(dt) {
      if (unique(dt$duration) == 0L) {
        dt[, grep("_prvl$", names(dt), value = TRUE) := NULL]
      } else {
        dt[, grep("_mrtl$|_size$", names(dt), value = TRUE) := NULL]
      }
    })

    dt <-
      lapply(dt,
             melt,
             c(setdiff(design$strata_for_outputs, "age"),
               "agegroup", "duration", "scenario", "mc"))
    dt <- rbindlist(dt)

    # clearly name incidence and prevalence without the need of duration (if > 1)
    if (max(dt$duration > 1L)) {
      # first sum to prvl including duration = 1
      tt <-
        dt[like(variable, "_prvl$"), .(sum(value), duration = 2L), by =
             c(setdiff(design$strata_for_outputs, "age"), "agegroup", "mc", "scenario", "variable")]
      dt[duration > 1L &
           like(variable, "_prvl$"), value := 0] # otherwise sum for agegroups would doublecount

      dt[tt, on = intersect(names(dt), names(tt)), value := i.V1]

      dt[duration == 1L & like(variable, "_prvl$"), variable :=
           gsub("_prvl$", "_incd", variable)] # rename to incd
    }

    if (max(dt$duration) > 2L) dt <- dt[duration <= 2L]
    dt[, duration := NULL]

    # upscale to the true Brazilian population
    dt[, value := round(value / pop_fraction)]

    # friendly scenario names
    replace_from_table(dt,
                       "scenario",
                       design$scenarios,
                       friendly_scenario_names)

    # to age groups
    # to_agegrp(dt, grp_width = 20L, max_age = design$ageH)
    # dt <- dt[, .(value = sum(value)), by =
    #            c(
    #              setdiff(design$strata_for_outputs, "age"),
    #              "agegroup",
    #              "mc",
    #              "scenario",
    #              "variable"
    #            )]

    # calculate rates
    pop <- dt[variable == "pop_size", ][, variable := NULL]
    setnames(pop, "value", "pop_size")
    tt <- dt[variable != "pop_size"]

    tt[pop, on = intersect(names(dt), names(pop)), value := per * x.value / i.pop_size]
    tt[, variable := paste0(variable, "_rate")]

    dt <- rbind(dt, tt)
    dt[pop, on = intersect(names(dt), names(pop)), pop_size := i.pop_size]

    invisible(dt)
  }

# Calculate quantiles TODO correct handling of rates
summarise_distr <-
  function(dt,
           design,
           probabilities,
           groupings = FALSE,
           rounding = FALSE
  ) {
    stopifnot(is.data.table(dt))
    past_keys <- key(dt)
    setkeyv(dt, "variable")
    if (groupings) { # TODO dynamic (not hardcoded)
      # Note: xps_out already in groupings
      out <- rbind(
        dt[, ifelse(grepl("_rate$", variable), weighted.mean(value, pop_size), sum(value)),
           by = c("scenario", "year", "variable", "mc")
           ][, fquantile_byid(V1, probabilities,
                              as.character(variable), rounding),
             by = c("scenario", "year")],
        dt[, ifelse(grepl("_rate$", variable), weighted.mean(value, pop_size), sum(value)),
           by =  c("scenario", "year", "sex", "variable", "mc")
           ][, fquantile_byid(V1, probabilities,
                              as.character(variable), rounding),
             by =  c("scenario", "year", "sex")],
        dt[, ifelse(grepl("_rate$", variable), weighted.mean(value, pop_size), sum(value)),
           by =  c("scenario", "year", "sex", "agegroup", "variable", "mc")
           ][, fquantile_byid(V1, probabilities,
                              as.character(variable), rounding),
             by =  c("scenario", "year", "sex", "agegroup")],
        fill = TRUE
      )
    } else {
      out <- # xps_out already in groupings. I don't use weighted mean cause xps_out have no pop_size
        dt[, ifelse(grepl("_rate$", variable), mean(value), sum(value)),
           by = c(setdiff(design$strata_for_outputs, "age"),
                  "agegroup",
                  "scenario", "variable", "mc")
           ][, fquantile_byid(V1, probabilities,
                              as.character(variable), rounding),
             by = c(setdiff(design$strata_for_outputs, "age"),
                    "agegroup",
                    "scenario")]
    }
    setkeyv(dt, past_keys)
    setnames(out, "V1", "variable")
    setnames(out, paste0("V", seq_along(pr) + 1L), scales::percent(pr))
    for (j in seq_len(ncol(out)))
      set(out, which(is.na(out[[j]])), j, "All")
    setkey(out, year, scenario)
    invisible(out)
  }


# Calculate cpp & dpp
calculate_cpp_dpp <-
  function(dt,
           design,
           probabilities,
           friendly_scenario_names,
           grouping = FALSE) {
    cpp_dpp_out <- dt[like(variable, "_prvl$|_mrtl$")]
    baseline    <-
      cpp_dpp_out[scenario == friendly_scenario_names[[1]]][, scenario := NULL]
    cpp_dpp_out  <-
      cpp_dpp_out[scenario != friendly_scenario_names[[1]]]

    cpp_dpp_out <-
      cpp_dpp_out[baseline, `:=` (value = i.value - x.value),
                  on = c("mc",
                         "variable",
                         "agegroup",
                         setdiff(design$strata_for_outputs, "age"))]

    # calculate cumulative cpp & dpp
    cpp_dpp_out[, variable := as.character(variable)]
    setkey(cpp_dpp_out, year, variable)
    cum_cpp_dpp_out <- copy(cpp_dpp_out)
    cum_cpp_dpp_out[, `:=` (value = cumsum(value),
                            pop_size = cumsum(pop_size),
                            variable = paste0("cumul_", variable)),
                    by = c("mc",
                    "variable", "scenario",
                    "agegroup",
                    setdiff(design$strata_for_outputs, c("age", "year")))]
    cpp_dpp_out <- rbind(cpp_dpp_out, cum_cpp_dpp_out)

    cpp_dpp_out <-
      summarise_distr(cpp_dpp_out, design, probabilities, grouping, TRUE)
    cpp_dpp_out[like(variable, "_prvl$"),
                variable := gsub("_prvl$", "_cpp", variable)]
    cpp_dpp_out[like(variable, "_mrtl$"),
                variable := gsub("_mrtl$", "_dpp", variable)]
  }

std_graph <-
  function(var,
           dt,
           subdir,
           y_axis_nam = var,
           filename = var) {
    gg <- ggplot(
      dt[variable == var],
      aes(
        x = year,
        y =    `50.0%`,
        ymin = `2.5%`,
        ymax = `97.5%`,
        col = scenario,
        fill = scenario
      )
    ) +
      geom_point(size = 1,
                 alpha = 5 / 5,
                 show.legend = FALSE) +
      geom_line(size = 1, alpha = 5 / 5) +
      geom_ribbon(alpha = 1 / 5,
                  linetype = 0,
                  show.legend = FALSE) +
      scale_x_continuous(name = "Year") +
      scale_y_continuous(name = y_axis_nam) +
      ggtitle(var)

    nam <- ""
    if (uniqueN(dt$sex) > 1L &&
        uniqueN(dt$agegroup) == 1L && unique(dt$agegroup) == "All") {
      gg <- gg + facet_grid(sex ~ .)
      nam <- "_sex"
    } else if (uniqueN(dt$sex) == 1L &&
               uniqueN(dt$agegroup) > 1L &&
               unique(dt$sex) == "All") {
      gg <- gg + facet_grid(. ~ agegroup)
      nam <- "_agegroup"
    } else if (uniqueN(dt$sex) > 1L &&
               uniqueN(dt$agegroup) > 1L) {
      gg <- gg + facet_grid(sex ~ agegroup)
      nam <- "_agegroup_sex"

    }
    cowplot::ggsave2(
      filename = paste0(filename, nam, ".png"),
      gg,
      height = 9,
      width = 16,
      units = "cm",
      path = output_dir(subdir)
    )
  }

# Aux function to get systems max memory
maxRAM <- function() {
  # in bytes
  if (Sys.info()[1] == "Windows") {
    out <- memory.size(NA) * 1024 ^ 2
  } else {
    out <-
      as.numeric(system("awk '/MemFree/ {print $2}' /proc/meminfo", intern = TRUE)) * 1024
  }
  out
}

pretty_numbers <- function(x, digits) {
  prettyNum(
    signif(x, digits),
    preserve.width = "individual",
    nsmall = 0L,
    scientific = F,
    big.mark = ",",
    big.interval = 3,
    drop0trailing = F
  )
}

# Capitalise first letter of a string
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

# Create output dirs
lapply(lapply(
  c("Summaries", "Tables", "Graphs", "Exposures", "Validation"),
  output_dir
), dir.create, showWarnings = FALSE)

# input$clusternumber <- ifelse(Sys.info()[1] == "Windows", 1L, input$clusternumber)
options(future.globals.maxSize = maxRAM() * 0.8)
