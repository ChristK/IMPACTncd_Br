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

# Load files --------------------------------------------------------------
cat("Generating life table...\n\n")
# for non hispanic black and whites (& other)


# file name for projection file
nam1 <-
  paste0(
    "./Population_statistics/",
    "mortality_projections_",
    paste(sort(design$fatal_diseases), collapse = "_"),
    "_",
    design$ageL,
    "_",
    design$ageH,
    "_",
    design$init_year + 2000L,
    "_",
    design$sim_horizon,
    ".fst"
  )


# All cause mortality ------------------------
if (file.exists(nam1)) {
  lifetable_all <- read_fst(nam1, as.data.table = TRUE)
} else {
  dependencies("demography")

  deaths <-
      fread(
        "./Population_statistics/mortality_all_cause_estimates_2000_2016_age_sex_race_educ.csv",
        header = TRUE,
        skip = 0
      )
  deaths <- deaths[, .(deaths = sum(deaths), pop = sum(pop)), by = .(year, agegroup, sex)]
  deaths_chd <-
    fread(
      "./Population_statistics/mortality_chd_estimates_2000_2016_age_sex_race_educ.csv",
      header = TRUE,
      skip = 0
    )
  deaths_chd <- deaths_chd[, .(deaths = sum(deaths), pop = sum(pop)), by = .(year, agegroup, sex)]
  deaths_stroke <-
      fread(
        "./Population_statistics/mortality_stroke_estimates_2000_2016_age_sex_race_educ.csv",
        header = TRUE,
        skip = 0
      )
  deaths_stroke <- deaths_stroke[, .(deaths = sum(deaths), pop = sum(pop)), by = .(year, agegroup, sex)]

  strata <- c("year", "sex", "agegroup")
  deaths[deaths_chd, chd_deaths := i.deaths, on = strata]
  deaths[is.na(chd_deaths), .N]
  deaths[deaths_stroke, stroke_deaths := i.deaths, on = strata]
  deaths[is.na(stroke_deaths), .N]

  deaths[, c("Mx_all", "Mx_chd", "Mx_stroke", "Mx_nonmodelled") := 0]
  deaths[pop > 0, Mx_all         := deaths / pop]
  deaths[pop > 0, Mx_chd         := chd_deaths / pop]
  deaths[pop > 0, Mx_stroke      := stroke_deaths / pop]
  deaths[pop > 0, Mx_nonmodelled := (deaths - chd_deaths - stroke_deaths) / pop]
  # deaths[, plot(year, Mx_all, main = paste0(agegroup, sex), ylim = c(0, 0.13)),
  #      by = .(sex, agegroup)]

  deaths[Mx_all == 0, counts(agegroup)]
  deaths[Mx_chd == 0, counts(agegroup)]
  deaths[Mx_stroke == 0, counts(agegroup)]
  deaths[Mx_nonmodelled == 0, counts(agegroup)]

  hor <- design$sim_horizon # maximum simulation horizon
  rate <- vector("list", 0)
  pop <- vector("list", 0)


  for (k in unique(deaths$sex)) {
        # Decompose mortality
        x1 <- dcast(deaths[sex == k, ],
                    agegroup ~ year, value.var = "Mx_all")
        x1[, agegroup := NULL]

        x2 <- dcast(deaths[sex == k, ],
                    agegroup ~ year, value.var = "pop")
        x2[, agegroup := NULL]

        nam <- paste0("Brazil_", k, "_allcause")
        rate[[nam]] <- as.matrix(x1)
        pop[[nam]] <- as.matrix(x2)
      }


  # demog data doesn't work on lists of matrices
  xx <- demogdata(
    rate[[1]],
    pop[[1]],
    c(seq(2, 79, 5), 90),
    sort(unique(deaths$year)),
    "mortality",
    paste0("Brazil"),
    names(rate[1]),
    lambda = 0
  )

  # work around of above limitation
  xx$rate <- rate
  xx$pop  <- pop
  xx$name <- names(rate)

  xx <-
    smooth.demogdata(xx, age.grid = (design$ageL - 10L):(design$ageH + 5L), obs.var = "theoretical") # age 30:79 does not forecast

  # xxx <-
  #   xx$pop # to be used later to extrapolate population aged >84

  mort.fit <-
    coherentfdm(xx,
                10,
                10,
                method = "rapca",
                weight = T,
                # beta = 0.2,
                max.age = (design$ageH + 5L)) # weight is the most important arguement

  mortf   <- forecast(mort.fit, h = hor, max.d = 1, level = 99)
  mortf60 <- forecast(mort.fit, h = hor, max.d = 1, level = 60)
  mortf70 <- forecast(mort.fit, h = hor, max.d = 1, level = 70)
  mortf80 <- forecast(mort.fit, h = hor, max.d = 1, level = 80)
  mortf90 <- forecast(mort.fit, h = hor, max.d = 1, level = 90)

  # produce lui & uui
  mortf.1  <- mortf.99 <- mortf
  mortf.40 <- mortf.60 <- mortf60
  mortf.30 <- mortf.70 <- mortf70
  mortf.20 <- mortf.80 <- mortf80
  mortf.10 <- mortf.90 <- mortf90

  output    <- vector("list", length(mortf) - 2L)
  output.1  <- vector("list", length(mortf) - 2L)
  output.99 <- vector("list", length(mortf) - 2L)
  output.40 <- vector("list", length(mortf) - 2L)
  output.60 <- vector("list", length(mortf) - 2L)
  output.30 <- vector("list", length(mortf) - 2L)
  output.70 <- vector("list", length(mortf) - 2L)
  output.20 <- vector("list", length(mortf) - 2L)
  output.80 <- vector("list", length(mortf) - 2L)
  output.10 <- vector("list", length(mortf) - 2L)
  output.90 <- vector("list", length(mortf) - 2L)

  strata <- c("sex", "type")
  for (ii in 1:(length(mortf) - 2)) {
    mortf.1[[ii]]$rate[[1]]  <- mortf[[ii]]$rate$lower
    mortf.99[[ii]]$rate[[1]] <- mortf[[ii]]$rate$upper
    mortf.40[[ii]]$rate[[1]] <- mortf60[[ii]]$rate$lower
    mortf.60[[ii]]$rate[[1]] <- mortf60[[ii]]$rate$upper
    mortf.30[[ii]]$rate[[1]] <- mortf70[[ii]]$rate$lower
    mortf.70[[ii]]$rate[[1]] <- mortf70[[ii]]$rate$upper
    mortf.20[[ii]]$rate[[1]] <- mortf80[[ii]]$rate$lower
    mortf.80[[ii]]$rate[[1]] <- mortf80[[ii]]$rate$upper
    mortf.10[[ii]]$rate[[1]] <- mortf90[[ii]]$rate$lower
    mortf.90[[ii]]$rate[[1]] <- mortf90[[ii]]$rate$upper

    output[[ii]] <-
      as.data.table(mortf[[ii]]$rate[[1]],
                    keep.rownames = T)[, `:=`(type = names(mortf[ii]))]
    output[[ii]][, (strata) :=
                   tstrsplit(type, "_", fixed = TRUE, keep = 2:(length(strata) +1L))]
    output.1[[ii]] <-
      as.data.table(mortf.1[[ii]]$rate[[1]],
                    keep.rownames = T)[, `:=`(type = names(mortf.1[ii]))]
    output.1[[ii]][, (strata) :=
                       tstrsplit(type, "_", fixed = TRUE, keep = 2:(length(strata) +1L))]
    output.99[[ii]] <-
      as.data.table(mortf.99[[ii]]$rate[[1]],
                    keep.rownames = T)[, `:=`(type = names(mortf.99[ii]))]
    output.99[[ii]][, (strata) :=
                       tstrsplit(type, "_", fixed = TRUE, keep = 2:(length(strata) +1L))]
    output.40[[ii]] <-
      as.data.table(mortf.40[[ii]]$rate[[1]],
                    keep.rownames = T)[, `:=`(type = names(mortf.40[ii]))]
    output.40[[ii]][, (strata) :=
                     tstrsplit(type, "_", fixed = TRUE, keep = 2:(length(strata) +1L))]
    output.60[[ii]] <-
      as.data.table(mortf.60[[ii]]$rate[[1]],
                    keep.rownames = T)[, `:=`(type = names(mortf.60[ii]))]
    output.60[[ii]][, (strata) :=
                      tstrsplit(type, "_", fixed = TRUE, keep = 2:(length(strata) +1L))]
    output.30[[ii]] <-
      as.data.table(mortf.30[[ii]]$rate[[1]],
                    keep.rownames = T)[, `:=`(type = names(mortf.30[ii]))]
    output.30[[ii]][, (strata) :=
                      tstrsplit(type, "_", fixed = TRUE, keep = 2:(length(strata) +1L))]
    output.70[[ii]] <-
      as.data.table(mortf.70[[ii]]$rate[[1]],
                    keep.rownames = T)[, `:=`(type = names(mortf.70[ii]))]
    output.70[[ii]][, (strata) :=
                      tstrsplit(type, "_", fixed = TRUE, keep = 2:(length(strata) +1L))]
    output.20[[ii]] <-
      as.data.table(mortf.20[[ii]]$rate[[1]],
                    keep.rownames = T)[, `:=`(type = names(mortf.20[ii]))]
    output.20[[ii]][, (strata) :=
                      tstrsplit(type, "_", fixed = TRUE, keep = 2:(length(strata) +1L))]
    output.80[[ii]] <-
      as.data.table(mortf.80[[ii]]$rate[[1]],
                    keep.rownames = T)[, `:=`(type = names(mortf.80[ii]))]
    output.80[[ii]][, (strata) :=
                      tstrsplit(type, "_", fixed = TRUE, keep = 2:(length(strata) +1L))]
    output.10[[ii]] <-
      as.data.table(mortf.10[[ii]]$rate[[1]],
                    keep.rownames = T)[, `:=`(type = names(mortf.10[ii]))]
    output.10[[ii]][, (strata) :=
                      tstrsplit(type, "_", fixed = TRUE, keep = 2:(length(strata) +1L))]
    output.90[[ii]] <-
      as.data.table(mortf.90[[ii]]$rate[[1]],
                    keep.rownames = T)[, `:=`(type = names(mortf.90[ii]))]
    output.90[[ii]][, (strata) :=
                      tstrsplit(type, "_", fixed = TRUE, keep = 2:(length(strata) +1L))]
  }

  output    <- rbindlist(output)
  output.1  <- rbindlist(output.1)
  output.99 <- rbindlist(output.99)
  output.40 <- rbindlist(output.40)
  output.60 <- rbindlist(output.60)
  output.30 <- rbindlist(output.30)
  output.70 <- rbindlist(output.70)
  output.20 <- rbindlist(output.20)
  output.80 <- rbindlist(output.80)
  output.10 <- rbindlist(output.10)
  output.90 <- rbindlist(output.90)

  strata <- c("rn", "sex", "type")
  output <- melt(output,
                 id.vars = strata,
                 value.name = "mx_total")
  output.1 <-
    melt(output.1,
         id.vars = strata,
         value.name = "mx_total_1")
  output.99 <-
    melt(output.99,
         id.vars = strata,
         value.name = "mx_total_99")
  output.10 <-
    melt(output.10,
         id.vars = strata,
         value.name = "mx_total_10")
  output.20 <-
    melt(output.20,
         id.vars = strata,
         value.name = "mx_total_20")
  output.30 <-
    melt(output.30,
         id.vars = strata,
         value.name = "mx_total_30")
  output.40 <-
    melt(output.40,
         id.vars = strata,
         value.name = "mx_total_40")
  output.60 <-
    melt(output.60,
         id.vars = strata,
         value.name = "mx_total_60")
  output.70 <-
    melt(output.70,
         id.vars = strata,
         value.name = "mx_total_70")
  output.80 <-
    melt(output.80,
         id.vars = strata,
         value.name = "mx_total_80")
  output.90 <-
    melt(output.90,
         id.vars = strata,
         value.name = "mx_total_90")

  strata <- c("rn", "sex", "type", "variable")
  output[output.1,  mx_total_1  := i.mx_total_1,  on = strata]
  output[output.99, mx_total_99 := i.mx_total_99, on = strata]
  output[output.10, mx_total_10 := i.mx_total_10, on = strata]
  output[output.20, mx_total_20 := i.mx_total_20, on = strata]
  output[output.30, mx_total_30 := i.mx_total_30, on = strata]
  output[output.40, mx_total_40 := i.mx_total_40, on = strata]
  output[output.60, mx_total_60 := i.mx_total_60, on = strata]
  output[output.70, mx_total_70 := i.mx_total_70, on = strata]
  output[output.80, mx_total_80 := i.mx_total_80, on = strata]
  output[output.90, mx_total_90 := i.mx_total_90, on = strata]

  test <- copy(xx$rate)
  original <- vector("list", length(test))

  for (ii in 1:(length(test))) {
    original[[ii]] <-
      as.data.table(test[[ii]], keep.rownames = TRUE)[, `:=`(type = names(test[ii]))]
    original[[ii]][, c("sex", "type") := tstrsplit(type, "_", fixed = TRUE, keep =
                                                             2:3)]
  }
  original <- rbindlist(original)
  original <-
    melt(original,
         id.vars = c("rn", "sex", "type"),
         value.name = "mx_total")
  original[, paste0("mx_total_",
                    c(1, 99, 10, 20, 30, 40, 60, 70, 80, 90)) :=
             mx_total]

  lifetable_all <-
    rbind(original, output, use.names = TRUE, fill = TRUE)
  lifetable_all[, age := as.integer(rn)]
  lifetable_all[, year := as.integer(as.character(variable))]
  lifetable_all <-
    lifetable_all[between(age, design$ageL, design$ageH) &
                    between(year, design$init_year + 2000L,
                            (design$init_year + 2000L + design$sim_horizon))]

  # lifetable_all[, quan := runif(1)]
  # lifetable_all[year >= 2030, quan := runif(1)]
  # #lifetable_all[year >= 2050, quan := runif(1)]
  #
  # lifetable_all[, quan := roll_mean_left(quan, 9), by = .(age, sex, race)]
  #
  # ggplot(lifetable_all[
  #   (age %% 10) == 5 & between(age, 20, 79) & year > 1000, ],
  #   aes(
  #     y = mx_total,
  #     ymin = mx_total_1,
  #     ymax = mx_total_99,
  #     x = year
  #   )) +
  #   geom_pointrange() +
  #   facet_grid(age~sex + race3 + educ, scales = "free")




  # cause specific mortality
  rate <- vector("list", 0)
  pop <- vector("list", 0)

  for (k in unique(deaths$sex)) {
    # Decompose mortality
    x1 <- dcast(deaths[sex == k,],
                agegroup ~ year, value.var = "Mx_nonmodelled")
    x1[, agegroup := NULL]
    for (j in 1:ncol(x1))
      set(x1, which(x1[[j]] == 0), j, 1e-8)

    x2 <- dcast(deaths[sex == k,],
                agegroup ~ year, value.var = "pop")
    x2[, agegroup := NULL]

    nam <- paste0("Brazil_", k, "_nonmodelled")
    rate[[nam]] <- as.matrix(x1)
    pop[[nam]] <- as.matrix(x2)

    x1 <- dcast(deaths[sex == k,],
                agegroup ~ year, value.var = "Mx_chd")
    x1[, agegroup := NULL]
    for (j in 1:ncol(x1))
      set(x1, which(x1[[j]] == 0), j, 1e-8)

    nam <- paste0("Brazil_", k, "_chd")
    rate[[nam]] <- as.matrix(x1)
    pop[[nam]] <- as.matrix(x2)

    x1 <- dcast(deaths[sex == k,],
                agegroup ~ year, value.var = "Mx_stroke")
    x1[, agegroup := NULL]
    for (j in 1:ncol(x1))
      set(x1, which(x1[[j]] == 0), j, 1e-8)

    nam <- paste0("Brazil_", k, "_stroke")
    rate[[nam]] <- as.matrix(x1)
    pop[[nam]] <- as.matrix(x2)
  }

  # demog data doesn't work on lists of matrices
  xx <- demogdata(
    rate[[1]],
    pop[[1]],
    c(seq(2, 79, 5), 90),
    sort(unique(deaths$year)),
    "mortality",
    paste0("Brazil"),
    names(rate[1]),
    lambda = 0
  )

  # work around of above limitation
  xx$rate <- rate
  xx$pop  <- pop
  xx$name <- names(rate)

  xx <-
    smooth.demogdata(xx, age.grid = (design$ageL - 10L):(design$ageH + 5L), obs.var = "theoretical")

  mort.fit <-
    coherentfdm(
      xx,
      10,
      10,
      method = "rapca",
      weight = T,
      # beta = 0.2,
      max.age = (design$ageH + 5L)
    ) # weight is the most important arguement
  mortf   <- forecast(mort.fit, h = hor, max.d = 1, level = 99)
  mortf60 <- forecast(mort.fit, h = hor, max.d = 1, level = 60)
  mortf70 <- forecast(mort.fit, h = hor, max.d = 1, level = 70)
  mortf80 <- forecast(mort.fit, h = hor, max.d = 1, level = 80)
  mortf90 <- forecast(mort.fit, h = hor, max.d = 1, level = 90)

  # produce lui & uui
  mortf.1  <- mortf.99 <- mortf
  mortf.40 <- mortf.60 <- mortf60
  mortf.30 <- mortf.70 <- mortf70
  mortf.20 <- mortf.80 <- mortf80
  mortf.10 <- mortf.90 <- mortf90

  output    <- vector("list", length(mortf) - 2)
  output.1  <- vector("list", length(mortf) - 2)
  output.99 <- vector("list", length(mortf) - 2)
  output.40 <- vector("list", length(mortf) - 2)
  output.60 <- vector("list", length(mortf) - 2)
  output.30 <- vector("list", length(mortf) - 2)
  output.70 <- vector("list", length(mortf) - 2)
  output.20 <- vector("list", length(mortf) - 2)
  output.80 <- vector("list", length(mortf) - 2)
  output.10 <- vector("list", length(mortf) - 2)
  output.90 <- vector("list", length(mortf) - 2)

  for (ii in 1:(length(mortf) - 2)) {
    mortf.1[[ii]]$rate[[1]]  <- mortf[[ii]]$rate$lower
    mortf.99[[ii]]$rate[[1]] <- mortf[[ii]]$rate$upper
    mortf.40[[ii]]$rate[[1]] <- mortf60[[ii]]$rate$lower
    mortf.60[[ii]]$rate[[1]] <- mortf60[[ii]]$rate$upper
    mortf.30[[ii]]$rate[[1]] <- mortf70[[ii]]$rate$lower
    mortf.70[[ii]]$rate[[1]] <- mortf70[[ii]]$rate$upper
    mortf.20[[ii]]$rate[[1]] <- mortf80[[ii]]$rate$lower
    mortf.80[[ii]]$rate[[1]] <- mortf80[[ii]]$rate$upper
    mortf.10[[ii]]$rate[[1]] <- mortf90[[ii]]$rate$lower
    mortf.90[[ii]]$rate[[1]] <- mortf90[[ii]]$rate$upper

    strata <- c("sex", "type")
    output[[ii]] <-
      as.data.table(mortf[[ii]]$rate[[1]],
                    keep.rownames = T)[, `:=`(type = names(mortf[ii]))]
    output[[ii]][, (strata) :=
                   tstrsplit(type, "_", fixed = TRUE, keep = 2:(length(strata) + 1L))]
    output.1[[ii]] <-
      as.data.table(mortf.1[[ii]]$rate[[1]],
                    keep.rownames = T)[, `:=`(type = names(mortf.1[ii]))]
    output.1[[ii]][, (strata) :=
                     tstrsplit(type, "_", fixed = TRUE, keep = 2:(length(strata) + 1L))]
    output.99[[ii]] <-
      as.data.table(mortf.99[[ii]]$rate[[1]],
                    keep.rownames = T)[, `:=`(type = names(mortf.99[ii]))]
    output.99[[ii]][, (strata) :=
                      tstrsplit(type, "_", fixed = TRUE, keep = 2:(length(strata) + 1L))]
    output.40[[ii]] <-
      as.data.table(mortf.40[[ii]]$rate[[1]],
                    keep.rownames = T)[, `:=`(type = names(mortf.40[ii]))]
    output.40[[ii]][, (strata) :=
                      tstrsplit(type, "_", fixed = TRUE, keep = 2:(length(strata) + 1L))]
    output.60[[ii]] <-
      as.data.table(mortf.60[[ii]]$rate[[1]],
                    keep.rownames = T)[, `:=`(type = names(mortf.60[ii]))]
    output.60[[ii]][, (strata) :=
                      tstrsplit(type, "_", fixed = TRUE, keep = 2:(length(strata) + 1L))]
    output.30[[ii]] <-
      as.data.table(mortf.30[[ii]]$rate[[1]],
                    keep.rownames = T)[, `:=`(type = names(mortf.30[ii]))]
    output.30[[ii]][, (strata) :=
                      tstrsplit(type, "_", fixed = TRUE, keep = 2:(length(strata) + 1L))]
    output.70[[ii]] <-
      as.data.table(mortf.70[[ii]]$rate[[1]],
                    keep.rownames = T)[, `:=`(type = names(mortf.70[ii]))]
    output.70[[ii]][, (strata) :=
                      tstrsplit(type, "_", fixed = TRUE, keep = 2:(length(strata) + 1L))]
    output.20[[ii]] <-
      as.data.table(mortf.20[[ii]]$rate[[1]],
                    keep.rownames = T)[, `:=`(type = names(mortf.20[ii]))]
    output.20[[ii]][, (strata) :=
                      tstrsplit(type, "_", fixed = TRUE, keep = 2:(length(strata) + 1L))]
    output.80[[ii]] <-
      as.data.table(mortf.80[[ii]]$rate[[1]],
                    keep.rownames = T)[, `:=`(type = names(mortf.80[ii]))]
    output.80[[ii]][, (strata) :=
                      tstrsplit(type, "_", fixed = TRUE, keep = 2:(length(strata) + 1L))]
    output.10[[ii]] <-
      as.data.table(mortf.10[[ii]]$rate[[1]],
                    keep.rownames = T)[, `:=`(type = names(mortf.10[ii]))]
    output.10[[ii]][, (strata) :=
                      tstrsplit(type, "_", fixed = TRUE, keep = 2:(length(strata) + 1L))]
    output.90[[ii]] <-
      as.data.table(mortf.90[[ii]]$rate[[1]],
                    keep.rownames = T)[, `:=`(type = names(mortf.90[ii]))]
    output.90[[ii]][, (strata) :=
                      tstrsplit(type, "_", fixed = TRUE, keep = 2:(length(strata) + 1L))]
  }

  output    <- rbindlist(output)
  output.1  <- rbindlist(output.1)
  output.99 <- rbindlist(output.99)
  output.40 <- rbindlist(output.40)
  output.60 <- rbindlist(output.60)
  output.30 <- rbindlist(output.30)
  output.70 <- rbindlist(output.70)
  output.20 <- rbindlist(output.20)
  output.80 <- rbindlist(output.80)
  output.10 <- rbindlist(output.10)
  output.90 <- rbindlist(output.90)

  strata <- c("rn", "sex", "type")
  output <- melt(output,
                 id.vars = strata,
                 value.name = "mx")
  output.1 <-
    melt(output.1,
         id.vars = strata,
         value.name = "mx_1")
  output.99 <-
    melt(output.99,
         id.vars = strata,
         value.name = "mx_99")
  output.10 <-
    melt(output.10,
         id.vars = strata,
         value.name = "mx_10")
  output.20 <-
    melt(output.20,
         id.vars = strata,
         value.name = "mx_20")
  output.30 <-
    melt(output.30,
         id.vars = strata,
         value.name = "mx_30")
  output.40 <-
    melt(output.40,
         id.vars = strata,
         value.name = "mx_40")
  output.60 <-
    melt(output.60,
         id.vars = strata,
         value.name = "mx_60")
  output.70 <-
    melt(output.70,
         id.vars = strata,
         value.name = "mx_70")
  output.80 <-
    melt(output.80,
         id.vars = strata,
         value.name = "mx_80")
  output.90 <-
    melt(output.90,
         id.vars = strata,
         value.name = "mx_90")

  strata <- c("rn", "sex", "type", "variable")
  output[output.1, mx_1 := i.mx_1, on = strata]
  output[output.99, mx_99 := i.mx_99, on = strata]
  output[output.10, mx_10 := i.mx_10, on = strata]
  output[output.20, mx_20 := i.mx_20, on = strata]
  output[output.30, mx_30 := i.mx_30, on = strata]
  output[output.40, mx_40 := i.mx_40, on = strata]
  output[output.60, mx_60 := i.mx_60, on = strata]
  output[output.70, mx_70 := i.mx_70, on = strata]
  output[output.80, mx_80 := i.mx_80, on = strata]
  output[output.90, mx_90 := i.mx_90, on = strata]

  test <- copy(xx$rate)
  original <- vector("list", length(test))

  for (ii in 1:(length(test))) {
    original[[ii]] <-
      as.data.table(test[[ii]], keep.rownames = T)[, `:=`(type = names(test[ii]))]
    original[[ii]][, c("sex", "type") := tstrsplit(type, "_", fixed = TRUE, keep =
                                                             2:3)]
  }
  original <- rbindlist(original)
  original <-
    melt(
      original,
      id.vars = c("rn", "sex", "type"),
      value.name = "mx"
    )
  original[, paste0("mx_", c(1, 99, 10, 20, 30, 40, 60, 70, 80, 90)) := mx]

  tt <- rbind(original, output, use.names = TRUE, fill = TRUE)
  tt[, age := as.integer(rn)]
  tt[, year := as.integer(as.character(variable))]

  lifetable_all[tt[type == "chd", ], `:=` (
    mx_chd    = mx,
    mx_chd_1  = mx_1,
    mx_chd_99 = mx_99,
    mx_chd_10 = mx_10,
    mx_chd_20 = mx_20,
    mx_chd_30 = mx_30,
    mx_chd_40 = mx_40,
    mx_chd_60 = mx_60,
    mx_chd_70 = mx_70,
    mx_chd_80 = mx_80,
    mx_chd_90 = mx_90
  ),
  on = c("age", "sex", "year")]

  lifetable_all[tt[type == "stroke", ], `:=` (
    mx_stroke    = mx,
    mx_stroke_1  = mx_1,
    mx_stroke_99 = mx_99,
    mx_stroke_10 = mx_10,
    mx_stroke_20 = mx_20,
    mx_stroke_30 = mx_30,
    mx_stroke_40 = mx_40,
    mx_stroke_60 = mx_60,
    mx_stroke_70 = mx_70,
    mx_stroke_80 = mx_80,
    mx_stroke_90 = mx_90
  ),
  on = c("age", "sex", "year")]

  lifetable_all[, `:=` (
    rn = NULL,
    type = NULL,
    variable = NULL
  )]

  # # include population estimates
  # output <- vector("list", length(xxx))
  # for (ii in 1:length(xxx)) {
  #   output[[ii]] <-
  #     data.table(xxx[[ii]], keep.rownames = T)[, `:=`(type = names(xxx[ii]))]
  #   output[[ii]][, c("sex", "race", "type") := tstrsplit(type, "_", fixed = TRUE, keep =
  #                                                          2:(length(strata) +1L))]
  # }
  #
  # output <- rbindlist(output)
  #
  # setnames(output, "rn", "age")
  # output[, age := as.integer(age)]
  # output <-
  #   melt(
  #     output,
  #     id.vars = c("age", "sex", "race", "type"),
  #     variable.name = "year",
  #     value.name = "population",
  #     variable.factor = F
  #   )
  # output[, `:=` (
  #   year = as.integer(year),
  #   population = round(population)
  # )]
  #
  # lifetable_all[output, `:=` (
  #   population = i.population
  # ),
  # on = c("age", "sex", "race", "year")]

  # calculate qx from mx
  lifetable_all[, sex := ifelse(sex == 1, "male", "female")] # male/female hardcoded to demography:::lt()
  setkey(lifetable_all, age)
  lifetable_all[, `:=` (
    qx_total    = demography:::lt(mx_total   , min(age), 1, sex)$qx,
    qx_total_1  = demography:::lt(mx_total_1 , min(age), 1, sex)$qx,
    qx_total_99 = demography:::lt(mx_total_99, min(age), 1, sex)$qx,
    qx_total_10 = demography:::lt(mx_total_10, min(age), 1, sex)$qx,
    qx_total_20 = demography:::lt(mx_total_20, min(age), 1, sex)$qx,
    qx_total_30 = demography:::lt(mx_total_30, min(age), 1, sex)$qx,
    qx_total_40 = demography:::lt(mx_total_40, min(age), 1, sex)$qx,
    qx_total_60 = demography:::lt(mx_total_60, min(age), 1, sex)$qx,
    qx_total_70 = demography:::lt(mx_total_70, min(age), 1, sex)$qx,
    qx_total_80 = demography:::lt(mx_total_80, min(age), 1, sex)$qx,
    qx_total_90 = demography:::lt(mx_total_90, min(age), 1, sex)$qx,
    qx_chd    = demography:::lt(mx_chd   , min(age), 1, sex)$qx,
    qx_chd_1  = demography:::lt(mx_chd_1 , min(age), 1, sex)$qx,
    qx_chd_99 = demography:::lt(mx_chd_99, min(age), 1, sex)$qx,
    qx_chd_10 = demography:::lt(mx_chd_10, min(age), 1, sex)$qx,
    qx_chd_20 = demography:::lt(mx_chd_20, min(age), 1, sex)$qx,
    qx_chd_30 = demography:::lt(mx_chd_30, min(age), 1, sex)$qx,
    qx_chd_40 = demography:::lt(mx_chd_40, min(age), 1, sex)$qx,
    qx_chd_60 = demography:::lt(mx_chd_60, min(age), 1, sex)$qx,
    qx_chd_70 = demography:::lt(mx_chd_70, min(age), 1, sex)$qx,
    qx_chd_80 = demography:::lt(mx_chd_80, min(age), 1, sex)$qx,
    qx_chd_90 = demography:::lt(mx_chd_90, min(age), 1, sex)$qx,
    qx_stroke    = demography:::lt(mx_stroke   , min(age), 1, sex)$qx,
    qx_stroke_1  = demography:::lt(mx_stroke_1 , min(age), 1, sex)$qx,
    qx_stroke_99 = demography:::lt(mx_stroke_99, min(age), 1, sex)$qx,
    qx_stroke_10 = demography:::lt(mx_stroke_10, min(age), 1, sex)$qx,
    qx_stroke_20 = demography:::lt(mx_stroke_20, min(age), 1, sex)$qx,
    qx_stroke_30 = demography:::lt(mx_stroke_30, min(age), 1, sex)$qx,
    qx_stroke_40 = demography:::lt(mx_stroke_40, min(age), 1, sex)$qx,
    qx_stroke_60 = demography:::lt(mx_stroke_60, min(age), 1, sex)$qx,
    qx_stroke_70 = demography:::lt(mx_stroke_70, min(age), 1, sex)$qx,
    qx_stroke_80 = demography:::lt(mx_stroke_80, min(age), 1, sex)$qx,
    qx_stroke_90 = demography:::lt(mx_stroke_90, min(age), 1, sex)$qx
  ), by = .(year, sex)]

  replace_from_table(lifetable_all, "sex", c("male", "female"), c("men", "women"))

  setcolorder(lifetable_all,
              c("year", "age", "sex", "qx_total",
                paste0("qx_total_", c(1, 99, 10, 90, 20, 80, 30, 70, 40, 60)),
                "qx_chd",
                paste0("qx_chd_", c(1, 99, 10, 90, 20, 80, 30, 70, 40, 60)),
                "qx_stroke",
                paste0("qx_stroke_", c(1, 99, 10, 90, 20, 80, 30, 70, 40, 60)),
                "mx_total",
                paste0("mx_total_", c(1, 99, 10, 90, 20, 80, 30, 70, 40, 60)),
                "mx_chd",
                paste0("mx_chd_", c(1, 99, 10, 90, 20, 80, 30, 70, 40, 60)),
                "mx_stroke",
                paste0("mx_stroke_", c(1, 99, 10, 90, 20, 80, 30, 70, 40, 60))
              ))



  # deaths <- deaths[year %in% 2014:2015 & age < 85, ]
  # deaths[, `:=` (
  # 	year = year - 2014L,
  # 	sex = ifelse(sex == "Male", "men", "women"),
  # 	race = tolower(race)
  # )]
  # xx <- deaths[race == "other"]
  # xx[, race := "hispanic"]
  # deaths <- rbind(deaths, xx)
  #
  # lifetable_all[deaths, on = c("age", "sex", "race", "year"),
  # 							mx_stroke := Mx_stroke]

  setkey(lifetable_all, year, age, sex)

  write_fst(lifetable_all, nam1)

  # validation
dependencies(c("ggplot2", "cowplot"))
theme_set(theme_cowplot())

gg <- ggplot(
  lifetable_all[age < 79],
  aes(
    x = year,
    y = qx_total,
    ymin = qx_total_10,
    ymax = qx_total_90,
    col = factor(age),
    fill = factor(age)
  )
) +
  geom_point(size = 1,
             alpha = 5 / 5,
             show.legend = F) +
  geom_line(size = 1, alpha = 5 / 5) +
  geom_ribbon(alpha = 1 / 5,
              linetype = 0,
              show.legend = F) +
  scale_x_continuous(name = "Year") +
  scale_y_continuous(name = "All-cause mortality") +
  facet_grid(sex ~ .)
cowplot::ggsave2(
  filename = "All_cause_mortality_projections_qx.png",
  gg, height = 9, width = 16, units = "in",
  path = "./Population_statistics/Graphs"
)

gg <- ggplot(
  lifetable_all[age < 79],
  aes(
    x = year,
    y = mx_total,
    ymin = mx_total_10,
    ymax = mx_total_90,
    col = factor(age),
    fill = factor(age)
  )
) +
  geom_point(size = 1,
             alpha = 5 / 5,
             show.legend = F) +
  geom_line(size = 1, alpha = 5 / 5) +
  geom_ribbon(alpha = 1 / 5,
              linetype = 0,
              show.legend = F) +
  scale_x_continuous(name = "Year") +
  scale_y_continuous(name = "All-cause mortality") +
  facet_grid(sex ~ .)
cowplot::ggsave2(
  filename = "All_cause_mortality_projections_mx.png",
  gg, height = 9, width = 16, units = "in",
  path = "./Population_statistics/Graphs"
)

gg <- ggplot(
  lifetable_all[age < 79],
  aes(
    x = year,
    y = qx_chd,
    ymin = qx_chd_10,
    ymax = qx_chd_90,
    col = factor(age),
    fill = factor(age)
  )
) +
  geom_point(size = 1,
             alpha = 5 / 5,
             show.legend = F) +
  geom_line(size = 1, alpha = 5 / 5) +
  geom_ribbon(alpha = 1 / 5,
              linetype = 0,
              show.legend = F) +
  scale_x_continuous(name = "Year") +
  scale_y_continuous(name = "CHD mortality") +
  facet_grid(sex ~ .)
cowplot::ggsave2(
  filename = "CHD_mortality_projections_qx.png",
  gg, height = 9, width = 16, units = "in",
  path = "./Population_statistics/Graphs"
)

gg <- ggplot(
  lifetable_all[age < 79],
  aes(
    x = year,
    y = mx_chd,
    ymin = mx_chd_10,
    ymax = mx_chd_90,
    col = factor(age),
    fill = factor(age)
  )
) +
  geom_point(size = 1,
             alpha = 5 / 5,
             show.legend = F) +
  geom_line(size = 1, alpha = 5 / 5) +
  geom_ribbon(alpha = 1 / 5,
              linetype = 0,
              show.legend = F) +
  scale_x_continuous(name = "Year") +
  scale_y_continuous(name = "CHD mortality") +
  facet_grid(sex ~ .)
cowplot::ggsave2(
  filename = "CHD_mortality_projections_mx.png",
  gg, height = 9, width = 16, units = "in",
  path = "./Population_statistics/Graphs"
)

gg <- ggplot(
  lifetable_all[age < 79],
  aes(
    x = year,
    y = qx_stroke,
    ymin = qx_stroke_10,
    ymax = qx_stroke_90,
    col = factor(age),
    fill = factor(age)
  )
) +
  geom_point(size = 1,
             alpha = 5 / 5,
             show.legend = F) +
  geom_line(size = 1, alpha = 5 / 5) +
  geom_ribbon(alpha = 1 / 5,
              linetype = 0,
              show.legend = F) +
  scale_x_continuous(name = "Year") +
  scale_y_continuous(name = "Stroke mortality") +
  facet_grid(sex ~ .)
cowplot::ggsave2(
  filename = "Stroke_mortality_projections_qx.png",
  gg, height = 9, width = 16, units = "in",
  path = "./Population_statistics/Graphs"
)

gg <- ggplot(
  lifetable_all[age < 79],
  aes(
    x = year,
    y = mx_stroke,
    ymin = mx_stroke_10,
    ymax = mx_stroke_90,
    col = factor(age),
    fill = factor(age)
  )
) +
  geom_point(size = 1,
             alpha = 5 / 5,
             show.legend = F) +
  geom_line(size = 1, alpha = 5 / 5) +
  geom_ribbon(alpha = 1 / 5,
              linetype = 0,
              show.legend = F) +
  scale_x_continuous(name = "Year") +
  scale_y_continuous(name = "Stroke mortality") +
  facet_grid(sex ~ .)
cowplot::ggsave2(
  filename = "Stroke_mortality_projections_mx.png",
  gg, height = 9, width = 16, units = "in",
  path = "./Population_statistics/Graphs"
)


  rm(
    rate,
    pop,
    hor,
    xx,
    x1,
    x2,
    output,
    deaths,
    deaths_chd,
    deaths_stroke,
    ii,
    mortf,
    mort.fit,
    k,
    tt
  )
}

# lifetable_all[, year := year - design$init_year + 2000L]
lifetable_all[, year := year - 2000L]

# Population distribution  ----------------------
# Match the sex and age structure of the initial year
population_projections <-
  na.omit(fread(
    "./Population_statistics/population_projections_2010_2060.csv",
    header = T,
    skip = 0
  ))

population_projections <- melt(population_projections, 1:2,
                               variable.name = "year",
                               value.name = "population")
population_projections[, year := as.integer(as.character(year))]
population_projections[, population := as.numeric(population)]

population_actual <-
  population_projections[year == design$init_year + 2000L &
                           between(age, design$ageL, design$ageH), ]

population_actual[, `:=` (year = NULL,
                          sex = factor(sex, levels = c("men", "women")))]

# Calculate the exact fraction of the mid design$init_year population this sample
# represents for ages between design$ageL and design$ageH
sum_pop <- population_actual[, sum(population)]
pop_fraction <- design$n / sum_pop
population_actual[, pct := round(population * pop_fraction)]

# Oversample population below ageL to approximate future population distribution
population_projections <-
  population_projections[age == design$ageL &
                           year  %in% (design$init_year + 2000L +
                                        seq_len(design$sim_horizon)), ]
population_projections[, age := age - (year - design$init_year - 2000L)]

population_projections[, `:=` (pct = round(population * pop_fraction),
                               year = NULL)]

# as the simulation progress these will become design$ageL year old and will enter the
# simulation. Not all of them will be needed as it depends on cvd lag.
# Because I used the number of design$ageL year old projections I disabled other
# mortality for these chaps until they reach this age

cat(
  paste0("Population fraction for ages ", design$ageL, " to ", design$ageH, " = ", pop_fraction, "\n"),
  file = output_dir("simulation parameters temp.txt"),
  append = TRUE
)

population_actual <- rbind(population_actual, population_projections)

rm(nam1, sum_pop)

# lifetable_comb[between(age, design$ageL, design$ageH) &
#                  type == "nonmodelled" ,
#                plot(year, qx, col = age, main = paste0(unique(type), .BY[[1]], .BY[[2]])),
#                keyby = .(race, sex)]
# lifetable_comb[between(age, design$ageL, design$ageH) &
#                  type == "chd" ,
#                plot(year, qx, col = age, main = paste0(unique(type), .BY[[1]], .BY[[2]])),
#                keyby = .(race, sex)]
# lifetable_comb[between(age, design$ageL, design$ageH) &
#                  type == "stroke" ,
#                plot(year, qx, col = age, main = paste0(unique(type), .BY[[1]], .BY[[2]])),
#                keyby = .(race, sex)]
