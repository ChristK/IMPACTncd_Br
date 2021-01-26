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


# CHD statistics ----------------------------------------------------------
if ("chd" %in% design$diseases) {
  if (file.exists("./CVD_statistics/chd_epi_l.qs")) {
    chd_epi_l <- qread("./CVD_statistics/chd_epi_l.qs")
  } else {
  xx <-  rbind(
    fread(paste0("./CVD_statistics/CHD_male_", design$init_year + 2000L ,".csv"),
          sep = ",", header = TRUE, stringsAsFactors = FALSE,
          skip = 3, nrows = design$ageH + 1L)[, `:=`(sex = "men")],
    fread(paste0("./CVD_statistics/CHD_female_", design$init_year + 2000L,".csv"),
          sep = ",", header = TRUE, stringsAsFactors = FALSE,
          skip = 3, nrows = design$ageH + 1L)[, `:=`(sex = "women")]
  )[, `:=` (sex = factor(sex), age = as.integer(Age))]

  chd_epi_l <- vector("list", 0)
  chd_epi_l$incidence <- setnames(copy(xx[, c(14, 15, 6), with = FALSE]),
                                "Incidence (rates)", "incidence")[
                                  , incidence := as.numeric(incidence)]

  chd_epi_l$prevalence <- setnames(copy(xx[, c(14, 15, 7), with = FALSE]),
                                 "Prevalence (rates)", "prevalence")[
                                   , prevalence := as.numeric(prevalence)]

  chd_epi_l$fatality <- setnames(copy(xx[, c(14, 15, 9), with = FALSE]),
                               "Case fatality (rates)", "fatality")[
                                 , fatality := as.numeric(fatality)]

  # chd_epi_l$duration <- setnames(copy(xx[, c(14, 15, 10), with = FALSE]),
  #                                "Duration (years)", "duration")[
  #                                  , duration := as.numeric(duration)]

  chd_epi_l$incidence[between(age, design$ageL, design$ageH),
                      c("shape1", "shape2") :=
                        as.list(get.beta.par(p = c(0.025, 0.5, 0.975), q = c(incidence * 0.5, incidence, incidence * 2),
                                     fit.weights = c(1, 2, 1), show.output = FALSE, plot = FALSE)),
                      by = .(age, sex)]
  chd_epi_l$prevalence[between(age, design$ageL, design$ageH),
                      c("shape1", "shape2") :=
                        as.list(get.beta.par(p = c(0.025, 0.5, 0.975), q = c(prevalence * 0.5, prevalence, prevalence * 2),
                                             fit.weights = c(1, 2, 1), show.output = FALSE, plot = FALSE)),
                      by = .(age, sex)]
  xx  <- fread("./CVD_statistics/chd_duration_FROM_US.csv")
  xxx <- data.table(age = 20:design$ageH)[, agegroup := agegroup_fn(age)]
  xx  <- xx[xxx, on = c(agegroup5 = "agegroup")]
  setnames(xx, "agegroup5", "agegroup")
  rm(xxx)

  chd_epi_l$duration <- copy(xx)

  # Estimate 30-day case fatality rate after a MI (corrected for angina cases).
  # Same as used for IMPACTncd for England from Smolina and BHF statistics
  # could not find Brazil data on 30-day case fatality for first MI by age that
  # also include not hospitalised fatal MIs

  rate <- c(0.1, 0.1, 0.1, 0.02, 0.02, 0.03, 0.03, 0.05, 0.05, 0.14, 0.14, 0.25,
            0.08, 0.08, 0.08, 0.01, 0.01, 0.03, 0.03, 0.05, 0.05, 0.16, 0.16, 0.3) # 30-day case fatality rate
  age <- rep(seq(30, 85, 5), 2)# data by age group for men and women
  # plot(age, rate, xlim = c(0, 100), ylim = c(0, 0.30))
  lm_fit <- lm(rate ~ poly(age, 2) )
  # lines(0:100, predict(lm_fit, newdata = data.frame(age = 0:100)), col = "red")
  chd_epi_l$day30_fatality <- data.table(age = (design$ageL:design$ageH),
  day30_fatality = predict(lm_fit, newdata = data.frame(age = design$ageL:design$ageH)))

  qsave(chd_epi_l, "./CVD_statistics/chd_epi_l.qs")
  }

}

# Stroke statistics -------------------------------------------------------
# Do I have to separate between ischaemic and haemorrhagic? The risk factors seems more ore less the same.
if ("stroke" %in% design$diseases) {
  if (file.exists("./CVD_statistics/stroke_epi_l.qs"))
  {
    stroke_epi_l <- qread("./CVD_statistics/stroke_epi_l.qs")
  } else {
  xx <- rbind(
    fread(paste0("./CVD_statistics/stroke_male_" ,design$init_year + 2000L,".csv"),
          sep = ",", header = T, stringsAsFactors = F,
          skip = 3, nrows = design$ageH + 1L)[, `:=`(sex = "men")],
    fread(paste0("./CVD_statistics/stroke_female_" ,design$init_year + 2000L,".csv"),
          sep = ",", header = T, stringsAsFactors = F,
          skip = 3, nrows = design$ageH + 1L)[, `:=`(sex = "women")]
  )[, `:=` (sex = factor(sex), age = as.integer(Age))]

  stroke_epi_l <- vector("list", 0)
  stroke_epi_l$incidence <- setnames(copy(xx[, c(14, 15, 6), with = F]),
                          "Incidence (rates)", "incidence")[
                            , incidence := as.numeric(incidence)]

  stroke_epi_l$prevalence <- setnames(copy(xx[, c(14, 15, 7), with = F]),
                           "Prevalence (rates)", "prevalence")[
                             , prevalence := as.numeric(prevalence)]

  stroke_epi_l$fatality <- setnames(copy(xx[, c(14, 15, 9), with = F]),
                          "Case fatality (rates)", "fatality")[
                            , fatality := as.numeric(fatality)]

  stroke_epi_l$incidence[between(age, design$ageL, design$ageH),
                      c("shape1", "shape2") :=
                        as.list(
                          get.beta.par(p = c(0.025, 0.5, 0.975),
                                             q = c(incidence * 0.5, incidence, incidence * 2),
                                             fit.weights = c(1, 2, 1), show.output = F, plot = F,
                                             tol = 0.005)),
                      by = .(age, sex)]
  stroke_epi_l$prevalence[between(age, design$ageL, design$ageH),
                       c("shape1", "shape2") :=
                         as.list(
                           get.beta.par(p = c(0.025, 0.5, 0.975),
                                              q = c(prevalence * 0.5, prevalence, prevalence * 2),
                                              fit.weights = c(1, 2, 1), show.output = F, plot = F)),
                       by = .(age, sex)]

  stroke_epi_l$duration <- fread("./CVD_statistics/stroke_duration_FROM_US.csv")

  # Estimate 30-day case fatality rate after a stroke.
  # from hollander_incidence_2003 table3
  # could not find USA data on 30-day case fatality for first stroke by age that
  # also include not hospitalised fatal strokes
  rate <- c(0.09, 0.139, 0.314, 0.521) # 30-day case fatality rate
  age <- c(60, 70, 80, 95)# data by age group. here is the
  # plot(age, rate, xlim = c(0, 100), ylim = c(0, 0.55))
  lm_fit <- lm(log(rate) ~ age)
  # exp(predict(lm_fit))
  # lines(age, exp(predict(lm_fit)))
  # lines(0:100, exp(predict(lm_fit, newdata = data.frame(age = 0:100))), col = "red")
  stroke_epi_l$day30_fatality <- data.table(age = (design$ageL:design$ageH),
                 day30_fatality = exp(predict(lm_fit, newdata = data.frame(age = design$ageL:design$ageH))))

  qsave(stroke_epi_l, "./CVD_statistics/stroke_epi_l.qs")
  }
}


if (exists("nam")) rm(nam)
if (exists("xx")) rm(xx)
if (exists("rate")) rm(rate)
if (exists("age")) rm(age)
if (exists("lm_fit")) rm(lm_fit)
