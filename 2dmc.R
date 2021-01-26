#cmpfile("./2dmc.R")
## US Sodium Policy model: A decision support tool for primary prevention of NCDs
## Copyright (C) 2016  Chris Kypridemos

## US Sodium Policy model is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program; if not, see <http://www.gnu.org/licenses/>
## or write to the Free Software Foundation, Inc., 51 Franklin Street,
## Fifth Floor, Boston, MA 02110-1301  USA.


cat("Sample RR values for 2d Monte Carlo\n\n")

cvd_lag_mc <- fixed_mc$cvd_lag_l[[mc_iter]]
salt_optim_mc <- fixed_mc$salt_optim_l[[mc_iter]]
death_sbp_rr_mc <- fixed_mc$death_sbp_rr_l[[mc_iter]]



## Stochastic lifetable. No univariate distribution was a good fit to the quantiles
## I will sample a year and a quantile separately for cvd and noncvd mortality.
## So I allow the trajectory to change once during the simulation

tt <- data.table( # all projections start from the median
  year = design$init_year:(design$init_year + design$sim_horizon),
  breakpoint_quant_noncvd = 0.5,
  breakpoint_quant_cvd = 0.5
)
tt[year >= fixed_mc$breakpoint_year_noncvd_l[[mc_iter]],
              breakpoint_quant_noncvd := fixed_mc$breakpoint_quant_noncvd_l[[mc_iter]]]
tt[year >= fixed_mc$breakpoint_year_cvd_l[[mc_iter]],
              breakpoint_quant_cvd := fixed_mc$breakpoint_quant_cvd_l[[mc_iter]]]
tt[, `:=` (
  breakpoint_quant_noncvd = roll_mean_left(breakpoint_quant_noncvd, 7),
  breakpoint_quant_cvd    = roll_mean_left(breakpoint_quant_cvd, 7)
)]
# tt[, plot(year, breakpoint_quant_cvd)]

lifetable_mc <- merge(lifetable_all, tt, by = "year", all.x = T) # need to copy lifetable_all within each cluster to avoid spill over between cores in parallel

# use mx insted of qx. Gives better validation graphs
lifetable_mc[, grep("qx_", names(lifetable_mc), value = T) := NULL]
setnames(lifetable_mc, names(lifetable_mc), gsub("mx_", "qx_", names(lifetable_mc)))

lifetable_mc[breakpoint_quant_noncvd < 0.1, qx_mc :=
                qunif(head(normalise(c(breakpoint_quant_noncvd, 0, 0.1)), -2),
                      qx_total_1, qx_total_10)]
lifetable_mc[between(breakpoint_quant_noncvd, 0.1, 0.2), qx_mc :=
                qunif(head(normalise(c(breakpoint_quant_noncvd, 0.1, 0.2)), -2),
                      qx_total_10, qx_total_20)]
lifetable_mc[between(breakpoint_quant_noncvd, 0.2, 0.3), qx_mc :=
                qunif(head(normalise(c(breakpoint_quant_noncvd, 0.2, 0.3)), -2),
                      qx_total_20, qx_total_30)]
lifetable_mc[between(breakpoint_quant_noncvd, 0.3, 0.4), qx_mc :=
                qunif(head(normalise(c(breakpoint_quant_noncvd, 0.3, 0.4)), -2),
                      qx_total_30, qx_total_40)]
lifetable_mc[between(breakpoint_quant_noncvd, 0.4, 0.5), qx_mc :=
                qunif(head(normalise(c(breakpoint_quant_noncvd, 0.4, 0.5)), -2),
                      qx_total_40, qx_total)]
lifetable_mc[between(breakpoint_quant_noncvd, 0.5, 0.6), qx_mc :=
                qunif(head(normalise(c(breakpoint_quant_noncvd, 0.5, 0.6)), -2),
                      qx_total, qx_total_60)]
lifetable_mc[between(breakpoint_quant_noncvd, 0.6, 0.7), qx_mc :=
                qunif(head(normalise(c(breakpoint_quant_noncvd, 0.6, 0.7)), -2),
                      qx_total_60, qx_total_70)]
lifetable_mc[between(breakpoint_quant_noncvd, 0.7, 0.8), qx_mc :=
                qunif(head(normalise(c(breakpoint_quant_noncvd, 0.7, 0.8)), -2),
                      qx_total_70, qx_total_80)]
lifetable_mc[between(breakpoint_quant_noncvd, 0.8, 0.9), qx_mc :=
                qunif(head(normalise(c(breakpoint_quant_noncvd, 0.8, 0.9)), -2),
                      qx_total_80, qx_total_90)]
lifetable_mc[between(breakpoint_quant_noncvd, 0.9, 1), qx_mc :=
                qunif(head(normalise(c(breakpoint_quant_noncvd, 0.9, 1)), -2),
                      qx_total_90, qx_total_99)]
lifetable_mc[breakpoint_quant_noncvd == 0.5, qx_mc := qx_total]

lifetable_mc[breakpoint_quant_cvd < 0.1, `:=` (
  qx_chd_mc = qunif(head(normalise(c(breakpoint_quant_cvd, 0, 0.1)), -2),
                     qx_chd_1, qx_chd_10),
  qx_stroke_mc = qunif(head(normalise(c(breakpoint_quant_cvd, 0, 0.1)), -2),
                    qx_stroke_1, qx_stroke_10)
  )]
lifetable_mc[between(breakpoint_quant_cvd, 0.1, 0.2),  `:=` (
  qx_chd_mc = qunif(head(normalise(c(breakpoint_quant_cvd, 0.1, 0.2)), -2),
                    qx_chd_10, qx_chd_20),
  qx_stroke_mc = qunif(head(normalise(c(breakpoint_quant_cvd, 0.1, 0.2)), -2),
                       qx_stroke_10, qx_stroke_20)
)]
lifetable_mc[between(breakpoint_quant_cvd, 0.2, 0.3),  `:=` (
  qx_chd_mc = qunif(head(normalise(c(breakpoint_quant_cvd, 0.2, 0.3)), -2),
                    qx_chd_20, qx_chd_30),
  qx_stroke_mc = qunif(head(normalise(c(breakpoint_quant_cvd, 0.2, 0.3)), -2),
                       qx_stroke_20, qx_stroke_30)
)]
lifetable_mc[between(breakpoint_quant_cvd, 0.3, 0.4),  `:=` (
  qx_chd_mc = qunif(head(normalise(c(breakpoint_quant_cvd, 0.3, 0.4)), -2),
                    qx_chd_30, qx_chd_40),
  qx_stroke_mc = qunif(head(normalise(c(breakpoint_quant_cvd, 0.3, 0.4)), -2),
                       qx_stroke_30, qx_stroke_40)
)]
lifetable_mc[between(breakpoint_quant_cvd, 0.4, 0.5),  `:=` (
  qx_chd_mc = qunif(head(normalise(c(breakpoint_quant_cvd, 0.4, 0.5)), -2),
                    qx_chd_40, qx_chd),
  qx_stroke_mc = qunif(head(normalise(c(breakpoint_quant_cvd, 0.4, 0.5)), -2),
                       qx_stroke_40, qx_stroke)
)]
lifetable_mc[between(breakpoint_quant_cvd, 0.5, 0.6),  `:=` (
  qx_chd_mc = qunif(head(normalise(c(breakpoint_quant_cvd, 0.5, 0.6)), -2),
                    qx_chd, qx_chd_60),
  qx_stroke_mc = qunif(head(normalise(c(breakpoint_quant_cvd, 0.5, 0.6)), -2),
                       qx_stroke, qx_stroke_60)
)]
lifetable_mc[between(breakpoint_quant_cvd, 0.6, 0.7),  `:=` (
  qx_chd_mc = qunif(head(normalise(c(breakpoint_quant_cvd, 0.6, 0.7)), -2),
                    qx_chd_60, qx_chd_70),
  qx_stroke_mc = qunif(head(normalise(c(breakpoint_quant_cvd, 0.6, 0.7)), -2),
                       qx_stroke_60, qx_stroke_70)
)]
lifetable_mc[between(breakpoint_quant_cvd, 0.7, 0.8),  `:=` (
  qx_chd_mc = qunif(head(normalise(c(breakpoint_quant_cvd, 0.7, 0.8)), -2),
                    qx_chd_70, qx_chd_80),
  qx_stroke_mc = qunif(head(normalise(c(breakpoint_quant_cvd, 0.7, 0.8)), -2),
                       qx_stroke_70, qx_stroke_80)
)]
lifetable_mc[between(breakpoint_quant_cvd, 0.8, 0.9),  `:=` (
  qx_chd_mc = qunif(head(normalise(c(breakpoint_quant_cvd, 0.8, 0.9)), -2),
                    qx_chd_80, qx_chd_90),
  qx_stroke_mc = qunif(head(normalise(c(breakpoint_quant_cvd, 0.8, 0.9)), -2),
                       qx_stroke_80, qx_stroke_90)
)]
lifetable_mc[between(breakpoint_quant_cvd, 0.9, 1),  `:=` (
  qx_chd_mc = qunif(head(normalise(c(breakpoint_quant_cvd, 0.9, 1)), -2),
                    qx_chd_90, qx_chd_99),
  qx_stroke_mc = qunif(head(normalise(c(breakpoint_quant_cvd, 0.9, 1)), -2),
                       qx_stroke_90, qx_stroke_99)
)]
lifetable_mc[breakpoint_quant_cvd == 0.5, `:=` (
  qx_chd_mc = qx_chd,
  qx_stroke_mc = qx_stroke
  )]

lifetable_mc <-
  lifetable_mc[, .(year, age, sex, qx_mc, qx_chd_mc, qx_stroke_mc,
                   qx_noncvd_mc = bound(abs(qx_mc - qx_chd_mc - qx_stroke_mc)))]

# calculate proportional change of mortalty from init_year
lifetable_mc[lifetable_mc[year == design$init_year, ], `:=` (
chd_prop_change = qx_chd_mc / i.qx_chd_mc,
stroke_prop_change = qx_stroke_mc / i.qx_stroke_mc
), on = c("age", "sex")]

lifetable_mc[, sex := factor(sex)]

# lifetable_mc[age == 30 & sex == "men" & race == "white", plot(year, qx_mc)]
# lifetable_mc[age == 60 & sex == "men" & race == "white", plot(year, breakpoint_quant_noncvd)]
# lifetable_mc[age == 30 & sex == "men" & race == "white", plot(year, qx_chd_mc)]
# lifetable_mc[age == 30 & sex == "men" & race == "white", plot(year, chd_prop_change)]
# lifetable_mc[age == 30 & sex == "men" & race == "white", plot(year, qx_stroke_mc)]
# lifetable_mc[age == 30 & sex == "men" & race == "white", plot(year, stroke_prop_change)]
# lifetable_mc[age == 60 & sex == "men" & race == "white", plot(year, breakpoint_quant_cvd)]

if ("chd" %in% design$diseases) {
  sbp_rr_chd_mc     <- fixed_mc$chd_sbp_rr_l[.id == mc_iter, ]
  # quick fix to split decrease in incidence. RR comes from studies that outcome was fatal & non fatal CHD
  # I assume that half of the dcrease in observational studies is because of decrease in incidence
  # (and half because of decrease in case fatality which we currently ignore)
  sbp_rr_chd_mc[, rr := 1 + (rr - 1) * 0.5]
  chd_sbp_tmred_mc  <- fixed_mc$chd_sbp_tmred_l[[mc_iter]]
  chd_incidence_ratio_mc <- fixed_mc$chd_incidence_ratio_l[[mc_iter]]
  chd_epi_mc <- vector("list", 0)
  chd_epi_mc$incidence <- # chd_epi_l$incidence[, .(age, sex, incidence)]
    chd_epi_l$incidence[,
                        .(age, sex,
                          incidence = ifelse(
                            between(age, design$ageL, design$ageH),
                            qbeta(fixed_mc$chd_burden_l[[mc_iter]], shape1, shape2),
                            incidence
                          ))]
  chd_epi_mc$prevalence <- # chd_epi_l$prevalence[, .(age, sex, prevalence)]
    chd_epi_l$prevalence[,
                        .(age, sex,
                          prevalence = ifelse(
                            between(age, design$ageL, design$ageH),
                            qbeta(fixed_mc$chd_burden_l[[mc_iter]], shape1, shape2),
                            prevalence
                          ))]
}

if ("stroke" %in% design$diseases) {
  sbp_rr_stroke_mc     <- fixed_mc$stroke_sbp_rr_l[.id == mc_iter, ]
  # quick fix to split decrease in incidence. RR comes from studies that outcome was fatal & non fatal stroke
  # I assume that half of the dcrease in observational studies is because of decrease in incidence
  # (and half because of decrease in case fatality which we currently ignore)
  sbp_rr_stroke_mc[, rr := 1 + (rr - 1) * 0.5]
  stroke_sbp_tmred_mc  <- fixed_mc$stroke_sbp_tmred_l[[mc_iter]]
  stroke_incidence_ratio_mc <- fixed_mc$stroke_incidence_ratio_l[[mc_iter]]
  stroke_epi_mc <- vector("list", 0)
  stroke_epi_mc$incidence <- # stroke_epi_l$incidence[, .(age, sex, incidence)]
    stroke_epi_l$incidence[,
                        .(age, sex,
                          incidence = ifelse(
                            between(age, design$ageL, design$ageH),
                            qbeta(fixed_mc$stroke_burden_l[[mc_iter]], shape1, shape2),
                            incidence
                          ))]
  stroke_epi_mc$prevalence <- # stroke_epi_l$prevalence[, .(age, sex, prevalence)]
    stroke_epi_l$prevalence[,
                         .(age, sex,
                           prevalence = ifelse(
                             between(age, design$ageL, design$ageH),
                             qbeta(fixed_mc$stroke_burden_l[[mc_iter]], shape1, shape2),
                             prevalence
                           ))]
}


if (exists("tt")) rm(tt)


