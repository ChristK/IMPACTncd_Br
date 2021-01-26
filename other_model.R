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

cat("Estimating deaths from other causes...\n")
ptm <- proc.time()
set(POP, NULL, "other_mrtl_sc00", 0L)
set(POP, NULL, "qx", 0)

#cat("Inflate mortality for hypertensives...\n\n")
# RR (= 1.3) roughly based on stringhini_socioeconomic_2017 figure 4
set(POP, NULL, "death_sbp_rr", 1)
POP[xps_sbp_cvdlag > 140, death_sbp_rr := death_sbp_rr_mc] # TODO more elegant/accurate way

# If mortality forecasts do not include uncertainty then PAF should be per year.
# Otherwise double counts expected reductions in SBP already included in forecasts.
# For this model I include forecast uncertainty so above doesn't apply
deathpaf <-
  POP[year == design$init_year,
      .(paf = 1 - 1 / (sum(death_sbp_rr) / .N)),
      keyby = .(age, sex)
      ]
deathpaf[, paf := predict(loess(paf~age, span = 0.50)),
         keyby = .(sex)]


lifetable_mc[deathpaf, `:=` (qx_noncvd_mc_adj = qx_noncvd_mc * (1 - paf),
                             qx_mc_adj = qx_mc * (1 - paf)),
             on = c("age", "sex")]
lifetable_mc[is.na(qx_noncvd_mc_adj), qx_noncvd_mc_adj := qx_noncvd_mc]
lifetable_mc[is.na(qx_mc_adj), qx_mc_adj := qx_mc]

# combine noncvd and total mortality depending on age
lifetable_mc[between(age, design$ageL, design$ageH), qx_mc_adj := qx_noncvd_mc_adj]

POP[lifetable_mc,
		qx_mc_adj := i.qx_mc_adj,
    on = c("year", "age", "sex")]

POP[, qx := qx_mc_adj * death_sbp_rr]

POP[year < design$init_year, qx := 0] # already survived to year 0
POP[age < design$ageL, qx := 0] # because birth engine use pop number for design$ageL
POP[age > design$ageH, qx := 1] # Need to be 100 for life expectancy calculations


if (!design$kismet) POP[, `:=`(rn_other_death = dice(.N),
                              rn_resolve_death = dice(.N))]

setkey(POP, pid, year)
POP[, other_mrtl_sc00 := mortality_type_1(year,
                                          new_simulant,
                                          qx,
                                          rn_other_death,
                                          rep(1L, .N),
                                          design$init_year,
                                          1L)]

print(proc.time() - ptm)
