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

# User input --------------------------------------------------------------
clusternumber <- 10L # Number of cores to be used ONLY VALID FOR LINUX
ncpu <- 10L
diagnostics <- FALSE

# Preample ----------------------------------------------------------------
if (!require(BrazilSaltModelmisc)) {
  if (!require(remotes)) install.packages("remotes")
  remotes::install_local("./Rpackage/Brazil_salt_model_pkg/")
  library(BrazilSaltModelmisc)
}

output_dir <- function(x = character(0)) (paste0("./Output/", x))
synthpop_dir <- "./SynthPop"
primer_dir <- "./SynthPop/primer"

dependencies(c("data.table", "gamlss", "gamlss.tr", "ggplot2", "haven"), TRUE, FALSE, FALSE, FALSE)



# Define function to split agegroups and create groups
group_pop <- function(x) {
  setorder(x, sex, agegroup)
  x[, group := rleid(sex, agegroup)]
}

agegroup_fn <-
function(x, lagtime = 0) {
    breaks <- c(0, seq(5, 85, 5), Inf)
    labels <- c(
      "00-04",
      "05-09",
      "10-14",
      "15-19",
      "20-24",
      "25-29",
      "30-34",
      "35-39",
      "40-44",
      "45-49",
      "50-54",
      "55-59",
      "60-64",
      "65-69",
      "70-74",
      "75-79",
      "80-84",
      "85+"
    )
    if (is.numeric(x)) {
      agegroup = cut(
        x + lagtime,
        breaks = breaks,
        labels = labels,
        include.lowest = T,
        right = F,
        ordered_result = T
      )
      return(invisible(agegroup))
    } else {
      if (is.data.table(x)) {
        x[, agegroup := cut(
          age + lagtime,
          breaks = breaks,
          labels = labels,
          include.lowest = T,
          right = F,
          ordered_result = T
        )]
        group_pop(x)
        return(invisible(x))
      } else
        return(print("only datatables and vectors are eligible inputs"))
    }
  }

# National Health Survey (PNS 2013)
# Dear Chris,
# I hope you are doing well.
# I have sent you a Dropbox link to the National Health Survey of 2013 micradata with the adjustments that were needed for the model, according do what we have talked before. Please let me know if it works.
# Following in some information about the most important variables with which we will be working.
#
# Variables in the National Health Survey microdata (PNS2013)
# C006 = sex (1=MALE , 2= FEMALE)
# C008 = age in years
# C009 = Race
# 1 = white
# 2 = black
# 3 = asian
# 4 = brown
# 5= indigenous
# 9 = not declared
#
# Q002 = self-declared medical diagnosis of hypertension
# Q006 = currently taking medication for hypertension
# Q06301 = self-declared medical diagnosis of CVDs - myocardial infarction
# Q06302 = self-declared medical diagnosis of CVDs - angina pectoris
# Q06303 = self-declared medical diagnosis of CVDs - heart failure
# Q06304 = self-declared medical diagnosis of CVDs - others
# Q068 = self-declared medical diagnosis of stroke
# Q121 = self-declared medical diagnosis of gastric cancer
# W00407 = final systolic blood pressure
# W00408 = final diastolic blood pressure
# (W00401 & 402, 403 & 404, 405 & 406 correspond to the 1st, 2nd and 3rd systolic and diastolic measurements)
#
# EDUCATION
# VDD004 = highest level of education (years)
# 1 = 0-3 years
# 2 = 4-7 years
# 3 = 8 years (complete primary school)
# 4 = 9-10- years
# 5 = 11 years (complete high school)
# 6 = 12-14 years
# 7 = >=15 years (complete university course)
#
# SAMPLE WEIGHTS
# V00282= sample weight for each person, corrected by non-response and adjusted to population projections
# (the other sample weights are used especially when working with both household and individual questionnaires)
# V00283 = post strata 1 domain
# V00293 = post strata 2 domain


# use  "C:\Users\Rafael\OneDrive\PNS 2013\2017.03.23\Dados\DOMPNS2013.dta", clear
# keep if V0015==1
# recast double V0001 V0024 UPA_PNS V0006_PNS
# *gen double cod_dom = (((( v0001 * 100000000) + v0024)*10000000)+ upa_pns *100)+v0006_pns
# gen double cod_dom_alt = (((( V0024)*10000000)+ UPA_PNS) *100)+V0006_PNS
# saveold "C:\data\DOMPNS2013_area.dta", replace version(12)
#
# use "C:\Users\Rafael\OneDrive\PNS 2013\2017.03.23\Dados\PESPNS2013.dta", clear
# *svyset UPA_PNS [pweight=V0029], strata(V0024)
# svyset UPA_PNS [pweight=V0029], strata(V0024) poststrata(V00293) postweight(V00292) singleunit(centered)
# recast double V0001 V0024 UPA_PNS V0006_PNS
# *gen double cod_dom = (((( V0001 * 100000000) + V0024)*10000000)+ UPA_PNS *100)+V0006_PNS
# gen double cod_dom_alt = (((( V0024)*10000000)+ UPA_PNS) *100)+V0006_PNS
# joinby cod_dom using "C:\Users\Rafael\OneDrive\PNS 2013\2017.03.23\Dados\DOMPNS2013_jb.dta", unmatched(both) _merge(_merge)
# tab _merge
# drop if _merge != 3
# drop if C008<18
# drop if V0025!=1
#
# LINES 4-5 ARE THE INDIVIDUAL CODES GENERATED BY JOINING BY USING V0001, V0024, UPA_PNS, V0006_PNS
#
# LINES 9-10 ARE THE SETTINGS FOR THE SURVEY FUNCTION IN STATA: V0029 U=IS THE SAMPLE WEIGHT, STRATA IS V0024, POST STRATA IS V00293 AND POSTWEIGHT IS V00292

# Saturated fats = SATURADA
# Cholesterol = COLESTER
# MONO_FIN stands for monounsaturated fats
# POLI_FIN for polyunsaturated fats


# sbp ----
tt <- setDT(haven::read_dta(paste0(primer_dir, "/PESPNS2013.dta")))
# tt <- setDT(foreign::read.dta(paste0(primer_dir, "/PESPNS2013.dta")))
tt <- tt[, .(C006, C008, C009, Q002, Q006, Q06301, Q06302, Q06303, Q06304, Q068,
             Q121, W00407, W00401, W00403, W00405, VDD004, V00282)]
setnames(tt, c("sex", "age", "race", "htn_diag", "bp_med", "mi_diag", "ap_diag",
               "hf_diag", "other_cvd_diag", "stroke_diag", "gc_diag",
               "sbp", "sbp1", "sbp2", "sbp3", "educ", "wt"))
replace_from_table(tt, "sex", 1:2, c("men", "women"))
replace_from_table(tt, "race", c(1:5, 9), c("other", "black", "other", "other",
                                            "other", NA), "black_race")
replace_from_table(tt, "race", c(1:5, 9), c("white", "other", "other", "brown",
                                           "other", NA))
tt[, sbp := (sbp + sbp2 + sbp3) / 3] # discard 1st sbp measurement and average 2, 3, and 4
tt[, paste0("sbp", 1:3) := NULL]

dt <- na.omit(tt[between(sbp, 70, 240) & between(age, 20, 84), .(
  sbp, age, sex, wt)]
)
dt[, wt := .N * wt/sum(wt)]
agegroup_fn(dt)
dt[, sex := factor(sex)]

dt[, mean(sbp), age][, plot(age, V1)]
dt[between(sbp, 70, 200), plot(density(sbp, weights = wt))]
dt[between(sbp, 70, 200), plot(density(sbp^-3, weights = wt))]
#dt[, sbp := sbp^-5]
dt[, age := scale(age, 43.33728, 16.54644)]

set.seed(44)
lns <- sample(nrow(dt), nrow(dt) * 0.8)
dt_trn   <- dt[lns] # train dataset
dt_crv   <- dt[!lns]  # cross-validation dataset
dt_small <- dt_trn[sample(.N, 2e4)] # draw a subsample for faster execution

gen.trun(
  par = c(dt[, min(sbp)], dt[, max(sbp)]),
  family = "BCPEo",
  type = "both",
  name = "_sbp"
)
marg_distr <- fitDistPred(dt_trn$sbp, type = "realplus", weights = dt_trn$wt,
                          extra = "BCPEo_sbp",
                          try.gamlss = FALSE, trace = TRUE, newdata = dt_crv$sbp)


marg_distr$fits
params <- vector("list")
for (i in seq_along(marg_distr$parameters)) {
  nam <- marg_distr$parameters[[i]]
  params[[i]] <- get(nam, marg_distr)
  names(params)[i] <- nam
}
params$p <- seq(0.00001, 0.99999, 0.00001)
distr_nam <- marg_distr$family[[1]]
y <- do.call(paste0("q", distr_nam), params)
y_wt <- rep(1/length(y), length(y))
# validation plots
# see http://www.csss.washington.edu/files/working-papers/2002/wp27.pdf
reldist_diagnostics(dt[, sbp], y, dt[, wt/sum(wt)], y_wt,
                    main = expression(bold(SBP~(mmHG))))

setMKLthreads(20L)
sbp_model <- gamlss(
  sbp ~ pvc(age, by = sex), sigma.formula = ~pvc(age, by = sex),
  nu.formula = ~pvc(age, by = sex), tau.formula = ~pvc(age, by = sex),
  family = BCPEo_sbp,
  weights = dt$wt,
  data = dt,
  method = mixed(20, 100)
)

sbp_model$data <- dt_trn

saveRDS(sbp_model, "./Lifecourse_models/sbp_model.rds")
print("Model saved.")

# setMKLthreads(1L)
# t1 <- chooseDist(sbp_model, type = "realplus", parallel = "multicore", ncpus = 15)
# getOrder(t1,3)
# fm <- update(sbp_model, family = "WEI3")


if (diagnostics) {
  sbp_model <- readRDS("./Lifecourse_models/sbp_model.rds")

  wp(sbp_model)
  wp(sbp_model, xvar = age)
  wp(sbp_model, xvar = ~sex)
  plot(sbp_model)

  zz <- validate_gamlss(dt, sbp_model, 10, dt)
  zz[, weight := wt/sum(wt), by = type]
  reldist_diagnostics(zz[type == "Observed", sbp],
                      zz[type == "Modelled", sbp],
                      zz[type == "Observed", weight],
                      zz[type == "Modelled", weight],
                      main = expression(bold(sbp)),
                      discrete = FALSE)
  dependencies("ggplot2")

  zz[, weight := wt/sum(wt), by = .(type, sex)]
  ggplot(zz, aes(sbp, colour = type, weight = wt, linetype = type)) +
    stat_ecdf() +
    facet_wrap(.~sex, nrow = 3) + ggtitle("Sex")

  zz[, weight := wt/sum(wt), by = .(type, sex, agegroup)]
  ggplot(zz, aes(sbp, colour = type, weight = wt)) +
    stat_ecdf() +
    facet_grid(sex~agegroup) + ggtitle("Sex ~ Age group")
}

# create table with distribution parameters
newdata <- CJ(age_int = 20:84, sex = factor(c("men", "women")))
newdata[, age := scale(age_int, 43.33728, 16.54644)]
newdata[, c("mu", "sigma", "nu", "tau") := predictAll(sbp_model, .SD[, .(age, sex)], data = dt_trn)]
newdata[, age := age_int]
newdata[, age_int := NULL]
fwrite(newdata, "./Lifecourse_models/sbp_table.csv")

# Black race prediction ----
# we need this from age 0 to acount for future populations
dt <- na.omit(tt[between(age, 0, 84), .(
  black_race, age, sex, wt)]
)
dt[, wt := .N * wt/sum(wt)]
agegroup_fn(dt)
dt[, sex := factor(sex)]
dt[, black_race := factor(black_race, c("other", "black"))]
dt[, age := scale(age, 43.33728, 16.54644)]

set.seed(45)
lns <- sample(nrow(dt), nrow(dt) * 0.8)
dt_trn   <- dt[lns] # train dataset
dt_crv   <- dt[!lns]  # cross-validation dataset
dt_small <- dt_trn[sample(.N, 2e4)] # draw a subsample for faster execution


marg_distr <- fitDistPred(dt_trn$black_race, type = "binom", weights = dt_trn$wt,
                          try.gamlss = FALSE, trace = TRUE,
                          newdata = dt_crv$black_race, bd = rep(1L, nrow(dt_crv)))


marg_distr$fits
params <- vector("list")
for (i in seq_along(marg_distr$parameters)) {
  nam <- marg_distr$parameters[[i]]
  params[[i]] <- get(nam, marg_distr)
  names(params)[i] <- nam
}
params$p <- seq(0.00001, 0.99999, 0.00001)
distr_nam <- marg_distr$family[[1]]
y <- do.call(paste0("q", distr_nam), params)
y_wt <- rep(1/length(y), length(y))
# validation plots
# see http://www.csss.washington.edu/files/working-papers/2002/wp27.pdf
reldist_diagnostics(dt[, sbp], y, dt[, wt/sum(wt)], y_wt,
                    main = expression(bold(SBP~(mmHG))))

setMKLthreads(20L)
black_race_model <- gamlss(
  black_race ~ pvc(age, by = sex),
  family = BI,
  weights = dt_trn$wt,
  data = dt_trn,
  method = mixed(20, 100)
)

newdata <- CJ(age_int = 20:84, sex = factor(c("men", "women")))
newdata[, age := scale(age_int, 43.33728, 16.54644)]
newdata[, c("mu") := predictAll(black_race_model, .SD[, .(age, sex)], data = dt_trn)]
newdata[, age := age_int]
newdata[, age_int := NULL]
saveRDS(black_race, "./Lifecourse_models/black_race_model.rds")
fwrite(newdata, "./Lifecourse_models/black_race_table.csv")


# Salt exposure ----
tt <- setDT(haven::read_dta(paste0(primer_dir, "/Food recordatory Sodium exposure 2008.dta")))
replace_from_table(tt, "sex", 1:2, c("men", "women"))
setnames(tt, "age_years", "age")
setnames(tt, "repwt_0", "wt")
setnames(tt, names(tt), tolower(names(tt)))
tt[, `:=` (salt_final        = sodium_final * 2.5 / 1e3,
           salt_othersources = sodium_othersources * 2.5 / 1e3,
           salt_added        = sodium_added * 2.5 / 1e3)] # convert sodium to salt

dt <- na.omit(tt[between(age, 20, 84) & salt_final < 30, .(
  salt_final, salt_othersources, salt_added, age, sex, wt)]
)

dt[, wt := .N * wt/sum(wt)]
agegroup_fn(dt)
dt[, sex := factor(sex)]
dt[, mean(salt_added), age][, plot(age, V1)]
dt[, plot(density(salt_added, weights = wt))]
dt[, age := scale(age, 41.90532, 15.52681)]

cov.wt(dt[, .(salt_added, salt_othersources)], dt$wt, TRUE)$cor # cor ~ -0.132
plot(dt$salt_othersources, dt$salt_added)
cov.wt(dt[, .(perc_rank(salt_added, .N), perc_rank(salt_othersources, .N))], dt$wt, TRUE)$cor # cor ~ -0.132

cov.wt(dt[, .(salt_final, salt_added/salt_final)], dt$wt, TRUE)$cor # cor ~ -0.063
cov.wt(dt[, .(salt_final, salt_added/salt_othersources)], dt$wt, TRUE)$cor # cor ~ 0.036
plot(dt$salt_added/dt$salt_final, dt$salt_final)
dt[, plot(density(salt_added/salt_final, weights = wt))]
dt[, salt_added_to_final := salt_added/salt_final]
# Correlation between sources of sodium is small. I will use copula-like dirty approach.
# I will model them as independent.


set.seed(49)
lns <- sample(nrow(dt), nrow(dt) * 0.8)
dt_trn   <- dt[lns] # train dataset
dt_crv   <- dt[!lns]  # cross-validation dataset

gen.trun(
  par = c(dt[, min(salt_final)], dt[, max(salt_final)]),
  family = "BCPEo",
  type = "both",
  name = "_salt"
)
marg_distr <- fitDistPred(dt_trn$salt_final, type = "realplus", weights = dt_trn$wt,
                          extra = "BCPEo_salt",
                          try.gamlss = FALSE, trace = TRUE, newdata = dt_crv$salt_final)


marg_distr$fits
params <- vector("list")
for (i in seq_along(marg_distr$parameters)) {
  nam <- marg_distr$parameters[[i]]
  params[[i]] <- get(nam, marg_distr)
  names(params)[i] <- nam
}
params$p <- seq(0.00001, 0.99999, 0.00001)
distr_nam <- marg_distr$family[[1]]
y <- do.call(paste0("q", distr_nam), params)
y_wt <- rep(1/length(y), length(y))
# validation plots
# see http://www.csss.washington.edu/files/working-papers/2002/wp27.pdf
reldist_diagnostics(dt[, salt_final], y, dt[, wt/sum(wt)], y_wt,
                    main = expression(bold(Salt~(g/d))))

setMKLthreads(20L)
salt_model <- gamlss(
  # salt_final ~ pvc(age, by = sex), sigma.formula = ~pvc(age, by = sex),
  # nu.formula = ~pvc(age, by = sex), tau.formula = ~pvc(age, by = sex),
  salt_final ~ pb(age) * sex, sigma.formula = ~pb(age) * sex,
  nu.formula = ~pb(age) * sex, tau.formula = ~pb(age) * sex,
  family = BCPEo_salt,
  mu.start = 11.2, sigma.start = 0.41, nu.start = 0.311, tau.start = 1.7,
  weights = dt$wt,
  data = dt,
  method = RS()
)

salt_model$data <- dt

saveRDS(salt_model, "./Lifecourse_models/salt_model.rds")
print("Model saved.")

if (diagnostics) {
  salt_model <- readRDS("./Lifecourse_models/salt_model.rds")

  wp(salt_model)
  wp(salt_model, xvar = age)
  wp(salt_model, xvar = ~sex)
  plot(salt_model)

  zz <- validate_gamlss(dt, salt_model, 10, dt)
  zz[, weight := wt/sum(wt), by = type]
  reldist_diagnostics(zz[type == "Observed", salt_final],
                      zz[type == "Modelled", salt_final],
                      zz[type == "Observed", weight],
                      zz[type == "Modelled", weight],
                      main = expression(bold(salt)),
                      discrete = FALSE)
  dependencies("ggplot2")

  zz[, weight := wt/sum(wt), by = .(type, sex)]
  ggplot(zz, aes(salt_final, colour = type, weight = wt, linetype = type)) +
    stat_ecdf() +
    facet_wrap(.~sex, nrow = 3) + ggtitle("Sex")

  zz[, weight := wt/sum(wt), by = .(type, sex, agegroup)]
  ggplot(zz, aes(salt_final, colour = type, weight = wt)) +
    stat_ecdf() +
    facet_grid(sex~agegroup) + ggtitle("Sex ~ Age group")
}

# create table with distribution parameters
newdata <- CJ(age_int = 20:84, sex = factor(c("men", "women")))
newdata[, age := scale(age_int, 41.90532, 15.52681)]
newdata[, c("mu", "sigma", "nu", "tau") := predictAll(salt_model, .SD[, .(age, sex)], data = salt_model$data)]
newdata[, age := age_int]
newdata[, age_int := NULL]
fwrite(newdata, "./Lifecourse_models/salt_table.csv")


# Salt othersources ----
dt[salt_othersources == 0, salt_othersources := 0.001]
dt_trn[salt_othersources == 0, salt_othersources := 0.001]
dt_crv[salt_othersources == 0, salt_othersources := 0.001]

gen.trun(
  par = c(dt[, min(salt_othersources)], dt[, max(salt_othersources)]),
  family = "BCPE",
  type = "both",
  name = "_salt2"
)
marg_distr <- fitDistPred(dt_trn$salt_othersources, type = "realplus", weights = dt_trn$wt,
                          extra = "BCPE_salt2",
                          try.gamlss = FALSE, trace = TRUE, newdata = dt_crv$salt_othersources)


marg_distr$fits
params <- vector("list")
for (i in seq_along(marg_distr$parameters)) {
  nam <- marg_distr$parameters[[i]]
  params[[i]] <- get(nam, marg_distr)
  names(params)[i] <- nam
}
params$p <- seq(0.00001, 0.99999, 0.00001)
distr_nam <- marg_distr$family[[1]]
y <- do.call(paste0("q", distr_nam), params)
y_wt <- rep(1/length(y), length(y))
# validation plots
# see http://www.csss.washington.edu/files/working-papers/2002/wp27.pdf
reldist_diagnostics(dt[, salt_othersources], y, dt[, wt/sum(wt)], y_wt,
                    main = expression(bold(salt_othersources)))

setMKLthreads(20L)
salt_othersources_model <- gamlss(
  # salt_final ~ pvc(age, by = sex), sigma.formula = ~pvc(age, by = sex),
  # nu.formula = ~pvc(age, by = sex), tau.formula = ~pvc(age, by = sex),
  salt_othersources ~ pb(age) * sex, sigma.formula = ~pb(age) * sex,
  nu.formula = ~pb(age) * sex, tau.formula = ~pb(age) * sex,
  family = BCPE_salt2,
  mu.start = params$mu, sigma.start = params$sigma, nu.start = params$nu, tau.start = params$tau,
  weights = dt$wt,
  data = dt,
  method = RS(100)
)

salt_othersources_model$data <- dt

saveRDS(salt_othersources_model, "./Lifecourse_models/salt_othersources_model.rds")
print("Model saved.")

if (diagnostics) {
  salt_othersources_model <- readRDS("./Lifecourse_models/salt_othersources_model.rds")

  wp(salt_othersources_model)
  wp(salt_othersources_model, xvar = age)
  wp(salt_othersources_model, xvar = ~sex)
  plot(salt_othersources_model)

  zz <- validate_gamlss(dt, salt_othersources_model, 10, dt)
  zz[, weight := wt/sum(wt), by = type]
  reldist_diagnostics(zz[type == "Observed", salt_othersources],
                      zz[type == "Modelled", salt_othersources],
                      zz[type == "Observed", weight],
                      zz[type == "Modelled", weight],
                      main = expression(bold(salt_othersources)),
                      discrete = FALSE)
  dependencies("ggplot2")

  zz[, weight := wt/sum(wt), by = .(type, sex)]
  ggplot(zz, aes(salt_othersources, colour = type, weight = wt, linetype = type)) +
    stat_ecdf() +
    facet_wrap(.~sex, nrow = 3) + ggtitle("Sex")

  zz[, weight := wt/sum(wt), by = .(type, sex, agegroup)]
  ggplot(zz, aes(salt_othersources, colour = type, weight = wt)) +
    stat_ecdf() +
    facet_grid(sex~agegroup) + ggtitle("Sex ~ Age group")
}

# create table with distribution parameters
newdata <- CJ(age_int = 20:84, sex = factor(c("men", "women")))
newdata[, age := scale(age_int, 41.90532, 15.52681)]
newdata[, c("mu", "sigma", "nu", "tau") :=
          predictAll(salt_othersources_model, .SD[, .(age, sex)],
                     data = salt_othersources_model$data)]
newdata[, age := age_int]
newdata[, age_int := NULL]
fwrite(newdata, "./Lifecourse_models/salt_othersources_table.csv")

# Salt added ----
dt[salt_added == 0, salt_added := 0.001]
dt_trn[salt_added == 0, salt_added := 0.001]
dt_crv[salt_added == 0, salt_added := 0.001]

gen.trun(
  par = c(dt[, min(salt_added)], dt[, max(salt_added)]),
  family = "BCT",
  type = "both",
  name = "_salt2"
)
marg_distr <- fitDistPred(dt_trn$salt_added, type = "realplus", weights = dt_trn$wt,
                          extra = "BCT_salt2",
                          try.gamlss = FALSE, trace = TRUE, newdata = dt_crv$salt_added)


marg_distr$fits
params <- vector("list")
for (i in seq_along(marg_distr$parameters)) {
  nam <- marg_distr$parameters[[i]]
  params[[i]] <- get(nam, marg_distr)
  names(params)[i] <- nam
}
params$p <- seq(0.00001, 0.99999, 0.00001)
distr_nam <- marg_distr$family[[1]]
y <- do.call(paste0("q", distr_nam), params)
y_wt <- rep(1/length(y), length(y))
# validation plots
# see http://www.csss.washington.edu/files/working-papers/2002/wp27.pdf
reldist_diagnostics(dt[, salt_added], y, dt[, wt/sum(wt)], y_wt,
                    main = expression(bold(salt_added)))

setMKLthreads(20L)
salt_added_model <- gamlss(
  # salt_final ~ pvc(age, by = sex), sigma.formula = ~pvc(age, by = sex),
  # nu.formula = ~pvc(age, by = sex), tau.formula = ~pvc(age, by = sex),
  salt_added ~ pb(age) * sex, sigma.formula = ~pb(age) * sex,
  nu.formula = ~pb(age) * sex, tau.formula = ~pb(age) * sex,
  family = BCT_salt2,
  mu.start = params$mu, sigma.start = params$sigma, nu.start = params$nu, tau.start = params$tau,
  weights = dt$wt,
  data = dt,
  method = RS()
)

salt_added_model$data <- dt

saveRDS(salt_added_model, "./Lifecourse_models/salt_added_model.rds")
print("Model saved.")

if (diagnostics) {
  salt_added_model <- readRDS("./Lifecourse_models/salt_added_model.rds")

  wp(salt_added_model)
  wp(salt_added_model, xvar = age)
  wp(salt_added_model, xvar = ~sex)
  plot(salt_added_model)

  zz <- validate_gamlss(dt, salt_added_model, 10, dt)
  zz[, weight := wt/sum(wt), by = type]
  reldist_diagnostics(zz[type == "Observed", salt_added],
                      zz[type == "Modelled", salt_added],
                      zz[type == "Observed", weight],
                      zz[type == "Modelled", weight],
                      main = expression(bold(salt_added)),
                      discrete = FALSE)
  dependencies("ggplot2")

  zz[, weight := wt/sum(wt), by = .(type, sex)]
  ggplot(zz, aes(salt_added, colour = type, weight = wt, linetype = type)) +
    stat_ecdf() +
    facet_wrap(.~sex, nrow = 3) + ggtitle("Sex")

  zz[, weight := wt/sum(wt), by = .(type, sex, agegroup)]
  ggplot(zz, aes(salt_added, colour = type, weight = wt)) +
    stat_ecdf() +
    facet_grid(sex~agegroup) + ggtitle("Sex ~ Age group")
}

# create table with distribution parameters
newdata <- CJ(age_int = 20:84, sex = factor(c("men", "women")))
newdata[, age := scale(age_int, 41.90532, 15.52681)]
newdata[, c("mu", "sigma", "nu", "tau") :=
          predictAll(salt_added_model, .SD[, .(age, sex)],
                     data = salt_added_model$data)]
newdata[, age := age_int]
newdata[, age_int := NULL]
fwrite(newdata, "./Lifecourse_models/salt_added_table.csv")

# Correlated uniforns ----
# from https://stats.stackexchange.com/questions/31771/generate-three-correlated-uniformly-distributed-random-variables
# and https://stats.stackexchange.com/questions/66610/generate-pairs-of-random-numbers-uniformly-distributed-and-correlated/66617#66617

u1 = runif(300)
u2 = runif(300)
# z = ifelse(rbinom(300,1,.132),u1,u2) # for positive correlations
z = ifelse(rbinom(300,1,.70),(1-u1),u2) # for negative correlations

cor(cbind(u1,z))
plot(u1, z)
# or

## Initialization and parameters
set.seed(123)
r <- -0.132                            # Target (Spearman) correlation
n <- 5e6                            # Number of samples

## Functions
gen.gauss.cop <- function(r, n){
  rho <- 2 * sin(r * pi/6)        # Pearson correlation
  P <- toeplitz(c(1, rho))        # Correlation matrix
  d <- nrow(P)                    # Dimension
  ## Generate sample
  U <- pnorm(matrix(rnorm(n*d), ncol = d) %*% chol(P))
  return(U)
}

## Data generation and visualization
U <- gen.gauss.cop(r = r, n = n)
pairs(U, diag.panel = function(x){
  h <- hist(x, plot = FALSE)
  rect(head(h$breaks, -1), 0, tail(h$breaks, -1), h$counts/max(h$counts))})
