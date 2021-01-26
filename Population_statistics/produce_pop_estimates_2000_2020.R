options(scipen = 999)
if (!require(Kmisc)) {devtools::install_github("kevinushey/Kmisc"); library("Kmisc")}
library("data.table")
library("cowplot")
tt <- fread("~/pCloudDrive/My Datasets/Salt_Brazil/Population/Census/ipumsi_00001.csv")
names(tt)
setnames(tt, tolower(names(tt)))
names(tt)
tt[, lapply(.SD, class)]
tt[, c("country", "sample", "serial", "hhwt", "pernum", "resident", "indig", "school", "lit") := NULL]

calib <- fread("~/pCloudDrive/My Datasets/Salt_Brazil/Population/officiall_pop_projection.csv") # official pop projections
calib <- melt(calib, 1:2, variable.name = "year", value.name = "pop")
calib[, year := as.integer(as.character(year))]
# compress age >80
calib[age > 80, age := 80]
calib <- calib[, .(pop = sum(pop)), keyby = .(year, age, sex)]

qplot(x = year, y = pop, col=age, data = calib, facets = sex~.)
calib[, sum(pop), keyby = year][, qplot(year, V1)]

t1 <- calib[year == 2020 & sex == 1, pop] / calib[year == 2010 & sex == 1, shift(pop, 2)]
t2 <- calib[year == 2020 & sex == 2, pop] / calib[year == 2010 & sex == 2, shift(pop, 2)]
summary(c(t1, t2))

## Recode education categories ----
# 0  = less than 4 years
# 1  = 4 to 11 years (inclussive)
# 2  = 12 or more years
# 99 = missing
#
# ttt <- tt[year > 1990 , table(educbr, yrschool)]
# ttt <- as.matrix(ttt)[]
# View(as.data.frame.matrix(ttt)[])
tt[educbr < 2141, educ := 0L]
tt[between(educbr, 2141, 3910), educ := 1L]
tt[between(educbr, 4110, 4900), educ := 2L]
tt[age < 6, educ := 0L]
tt[is.na(educ), .N, keyby = year] # all missing cases before 1991
tt[is.na(educ), educ := 99L]

## Recode race categories ----
# 0 = non-white, non-brown, non-black
# 1 = brown/black
# 2 = white
# 99 = missing
tt[, counts(race)]
tt[is.na(race) | race == 99, .N, keyby = year]
tt[, race3 := 0L]
tt[race == 51 | between(race, 20, 24), race3 := 1L] # brown/black
tt[race == 10, race3 := 2L] # white
tt[is.na(race) | race == 99, race3 := 99L] # ignored
tt[, counts(race3)]

#create data table for 10 years and older from 1991
xx <- tt[year > 1990, .(pop = sum(perwt, na.rm = T)), keyby = .(year, age, sex, race3, educ)]
all.equal(tt[year > 1990, sum(perwt, na.rm = T)], xx[, sum(pop)])

## Treat missing values ----
# assumes missing race cases have similar age/sex distribution as the non missing
xx[, table(race3, educ)]
strata <- c("year", "age", "sex")
gg <- xx[race3 == 99, .(pop = sum(pop)), keyby = strata]
mm <- xx[race3 != 99, .(pop = sum(pop)), keyby = c(strata, "race3")]
mm[, per := pop/sum(pop), keyby = strata]
mm[gg, race_99 := i.pop, on = strata]
mm[race3 == 99 & per == 1, .N]
mm[!is.na(race_99), pop := pop + race_99 * per]
setkeyv(mm, c(strata, "race3"))
mm[, c("per","race_99") := NULL]
all.equal(xx[, sum(pop)], mm[, sum(pop)])
mm[, table(race3)]

# NOT NECESSARY as no more missing Educ cases
# assumes missing educ cases have similar age/sex/race distribution as the non missing
gg <- xx[educ != 99 & race3 != 99, pop, keyby = c(strata, "race3", "educ")]
gg[, per := pop/sum(pop), keyby = c(strata, "race3")]
gg[, sum(per), keyby = c(strata, "race3")][, table(V1)]
gg[mm, pop := per * i.pop, on = c(strata, "race3")]
gg[, per := NULL]
all.equal(xx[, sum(pop)], gg[, sum(pop)])

xx <- copy(gg)
rm(gg, mm, tt)
xx[, table(race3, educ)]



xx[, .(pop = sum(pop)), keyby = .(year, sex)
   ][, qplot(year, pop, col = factor(sex), geom = "line")]

xx[, .(pop = sum(pop)), keyby = .(year, race3)
   ][, qplot(year, pop, col = factor(race3), geom = "line")]

xx[, .(pop = sum(pop)), keyby = .(year, educ)
   ][, qplot(year, pop, col = factor(educ), geom = "line")]

xx[, .(pop = sum(pop)), keyby = .(year, age)
   ][, qplot(year, pop, col = factor(age), geom = "line")]

ggplot(xx[, .(pop = sum(pop)), keyby = .(year, sex, race3, educ)
          ], aes(x = year, y = pop, col = factor(sex))) +
  geom_line() +
  facet_grid(race3~educ)

# Hamilton-Perry method for intracensus years -----------
# see Swanson et al 2010
zz <- dcast(xx[year > 1999], age + sex + race3 + educ ~ year, value.var = "pop") # years appear as columns
setnames(zz, c("2000","2010"), c("y2000","y2010"))

ttt <- CJ(age = 0:100, sex = 1:2, race3 = 0:2, educ = 0:2) # all combinations of age/sex/race/educ
zz[ttt, on = c("age", "sex", "race3", "educ")
   ][is.na(y2000) & is.na(y2010), table(age)]

zz <- zz[ttt, on = c("age", "sex", "race3", "educ")]
zz[age > 80, age := 80]
zz[is.na(y2000), y2000 := 0]
zz[is.na(y2010), y2010 := 0]

zz <- zz[, .(y2000 = sum(y2000), y2010 = sum(y2010)),
         keyby = .(age, sex, race3, educ)]

# shifting columns
# data.table(a = 1:20, b = shift(1:20,10,0,"lag"))
zz[, .N, by = .(race3, educ, sex)][, unique(N)]

# Calculate CCR
# for ages <30 do not use the lagged vars because educ evolves
zzyoung <- zz[age < 40, ]

# limits and calibration
zzyoung[, CCR10_00   := y2010 / y2000]
zzyoung[CCR10_00 < 0.5, CCR10_00 := 0.5]
zzyoung[CCR10_00 > 2,   CCR10_00 := 2]
zzyoung[is.na(CCR10_00), CCR10_00 := 0]

zzyoung[, y2020 := y2010 * CCR10_00]
tt <- zzyoung[, .(pop = sum(y2020)), keyby = .(age, sex)]
tt[calib[year == 2020, ], on = c("age", "sex"), calibpop := i.pop]
tt[, calib_wt := calibpop/pop][]
zzyoung[tt, on = c("age", "sex"), CCR10_00 := CCR10_00 * i.calib_wt]
zzyoung[, y2020      := y2010 * CCR10_00]
tt <- zzyoung[, .(pop = sum(y2020)), keyby = .(age, sex)]
tt[calib[year == 2020, ], on = c("age", "sex"), calibpop := i.pop]
tt[, calib_wt := calibpop/pop][]


# zzyoung[, diff00_10  := (y2010 - y2000)/10]
# zzyoung[, diff10_20  := (y2020 - y2010)/10]
# for (ii in 1:9) zzyoung[,paste0("y200", ii) := y2000 + (diff00_10 * ii)]
# for (ii in 1:9) zzyoung[,paste0("y201", ii) := y2010 + (diff10_20 * ii)]
# for (j in 5:ncol(zzyoung)) set(zzyoung, which(zzyoung[[j]] == 1), j, 0)

zzyoung[y2000 == 0, counts(educ)]
zzyoung[, diff00_10  := y2010 / y2000]
zzyoung[, diff10_20  := y2020 / y2010]

zzyoung[is.na(diff00_10) | is.infinite(diff00_10), diff00_10 := 0]
zzyoung[is.na(diff10_20) | is.infinite(diff10_20), diff10_20 := 0]

for (ii in 1:9) zzyoung[,paste0("y200", ii) := y2000 * diff00_10 ^ (ii/10)]
for (ii in 1:9) zzyoung[,paste0("y201", ii) := y2010 * diff10_20 ^ (ii/10)]

zzyoung <- zzyoung[, c("age", "sex", "race3", "educ", paste0("y", 2000:2020))]
setnames(zzyoung, paste0("y", 2000:2020), paste0(2000:2020))
zzyoung <- melt(zzyoung, 1:4, variable.name = "year", value.name = "pop")
zzyoung[, year := as.integer(as.character(year))]

zzold   <- zz[age >= 20, ]
zzold[, lag_y2000 := shift(y2000, 10, NA, "lag"), by = .(sex, race3, educ)]
zzold[, CCR10_00  := y2010 / lag_y2000]

# limits and calibration
zzold[CCR10_00 < 0.7, CCR10_00 := 0.7]
zzold[CCR10_00 > 1.5,   CCR10_00 := 1.5]
zzold[is.na(CCR10_00), CCR10_00 := 0]
zzold[, CCR10_00 := shift(CCR10_00, 10, 0, "lag"), by = .(sex, race3, educ)]
zzold[, y2020 := y2010 * CCR10_00]

tt <- zzold[, .(pop = sum(y2020)), keyby = .(age, sex)]
tt[calib[year == 2020, ], on = c("age", "sex"), calibpop := i.pop]
tt[, calib_wt := calibpop/pop]
tt[is.infinite(calib_wt), calib_wt := 1][]
zzold[tt, on = c("age", "sex"), CCR10_00 := CCR10_00 * i.calib_wt]
zzold[, y2020      := y2010 * CCR10_00]
tt <- zzold[, .(pop = sum(y2020)), keyby = .(age, sex)]
tt[calib[year == 2020, ], on = c("age", "sex"), calibpop := i.pop]
tt[, calib_wt := calibpop/pop][]
zzold[, summary(CCR10_00)]

#calculate differnce
zzold[y2000 == 0, y2000 := 10] # to avoid divide by 0
zzold[, diff00_10 := y2010 / y2000]
zzold[, diff10_20 := y2020 / y2010]
zzold[is.infinite(diff00_10), ]

for (ii in 1:9) zzold[, paste0("y200", ii) := y2000 * diff00_10 ^ (ii/10)]

for (ii in 1:9) zzold[, paste0("y201", ii) := y2010 * diff10_20 ^ (ii/10)]
# View(zzold)

# remove 1 values
#for (j in 5:ncol(zzold)) set(zzold, which(zzold[[j]] == 1), j, 0)
zzold <- zzold[age >= 40, ]
zzold <- zzold[, c("age", "sex", "race3", "educ", paste0("y", 2000:2020))]
setnames(zzold, paste0("y", 2000:2020), paste0(2000:2020))
zzold <- melt(zzold, 1:4, variable.name = "year", value.name = "pop")
zzold[, year := as.integer(as.character(year))]
zzold[, summary(pop)]

zz <- rbind(zzyoung, zzold)

# calibration
tt <- zz[, .(pop = sum(pop)), keyby = .(year, age, sex)]
tt[calib[between(year, 2000, 2020) ], on = c("year", "age", "sex"), calibpop := i.pop]
tt[, calib_wt := calibpop/pop][, summary(calib_wt)]
zz[tt, on = c("year", "age", "sex"), pop := pop * i.calib_wt]
tt <- zz[, .(pop = sum(pop)), keyby = .(year, age, sex)]
tt[calib[between(year, 2000, 2020) ], on = c("year", "age", "sex"), calibpop := i.pop]
tt[, calib_wt := calibpop/pop][, summary(calib_wt)]


# Graphs --------------------------------------
zz[, .(pop = sum(pop)), keyby = .(year, sex)
   ][, qplot(year, pop, col = factor(sex), geom = "line")]
xx[year > 1999, .(pop = sum(pop)), keyby = .(year, sex)
   ][, qplot(year, pop, col = factor(sex), geom = "line")]
zz[, .(pop = sum(pop)), keyby = .(year, race3)
   ][, qplot(year, pop, col = factor(race3), geom = "line")]
xx[year > 1999, .(pop = sum(pop)), keyby = .(year, race3)
   ][, qplot(year, pop, col = factor(race3), geom = "line")]
zz[, .(pop = sum(pop)), keyby = .(year, educ)
   ][, qplot(year, pop, col = factor(educ), geom = "line")]
xx[year > 1999, .(pop = sum(pop)), keyby = .(year, educ)
   ][, qplot(year, pop, col = factor(educ), geom = "line")]
zz[, .(pop = sum(pop)), keyby = .(year, age)
   ][, qplot(year, log(pop), col = factor(age), geom = "line")]
xx[year > 1999, .(pop = sum(pop)), keyby = .(year, age)
   ][, qplot(year, log(pop), col = factor(age), geom = "line")]

ggplot(zz[, .(pop = sum(pop)), keyby = .(year, sex, race3, educ)
          ], aes(x = year, y = pop, col = factor(sex))) +
  geom_line() +
  facet_grid(race3~educ)

# zz[, pop := round(pop)]


fwrite(zz,"~/pCloudDrive/My Models/Brazil_salt_model/Population_statistics/pop_estimates_2000_2020_age_sex_race_educ.csv")


