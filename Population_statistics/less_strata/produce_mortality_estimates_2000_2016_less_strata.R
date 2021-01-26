options(scipen = 999)
if (!require(Kmisc)) {devtools::install_github("kevinushey/Kmisc"); library("Kmisc")}
library("data.table")
library("cowplot")

ll <- vector("list", 4)
names(ll) <- c("mort_total", "mort_CHD", "mort_gastric", "mort_stroke")
ll$mort_total <- fread("~/pCloudDrive/My Datasets/Salt_Brazil/Mortality/Brazil total deaths Brazil.csv")[year > 1999, ]
ll$mort_CHD   <- fread("~/pCloudDrive/My Datasets/Salt_Brazil/Mortality/Brazil CHD deaths Brazil.csv")[year > 1999, ]
ll$mort_gastric <- fread("~/pCloudDrive/My Datasets/Salt_Brazil/Mortality/Brazil GCa deaths Brazil.csv")[year > 1999, ]
ll$mort_stroke <- fread("~/pCloudDrive/My Datasets/Salt_Brazil/Mortality/Brazil stroke deaths Brazil.csv")[year > 1999, ]

nam <- c("sex", "race", "education", "age_group")
check_cols <- function(x) unlist(unique(lapply(ll, function(dt) unique(dt[, get(x)]))))
lapply(nam, check_cols)
lapply(ll, function(x) x[, summary(deaths)])

treat <- function(dt, colname, oldvals, newvals, newcolname = NULL, erase = TRUE) {
  if (identical(colname, newcolname)) stop("colname and newcolname cannot be the same")
  tt <- data.table("V1" = oldvals, "V2" = newvals)
  if (is.null(newcolname)) {
    setnames(tt, "V1", ".")
    setnames(dt, colname, ".")
    dt[tt, on = ".", (colname) := i.V2]
    dt[, "." := NULL]
  } else {
    setnames(tt, "V1", colname)
    dt[tt, on = colname, (newcolname) := i.V2]
  }
  if (!is.null(newcolname) && erase) dt[, (colname) := NULL]
  dt
}

lapply(ll, treat, "sex", c("Man", "Woman"), 1:2)
lapply(ll, treat, "race", c("White", "Black", "Asian", "Brown", "Indigenous", "Ignored"),
       c(2L, 1L, 1L, 1L, 1L, 99L), "race3")
lapply(ll, treat, "education", c("None", "1 to 3 years", "4 to 7 years", "8 to 11 years", ">=12 years", "Ignored"),
       c(1L, 1L, 2L, 2L, 2L, 99L), "educ")
lapply(ll, treat, "age_group", c("<1 year", "1 to 4 years", "5 to 9 years",
                                 "10 to 14 years", "15 to 19 years", "20 to 24 years",
                                 "25 to 29 years", "30 to 34 years", "35 to 39 years",
                                 "40 to 44 years", "45 to 49 years", "50 to 54 years",
                                 "55 to 59 years", "60 to 64 years", "65 to 69 years",
                                 "70 to 74 years", "75 to 79 years", ">=80 years",
                                 "Age ignored"),
       c("00-04", "00-04","05-09", "10-14", "15-19", "20-24", "25-29", "30-34", "35-39",
         "40-44", "45-49", "50-54", "55-59", "60-64", "65-69", "70-74", "75-79",
         "80+", "99"), "agegroup")

ll <- lapply(ll, function(x) x[, .(deaths = sum(deaths, na.rm = TRUE)), keyby = .(year, agegroup, sex, race3, educ)])
names(ll) <- c("mort_total", "mort_CHD", "mort_gastric", "mort_stroke")
lapply(ll, function(x) x[, summary(deaths)])
lapply(c("sex", "race3", "educ", "agegroup"), check_cols)


# Redistribute missing cases ----------------------------------------------
ggplot(ll$mort_total[, .(deaths = sum(deaths)), keyby = .(year, agegroup)],
       aes(x = year, y = log(deaths), col = agegroup)) +
  geom_line()

lapply(ll, function(x) signif(x[agegroup == 99, sum(deaths)/sum(x$deaths)], 2))
lapply(ll, function(x) signif(x[sex == 99, sum(deaths)/sum(x$deaths)], 2))
lapply(ll, function(x) signif(x[race3 == 99, sum(deaths)/sum(x$deaths)], 2))
lapply(ll, function(x) signif(x[educ == 99, sum(deaths)/sum(x$deaths)], 2))

# function that calculates weights for the known cases to redistribute unknown cases assuming similar characteristics
redistribute_99 <- function(x) {
  dt <- copy(x) # for safety
  dt[agegroup == "00-04", educ := 1]
  dt[agegroup %in% c("05-09", "10-14") & educ == 2, educ := 99]

  dt <-
    dt[, .("deaths" = as.numeric(sum(deaths))), keyby = c("year", "sex", "agegroup", "race3", "educ")]

  all_deaths <- as.numeric(sum(dt$deaths))

  out <-
    CJ(
      year = 2000:2016,
      agegroup = c(
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
        "80+"
      ),
      sex = c(1:2),
      race3 = c(1:2),
      educ = c(1:2)
    )
  # year
  all_deaths_nomissing <-
    dt[year != 99, as.numeric(sum(deaths)), keyby = NULL]
  if (any(all_deaths_nomissing$V1 == 0)) print("0 in year")
  t1 <-
    dt[year != 99, sum(deaths) / all_deaths_nomissing$V1, keyby = "year"]
  if (any(is.na(t1$V1))) print("NA in year")
  out[t1, wt_year := i.V1, on = "year"]
  stopifnot(all_deaths == sum(dt$deaths))

  # year + sex
  all_deaths_nomissing <-
    dt[year != 99 & sex != 99, as.numeric(sum(deaths)), keyby = "year"]
  if (any(all_deaths_nomissing$V1 == 0)) print("0 in sex")
  t2 <-
    dt[year != 99 &
         sex != 99, sum(deaths) / all_deaths_nomissing[year == .BY[[1]], V1], keyby = c("year", "sex")]
  if (any(is.na(t2$V1))) print("NA in sex")
  out[t2, wt_sex := i.V1, on = c("year", "sex")]
  # year + sex + agegroup
  all_deaths_nomissing <-
    dt[sex != 99 &
         year != 99 &
         agegroup != 99, as.numeric(sum(deaths)), keyby = c("year", "sex")]
  if (any(all_deaths_nomissing$V1 == 0)) print("0 in age group")
  t3 <- dt[sex != 99 & year != 99 & agegroup != 99,
           sum(deaths) / all_deaths_nomissing[year == .BY[[1]] &
                                                sex ==  .BY[[2]], V1], keyby = c("year", "sex", "agegroup")]
  if (any(is.na(t3$V1))) print("NA in age group")
  out[t3, wt_agegroup := i.V1, on = c("year", "sex", "agegroup")]
  stopifnot(all_deaths == sum(dt$deaths))

  # year + sex + agegroup + race3
  all_deaths_nomissing <-
    dt[sex != 99 & year != 99 & agegroup != 99 & race3 != 99,
       as.numeric(sum(deaths)),
       keyby = c("year", "sex", "agegroup")]
  if (any(all_deaths_nomissing$V1 == 0)) {
    print("0 in race3")
    tt <- dt[sex != 99 & year != 99 & agegroup != 99 & race3 == 99,
             as.numeric(sum(deaths)),
             keyby = c("year", "sex", "agegroup")
             ][all_deaths_nomissing[V1 == 0], on = c("year", "sex", "agegroup")
               ][V1 > 0]
    tt[, i.V1 := NULL]
    if (nrow(tt) > 0) {
      print("Race3 has strata with exclussively missing cases. I evenly redistribute them")
      dt[tt, on = c("year", "sex", "agegroup"), trgt := V1]
      dt[trgt > 0 & race3 != 99 & educ != 99, deaths := trgt / .N,
         by = c("year", "sex", "agegroup")]
      dt[trgt > 0 & (race3 == 99 | educ == 99), deaths := 0]
      dt[, trgt := NULL]
    }
    all_deaths_nomissing <-
      dt[sex != 99 & year != 99 & agegroup != 99 & race3 != 99,
         as.numeric(sum(deaths)),
         keyby = c("year", "sex", "agegroup")]
  }


  t4 <- dt[sex != 99 & year != 99 & agegroup != 99 & race3 != 99,
           sum(deaths) / all_deaths_nomissing[year == .BY[[1]] &
                                                sex ==  .BY[[2]] & agegroup ==  .BY[[3]], V1],
           keyby = c("year", "sex", "agegroup", "race3")]
  if (any(is.na(t4$V1))) {
    print("NA in race3 replaced with 0")
    t4[is.na(V1), V1 := 0]
  }
  out[t4, wt_race3 := i.V1, on = c("year", "sex", "agegroup", "race3")]
  stopifnot(all_deaths == sum(dt$deaths))

  # year + sex + agegroup + race3 + educ
  all_deaths_nomissing <-
    dt[sex != 99 &
         year != 99 & agegroup != 99 & race3 != 99 & educ != 99,
       as.numeric(sum(deaths)),
       keyby = c("year", "sex", "agegroup", "race3")]
  if (any(all_deaths_nomissing$V1 == 0)) {
    print("0 in educ")
    tt <- dt[sex != 99 & year != 99 & agegroup != 99 & race3 != 99 & educ == 99,
             as.numeric(sum(deaths)),
             keyby = c("year", "sex", "agegroup", "race3")
             ][all_deaths_nomissing[V1 == 0], on = c("year", "sex", "agegroup", "race3")
               ][V1 > 0]

    if (nrow(tt) > 0) {
      print("Educ has strata with exclussively missing cases. I evenly redistribute them")
      dt[tt, on = c("year", "sex", "agegroup", "race3"), trgt := V1]
      dt[trgt > 0 & educ != 99, deaths := trgt / .N,
         by = c("year", "sex", "agegroup", "race3")]
      dt[trgt > 0 & (race3 == 99 | educ == 99), deaths := 0]
      dt[, trgt := NULL]
    }

    all_deaths_nomissing <-
      dt[sex != 99 &
           year != 99 & agegroup != 99 & race3 != 99 & educ != 99,
         as.numeric(sum(deaths)),
         keyby = c("year", "sex", "agegroup", "race3")]
  }

  t5 <-
    dt[sex != 99 &
         year != 99 & agegroup != 99 & race3 != 99 & educ != 99,
       sum(deaths) / all_deaths_nomissing[year == .BY[[1]] &
                                            sex ==  .BY[[2]] & agegroup ==  .BY[[3]] & race3 ==  .BY[[4]], V1],
       keyby = c("year", "sex", "agegroup", "race3", "educ")]
  if (any(is.na(t5$V1))) {
    print("NA in educ replaced with 0")
    t5[is.na(V1), V1 := 0]
  }
  out[t5, wt_educ := i.V1, on = c("year", "sex", "agegroup", "race3", "educ")]
  stopifnot(all_deaths == sum(dt$deaths))


  out[, deaths := all_deaths * wt_year * wt_agegroup * wt_sex * wt_race3 * wt_educ]
  out[is.na(deaths), deaths := 0]

  print(paste0("Original ", all_deaths, ", imputed ", sum(out$deaths),
               ", difference ~",round(all_deaths - sum(out$deaths))))

  out[, c("wt_year", "wt_sex", "wt_agegroup", "wt_race3", "wt_educ") := NULL]
  return(out)
}

lll <- lapply(ll, redistribute_99)

ggplot(ll$mort_total[, .(deaths = sum(deaths)), keyby = .(year, agegroup)],
       aes(x = year, y = log(deaths), col = agegroup)) +
  geom_line()
ggplot(lll$mort_total[, .(deaths = sum(deaths)), keyby = .(year, agegroup)],
       aes(x = year, y = log(deaths), col = agegroup)) +
  geom_line()

ggplot(ll$mort_total[, .(deaths = sum(deaths)), keyby = .(year, sex)],
       aes(x = year, y = log(deaths), col = factor(sex))) +
  geom_line()
ggplot(lll$mort_total[, .(deaths = sum(deaths)), keyby = .(year, sex)],
       aes(x = year, y = log(deaths), col = factor(sex))) +
  geom_line()

ggplot(ll$mort_total[, .(deaths = sum(deaths)), keyby = .(year, race3)],
       aes(x = year, y = log(deaths), col = factor(race3))) +
  geom_line()
ggplot(lll$mort_total[, .(deaths = sum(deaths)), keyby = .(year, race3)],
       aes(x = year, y = log(deaths), col = factor(race3))) +
  geom_line()

ggplot(ll$mort_total[, .(deaths = sum(deaths)), keyby = .(year, educ)],
       aes(x = year, y = log(deaths), col = factor(educ))) +
  geom_line()
ggplot(lll$mort_total[, .(deaths = sum(deaths)), keyby = .(year, educ)],
       aes(x = year, y = log(deaths), col = factor(educ))) +
  geom_line()


# Rates -------------------------------------------------------------------

pop <- fread("~/pCloudDrive/My Models/Brazil_salt_model/Population_statistics/less_strata/pop_estimates_2000_2020_age_sex_race_educ_less_strata.csv")

age_to_agegroup <- function(dt) {
  stopifnot(is.data.table(dt))
  aux <- data.table(age = 0:80,
                    agegroup = c(
                      rep(c(
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
                      "75-79"), each = 5),
                      "80+"
                    ))
  dt[aux, on = "age", agegroup := i.agegroup]
  dt[, .("pop" = sum(pop)), keyby = c("year", "agegroup", "sex", "race3", "educ")]
}
pop <- age_to_agegroup(pop)

lapply(lll, function(x) {
  x[pop, on = c("year", "agegroup", "sex", "race3", "educ"), pop := i.pop]
  x[, mortality_rate := deaths / pop]
  x[is.na(mortality_rate), mortality_rate := 0]
})

# Graphs ------------------------------------------------------------------
tt <- copy(pop)
tt[, agegroup := factor(agegroup)]
tt[, sex := factor(sex, 1:2, c("men", "women"))]
tt[, race := factor(race3, 1:2, c("other", "white"))]
tt[, education := factor(educ, 1:2, c("<4 years", "4+ years"))]
gg <- ggplot(tt, aes(x = year, y = pop/1e6, col = education)) +
  geom_path() +
  facet_grid(sex + race ~ agegroup) +
  ggtitle("Population size", subtitle = "Brazil 2000-2020") +
  labs(x = "Year", y = "Population (millions)")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggsave("~/pCloudDrive/My Models/Brazil_salt_model/Population_statistics/less_strata/population.tiff", gg, width = 16, height = 9, units = "cm", scale = 3, compression = "lzw")

tt <- copy(lll$mort_total)
tt[, sex := factor(sex, 1:2, c("men", "women"))]
tt[, race := factor(race3, 1:2, c("other", "white"))]
tt[, education := factor(educ, 1:2, c("<4 years", "4+ years"))]

gg <- ggplot(tt, aes(x = year, y = mortality_rate * 1e5, col = education)) +
  geom_line() +
  facet_grid(sex + race ~ agegroup) +
  ggtitle("All-cause mortality", subtitle = "Brazil 2000-2016") +
  labs(x = "Year", y = "Mortality per 100,000")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggsave("~/pCloudDrive/My Models/Brazil_salt_model/Population_statistics/less_strata/mortality_all-cause.tiff", gg, width = 16, height = 9, units = "cm", scale = 3, compression = "lzw")

tt <- copy(lll$mort_CHD)
tt[, sex := factor(sex, 1:2, c("men", "women"))]
tt[, race := factor(race3, 1:2, c("other", "white"))]
tt[, education := factor(educ, 1:2, c("<4 years", "4+ years"))]

gg <- ggplot(tt, aes(x = year, y = mortality_rate * 1e5, col = education)) +
  geom_line() +
  facet_grid(sex + race ~ agegroup) +
  ggtitle("CHD mortality", subtitle = "Brazil 2000-2016") +
  labs(x = "Year", y = "Mortality per 100,000") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave("~/pCloudDrive/My Models/Brazil_salt_model/Population_statistics/less_strata/mortality_chd.tiff", gg, width = 16, height = 9, units = "cm", scale = 3, compression = "lzw")

tt <- copy(lll$mort_stroke)
tt[, sex := factor(sex, 1:2, c("men", "women"))]
tt[, race := factor(race3, 1:2, c("other", "white"))]
tt[, education := factor(educ, 1:2, c("<4 years", "4+ years"))]

gg <- ggplot(tt, aes(x = year, y = mortality_rate * 1e5, col = education)) +
  geom_line() +
  facet_grid(sex + race ~ agegroup) +
  ggtitle("Stroke mortality", subtitle = "Brazil 2000-2016") +
  labs(x = "Year", y = "Mortality per 100,000")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggsave("~/pCloudDrive/My Models/Brazil_salt_model/Population_statistics/less_strata/mortality_stroke.tiff", gg, width = 16, height = 9, units = "cm", scale = 3, compression = "lzw")

tt <- copy(lll$mort_gastric)
tt[, sex := factor(sex, 1:2, c("men", "women"))]
tt[, race := factor(race3, 1:2, c("other", "white"))]
tt[, education := factor(educ, 1:2, c("<4 years", "4+ years"))]

gg <- ggplot(tt, aes(x = year, y = mortality_rate * 1e5, col = education)) +
  geom_line() +
  facet_grid(sex + race ~ agegroup) +
  ggtitle("Gastric cancer mortality", subtitle = "Brazil 2000-2016") +
  labs(x = "Year", y = "Mortality per 100,000")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggsave("~/pCloudDrive/My Models/Brazil_salt_model/Population_statistics/less_strata/mortality_gastric.tiff", gg, width = 16, height = 9, units = "cm", scale = 3, compression = "lzw")

tt <- copy(lll$mort_total)
tt[, sex := factor(sex, 1:2, c("men", "women"))]
tt[, race := factor(race3, 1:2, c("other", "white"))]
tt[, education := factor(educ, 1:2, c("<4 years", "4+ years"))]

gg <- ggplot(tt[!agegroup %in% c("00-04", "05-09", "10-14", "15-19", "20-24", "25-29"),
                sum(deaths)/sum(pop), keyby = .(year, education)],
             aes(x = year, y = V1 * 1e5, col = education)) +
  geom_line() +
  ggtitle("All-cause mortality", subtitle = "Brazil 2000-2016, ages 30+") +
  labs(x = "Year", y = "Mortality per 100,000")
ggsave("~/pCloudDrive/My Models/Brazil_salt_model/Population_statistics/less_strata/mortality_all-cause_educ.tiff", gg, width = 16, height = 9, units = "cm", scale = 3, compression = "lzw")

gg <- ggplot(tt[!agegroup %in% c("00-04", "05-09", "10-14", "15-19", "20-24", "25-29"),
                sum(deaths)/sum(pop), keyby = .(year, race)],
             aes(x = year, y = V1 * 1e5, col = race)) +
  geom_line() +
  ggtitle("All-cause mortality", subtitle = "Brazil 2000-2016, ages 30+") +
  labs(x = "Year", y = "Mortality per 100,000")
ggsave("~/pCloudDrive/My Models/Brazil_salt_model/Population_statistics/less_strata/mortality_all-cause_race.tiff", gg, width = 16, height = 9, units = "cm", scale = 3, compression = "lzw")

gg <- ggplot(tt[!agegroup %in% c("00-04", "05-09", "10-14", "15-19", "20-24", "25-29"),
                sum(deaths)/sum(pop), keyby = .(year, sex)],
             aes(x = year, y = V1 * 1e5, col = sex)) +
  geom_line() +
  ggtitle("All-cause mortality", subtitle = "Brazil 2000-2016, ages 30+") +
  labs(x = "Year", y = "Mortality per 100,000")
ggsave("~/pCloudDrive/My Models/Brazil_salt_model/Population_statistics/less_strata/mortality_all-cause_sex.tiff", gg, width = 16, height = 9, units = "cm", scale = 3, compression = "lzw")

gg <- ggplot(tt[!agegroup %in% c("00-04", "05-09", "10-14", "15-19", "20-24", "25-29"),
                sum(deaths)/sum(pop), keyby = .(year, agegroup)],
             aes(x = year, y = V1 * 1e5, col = agegroup)) +
  geom_line() +
  ggtitle("All-cause mortality", subtitle = "Brazil 2000-2016, ages 30+") +
  labs(x = "Year", y = "Mortality per 100,000")
ggsave("~/pCloudDrive/My Models/Brazil_salt_model/Population_statistics/less_strata/mortality_all-cause_agegroup.tiff", gg, width = 16, height = 9, units = "cm", scale = 3, compression = "lzw")

gg <- ggplot(tt, aes(x = year, y = deaths/1e3, col = education)) +
  geom_path() +
  facet_grid(sex + race ~ agegroup) +
  ggtitle("All-cause deaths", subtitle = "Brazil 2000-2016") +
  labs(x = "Year", y = "Deaths (in thousands)")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggsave("~/pCloudDrive/My Models/Brazil_salt_model/Population_statistics/less_strata/deaths.tiff", gg, width = 16, height = 9, units = "cm", scale = 3, compression = "lzw")

# Write files -------------------------------------------------------------
fwrite(lll$mort_total,"~/pCloudDrive/My Models/Brazil_salt_model/Population_statistics/less_strata/mortality_all_cause_estimates_2000_2016_age_sex_race_educ_less_strata.csv")
fwrite(lll$mort_CHD,"~/pCloudDrive/My Models/Brazil_salt_model/Population_statistics/less_strata/mortality_chd_estimates_2000_2016_age_sex_race_educ_less_strata.csv")
fwrite(lll$mort_stroke,"~/pCloudDrive/My Models/Brazil_salt_model/Population_statistics/less_strata/mortality_stroke_estimates_2000_2016_age_sex_race_educ_less_strata.csv")
fwrite(lll$mort_gastric,"~/pCloudDrive/My Models/Brazil_salt_model/Population_statistics/less_strata/mortality_gastric_cancer_estimates_2000_2016_age_sex_race_educ_less_strata.csv")
