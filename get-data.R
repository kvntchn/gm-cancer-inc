# Get data for cancer incidence analyses ####
# Kevin Chen
# December 15, 2020

library(here)
# here::set_here()

# rm(list = ls())
# rm(list = ls()[-grep("outcome.selected", ls())])
# rm(list = ls()[-grep('cohort', ls())])

library(survival)

if (!("exposure.lag" %in% ls())) {
	exposure.lag <- 21
	# rm(cohort_analytic)
}

if (!("get.cohort2" %in% ls())) {
	get.cohort2 <- T
	# rm(cohort_analytic)
}

if (!("yout.which" %in% ls())) {
	yout.which <- "yout"
	# rm(cohort_analytic)
}

year.max <- 2015

# Source 00-hello.R
if (!('cohort' %in% ls(envir = .GlobalEnv))) {
	source(here::here('../gm-wrangling/wrangling', '00-hello.R'))
}

# Getting posterior distribution for imputation ####
if (!"race.imputed.melt" %in% ls()) {
message("Loading posterior distribution for imputation...")
	# race.imputed.melt <- box_read(ifelse(exposure.lag == 21, 780327112319, NaN))
	box_load(ifelse(exposure.lag == 21, 783387319906, NaN))
}

# Get cohort analytic ####
if (!('cohort_analytic' %in% ls())) {
	outcome.type <- 'incidence'
	cohort_analytic <- get.cohort_analytic(
		outcome_type = outcome.type,
		exposure.lag = 1,
		deathage.max = if ("deathage.max" %in% ls(envir = .GlobalEnv)) {deathage.max} else {NULL},
		year.max = year.max,
		hire.year.min = -Inf,
		hire.year.max = Inf,
		use_seer = ifelse("use_seer" %in% ls(envir = .GlobalEnv), use_seer, T)
	)
	setorder(cohort_analytic, studyno, year)
	cohort_analytic[, `:=`(yin.gm = date.to.gm(yin))]

	# Keep only people who appear in the exposure data
	if (!"restrict_by_exposure" %in% ls(envir = .GlobalEnv)) {
		cohort_analytic <- cohort_analytic[studyno %in% unique(exposure$studyno)]
	} else {
		if (restrict_by_exposure) {
			cohort_analytic <- cohort_analytic[studyno %in% unique(exposure$studyno)]
		}
	}

	# PICK YOUT ####
	cohort_analytic[, jobloss.date := get(yout.which)]

	# Time off
	cohort_analytic[,`:=`(off = off.gan + off.san + off.han)]
	cohort_analytic[is.na(off), off := 0]
	cohort_analytic[, cum_off := cumsum(off), by = .(studyno)]

	exposure.names <- c("straight", "soluble", "synthetic", "bio", "cl",
											"ea", "tea", "trz", "s", "no2")

	# Exposure after leaving work is 0
	cohort_analytic[year > (year(jobloss.date) + 1), (exposure.names):=0]
	# NA fill
	cohort_analytic[year <= (year(jobloss.date) + 1), (exposure.names):=lapply(exposure.names, function (x) {
		x <- get(x)
		return(zoo::na.locf(x))
	}), by = .(studyno)]

	# table(cohort_analytic[,.(N = length(table(jobloss.date))), by = .(studyno)]$N)

	if (get.cohort2 & !("cohort2" %in% ls())) {
		cohort2 <- data.table::copy(cohort_analytic)
	}

	# Lag exposure appropriately
	cohort_analytic[, (exposure.names):=lapply(exposure.names, function (x) {
		x <- get(x)
		return(shift(x, get("exposure.lag", envir = .GlobalEnv) - 1, 0))
	}), by = .(studyno)]

	cohort_analytic[, (paste0("cum_", exposure.names)):=lapply(
		exposure.names, function(x) {
			x <- get(x)
			return(cumsum(x))
		}), by = .(studyno)]

	# Which columns ####
	col.names <- names(cohort_analytic[, c(
		"studyno",
		"age.year1",
		"age.year2",
		"year1",
		"year2",
		grep("canc\\_", names(cohort_analytic), value = T),
		paste0(c("cum_", ""), rep(c("straight", "soluble", "synthetic"), each = 2)),
		paste0(c("cum_", ""), rep(c("bio", "cl", "ea", "tea", "trz", "s", "no2"), each = 2)),
		"year",
		"yin.gm",
		"yin",
		"yrin",
		"yrin16",
		"race",
		"finrace",
		"plant",
		grep("ddiag", names(cohort_analytic), value = T),
		"yod",
		"yoc",
		"yob",
		"yob.gm",
		"sex",
		"dateout.date",
		"employment_end.date",
		"employment_end.date.legacy",
		"off", "cum_off",
		"yout", "yout95", "yout15", "yout16",
		"yout_recode",
		"jobloss.date",
		"All causes",
		"Chronic obstructive pulmonary disease",
		"All external causes",
		"nohist", "wh", "immortal", "right.censored",
		"possdiscr_new", "flag77", "oddend",
		"status15",
		"cancinccoh", "cancinccoh15_new", "alive85"), with = F])

	# Drop unnecessary data ####
	cohort_analytic <- cohort_analytic[
		wh == 1 & nohist == 0,
			# cancinccoh15_new == 1 &
			possdiscr_new == 0,
		# & immortal == 0 & right.censored == 0
		,
		col.names, with = F]
}

Sys.sleep(0)
# Get cohort HWSE 2 ####
if (get.cohort2 | !("cohort2" %in% ls())) {

	cohort2 <- cohort2[
		wh == 1 & nohist == 0,
			# cancinccoh15_new == 1 &
			possdiscr_new == 0,
		# & immortal == 0 & right.censored == 0
		,
		col.names, with = F]

	cohort2[, (paste0("cum_", exposure.names)):=lapply(
		exposure.names, function(x) {
			x <- get(x)
			return(cumsum(x))
		}), by = .(studyno)]

	cohort2.og <- data.table::copy(cohort2)

	# # Censor those still at work in 1995
	# cohort2 <- cohort2[year < 1995 | year(jobloss.date) < 1995]
	# # Stop FU in 1994
	# cohort2 <- cohort2[year < 1995,]

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Mortality cohort for positive/negative controls ####
	mortality.cohort2 <- data.table::copy(cohort2)

	# # Censor at 80
	# mortality.cohort2[, `:=`(yoc = as.Date(apply(data.frame(
	# 	yoc,
	# 	yob + years(80)
	# ), 1, min, na.rm = T)
	# ))]
	# # Anybody enter cohort too late?
	# mortality.cohort2 <- mortality.cohort2[yin + years(3) < as.Date(apply(data.frame(
	# 	yoc,
	# 	yob + years(80)
	# ), 1, min, na.rm = T)
	# )]
	# mortality.cohort2 <- mortality.cohort2[year <= year(yoc) | is.na(yoc)]
	# mortality.cohort2[year == year(yoc), `:=`(age.year2 = time_length(difftime(as.Date(
	# 	apply(data.frame(
	# 		as.Date(paste0(year.max + 1, "-01-01")),
	# 		yod + days(1),
	# 		yoc + days(1)
	# 	), 1, min, na.rm = T)
	# ), yob), 'day'))]

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	# # Censor at 80 or first cancer
	# cohort2[, `:=`(yoc = as.Date(apply(data.frame(
	# 	yoc,
	# 	yob + years(80),
	# 	ddiag_first
	# ), 1, min, na.rm = T)
	# ))]
	# # Anybody enter cohort too late?
	# cohort2 <- cohort2[yin + years(3) < as.Date(apply(data.frame(
	# 																					yoc,
	# 																					yob + years(80),
	# 																					ddiag_first
	# 																					), 1, min, na.rm = T)
	# 																					)]
	# cohort2 <- cohort2[year <= year(yoc) | is.na(yoc)]
	# cohort2[year == year(yoc), `:=`(age.year2 = time_length(difftime(as.Date(
	# 	apply(data.frame(
	# 		as.Date(paste0(year.max + 1, "-01-01")),
	# 		yod + days(1),
	# 		yoc + days(1)
	# 	), 1, min, na.rm = T)
	# ), yob), 'day'))]

} # End get HWSE 2 analytic data script

# # Clean up environment
# rm(list = grep('jobhist|exposure$', ls(), value = T))
# rm(list = ls()[!grepl("drive\\_D|end\\.year|cohort\\_analytic|cohort2|exposure\\.lag|^incidence", ls())])

incidence.key <- data.table::fread(here::here("../gm-wrangling/cancer incidence", 'cancer-key.tsv'))
# Expand incidence key ####
incidence.key[, `:=`(
	var.name = paste0("canc_", code),
	date.name = paste0("ddiag_", code))]
incidence.key <- rbindlist(
	list(incidence.key,
			 data.table(
			 	code = c("first", "copd", "external"),
			 	description = c(
			 		"All cancers",
			 		"Chronic obstructive pulmonary disease",
			 		"All external causes"),
			 	var.name = c(
			 		"canc_first",
			 		"Chronic obstructive pulmonary disease",
			 		"All external causes"),
			 	date.name = c(
			 		"ddiag_first",
			 		"yod", "yod")
			 )))

# fwrite(incidence.key, here::here("resources", 'cancer-key-expanded.csv'))

# Cases per outcome ####
get.nevent <- function(outcomes = 1:nrow(incidence.key),
											 cohort_name = "cohort_analytic") {
	data.table(incidence.key[outcomes, 2],
						 cases = sapply(outcomes, function(i = outcomes.which[1]) {
						 	# Get data ####
						 	var.name <- unlist(incidence.key[i, 3])
						 	sum(get(cohort_name)[year >= 1973,
						 											 var.name, with = F] == 1)
						 }))
}

nevent <- get.nevent()

# saveRDS(nevent,
# 				here::here("cancer incidence/resources",
# 									 "nevent.rds"))
# nevent <- readRDS(here::here("cancer incidence/resources",
# 														 "nevent.rds"))

