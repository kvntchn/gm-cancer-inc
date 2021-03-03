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

year.max <- 2015
employment_status.lag <- 0

# Source 00-hello.R
if (!('cohort' %in% ls(envir = .GlobalEnv))) {
	source(here::here('../gm-wrangling/wrangling', '00-hello.R'))
}

# Get cohort analytic ####
if (!('cohort_analytic' %in% ls())) {
	outcome.type <- 'incidence'
	cohort_analytic <- get.cohort_analytic(
		outcome_type = outcome.type,
		exposure.lag = exposure.lag,
		deathage.max = NULL,
		year.max = year.max,
		hire.year.min = -Inf,
		use_seer = T
	)
	setorder(cohort_analytic, studyno, year)
	cohort_analytic[, `:=`(yin.gm = date.to.gm(yin))]

	# Keep only people who appear in the exposure data
	cohort_analytic <- cohort_analytic[studyno %in% unique(exposure$studyno)]

	# PICK YOUT ####
	cohort_analytic[, jobloss.date := get(yout.which)]

	# Exposure after leaving work is 0
	cohort_analytic[year > (year(jobloss.date) + exposure.lag), `:=`(
		straight = 0,
		soluble = 0,
		synthetic = 0,
		bio = 0,
		cl = 0,
		ea = 0,
		tea = 0,
		trz = 0,
		s = 0,
		no2 = 0)]
	# NA fill
	cohort_analytic[year <= (year(jobloss.date) + exposure.lag), `:=`(
		straight = zoo::na.locf(straight),
		soluble = zoo::na.locf(soluble),
		synthetic = zoo::na.locf(synthetic),
		bio = zoo::na.locf(bio),
		cl = zoo::na.locf(cl),
		ea = zoo::na.locf(ea),
		tea = zoo::na.locf(tea),
		trz = zoo::na.locf(trz),
		s = zoo::na.locf(s),
		no2 = zoo::na.locf(no2)
	), by = .(studyno)]
	table(cohort_analytic[,.(N = length(table(jobloss.date))), by = .(studyno)]$N)
	cohort_analytic[, `:=`(
		cum_straight = cumsum(straight),
		cum_soluble = cumsum(soluble),
		cum_synthetic = cumsum(synthetic),
		cum_bio = cumsum(bio),
		cum_cl = cumsum(cl),
		cum_ea = cumsum(ea),
		cum_tea = cumsum(tea),
		cum_trz = cumsum(trz),
		cum_s = cumsum(s),
		cum_no2 = cumsum(no2)
	), by = .(studyno)]

	# Time off
	cohort_analytic[,`:=`(
		off = off.gan + off.san + off.han
	)]
	cohort_analytic[is.na(off), off := 0]
	cohort_analytic[, cum_off := cumsum(off), by = .(studyno)]

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
		"sex",
		"dateout.date",
		"employment_end.date",
		"employment_end.date.legacy",
		"off", "cum_off",
		"yout",
		"yout_recode",
		"jobloss.date",
		"All causes",
		"Chronic obstructive pulmonary disease",
		"All external causes",
		"nohist", "wh", "immortal", "right.censored",
		"possdiscr_new", "flag77", "oddend",
		"status15", "cancinccoh15_new"), with = F])

	# Drop unnecessary data ####
	cohort_analytic <- cohort_analytic[
		wh == 1 & nohist == 0 &
			# cancinccoh15_new == 1 &
			possdiscr_new == 0,
		# & immortal == 0 & right.censored == 0,
		col.names, with = F]
}

Sys.sleep(0)
# Get cohort HWSE 2 ####
if (get.cohort2 & !("cohort2" %in% ls())) {

		cohort2 <- as.data.table(as.data.frame(cohort_analytic))
		setorder(cohort2, studyno, year)
		cohort2.og <- as.data.table(as.data.frame(cohort2))

	if (exposure.lag + 1994 > 2015) {stop("Needlessly large exposure lag.")
	} else {
		cohort2[, `:=`(
			straight = shift(straight, 1 - exposure.lag, fill = 0),
			soluble = shift(soluble, 1 - exposure.lag, fill = 0),
			synthetic = shift(synthetic, 1 - exposure.lag, fill = 0),
			bio = shift(bio, 1 - exposure.lag, fill = 0),
			cl = shift(cl, 1 - exposure.lag, fill = 0),
			ea = shift(ea, 1 - exposure.lag, fill = 0),
			tea = shift(tea, 1 - exposure.lag, fill = 0),
			trz = shift(trz, 1 - exposure.lag, fill = 0),
			s = shift(s, 1 - exposure.lag, fill = 0),
			no2 = shift(no2, 1 - exposure.lag, fill = 0)),
			by = .(studyno)]
	}

	cohort2[, `:=`(
		cum_straight = cumsum(straight),
		cum_soluble = cumsum(soluble),
		cum_synthetic = cumsum(synthetic),
		cum_bio = cumsum(bio),
		cum_cl = cumsum(cl),
		cum_ea = cumsum(ea),
		cum_tea = cumsum(tea),
		cum_trz = cumsum(trz),
		cum_s = cumsum(s),
		cum_no2 = cumsum(no2)
	), by = .(studyno)]

	# # Censor those still at work in 1995
	# cohort2 <- cohort2[year < 1995 | year(jobloss.date) < 1995]
	# # Stop FU in 1994
	# cohort2 <- cohort2[year < 1995,]

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Mortality cohort for positive/negative controls ####
	mortality.cohort2 <- as.data.table(as.data.frame(cohort2))

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

