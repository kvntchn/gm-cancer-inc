# Standardized incidence ratios ####
# Kevin Chen
# March 19, 2020

# rm(list = ls())
# rm(list = ls()[-grep("outcome.selected", ls())])
# rm(list = ls()[-grep('cohort', ls())])

if (!grepl("gm", getwd(), ignore.case = T)) {
	if ("here" %in% .packages()) {
		detach("package:here", unload = T)
	}
	setwd('gm')
	library(here)
}

nevent <- readRDS(here::here("cancer incidence/resources",
														 "nevent.rds"))
incidence.key <- fread(here::here("cancer incidence", 'cancer-key.tsv'))

cancers.which <- c(3:8, 10:12, 15, 18, 19, 21, 25, 31)

# Read cleaned reference rates
# Note rates are per 100,000 (hundred thousand)
ref_rates <- rbindlist(apply(data.frame(
	race = c("white", "white", "black", "black"),
	gender = rep(c("male", "female"), 2)
), 1, function(x){
	race <- unlist(x[1]); gender <- unlist(x[2])
	fread(here::here("seer rates/rates cleaned", paste0(race, "_", gender, ".csv")))
}))

ref_rates_period <- ref_rates[is.na(as.numeric(year)) & !(
	diag_age %in% c("00", "All ages")), ]
ref_rates_period
# Deleting unnecessary rows
ref_rates <- ref_rates[!is.na(as.numeric(year)) & !(
	diag_age %in% c("00", "All ages")), ]

# Columns for easier age cutting
ref_rates[, `:=`(
	age.lower = as.numeric(ifelse(
		grepl("85", diag_age), 85, substr(diag_age, 1, unlist(gregexpr("-", diag_age)) - 1))),
	age.upper = as.numeric(ifelse(
		grepl("85", diag_age), "Inf",
		substr(diag_age, unlist(gregexpr("-", diag_age)) + 1, nchar(diag_age)))),
	year = as.integer(year)
)]
ref_rates_period[, `:=`(
	age.lower = as.numeric(ifelse(
		grepl("85", diag_age), 85, substr(diag_age, 1, unlist(gregexpr("-", diag_age)) - 1))),
	age.upper = as.numeric(ifelse(
		grepl("85", diag_age), "Inf",
		substr(diag_age, unlist(gregexpr("-", diag_age)) + 1, nchar(diag_age)))),
	year.lower = as.numeric(substr(year, 1, unlist(gregexpr("-", year)) - 1)),
	year.upper = as.numeric(substr(year, unlist(gregexpr("-", year)) + 1, nchar(year)))
)]
# make equivalent levels
ref_rates[, age := paste0(
	"[", age.lower,
	",", age.upper + 1, ")"
)]
ref_rates
ref_rates_period[, `:=`(
	age = paste0(
		"[", age.lower,
		",", age.upper + 1, ")")
)]
ref_rates_period

ref_rates <- ref_rates[age.lower >= 15]
ref_rates_period <- ref_rates_period[age.lower >= 15]

# Make long
# Calendar-averaged rates overlap...
year.lowerbounds <- c(sort(unique(ref_rates_period$year.lower)), 2017)
ref_rates_period[, `:=`(
	next.year.lower = sapply(year.lower, function(x) {
		year.lowerbounds[x < year.lowerbounds][1] - 1
	}))]
# Long-ify
ref_rates_period_long <- ref_rates_period[!is.na(rate), -c("year")][,.(
	year = seq(year.lower, next.year.lower),
	year.lower = rep(year.lower,
									 length(seq(year.lower, next.year.lower))),
	year.upper = rep(year.upper,
									 length(seq(year.lower, next.year.lower)))
), by = .(code, diag_age, rate, lower95, upper95, Race, Sex, age)]

# Fill NA rates
ref_rates <- rbindlist(list(
	ref_rates[!is.na(rate)],
	merge(ref_rates[is.na(rate), -c("rate", "lower95", "upper95")],
				ref_rates_period_long[!is.na(rate), -c("year.lower", "year.upper")],
				on = .(code, diag_age, rate, lower95, upper95, Race, Sex, age),
				all.x = T)
), use.names = T)

# Age cutoffs
age.cutpoints <- ref_rates[,.(lower = age.lower, upper = age.upper)]
age.cutpoints <- age.cutpoints[!duplicated(age.cutpoints)]

# SIR function ####
sir <- function(
	outcome = 1,
	cohort_name = "cohort_analytic",
	ref_name = "ref_rates",
	detail = T,
	significance.level = 0.95,
	verbose = T,
	return_output = F) {
	# Outomce of interest
	present_code <- incidence.key[outcome, code]
	description <- incidence.key[outcome, description]

	# Get analytic data
	obs <- as.data.table(as.data.frame(get(
		paste0(present_code, ".dat", ifelse(cohort_name == "cohort2", "2", ""))
	)))
	obs[, Race := factor(Race, levels = c("White", "Not white"),
											 labels = c("White", "Black"))]
	ref <- as.data.table(as.data.frame(get(ref_name)))

	# Gender exclusion, if needed
	if (grepl("Prosta|Male", description)) {
		obs <- obs[Sex == "M"]
		ref <- ref[Sex == "M"]
	}

	if (grepl("Breast|Fema", description)) {
		obs <- obs[Sex == "F"]
		ref <- ref[Sex == "F"]
	}

	# Calculate person-years in each row
	setorder(obs, studyno, year)
	obs[,`:=`(
		I = 1:.N,
		N = .N
	), by = .(studyno)]
	obs[I != N & year != year(yin) + 3, py := 1]
	obs[I == 1 & year == year(yin) + 3,
			py := time_length(difftime(
				as.Date(paste0(year + 1, "-01-01")),
				yin + years(3)
			), 'year')]
	obs[I == N, py := {
		time_length(difftime(
			as.Date(apply(data.frame(
				as.Date(paste0(end.year + 1, "-01-01")),
				yod + days(1),
				yoi + days(1),
				yoc + days(1)
			), 1, min, na.rm = T)),
			as.Date(apply(data.frame(
				paste0(year, "-01-01"),
				yin + years(3)
			), 1, max, na.rm = T))
		), 'year')
	}]

	# Cut age appropriately
	obs[, age := cut(age.year2 / 365, c(age.cutpoints$lower, Inf), right = F)]


	obs <- obs[,
						 .(
						 	count = sum(status),
						 	# Rate per 100,000
						 	rate = sum(status) / sum(py) * 100000,
						 	risk = sum(status) / length(table(studyno)) * 100000,
						 	py = sum(py),
						 	N = length(table(studyno))
						 ),
						 by = .(Sex,
						 			 Race,
						 			 year,
						 			 age)]

	exp <- ref[code == paste(present_code)]
	names(exp)[names(exp) == "rate"] <- "ref_rate"
	exp <- exp[obs, .(
		Sex,
		Race,
		year,
		age,
		ref_rate,
		obs.py = py,
		# incidence rate or cumulative incidence?
		exp.count = ref_rate * py / 100000,
		exp.count_byrisk = ref_rate * N / 100000,
		count = count,
		N
	),
	on = .(Race, Sex, age, year)]
	exp
	excluded <- exp[is.na(exp.count)]
	exp <- exp[!is.na(exp.count)]

	sir <- sum(exp$count) / sum(exp$exp.count)

	# Assuming A_i ~ Pois(A_i/(T_i^2))
	# Whether you assume population estimates makes a small
	# difference in the point estimates of the variance
	lsir.var <-	1 / sum(exp$count) + 1 / sum(exp$exp.count)
	lsir.se  <- sqrt(lsir.var)
	sir.lower2.5 <-
		exp(log(sir) + qnorm((1 - significance.level) / 2) * lsir.se)
	sir.upper2.5 <-
		exp(log(sir) - qnorm((1 - significance.level) / 2) * lsir.se)

	# smd <- sum(exp$count) - sum(exp$exp.count)
	#
	# smd.var <- sum(exp$count) / sum(exp$obs.py) ^ 2 +
	# 	sum(exp$exp.count) / sum(exp$obs.py) ^ 2
	#
	# smd.se  <- sqrt(smd.var)
	# smd.lower2.5 <-
	# 	smd + qnorm((1 - significance.level) / 2) * smd.se
	# smd.upper2.5 <-
	# 	smd - qnorm((1 - significance.level) / 2) * smd.se

	N <- sum(exp$count)

	if (verbose) {
		cat(paste0(rep("_", 85), collapse = ""))
		cat(paste0("\n", "SIR for '", tolower(description), "'",
							 " (n = ", N,")"))
		cat(paste0("\n", round(sir, 1), " (", round(sir.lower2.5, 1), ", ", round(sir.upper2.5, 1), ")"))
		if (sum(exp$count) == 0) {
			cat(paste0("\n", "No observed cases of '", description, "'"))
		}
		if (nrow(excluded) > 0) {
			cat(paste0("\n", "Reference rates not available for ", nrow(excluded), " strat",
								 ifelse(nrow(excluded) == 1, "um", "a"),
								 " (", sum(excluded$count), " cases, ",
								 round(sum(excluded$obs.py)), " person-years", ")"))
		}
		cat("\n")}

	output <- list(
		observed = obs,
		expected = exp,
		excluded = excluded,

		sir = sir,
		lsir.se = lsir.se,
		sir.lower2.5 = sir.lower2.5,
		sir.upper2.5 = sir.upper2.5,

		# smd = smd,
		# smd.se = smd.se,
		# smd.lower2.5 = smd.lower2.5,
		# smd.upper2.5 = smd.upper2.5,

		N = N
	)

	if (!detail) {
		output[-(1:3)]
	}
	if (return_output) {
		return(output)
	}
}

# Run SIR ####
# get.coxph(run_model = F, messy_sol = NA, time_scale = "age")
# sir.list <- invisible(lapply(cancers.which, sir, return_output = T))
#
# saveRDS(sir.list,
# 		 file = here::here("reports/cancer incidence/resources",
# 		 									"sir.rds"))
