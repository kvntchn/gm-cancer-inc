# Cancer incidence helper functions ####
# Kevin Chen
# February 28, 2020

library(here)
# here::set_here()

# rm(list = ls(envir = .GlobalEnv))
# rm(list = ls(envir = .GlobalEnv)[-grep("outcome.selected", ls(envir = .GlobalEnv))])
# rm(list = ls(envir = .GlobalEnv)[-grep('cohort', ls(envir = .GlobalEnv))])

library(survival)

if (!("exposure.lag" %in% ls(envir = .GlobalEnv))) {
	exposure.lag <- 21
	# rm(cohort_analytic)
}

if (!("get.cohort2" %in% ls(envir = .GlobalEnv))) {
	get.cohort2 <- T
	# rm(cohort_analytic)
}

# # Posterior probabilities for mi race ####
# if (!"race.imputed.melt" %in% ls(envir = .GlobalEnv)) {
# 	race.imputed.melt <- box_read(ifelse(exposure.lag == 21, 780327112319, NaN))
# }

get.mi_race <- function(dat, M, posterior_probs, threshold = 0.5) {
	lapply(1:M, function(m) {

		b <- sample(posterior_probs$variable, 1)

		to_impute <- posterior_probs[variable == b, .(
			studyno,
			year,
			Race.impute = factor(as.numeric(plogis(value) >= threshold),
													 levels = c(1, 0),
													 labels = c("Black", "White")),
			b = b)]

		return(to_impute)
	}) # end lapply
}

# Get Cox PH ####
get.coxph <- function(
	outcomes = outcomes.which, #c(32, 33),
	cohort_name = NULL,
	run_model = F,
	save_dat = T,
	messy_sol = 0.05,
	spline_year = F,
	year.df = 0,
	spline_yin = F,
	yin.df = 0,
	time_scale = "age",
	hwse2 = F,
	hwse3 = F,
	start.year = NULL,
	additional.lag = 0,
	directory.name = NULL,
	employment_status.lag = 0,
	year.max = 2015,
	start_m = 1,
	mi = 0,
	age_under = Inf) {

	# start time
	start <- Sys.time()

	options(warn = 2)

	# invisible(sapply(outcomes, function(i = outcomes[1]) # sapply over outcomes
	for (i in outcomes) {
		# Get data ####

		# Two-to-three letter code indicating cancer type
		code <- unlist(incidence.key[i, 1])
		# Cancer type
		description <- unlist(incidence.key[i, 2])

		var.name <- unlist(incidence.key[i, 3])
		date.name <- unlist(incidence.key[i, 4])

		# If looking at employment status, either as an outcome or as an exposure
		if (is.null(cohort_name)) {
			if (hwse3 | hwse2) {
				cohort_name <- "cohort2"
				if (i %in% which(grepl("copd|external", incidence.key$code))) {
					cohort_name <- "mortality.cohort2"
				}
			} else {
				cohort_name <- "cohort_analytic"
			}}

		dat <- data.table::copy(get(cohort_name))

		# Additional lag ####
		setorder(dat, studyno, year)
		dat[,`:=`(I = 1:.N,
							N = .N), by = .(studyno)]
		if (additional.lag > 0) {
			dat[,`:=`(
				straight = shift(straight, additional.lag, fill = 0),
				cum_straight = shift(cum_straight, additional.lag, fill = 0),
				soluble = shift(soluble, additional.lag, fill = 0),
				cum_soluble = shift(cum_soluble, additional.lag, fill = 0),
				synthetic = shift(synthetic, additional.lag, fill = 0),
				cum_synthetic = shift(cum_synthetic, additional.lag, fill = 0)
				# cum_bio = shift(cum_bio, additional.lag, fill = 0),
				# cum_cl = shift(cum_cl, additional.lag, fill = 0),
				# cum_ea = shift(cum_ea, additional.lag, fill = 0),
				# cum_tea = shift(cum_tea, additional.lag, fill = 0),
				# cum_trz = shift(cum_trz, additional.lag, fill = 0),
				# cum_s = shift(cum_s, additional.lag, fill = 0),
				# cum_no2 = shift(cum_no2, additional.lag, fill = 0),
				# bio = shift(bio, additional.lag, fill = 0),
				# cl = shift(cl, additional.lag, fill = 0),
				# ea = shift(ea, additional.lag, fill = 0),
				# tea = shift(tea, additional.lag, fill = 0),
				# trz = shift(trz, additional.lag, fill = 0),
				# s = shift(s, additional.lag, fill = 0),
				# no2 = shift(no2, additional.lag, fill = 0)
			), by = .(studyno)]
		}
		if (additional.lag < 0) {
			dat[,`:=`(
				straight = shift(straight, additional.lag, fill = straight[I == N]),
				soluble = shift(soluble, additional.lag, fill = soluble[I == N]),
				synthetic = shift(synthetic, additional.lag, fill = synthetic[I == N])
				# bio = shift(bio, additional.lag, fill = bio[I == N]),
				# cl = shift(cl, additional.lag, fill = cl[I == N]),
				# ea = shift(ea, additional.lag, fill = ea[I == N]),
				# tea = shift(tea, additional.lag, fill = tea[I == N]),
				# trz = shift(trz, additional.lag, fill = trz[I == N]),
				# s = shift(s, additional.lag, fill = s[I == N]),
				# no2 = shift(no2, additional.lag, fill = no2[I == N])
			), by = .(studyno)]
			# Exposure after leaving work is 0
			dat[year > (year(jobloss.date) + exposure.lag + additional.lag), `:=`(
				straight = 0,
				soluble = 0,
				synthetic = 0
				# bio = 0,
				# cl = 0,
				# ea = 0,
				# tea = 0,
				# trz = 0,
				# s = 0,
				# no2 = 0
			)]
			dat[,`:=`(
				cum_straight = cumsum(straight),
				cum_soluble = cumsum(soluble),
				cum_synthetic = cumsum(synthetic)
				# cum_bio = cumsum(bio),
				# cum_cl = cumsum(cl),
				# cum_ea = cumsum(ea),
				# cum_tea = cumsum(tea),
				# cum_trz = cumsum(trz),
				# cum_s = cumsum(s),
				# cum_no2 = cumsum(no2),
			), by = .(studyno)]
		}

		# # New category: Alkanolamines
		# dat[,`:=`(
		# 	ohnh = ea + tea,
		# 	cum_ohnh = cum_ea + cum_tea
		# )]


		# When to start FU ####
		# HWSE 3 is for the MWF-leaving work association: no bounds to FU
		# If looking at mortality outcome, no need to start FU later
		if (is.null(start.year)) {
			if (#F # No more FU starting in 1941
				hwse3 | i %in% which(grepl("copd|external", incidence.key$code))
			) {
				dat$start.year <- 1941} else if (
					# If looking at cancer incidence, should use 1973 or 1985
					min(dat$ddiag_first, na.rm = T) < as.Date("1985-01-01")) {
					dat$start.year <- 1973
					dat[plant == 3, start.year := 1985]
				} else {
					dat$start.year <- 1985
				}
		}

		dat <- dat[year >= start.year]

		if (!hwse3 & nrow(dat[
			plant == 3 &
			!i %in% which(grepl("copd|external", incidence.key$code)) &
			year < 1985]) > 0) {
			warning("Cancer incidence FU for plant 3 before 1985?")}

		# Lag leaving work ####
		dat[,`:=`(
			employment_status_lag = employment_status.lag,
			unlagged_jobloss.date = jobloss.date,
			unlagged_age.leavework = time_length(difftime(jobloss.date, yob), 'year'))]

		# Censor when (lagged) employment status not known ####
		if (hwse2 | hwse3) {
			if (employment_status.lag > 0) {
				dat[,`:=`(jobloss.date = jobloss.date[1] + years(employment_status.lag)
				), by = .(studyno)]
				dat <- dat[jobloss.date <= as.Date(paste0(
					1994 + employment_status.lag, "-12-31")) |
						year <= 1994 + employment_status.lag]
			} else {
				dat <- dat[jobloss.date <= as.Date("1994-12-31") | year <= 1994 ]
				dat[,`:=`(employment_status_lag = 0)]
			}}

		if (!hwse3) {
			dat[,`:=`(status = get(var.name),
								yoi = get(date.name))]
		} else {
			code <- ""
			dat[year < year(jobloss.date), status := 0]
			dat[year == year(jobloss.date), status := 1]
			dat[year > year(jobloss.date), status := 2]
			dat$yoi <- dat$jobloss.date
		}

		# If incidence is after date of death, make it the date of death ####
		dat[yoi > yod, yoi := yod]

		# if (hwse2) {
		# 	# If job loss in the same year as incidence, and yoi was arbitrary ####
		# 	dat[year(jobloss.date) == year(yoi) & month(yoi) == 7 & day(yoi) == 1,
		# 			yoi := jobloss.date]
		# }

		# Drop unnecessary data ####
		dat <- dat[, -grep(
			"canc_|ddiag|racon", names(dat)[
				-which(names(dat) %in% c(
					"ddiag_first", "canc_first", var.name))],
			value = T), with = F]

		dat <- dat[year <= as.integer(apply(data.frame(
			year(yoi),
			year(yod),
			year(yoc),
			year.max), 1, min, na.rm = T)),]

		# If experienced cancer after censor date, make sure status is 0
		dat[yoi >= yoc, status := 0]
		# Check age and year index at last row
		dat[year == year(yoc) & yoi >= yoc, `:=`(
			age.year2 = time_length(difftime(
				as.Date(apply(
					data.frame(as.Date(paste0(year.max + 1, "-01-01")),
										 yod + days(1),
										 yoc + days(1)), 1, min, na.rm = T
				)), yob), 'day'),
			year2 = as.Date(apply(data.frame(
				as.Date(paste0(year, "-12-31")),
				as.Date(paste0(year.max, "-12-31")),
				yod,
				yoc), 1, min, na.rm = T
			))
		)]

		# Order data
		setorder(dat, studyno, year)
		dat[, `:=`(I = 1:.N,
							 N = .N), by = .(studyno)]

		# Fix age.year1 for those hired before start.year - 3
		dat[year(yin) + 3 < start.year & I == 1, `:=`(
			age.year1 = time_length(
				difftime(as.Date(paste0(start.year, "-01-01")), yob), 'days'),
			year1 = as.Date(paste0(start.year - 1, "-12-31"))
		)]

		# Must be at-risk for incidence (or jobloss)
		dat <- dat[yoi >= yin + years(3) | is.na(yoi)]
		dat <- dat[status != 2]

		# Fix age.year2 at time of incidence (or jobloss)
		dat[status == 1, `:=`(age.year2 = floor(time_length(difftime(
			as.Date(apply(
				data.frame(as.Date(paste0(
					year.max + 1, "-01-01"
				)),
				yod + days(1),
				yoc + days(1),
				yoi + days(1)),
				1,
				min,
				na.rm = T
			)),
			yob
		), 'day')),
		year2 = as.Date(apply(
			data.frame(as.Date(paste0(
				year.max, "-12-31"
			)),
			yod,
			yoc,
			yoi),
			1,
			min,
			na.rm = T
		)))]


		# Don't confuse yoi with jobloss.date in the case of HWSE 3
		if (hwse3) {
			dat$yoi <- NA
		}

		# Duration of employment ####
		dat[, `:=`(
			employment.years = time_length(difftime(as.Date(
				apply(data.frame(
					as.Date(paste0(year + 1, "-01-01")),
					as.Date("1995-01-01")
				), 1, min, na.rm = T)
			),
			yin), 'year'))]

		# Age at leaving work (before employment status lag)
		dat[,`:=`(age.leavework = time_length(difftime(jobloss.date, yob), 'year'))]

		# Define quantiles ####
		covariate.breaks <- apply(dat[status == 1, .(
			year,
			yin = yin.gm,
			employment.years,
			age = age.year2 / 365)], 2, function(x) {
				if (length(x) > 100) {
					breaks <- quantile(x, seq(0, 1, 1 / 6))
				} else if (length(x) > 80) {
					breaks <- quantile(x, seq(0, 1, 1 / 5))
				} else if (length(x) > 60) {
					breaks <- quantile(x, seq(0, 1, 1 / 4))
				} else if (length(x) > 40) {
					breaks <- quantile(x, seq(0, 1, 1 / 3))
				} else {
					breaks <- quantile(x, seq(0, 1, 1 / 2))
				}
				breaks[-length(breaks)] <- floor(breaks[-length(breaks)])
				breaks[length(breaks)] <- ceiling(breaks[length(breaks)])
				breaks[c(1, length(breaks))] <- c(-Inf, Inf)
				breaks
			})

		covariate.breaks <- as.data.table(covariate.breaks)

		# # Pool last two calendar year levels
		# covariate.breaks[is.finite(year), year := {
		# 	if (length(year) > 1) {
		# 		c(year[-length(year)], NA)
		# 	} else {
		# 		year
		# 	}
		# }]
		# # Restrict lower bound of last Year category?
		# covariate.breaks[year > 2010 & is.finite(year), year := NA]

		mwf.breaks <- cbind(apply(dat[status == 1, .(
			straight,
			soluble,
			synthetic,
			off,
			cum_straight,
			# cum_soluble,
			cum_synthetic,
			cum_off)], 2, function(x) {
				x <- x[x > 0]
				if (length(x) > 40) {
					if (length(x) > 60) {
						probs <- seq(0, 1, 1 / 3)
					} else {
						probs <- seq(0, 1, 1 / 2)
					}
					breaks <- quantile(x, probs)
					breaks[c(1, length(probs))] <- c(0, Inf)
					breaks <- c(-Inf, breaks)
				} else {
					breaks <- c(-Inf, 0, Inf)
				}
				if (length(breaks) < 5) {
					breaks <- c(breaks, rep(NA, 5 - length(breaks)))
				}
				names(breaks) <- NULL
				breaks
			}),
			apply(dat[status == 1, .(cum_soluble)], 2, function(x) {
				x <- x[x > 0]
				if (length(x) > 40) {
					if (length(x) > 60) {
						probs <- seq(0, 1, 1 / 3)
					} else {
						probs <- seq(0, 1, 1 / 2)
					}
					breaks <- quantile(x[x > 0], probs)
					breaks[c(1, length(probs))] <- c(0, Inf)
					breaks <- c(-Inf, breaks)
				} else {
					breaks <- c(-Inf, 0, Inf)
				}
				if (length(breaks) < 5) {
					breaks <- c(breaks, rep(NA, 5 - length(breaks)))
				}
				breaks
			}),
			apply(dat[status == 1, .(cum_soluble5 = cum_soluble,
															 soluble5 = soluble)], 2, function(x) {
															 	x <- x[x > 0.05]
															 	if (length(x) > 40) {
															 		if (length(x) > 60) {
															 			probs <- seq(0, 1, 1 / 3)
															 		} else {
															 			probs <- seq(0, 1, 1 / 2)
															 		}
															 		breaks <- quantile(x[x > 0], probs)
															 		breaks[c(1, length(probs))] <- c(0.05, Inf)
															 		breaks <- c(-Inf, breaks)
															 	} else {
															 		breaks <- c(-Inf, 0.05, Inf)
															 	}
															 	if (length(breaks) < 5) {
															 		breaks <- c(breaks, rep(NA, 5 - length(breaks)))
															 	}
															 	breaks
															 }),
			apply(dat[status == 1, .(cum_soluble127 = cum_soluble,
															 soluble127 = soluble)], 2, function(x) {
															 	x <- x[x > 1.2742]
															 	if (length(x) > 40) {
															 		if (length(x) > 60) {
															 			probs <- seq(0, 1, 1 / 3)
															 		} else {
															 			probs <- seq(0, 1, 1 / 2)
															 		}
															 		breaks <- quantile(x[x > 0], probs)
															 		breaks[c(1, length(probs))] <- c(1.2742, Inf)
															 		breaks <- c(-Inf, breaks)
															 	} else {
															 		breaks <- c(-Inf, 1.2742, Inf)
															 	}
															 	if (length(breaks) < 5) {
															 		breaks <- c(breaks, rep(NA, 5 - length(breaks)))
															 	}
															 	breaks
															 }),
			apply(dat[status == 1, .(cum_soluble11 = cum_soluble,
															 soluble11 = soluble)], 2, function(x) {
															 	x <- x[x > 0.11]
															 	if (length(x) > 40) {
															 		if (length(x) > 60) {
															 			probs <- seq(0, 1, 1 / 3)
															 		} else {
															 			probs <- seq(0, 1, 1 / 2)
															 		}
															 		breaks <- quantile(x[x > 0], probs)
															 		breaks[c(1, length(probs))] <- c(0.11, Inf)
															 		breaks <- c(-Inf, breaks)
															 	} else {
															 		breaks <- c(-Inf, 0.11, Inf)
															 	}
															 	if (length(breaks) < 5) {
															 		breaks <- c(breaks, rep(NA, 5 - length(breaks)))
															 	}
															 	breaks
															 })
		)

		mwf.breaks <- as.data.table(mwf.breaks)

		# component.breaks <- apply(dat[status == 1, .(
		# 	bio,
		# 	cl,
		# 	ea,
		# 	tea,
		# 	ohnh,
		# 	trz,
		# 	s,
		# 	no2,
		# 	cum_bio,
		# 	cum_cl,
		# 	cum_ea,
		# 	cum_tea,
		# 	cum_ohnh,
		# 	cum_trz,
		# 	cum_s,
		# 	cum_no2)], 2, function(x) {
		# 		x <- x[x > 0]
		# 		if (length(x) > 40) {
		# 			if (length(x) > 60) {
		# 				probs <- seq(0, 1, 1 / 3)
		# 			} else {
		# 				probs <- seq(0, 1, 1 / 2)
		# 			}
		# 			breaks <- quantile(x, probs)
		# 			breaks[c(1, length(probs))] <- c(0, Inf)
		# 			breaks <- c(-Inf, breaks)
		# 		} else {
		# 			breaks <- c(-Inf, 0, Inf)
		# 		}
		# 		if (length(breaks) < 5) {
		# 			breaks <- c(breaks, rep(NA, 5 - length(breaks)))
		# 		}
		# 		names(breaks) <- NULL
		# 		breaks
		# 	})
		#
		# component.breaks <- as.data.table(component.breaks)

		# Make categorical variables ####
		get.cut <- function(x, breaks, y = NULL, include.lowest = T, dig.lab = 4) {
			if (is.null(y)) {y <- deparse(substitute(x))}
			cut(
				x,
				unlist(unique(na.exclude(breaks[, y, with = F]))),
				include.lowest = include.lowest,
				dig.lab = dig.lab
			)
		}
		dat[, `:=`(
			Year = get.cut(year, covariate.breaks),
			`Year of hire` = get.cut(yin.gm, covariate.breaks, "yin"),
			`Duration of employment` = get.cut(employment.years, covariate.breaks),
			Age = get.cut(age.year2 / 365, covariate.breaks, "age"),
			Straight = get.cut(straight, mwf.breaks, dig.lab = 3),
			Soluble = get.cut(soluble, mwf.breaks, "soluble5", dig.lab = 3),
			Synthetic =  get.cut(synthetic, mwf.breaks, dig.lab = 3),
			`Time off` =  get.cut(off, mwf.breaks, dig.lab = 3),
			`Cumulative straight` = get.cut(cum_straight, mwf.breaks, dig.lab = 3),
			`Cumulative soluble` = get.cut(cum_soluble, mwf.breaks, dig.lab = 3),
			`Cumulative soluble 5` = get.cut(cum_soluble, mwf.breaks, "cum_soluble5", dig.lab = 3),
			`Cumulative soluble 127` = get.cut(cum_soluble, mwf.breaks, "cum_soluble127", dig.lab = 3),
			`Cumulative soluble 10` = get.cut(cum_soluble, mwf.breaks, "cum_soluble11", dig.lab = 3),
			`Cumulative synthetic` = get.cut(cum_synthetic, mwf.breaks, dig.lab = 3),
			`Cumulative time off` =  get.cut(cum_off, mwf.breaks, dig.lab = 3)
			# Biocide = get.cut(bio, component.breaks, dig.lab = 3),
			# Chlorine = get.cut(cl, component.breaks, dig.lab = 3),
			# Ethanolamine = get.cut(ea, component.breaks, dig.lab = 3),
			# Triethanolamine = get.cut(tea, component.breaks, dig.lab = 3),
			# Alkanolamine = get.cut(ohnh, component.breaks, dig.lab = 3),
			# Triazine = get.cut(trz, component.breaks, dig.lab = 3),
			# Sulfur = get.cut(s, component.breaks, dig.lab = 3),
			# Nitrites = get.cut(no2, component.breaks, dig.lab = 3),
			# `Cumulative biocide` = get.cut(cum_bio, component.breaks, dig.lab = 3),
			# `Cumulative chlorine` = get.cut(cum_cl, component.breaks, dig.lab = 3),
			# `Cumulative ethanolamine` = get.cut(cum_ea, component.breaks, dig.lab = 3),
			# `Cumulative triethanolamine` = get.cut(cum_tea, component.breaks, dig.lab = 3),
			# `Cumulative alkanolamine` = get.cut(cum_ohnh, component.breaks, dig.lab = 3),
			# `Cumulative triazine` = get.cut(cum_trz, component.breaks, dig.lab = 3),
			# `Cumulative sulfur` = get.cut(cum_s, component.breaks, dig.lab = 3),
			# `Cumulative nitrosamine` = get.cut(cum_no2, component.breaks, dig.lab = 3)
		)]

		dat[,`:=`(
			Race = {
				race <- relevel(factor(race), "White")
				levels(race)[grep("Not white", levels(race))] <- "Black"
				race},
			Sex = factor(sex, levels = c("M", "F"), labels = c("Male", "Female")),
			Plant = relevel(factor(plant), "1")
		)]

		dat[, `:=`(
			year1.date = year1,
			year2.date = year2,
			year1 = as.numeric(year1, origin = as.Date("1970-01-01")),
			year2 = as.numeric(year2, origin = as.Date("1970-01-01"))
		)]

		# dat <- dat[studyno %in% sample(studyno, 5000)]
		setorder(dat, studyno, age.year2)
		# Save data ####
		if (save_dat) {
			assign(
				paste0(code, ifelse(!hwse3, ".", ""), "dat", ifelse(hwse2, 2, ifelse(hwse3, 3, ""))),
				dat,
				inherits = T,
				envir = .GlobalEnv)
			Sys.sleep(0)
		}

		# Restrict age ####
		if (is.finite(age_under)) {
			setorder(dat, studyno, age.year2)
			dat[,I := 1:.N, by = .(studyno)]
			include.who <- dat[age.year1/365 < age_under & I == 1, studyno]
			dat <- dat[studyno %in% include.who]
		}

		# Get data to impute
		if (mi > 0) {
			if (!paste0(gsub(" ", "_", code), ".mi_race.M", mi) %in% ls(envir = .GlobalEnv)) {
				to_impute <- get.mi_race(dat, mi, race.imputed.melt)
				assign(paste0(
					gsub(" ", "_", code),
					".mi_race.M", mi),
					to_impute, inherits = T)
			} else {
				to_impute <- get(paste0(gsub(" ", "_", code),".mi_race.M", mi))
			}
		}

		# Multiple imputation loop ####
		# Loop over mi ####
		lapply(if (mi == 0) {0} else {start_m:mi}, function(m = 1) {

			if (mi > 0) {

				to_impute[[m]] <- to_impute[[m]][,.(
					studyno, year, Race.impute, finrace = 9)]

				table(dat$Race, useNA = "always")
				nrow(dat)
				dat[finrace %in% c(0, 9), finrace := 9]
				dat$Race.impute <- NA
				dat <- merge(
					dat[,-"Race.impute"], to_impute[[m]],
					by = c("studyno", "year", "finrace"), all.x = T)
				nrow(dat)
				nrow(to_impute[[m]])
				dat[
					finrace %in% c(0, 9),
					Race := Race.impute]
				table(dat[immortal != 1 & right.censored != 1]$Race, useNA = "always")
			}

			# Run model ####
			if (run_model) {

				dat <- dat[immortal != 1 & right.censored != 1]
				Sys.sleep(0)
				tmp.coxph <- coxph(as.formula(
					paste(
						ifelse(
							grepl("age", time_scale),
							"Surv(age.year1, age.year2, status) ~",
							"Surv(year1, year2, status) ~"
						),
						"`Cumulative straight` +",
						ifelse(
							!is.finite(messy_sol),
							"`Cumulative soluble` +",
							paste0(
								"`Cumulative soluble",
								paste0(" ", messy_sol %/% 0.01),
								"` +"
							)
						),
						"`Cumulative synthetic` +",
						ifelse(
							grepl("age", time_scale),
							ifelse(!spline_year, paste0("Year +", "pspline(year, df = ", year.df, ") +")),
							"Age +"
						),
						ifelse(!spline_yin, "`Year of hire` + ", paste0("pspline(yin.gm, df = ", yin.df, ") +")),
						"Race + Plant",
						ifelse(!grepl('breast|male|prostate', description, ignore.case = T), "+ Sex", "")
					)
				),
				data = {
					if (grepl('breast|male|prostate', description, ignore.case = T)) {
						dat[Sex == ifelse(grepl('breast|female', description, ignore.case = T), "Female", "Male")]
					} else {
						dat
					}
				},
				method = "efron")

				# Save model  ####
				if (is.null(directory.name)) {
					directory.name <- to_drive_D(here::here(
						paste0("resources/Lag ",
									 exposure.lag + additional.lag,
									 ifelse(!grepl("age", time_scale),
									 			 "/indexed by calendar",
									 			 "/indexed by age"),
									 ifelse(mi > 0, "/mi race", ""),
						)
					))
				}

				if (mi > 0) {
					directory.name <- paste0(directory.name, "/", code)
				}

				dir.create(directory.name, showWarnings = F, recursive = T)

				final_name <- paste0(
					directory.name, "/",
					code,
					ifelse(
						is.finite(messy_sol),
						paste0("_sol", messy_sol %/% .01),
						""),
					ifelse(spline_year,
								 ifelse(spline_yin, paste0("_splined"), "_splinedyear"),
								 ifelse(spline_yin, paste0("_splinedyin"), "")),
					ifelse(mi > 0, paste0(".m", m), ""),
					".coxph.rds"
				)

				message(paste0("\nSaving model to:\n",
											 final_name))
				saveRDS(tmp.coxph,
								file = final_name)

				# if (nrow(dat[age.year2 <= age.year1]) > 0) {
				# 	nrow(dat[age.year2 <= age.year1])
				# }

				# Print to console ####
				cat(paste0("\n", paste0(rep("_", 80), collapse = ""), "\n"))
				cat(unlist(description))
				cat("\n")
				if (mi > 0) {cat(paste0("m = ", m, "\n"))}
				cat(paste0("Exposure lagged ", exposure.lag + additional.lag, " years\n"))
				cat(paste0(tmp.coxph$nevent, " cases"))
				# Time elapsed
				since_start <- time_length(difftime(Sys.time(), start), "minutes")
				if (since_start <= 90) {
					cat(paste0("\n", round(since_start, 2), " minutes since get.coxph() was called."))
				} else {
					cat(paste0("\n", round(since_start/60, 3), " hours since get.coxph() was called."))
				}
				cat(paste0("\n", paste0(rep("_", 80), collapse = ""), "\n"))
				print(summary(tmp.coxph)$coefficients)

			}
		}) # End MI loop
	}
	# )) # end sapply over outcomes
	options(warn = 0)
}

# # Try out get.coxph() ####
# get.coxph(outcomes = c(32:33))
#
# first.dat[,.(
# 	final.year = as.Date(year2[.N], as.Date("1970-01-01")),
# 	status = status[.N],
# 	plant = plant[.N],
# 	yoi = yoi[.N],
# 	yod = yod[.N],
# 	year = year[.N],
# 	yoc = yoc[.N]
# ), by = .(studyno)][,.(plant, year)] %>% table
# first.dat[,.(
# 	final.year = as.Date(year2[.N], as.Date("1970-01-01")),
# 	status = status[.N],
# 	plant = plant[.N],
# 	yoi = yoi[.N],
# 	yod = yod[.N],
# 	year = year[.N],
# 	yoc = yoc[.N]
# ), by = .(studyno)][final.year != yoi & !is.na(yoi)]
# sum(first.dat$status)
# sum(first.dat[status == 1 & yoi < as.Date("2016-01-01"), 1, by = .(studyno)]$V1)
#
#
# copd.dat[,.(
# 	final.year = as.Date(year2[.N], as.Date("1970-01-01")),
# 	status = status[.N],
# 	plant = plant[.N],
# 	yoi = yoi[.N],
# 	yod = yod[.N],
# 	year = year[.N],
# 	yoc = yoc[.N]
# ), by = .(studyno)][,.(plant, year)] %>% table
# copd.dat[yoi != yod]
# sum(copd.dat$status)
# sum(copd.dat[status == 1 & yoi < as.Date("2016-01-01"), 1, by = .(studyno)]$V1)

# HWSE 2 - Leaving work and outcome ####
get.hwse2.coxph <- function(
	outcomes =  outcomes.which,
	cohort_name = NULL,
	run_model = F,
	new_dat = T,
	save_dat = T,
	messy_sol = 0.05,
	spline_year = F,
	year.df = 0,
	spline_yin = F,
	yin.df = 0,
	time_scale = "age",
	additional.lag = 0,
	employment_status.lag = 0,
	year.max = Inf,
	directory.name = NULL,
	start_m = 1,
	mi = 0,
	age_under = Inf) {

	# start time
	start <- Sys.time()

	options(warn = 2)

	# invisible(sapply(outcomes, function(i = outcomes.which[1])
	for (i in outcomes) {

		code <- unlist(incidence.key[i, 1])
		description <- unlist(incidence.key[i, 2])
		var.name <- unlist(incidence.key[i, 3])

		if (!(paste0(code, ".dat2") %in% ls(envir = .GlobalEnv)) | new_dat) {

			get.coxph(
				cohort_name = cohort_name,
				run_model = F,
				outcomes = i,
				hwse2 = T,
				additional.lag = additional.lag,
				employment_status.lag = employment_status.lag,
				year.max = year.max
			)

			dat <- data.table::copy(get(paste0(code, ".dat2"), envir = .GlobalEnv))

			# Special attention to (lagged) leaving work ####
			dat[year >= year(jobloss.date), `Employment status` := 1]
			dat[year < year(jobloss.date) | yoi <= jobloss.date, `Employment status` := 0]

			{
				# Years other than year of leaving work
				dat1 <- dat[year != year(jobloss.date) |
											year(jobloss.date) >= (1995 + employment_status.lag) |
											yoc <= jobloss.date |
											yod <= jobloss.date |
											yoi <= jobloss.date]
				dat1[year == year(jobloss.date) & (
					year(jobloss.date) >= (1995 + employment_status.lag) |
						yoc <= jobloss.date |
						yod <= jobloss.date |
						yoi <= jobloss.date), `:=`(
							age.year2 = time_length(difftime(as.Date(
								apply(
									data.frame(as.Date(paste0(
										year + 1, "-01-01"
									)),
									yoc + days(1),
									yod + days(1),
									yoi + days(1)),
									1,
									min,
									na.rm = T
								)
							), yob), 'day'),
							year2 = as.Date(apply(
								data.frame(as.Date(paste0(
									year, "-12-31"
								)),
								yoc, yod, yoi),
								1,
								min,
								na.rm = T
							)),
							status = ifelse(status == 1 &
																(yoi <= yoc | is.na(yoc)) &
																(yoi <= yod | is.na(yod)) &
																(yoi <= yoi | is.na(yoi)), 1, 0),
							`Employment status` = 0
						)]

				# Year of leaving work
				dat2 <- dat[year == year(jobloss.date) &
											year(jobloss.date) < (1995  + employment_status.lag) &
											(yoc > jobloss.date | is.na(yoc)) &
											(yod > jobloss.date | is.na(yod)) &
											(yoi > jobloss.date | is.na(yoi))]
				dat2[, `:=`(
					age.year2 = time_length(difftime(as.Date(
						apply(
							data.frame(as.Date(paste0(
								year + 1, "-01-01"
							)),
							jobloss.date + days(1)),
							1,
							min,
							na.rm = T
						)
					), yob), 'day'),
					year2 = as.Date(apply(
						data.frame(as.Date(paste0(
							year, "-12-31"
						)),
						jobloss.date),
						1,
						min,
						na.rm = T
					)),
					status = ifelse(status == 1 & yoi <= jobloss.date, 1, 0),
					`Employment status` = 0
				)]
				# Year of leaving jobloss end of year
				dat3 <-	dat[year == year(jobloss.date) &
											year(jobloss.date) < (1995  + employment_status.lag) &
											(yoc > jobloss.date | is.na(yoc)) &
											(yod > jobloss.date | is.na(yod)) &
											(yoi > jobloss.date | is.na(yoi)) &
											!(month(jobloss.date) == ifelse(grepl("age", time_scale), 12, 1) &
													day(jobloss.date) == ifelse(grepl("age", time_scale), 31, 1)) & (
														jobloss.date < yoi | is.na(yoi))]
				dat3[, `:=`(
					age.year1 = floor(time_length(difftime(
						as.Date(apply(
							data.frame(as.Date(paste0(
								year + 1, "-01-01"
							)),
							jobloss.date + days(1)),
							1,
							min,
							na.rm = T
						)), yob
					), 'day')),
					year1 = as.Date(apply(
						data.frame(as.Date(paste0(
							year, "-12-31"
						)),
						jobloss.date),
						1,
						min,
						na.rm = T
					)),
					status = ifelse(status == 1 & yoi > jobloss.date, 1, 0),
					`Employment status` = 1
				)]

				dat1[, `:=`(
					year1 = as.numeric(year1, origin = as.Date("1970-01-01")),
					year2 = as.numeric(year2, origin = as.Date("1970-01-01"))
				)]
				dat2[, `:=`(
					year1 = as.numeric(year1, origin = as.Date("1970-01-01")),
					year2 = as.numeric(year2, origin = as.Date("1970-01-01"))
				)]
				dat3[, `:=`(
					year1 = as.numeric(year1, origin = as.Date("1970-01-01")),
					year2 = as.numeric(year2, origin = as.Date("1970-01-01"))
				)]

				dat <- rbindlist(list(dat1, dat2, dat3))
				dat[, `:=`(
					year1.date = as.Date(year1, origin = as.Date("1970-01-01")),
					year2.date = as.Date(year2, origin = as.Date("1970-01-01"))
				)]
				# If death/incidence date is date of leaving work, make it off the job ####
				dat[yoi == jobloss.date & year == year(jobloss.date), `Employment status` := 1]
				dat[year(jobloss.date) >= 1995 + employment_status.lag | jobloss.date > yoc, `Employment status` := 0]

				dat[year > year(jobloss.date),
						`Employment status` := 1]

				setorder(dat, studyno, age.year2)
			}

			# Indicator for age at leaving work ####
			dat[,`:=`(
				age.leavework = time_length(difftime(jobloss.date, yob), 'year')
			)]
			dat[, `:=`(
				`Employment status` = factor(
					`Employment status`, 0:1, c(
						ifelse(employment_status.lag == 0, "Still at work",
									 paste0("Not yet ", employment_status.lag, " years since leaving work")),
						ifelse(employment_status.lag == 0, "Left work (any age)",
									 paste0("Left work at least ", employment_status.lag, " years ago"))
					)),
				`Employment status 50` = {
					tmp <- ifelse(age.leavework < 51,	2, 3)
					tmp[`Employment status` == 0] <- 1
					factor(tmp, levels = 1:3,
								 labels = c(
								 	ifelse(employment_status.lag == 0, "Still at work",
								 				 paste0("Not yet ", employment_status.lag, " years since leaving work")),
								 	ifelse(employment_status.lag == 0, "Left work (Age 50 or younger)",
								 				 paste0("Left work at least ", employment_status.lag, " years ago (Age 50 or younger)")),
								 	ifelse(employment_status.lag == 0, "Left work (Age 51 or older)",
								 				 paste0("Left work at least ", employment_status.lag, " years ago (Age 51 or older)"))))},
				`Employment status 55` = {
					tmp <- ifelse(age.leavework < 56,	2, 3)
					tmp[`Employment status` == 0] <- 1
					factor(tmp, levels = 1:3,
								 labels = c(
								 	ifelse(employment_status.lag == 0, "Still at work",
								 				 paste0("Not yet ", employment_status.lag, " years since leaving work")),
								 	ifelse(employment_status.lag == 0, "Left work (Age 55 or younger)",
								 				 paste0("Left work at least ", employment_status.lag, " years ago (Age 55 or younger)")),
								 	ifelse(employment_status.lag == 0, "Left work (Age 56 or older)",
								 				 paste0("Left work at least ", employment_status.lag, " years ago (Age 56 or older)"))
								 ))},
				`Employment status 60` = {
					tmp <- ifelse(age.leavework < 61,	2, 3)
					tmp[`Employment status` == 0] <- 1
					factor(tmp, levels = 1:3,
								 labels = c(
								 	ifelse(employment_status.lag == 0, "Still at work",
								 				 paste0("Not yet ", employment_status.lag, " years since leaving work")),
								 	ifelse(employment_status.lag == 0, "Left work (Age 60 or younger)",
								 				 paste0("Left work at least ", employment_status.lag, " years ago (Age 60 or younger)")),
								 	ifelse(employment_status.lag == 0, "Left work (Age 61 or older)",
								 				 paste0("Left work at least ", employment_status.lag, " years ago (Age 61 or older)"))
								 ))}
			)]

			setorder(dat, studyno, year, year2)

			# dat <- dat[!(year == year(jobloss.date) &
			# 						 	yoi < jobloss.date &
			# 						 	`Employment status` == 1),]

			# Check time 2 > time 1
			# dat[year2 <= year1, .(year1, year2, yob, yoc, jobloss.date, yin)]
		} else {
			dat <- get(paste0(code, ".dat2"))
		}

		# Apply year.max
		if (is.finite(year.max)) {
			dat <- dat[year <= year.max]
		}

		# dat <- dat[studyno %in% sample(unique(studyno), 8000)]
		# Save data ####
		if (save_dat) {
			assign(paste0(code, ".dat2"),
						 dat,
						 inherits = T,
						 envir = .GlobalEnv)
			Sys.sleep(0)
		}

		# Restrict age ####
		if (is.finite(age_under)) {
			setorder(dat, studyno, age.year2)
			dat[,I := 1:.N, by = .(studyno)]
			include.who <- dat[age.year1/365 < age_under & I == 1, studyno]
			dat <- dat[studyno %in% include.who]
		}

		basic_formula <- lapply(c(
			"binary",
			"age 50",
			"age 55",
			"age 60"), function(x = "age 50") {
				paste(
					ifelse(
						grepl("age", time_scale),
						"Surv(age.year1, age.year2, status) ~",
						"Surv(year1, year2, status) ~"
					),
					paste0("`Employment status", gsub("[a-z]", "", x), "`"),
					"+ `Duration of employment`",
					"+ `Cumulative straight`",
					ifelse(
						!is.finite(messy_sol),
						"+ `Cumulative soluble`",
						paste0("+ `Cumulative soluble", paste0(" ", messy_sol %/% 0.01), "`")),
					"+ `Cumulative synthetic`",
					ifelse(
						grepl("age", time_scale),
						ifelse(!spline_year, "+ Year",
									 paste0("+ pspline(year, df = ", year.df, ")")),
						"+ Age"
					),
					ifelse(
						!spline_yin,
						"+ `Year of hire`",
						paste0("+ pspline(yin.gm, df = ", yin.df, ")")
					),
					if (length(table(dat$Race)) > 1) {"+ Race"},
					if (length(table(dat$Plant)) > 1) {"+ Plant"},
					ifelse(!grepl('breast|male|prosate', description, ignore.case = T), "+ Sex", "")
				)
			})
		names(basic_formula) <- c("binary", paste0("age", seq(50, 60, 5)))

		dat <- dat[immortal != 1 & right.censored != 1]
		Sys.sleep(0)

		# Directory names model  ####
		if (is.null(directory.name)) {
			directory.name <- to_drive_D(gsub("//", "/", here::here(
				paste('./resources/hwse 2',
							ifelse(is.finite(year.max), paste0("FU through ", year.max), ""),
							paste0("lag ", 1 + additional.lag),
							ifelse(employment_status.lag != 0,
										 paste0("Employment status lagged ", employment_status.lag, " years"),
										 ""
							),
							ifelse(!grepl("age", time_scale),
										 "indexed by calendar",
										 "indexed by age"
							), sep = "/"))))
		}

		if (length(directory.name) <= 1) {
			directory.name <- paste(
				directory.name, c("Binary", "Age 50", "Age 55", "Age 60"), sep = "/")
		}

		if (mi > 0) {
			directory.name <- paste0(directory.name, "/mi race")
		}


		# # Test
		# cat("\n")
		# cat(paste0(gsub("//", "/", paste0(directory.name,
		# 																	ifelse(mi > 0, paste0("/", code), ""))), collapse = "\n"))

		# Make directories, if they don't exist
		sapply(gsub("//", "/", paste0(directory.name, ifelse(mi > 0, paste0("/", code), ""))), dir.create, showWarnings = F, recursive = T)

		# Get data to impute
		if (mi > 0) {
			if (!paste0(gsub(" ", "_", code), ".mi_race.M", mi, ".hwse2") %in% ls(envir = .GlobalEnv)) {
				to_impute <- get.mi_race(dat, mi, race.imputed.melt)
				assign(paste0(
					gsub(" ", "_", code),
					".mi_race.M", mi, ".hwse2"),
					to_impute, inherits = T)
			} else {
				to_impute <- get(paste0(gsub(" ", "_", code), ".mi_race.M", mi, ".hwse2"))
			}
		}

		# Run model conditional
		if (run_model) {
			# Multiple imputation loop ####
			lapply(if (mi == 0) {0} else {start_m:mi}, function(m = 1) {

				if (mi > 0) {

					to_impute[[m]] <- to_impute[[m]][,.(
						studyno, year, Race.impute, finrace = 9)]

					table(dat$Race, useNA = "always")
					nrow(dat)
					dat[finrace %in% c(0, 9), finrace := 9]
					dat$Race.impute <- NA
					dat <- merge(
						dat[,-"Race.impute"], to_impute[[m]],
						by = c("studyno", "year", "finrace"), all.x = T)
					nrow(dat)
					nrow(to_impute[[m]])
					dat[
						finrace %in% c(0, 9),
						Race := Race.impute]
					table(dat[immortal != 1 & right.censored != 1]$Race, useNA = "always")
				}


				# Fit model ####
				invisible(sapply(1:length(basic_formula), function(j = 1) {
					tmp.coxph <- coxph(
						as.formula(basic_formula[[j]]),
						data = {
							if (!grepl('breast|male|prosate', description, ignore.case = T)) {
								dat
							} else {
								dat[Sex == ifelse(grepl('breast|female', description, ignore.case = T), "Female", "Male")]
							}
						},
						method = "efron")

					# Save model  ####
					final_name <- paste(
						gsub("//", "/", paste0(directory.name[j], ifelse(mi > 0, paste0("/", code), ""))),
						paste0(code,
									 ifelse(
									 	is.finite(messy_sol),
									 	paste0("_sol", messy_sol %/% .01),
									 	""),
									 ifelse(spline_year,
									 			 ifelse(spline_yin, paste0("_splined"), "_splinedyear"),
									 			 ifelse(spline_yin, paste0("_splinedyin"), "")),
									 ifelse(mi > 0, paste0(".m", m), ""),
									 ".coxph.rds"),
						sep = "/"
					)

					message(paste0("\nSaving model to:\n",
												 final_name))
					saveRDS(tmp.coxph,
									file = final_name)

					# Print to console ####
					cat(paste0("\n", paste0(rep("_", 80), collapse = ""), "\n"))
					cat(unlist(description))
					cat("\n")
					if (mi > 0) {cat(paste0("m = ", m, "\n"))}
					cat(paste0(tmp.coxph$nevent, " cases"))
					# Time elapsed
					since_start <- time_length(difftime(Sys.time(), start), "minutes")
					if (since_start <= 90) {
						cat(paste0("\n", round(since_start, 2), " minutes since get.hwse2.coxph() was called."))
					} else {
						cat(paste0("\n", round(since_start/60, 3), " hours since get.hwse2.coxph() was called."))
					}
					cat(paste0("\n", paste0(rep("_", 80), collapse = ""), "\n"))
					print(summary(tmp.coxph)$coefficients)
				}))
			}) # End MI loop
		} # Run model conditional
	} # End loop over outcomes
	# )) # End sapply over outcomes
	options(warn = 0)
}

# # Try out get.hwse2.coxph() ####
# get.hwse2.coxph(outcomes = c(32:33))
#
# first.dat2[,.(
# 	final.year = as.Date(year2[.N], as.Date("1970-01-01")),
# 	status = status[.N],
# 	plant = plant[.N],
# 	yoi = yoi[.N],
# 	yod = yod[.N],
# 	year = year[.N],
# 	yoc = yoc[.N]
# ), by = .(studyno)][,.(plant, year)] %>% table
# first.dat2[,.(
# 	final.year = as.Date(year2[.N], as.Date("1970-01-01")),
# 	status = status[.N],
# 	plant = plant[.N],
# 	yoi = yoi[.N],
# 	yod = yod[.N],
# 	year = year[.N],
# 	yoc = yoc[.N]
# ), by = .(studyno)][final.year != yoi & !is.na(yoi) & (is.na(yoc) | yoc >= yoi)]
# first.dat2[,.(
# 	final.year = as.Date(year2[.N], as.Date("1970-01-01")),
# 	status = status[.N],
# 	plant = plant[.N],
# 	yoi = yoi[.N],
# 	yod = yod[.N],
# 	year = year[.N],
# 	yoc = yoc[.N]
# ), by = .(studyno)][yoc == yoi]
# first.dat2[year == year(jobloss.date),.(
# 	studyno, year, year1.date, year2.date, status, plant,
# 	`Employment status` = as.numeric(`Employment status`),
# 	yoi, yod, yoc,
# 	jobloss.date, unlagged_jobloss.date)][,.(
# 		N = .N,
# 		year2.date = year2.date[1],
# 		jobloss.date = jobloss.date[1],
# 		unlagged_jobloss.date = unlagged_jobloss.date[1],
# 		yoi = yoi[1],
# 		yod = yod[1],
# 		yoc = yoc[1]), by = .(studyno)][year2.date != jobloss.date]
# sum(first.dat2$status)
# sum(first.dat2[status == 1 & yoi < as.Date("2016-01-01"), 1, by = .(studyno)]$V1)

# HWSE Condition 3 - Exposure to leaving work ####
get.hwse3.coxph <- function(
	cohort_name =  NULL,
	run_model = F,
	new_dat = T,
	save_dat = T,
	messy_sol = 0.05,
	spline_year = F,
	year.df = 0,
	include_yin = F,
	spline_yin = F,
	yin.df = 0,
	time_scale = "age",
	additional.lag = 0,
	employment_status.lag = 0,
	directory.name = NULL,
	year.max = Inf,
	mi = 0,
	start.year = NULL,
	age_under = Inf) {

	# start time
	start <- Sys.time()

	options(warn = 2)

	# Clean up environment
	rm(list = grep("*\\.dat$", ls(envir = .GlobalEnv), value = T), envir = .GlobalEnv)
	rm(list = grep("tmp.coxph", ls(envir = .GlobalEnv), value = T), envir = .GlobalEnv)

	# Get data ####
	if (!("dat3" %in% ls(envir = .GlobalEnv)) | new_dat) {

		get.coxph(
			cohort_name = cohort_name,
			run_model = F,
			outcomes = 1,
			hwse2 = F,
			hwse3 = T,
			additional.lag = additional.lag,
			employment_status.lag = employment_status.lag,
			year.max = year.max,
			start.year = start.year
		)

		dat <- get("dat3", envir = .GlobalEnv)

		# Indicator for age at leaving work ####
		dat[,`:=`(
			age.leavework = time_length(difftime(jobloss.date, yob), 'year')
		)]
		dat$`Employment status` <- dat$status
		dat[, `:=`(
			`Employment status` = factor(
				`Employment status`, 0:1, c(
					ifelse(employment_status.lag == 0, "Still at work",
								 paste0("Not yet ", employment_status.lag, " years since leaving work")),
					ifelse(employment_status.lag == 0, "Left work (any age)",
								 paste0("Left work at least ", employment_status.lag, " years ago"))
				)),
			`Employment status 50` = {
				tmp <- ifelse(age.leavework < 51,	2, 3)
				tmp[`Employment status` == 0] <- 1
				factor(tmp, levels = 1:3,
							 labels = c(
							 	ifelse(employment_status.lag == 0, "Still at work",
							 				 paste0("Not yet ", employment_status.lag, " years since leaving work")),
							 	ifelse(employment_status.lag == 0, "Left work (Age 50 or younger)",
							 				 paste0("Left work at least ", employment_status.lag, " years ago (Age 50 or younger)")),
							 	ifelse(employment_status.lag == 0, "Left work (Age 51 or older)",
							 				 paste0("Left work at least ", employment_status.lag, " years ago (Age 51 or older)"))))},
			`Employment status 55` = {
				tmp <- ifelse(age.leavework < 56,	2, 3)
				tmp[`Employment status` == 0] <- 1
				factor(tmp, levels = 1:3,
							 labels = c(
							 	ifelse(employment_status.lag == 0, "Still at work",
							 				 paste0("Not yet ", employment_status.lag, " years since leaving work")),
							 	ifelse(employment_status.lag == 0, "Left work (Age 55 or younger)",
							 				 paste0("Left work at least ", employment_status.lag, " years ago (Age 55 or younger)")),
							 	ifelse(employment_status.lag == 0, "Left work (Age 56 or older)",
							 				 paste0("Left work at least ", employment_status.lag, " years ago (Age 56 or older)"))
							 ))},
			`Employment status 60` = {
				tmp <- ifelse(age.leavework < 61,	2, 3)
				tmp[`Employment status` == 0] <- 1
				factor(tmp, levels = 1:3,
							 labels = c(
							 	ifelse(employment_status.lag == 0, "Still at work",
							 				 paste0("Not yet ", employment_status.lag, " years since leaving work")),
							 	ifelse(employment_status.lag == 0, "Left work (Age 60 or younger)",
							 				 paste0("Left work at least ", employment_status.lag, " years ago (Age 60 or younger)")),
							 	ifelse(employment_status.lag == 0, "Left work (Age 61 or older)",
							 				 paste0("Left work at least ", employment_status.lag, " years ago (Age 61 or older)"))
							 ))}
		)]

		setorder(dat, studyno, year, year2)

		# dat <- dat[!(year == year(jobloss.date) &
		# 						 	yoi < jobloss.date &
		# 						 	`Employment status` == 1),]

		# Check time 2 > time 1
		# dat[year2 <= year1, .(year1, year2, yob, yoc, jobloss.date, yin)]
	} else {dat <- get("dat3")}

	# dat$studyno %>% n_distinct

	# Apply year.max
	if (is.finite(year.max)) {
		dat <- dat[year <= year.max]
	}

	# Save data ####
	if (save_dat) {
		assign("dat3",
					 dat,
					 inherits = T,
					 envir = .GlobalEnv)
		Sys.sleep(0)
	}

	# Restrict age ####
	if (is.finite(age_under)) {
		setorder(dat, studyno, age.year2)
		dat[,I := 1:.N, by = .(studyno)]
		include.who <- dat[age.year1/365 < age_under & I == 1, studyno]
		dat <- dat[studyno %in% include.who]
	}

	# Get data to impute
	if (mi > 0) {
		if (!paste0("mi_race.M", mi, ".hwse3") %in% ls(envir = .GlobalEnv)) {
			to_impute <- get.mi_race(dat, mi, race.imputed.melt)
			assign(paste0("mi_race.M", mi, ".hwse3"), to_impute, inherits = T)
		} else {
			to_impute <- get(paste0("mi_race.M", mi, ".hwse3"))
		}
	}

	if (run_model) {
		# Multiple imputation loop ####
		lapply(if (mi == 0) {0} else {1:mi}, function(m = 1) {

			if (mi > 0) {

				to_impute[[m]] <- to_impute[[m]][,.(
					studyno, year, Race.impute, finrace = 9)]

				table(dat$Race, useNA = "always")
				nrow(dat)
				dat[finrace %in% c(0, 9), finrace := 9]
				dat$Race.impute <- NA
				dat <- merge(
					dat[,-"Race.impute"], to_impute[[m]],
					by = c("studyno", "year", "finrace"), all.x = T)
				nrow(dat)
				nrow(to_impute[[m]])
				dat[
					finrace %in% c(0, 9),
					Race := Race.impute]
				table(dat[immortal != 1 & right.censored != 1]$Race, useNA = "always")
			}

			dat <- dat[immortal != 1 & right.censored != 1]
			Sys.sleep(0)
			# Fit model ####
			tmp.coxph <- coxph(as.formula(
				paste(
					ifelse(
						grepl("age", time_scale),
						"Surv(age.year1, age.year2, status) ~",
						"Surv(year1, year2, status) ~"
					),
					"`Cumulative straight` +",
					ifelse(
						!is.finite(messy_sol),
						"`Cumulative soluble` +",
						paste0("`Cumulative soluble", paste0(" ", messy_sol %/% 0.01), "` +")),
					"`Cumulative synthetic` +",
					ifelse(
						grepl("age", time_scale),
						ifelse(!spline_year, "Year +", paste0("pspline(year, df = ", year.df, ") +")),
						"Age +"
					),
					ifelse(include_yin, ifelse(
						!spline_yin,
						"`Year of hire` + ", paste0("pspline(yin.gm, df = ", yin.df, ") +")), ""),
					"Race + Plant + Sex"
				)
			),
			data = dat,
			method = "efron")

			if (is.null(directory.name)) {
				directory.name <- to_drive_D(here::here(paste0(
					'./resources/hwse 3',
					ifelse(is.finite(year.max), paste0("/FU through ", year.max), ""),
					"/Lag ", 1 + additional.lag,
					ifelse(employment_status.lag != 0,
								 paste0("/Employment status lagged ", employment_status.lag, " years"),
								 ""),
					ifelse( !grepl("age", time_scale),
									"/indexed by calendar",
									"/indexed by age")
				)))
			}

			if (mi > 0) {
				directory.name <- paste0(directory.name, ifelse(mi > 0, "/mi race"))
			}

			dir.create(directory.name, showWarnings = F, recursive = T)

			# Save model  ####
			saveRDS(tmp.coxph,
							file = paste(
								directory.name,
								gsub("^_", "", paste0(
									ifelse(
										is.finite(messy_sol),
										paste0("_sol", messy_sol %/% .01),
										""),
									ifelse(spline_year,
												 ifelse(spline_yin, paste0("_splined"), "_splinedyear"),
												 ifelse(spline_yin, paste0("_splinedyin"), "")),
									ifelse(mi > 0, paste0(".m", m), ""),
									".coxph.rds")),
								sep = "/"
							))

			# if (nrow(dat[age.year2 <= age.year1]) > 0) {
			# 	nrow(dat[age.year2 <= age.year1])
			# }

			# Print to console ####
			cat(paste0("\n", paste0(rep("_", 80), collapse = ""), "\n"))
			print("HWSE Condition 3: previous MWF exposure as a predictor of leaving work")
			cat("\n")
			if (mi > 0) {cat(paste0("m = ", m, "\n"))}
			cat(paste0(tmp.coxph$nevent, " cases"))
			# Time elapsed
			since_start <- time_length(difftime(Sys.time(), start), "minutes")
			if (since_start <= 90) {
				cat(paste0("\n", round(since_start, 2), " minutes since get.hwse2.coxph() was called."))
			} else {
				cat(paste0("\n", round(since_start/60, 3), " hours since get.hwse2.coxph() was called."))
			}
			cat(paste0("\n", paste0(rep("_", 80), collapse = ""), "\n"))
			print(summary(tmp.coxph)$coefficients)

		}) # End mi loop
	}
	options(warn = 0)
}

# Get Coef ####
get.coef <- function(
	outcomes = outcomes.which,
	cohort_name = NULL,
	analytic.name = NULL,
	new_dat = T,
	messy_sol = 0.05,
	spline_year = F,
	spline_yin = F,
	time_scale = "age",
	employment.which = "Binary",
	hwse2 = F,
	hwse3 = F,
	additional.lag = 0,
	employment_status.lag = 0,
	mod.name = NULL,
	mod.directory = NULL,
	directory = NULL,
	file.prefix = NULL,
	year.max = Inf,
	start.year = NULL,
	mi = 0,
	age_under = Inf
) {

	# start time
	start <- Sys.time()

	if (hwse3) {
		outcomes <- 1
	}

	invisible(sapply(outcomes, function(i = outcomes[1]) {

		code <- unlist(incidence.key[i, 1])
		description <- unlist(incidence.key[i, 2])
		var.name <- unlist(incidence.key[i, 3])

		if (hwse3) {
			code <- ""
			description <- "Leaving work"
		}

		# Get data ####
		if (is.null(analytic.name)) {
			if (!hwse2 & !hwse3) {
				if (!(paste0(code, ".dat") %in% ls(envir = .GlobalEnv)) | new_dat |
						!(paste0(gsub(" ", "_", code), ".mi_race.M", mi) %in% ls(envir = .GlobalEnv))) {
					get.coxph(cohort_name = cohort_name, run_model = F, outcomes = i, additional.lag = additional.lag)
				}
				dat.og <- data.table::copy(get(paste0(code, ".dat")))
			} else  if (hwse2) {
				if (!(paste0(code, ".dat2") %in% ls(envir = .GlobalEnv)) | new_dat |
						!(paste0(gsub(" ", "_", code), ".mi_race.M", mi, ".hwse2") %in% ls(envir = .GlobalEnv))) {
					get.hwse2.coxph(
						cohort_name = cohort_name, run_model = F, outcomes = i, additional.lag = additional.lag, employment_status.lag = employment_status.lag, year.max = year.max)}
				dat.og <- data.table::copy(get(paste0(code, ".dat2")))
			} else if (hwse3) {
				if (!("dat3" %in% ls(envir = .GlobalEnv)) | new_dat |
						!(paste0(gsub(" ", "_", code), ".mi_race.M", mi, ".hwse3") %in% ls(envir = .GlobalEnv))) {
					get.hwse3.coxph(cohort_name = cohort_name, run_model = F, additional.lag = additional.lag,
													employment_status.lag = employment_status.lag,
													year.max = year.max, start.year = start.year)
				}
				dat.og <- data.table::copy(get("dat3"))
			}

			if (!grepl('breast|male|prosate', description, ignore.case = T)) {
				dat.og <- dat.og
			} else if (grepl("prostate|^male", description, ignore.case = T)) {
				dat.og <- dat.og[Sex == "Male"]
			} else {
				dat.og <- dat.og[Sex == "Female"]
			}
		} else {dat.og <- get(analytic.name)}

		# Restrict age ####
		if (is.finite(age_under)) {
			setorder(dat.og, studyno, age.year2)
			dat.og[,I := 1:.N, by = .(studyno)]
			include.who <- dat.og[age.year1/365 < age_under & I == 1, studyno]
			dat.og <- dat.og[studyno %in% include.who]
		}

		dat <- data.table::copy(dat.og)

		# Get MI race data
		if (mi > 0) {
			to_impute <- get(
				gsub("^\\.", "", paste0(gsub(" ", "_", code), ".mi_race.M", mi, ifelse(!hwse2 & !hwse3, "", ifelse(!hwse2, ".hwse3", ".hwse2")))
				))}

		# Loop over mi ####
		if (mi > 0) {
			message(paste0("Combining M = ", mi, " model estimates for: ", description,
										 ifelse(hwse2, paste0(" and leaving work (", employment.which, ")"), "")))
			pb <- txtProgressBar(min = 0, max = mi, style = 3)
		}
		tmp.coef <- lapply(if (mi == 0) {0} else {1:mi}, function(m = 0) {

			if (is.null(mod.name)) {
				mod.name <- paste0(
					code,
					ifelse(is.finite(messy_sol), # & !(hwse2 | hwse3),
								 paste0("_sol", messy_sol %/% .01), ""),
					ifelse(spline_year,
								 ifelse(spline_yin,
								 			 paste0("_splined"),
								 			 "_splinedyear"),
								 ifelse(spline_yin,
								 			 paste0("_splinedyin"),
								 			 "")),
					ifelse(mi > 0, paste0(".m", m), ""),
					".coxph", '.rds')

				mod.name <- gsub("^_", "", mod.name)
			}

			if (is.null(mod.directory)) {
				mod.directory <- to_drive_D(here::here(
					paste(
						"./resources",
						ifelse(!(hwse2 | hwse3),
									 paste0("lag ", exposure.lag + additional.lag),
									 paste0(
									 	ifelse(hwse2, "hwse 2", "hwse 3"),
									 	ifelse(is.finite(year.max), paste0("/FU through ", year.max), ""),
									 	"/Lag ", 1 + additional.lag,
									 	ifelse(employment_status.lag != 0,
									 				 paste0("/Employment status lagged ", employment_status.lag, " years"),
									 				 ""))),
						ifelse(!grepl("age", time_scale),
									 paste("indexed by calendar",
									 			ifelse(hwse2, employment.which, ""),
									 			sep = "/"),
									 paste("indexed by age",
									 			ifelse(hwse2, employment.which, ""), sep = "/")),
						sep = "/")))
			}

			if (mi > 0) {
				mod.directory <- paste(mod.directory, ifelse(mi > 0, paste0("mi race/", code), ""),
															 sep = "/")
			}

			# Get model ####
			# message(gsub("/./", "/", paste0("Reading file: ", gsub("//", "/", paste0(mod.directory, "/", mod.name)))))
			if (mi > 0) {setTxtProgressBar(pb, m)}
			tmp.coxph <- readRDS(gsub("//", "/", gsub("//", "/", paste0(mod.directory, "/", mod.name))))
			# tmp.coxph

			library(mgcv)
			if (!"gam" %in% class(tmp.coxph)) {
				tmp.coef <- summary(tmp.coxph)$coefficients
			} else {
				tmp.coef <- rbind(summary(tmp.coxph)$p.table,
													summary(tmp.coxph)$s.table)
			}
			tmp.coef <- cbind(rownames(tmp.coef), tmp.coef)
			if ("coxph" %in% class(tmp.coxph)) {
				colnames(tmp.coef)[grep(
					"^$|Covariate|^coef|^Estimate|^se|^SE|^z|^t|^Pr|^p",
					colnames(tmp.coef))] <- c("Covariate", "Estimate", "SE", "t", "p")
			} else {
				colnames(tmp.coef)[1:5] <- c("Covariate", "Estimate", "SE", "t", "p")
			}
			tmp.coef <- as.data.table(tmp.coef)

			# Clean covariate name
			covariates <- names(attr(tmp.coxph$terms, "dataClasses"))
			covariates <- covariates[!grepl("Surv|strata\\(.*\\)|status", covariates)]
			quote.which <- which(grepl(" ", covariates) & !grepl("pspline|", covariates))
			covariates[quote.which] <- paste0("`", covariates[quote.which], "`")
			if ("gam" %in% class(tmp.coxph)) {
				covariates <- c("(Intercept)", covariates)
			}

			covariate.levels <- tmp.coxph$xlevels

			tmp.coef[!grepl(ifelse("gam" %in% class(tmp.coxph),
														 "^s\\(", "pspline\\("), Covariate),
							 Covariate := {
							 	covariates <- sapply(covariate.levels, function(x) {length(x)})
							 	unlist(sapply(1:length(covariates), function(i) {
							 		rep(names(covariates)[i], covariates[i] - 1)}))
							 }]

			tmp.coef[grepl(
				ifelse("gam" %in% class(tmp.coxph), "^s\\(", "pspline\\("), Covariate),`:=`(
					Covariate = as.vector(sapply(covariates[grepl(ifelse("gam" %in% class(tmp.coxph),
																															 "^s\\(", "pspline\\("), covariates)], rep, 2))
				)]

			# Clean levels
			tmp.coef[grep("soluble_", Covariate), Covariate := "Cumulative_soluble_5"]
			tmp.coef[
				gsub("`", "", Covariate) %in% names(covariate.levels), `:=`(
					level = unlist(
						sapply(names(covariate.levels)[
							names(covariate.levels) %in% gsub("`", "", tmp.coef$Covariate)],
							function(x) {
								unlist(covariate.levels[names(covariate.levels) == x])[-1]
							}))
				)]

			# Confidence intervals
			tmp.coef <- tmp.coef[, .(
				Covariate = gsub("`", "", Covariate),
				level,
				HR = exp(as.numeric(Estimate)),
				lower.ci = exp(as.numeric(Estimate) - qnorm(1 - 0.05/2) * as.numeric(SE)),
				upper.ci = exp(as.numeric(Estimate) + qnorm(1 - 0.05/2) * as.numeric(SE)),
				SE = as.numeric(SE),
				p = as.numeric(`p`),
				lower = substr(level,
											 unlist(gregexpr(
											 	"\\(", level
											 )) + 1,
											 unlist(gregexpr(",", level)) - 1),
				upper = substr(level,
											 unlist(gregexpr(",", level)) + 1,
											 unlist(gregexpr(
											 	"\\]", level
											 )) - 1)
			)]

			# Need for getting spline df
			full.covariates <- unique(tmp.coef$Covariate)

			# Set aside splined stuff
			if (spline_year | spline_yin) {
				tmp.spline.coef <- tmp.coef[grepl("pspline|^s\\(", Covariate)]
				tmp.spline.coef[, `:=`(
					Covariate = substr(
						Covariate, 1, gregexpr("\\)", Covariate))
				)]
				tmp.coef <- tmp.coef[!grepl("pspline|^s\\(", Covariate)]
			}

			covariates <- unique(tmp.coef$Covariate)

			# Pretty covariate/level ####
			tmp.coef <- rbindlist(lapply(
				covariates[which(covariates != "(Intercept)")],
				function(x = covariates[4]) {
					# Get relevant rows
					dat <- tmp.coef[grep(paste0(x, "$"), tmp.coef$Covariate), ]

					# Make pretty level
					lower <- dat$lower
					upper <- dat$upper
					if (!hwse2) {
						# Pretty levels Exposure-outcome ####
						if (grepl("Cumu", x)) {
							level <- paste0("$>", lower, "$ to $", upper, "$")
							level[grep("Inf", upper)] <-
								paste0("$>", lower[grep("Inf", upper)], "$")
						} else if (grepl("Year|Age", x)) {
							level <- paste0("$", as.numeric(lower) + 1, "$ to $", upper, "$")
							level[grep("Inf", upper)] <-
								paste0("$>", lower[grep("Inf", upper)], "$")
						} else if (grepl("Race", x)) {
							level <- levels(dat.og$Race)[-1]
						} else if (grepl("Plant", x)) {
							level <- levels(dat.og$Plant)[-1]
						} else if (grepl("Sex", x)) {
							level <- levels(dat.og$Sex)[-1]
						}

						dat$level <- level

						# Get referent level
						if (grepl("Cumu", x)) {
							ref.lower <- -Inf
							ref.upper <- min(as.numeric(lower))
							ref.level <-
								ifelse(ref.upper == 0,
											 "$0$",
											 paste0("$", 0, "$ to $", ref.upper, "$"))
						} else if (grepl("Year|Age", x)) {
							ref.lower <- ifelse(grepl("Year$", x), 1973,
																	ifelse(grepl("hire$", x), 1938,
																				 floor(min(
																				 	dat.og$age.year2
																				 ) / 365)))
							ref.upper <- min(as.numeric(lower))
							ref.level <- paste0("$", ref.lower, "$ to $", ref.upper, "$")
						} else if (grepl("Race", x)) {
							ref.upper <- NA
							ref.lower <- NA
							ref.level = "White"
						} else if (grepl("Plant", x)) {
							ref.upper <- NA
							ref.lower <- NA
							ref.level <- "1"
						} else if (grepl("Sex", x)) {
							ref.upper <- NA
							ref.lower <- NA
							ref.level <- "Male"
							level <- levels(dat.og$Sex)[-1]
						}
					}

					if (hwse2) {
						# Pretty levels HWSE 2 ####
						if (grepl("Employ", x)) {
							if (grepl("Binary", employment.which, ignore.case = T)) {
								level <- c("Not employed")}
							if (grepl("50", employment.which, ignore.case = T)) {
								level <- c("Left work (50 or younger)",
													 "Left work (51 or older)")}
							if (grepl("55", employment.which, ignore.case = T)) {
								level <- c("Left work (55 or younger)",
													 "Left work (56 or older)")}
							if (grepl("60", employment.which, ignore.case = T)) {
								level <- c("Left work (60 or younger)",
													 "Left work (61 or older)")}
						} else  if (grepl("Duration", x)) {
							level <- paste0("$>", lower, "$ to $", upper, "$")
							level[grep("Inf", upper)] <-
								paste0("$>", lower[grep("Inf", upper)], "$")
						} else if (grepl("Cumu", x)) {
							level <- paste0("$>", lower, "$ to $", upper, "$")
							level[grep("Inf", upper)] <-
								paste0("$>", lower[grep("Inf", upper)], "$")
						} else if (grepl("Year|Age", x)) {
							level <- paste0("$", as.numeric(lower) + 1, "$ to $", upper, "$")
							level[grep("Inf", upper)] <-
								paste0("$>", lower[grep("Inf", upper)], "$")
						} else if (grepl("Race", x)) {
							level <- levels(dat.og$Race)[-1]
						} else  if (grepl("Plant", x)) {
							level <- levels(dat.og$Plant)[-1]
						} else if (grepl("Sex", x)) {
							level <- levels(dat.og$Sex)[-1]
						}

						if (length(level) != length(dat$level)) {
							stop(paste0("\nIncorrect number of levels for covariate: ", x,
													ifelse(mi > 0, paste0("\n m = ", m), "")))
						}
						dat$level <- level

						# Get referent level
						if (grepl("Employ", x)) {
							ref.upper <- ref.lower <- NA
							ref.level <- "Still employed"
						} else  if (grepl("Duration", x)) {
							ref.lower <- -Inf
							ref.upper <- min(as.numeric(lower))
							ref.level <-
								ifelse(ref.upper == 0,
											 "$0$",
											 paste0("$", 0, "$ to $", ref.upper, "$"))
						} else if (grepl("Cumu", x)) {
							ref.lower <- -Inf
							ref.upper <- min(as.numeric(lower))
							ref.level <-
								ifelse(ref.upper == 0,
											 "$0$",
											 paste0("$", 0, "$ to $", ref.upper, "$"))
						} else if (grepl("Year|Age", x)) {
							ref.lower <- ifelse(grepl("Year$", x),
																	1973,
																	ifelse(grepl("hire$", x), 1938,
																				 floor(
																				 	min(dat.og$age.year2) / 365
																				 )))
							ref.upper <- min(as.numeric(lower))
							ref.level <-
								paste0("$", ref.lower, "$ to $", ref.upper, "$")
						} else if (grepl("Race", x)) {
							ref.upper <- NA
							ref.lower <- NA
							ref.level = "White"
						} else if (grepl("Plant", x)) {
							ref.upper <- NA
							ref.lower <- NA
							ref.level <- "1"
						} else if (grepl("Sex", x)) {
							ref.upper <- NA
							ref.lower <- NA
							ref.level <- "Male"
							level <- levels(dat.og$Sex)[-1]
						}
					}

					# Get case count
					if (!x %in% names(dat.og)) {
						x <- gsub("`", "", x)
					}
					n <- as.vector(table(dat.og[status == 1, x, with = F]))
					# MI race?
					if (grepl("Race", x) & mi > 0) {
						race.mi <- merge(dat.og[status == 1, .(studyno, year, finrace, Race)], to_impute[[m]][,.(studyno, year, Race.impute)],
														 on = c("studyno", "Year"),
														 all.x = T)
						race.mi[finrace %in% c(0, 9), `:=`(Race = Race.impute)]
						n <- as.vector(table(race.mi[, x, with = F]))
					}
					if (length(n) - 1 != length(dat$level)) {
						stop(paste0("\nIncorrect number of levels for covariate: ", x,
												ifelse(mi > 0, paste0("\n m = ", m), "")))
					}
					dat$n <- n[-1]

					# Make data.table
					rbindlist(list(
						data.table(
							Covariate = paste(x),
							level = ref.level,
							n = n[1],
							HR = 1,
							lower.ci = NA,
							upper.ci = NA,
							SE = NA,
							p = NA,
							lower = ref.lower,
							upper = ref.upper
						),
						dat
					),
					use.names = T)

				}))

			if (!"gam" %in% class(tmp.coxph)) {
				tmp.coef$events <- tmp.coxph$nevent
			} else {
				tmp.coef$events <- sum(tmp.coxph$y)
			}

			tmp.coef$df <- NA
			# Add back splined stuff ####
			if (spline_year | spline_yin) {
				tmp.coef <- rbindlist(list(
					tmp.coef,
					tmp.spline.coef[, .(
						Covariate = gsub("yin.gm", "year of hire", gsub(
							"year", "calendar year",
							paste0("P-spline of ", gsub("pspline\\(|, df.*$", "", Covariate)))),
						n = NA,
						p,
						events = tmp.coef$events[1],
						df = {if ("gam" %in% class(tmp.coxph)) {
							round(sapply(1:length(tmp.coxph$smooth), function(i) {tmp.coxph$smooth[[i]]$df}), 2)
						} else {round(tmp.coxph$df[sapply(
							Covariate,
							match, gsub("`", "", full.covariates))],
							2)}}
					)]
				), use.names = T, fill = T)

			}

			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			return(tmp.coef)
		}) # End MI loop

		# tmp.coef

		if (mi > 0) {
			names.og <- sapply(tmp.coef, names)
			nrow.og <- sapply(tmp.coef, nrow)
			tmp.coef <- rbindlist(tmp.coef)[,.(
				n = mean(n),
				HR = mean(HR),
				lower = unique(lower),
				upper = unique(upper),
				SE = sqrt(mean(SE^2) + (1 + 1/mi) * sum((log(HR) - mean(log(HR)))^2)/(mi - 1)),
				events = unique(events),
				df = mean(df)
			), by = .(Covariate, level)]
			tmp.coef[,`:=`(
				lower.ci = exp(log(HR) - SE * qnorm(1 - 0.05 / 2)),
				upper.ci = exp(log(HR) + SE * qnorm(1 - 0.05 / 2)),
				p = (1 - pnorm(abs(log(HR)), sd = SE)) * 2
			)]
			tmp.coef <- tmp.coef[,names.og[,1], with = F]
		} else {
			tmp.coef <- tmp.coef[[1]]
		}
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		# Cleaned table ####
		coef.tab <- tmp.coef[, .(
			Covariate = {
				tmp <- Covariate
				tmp[duplicated(Covariate)] <- NA
				tmp
			},
			level,
			n,
			HR = paste0("$", formatC(
				round(HR, 2), format = "f", digits = 2
			), "$"),
			`(95% CI)` = {
				lower <- formatC(round(lower.ci, 2),
												 format = "f",
												 digits = 2)
				upper <- formatC(round(upper.ci, 2),
												 format = "f",
												 digits = 2)
				ci <- paste0("$(", lower, ",\\,", upper, ")$")
				ci[grepl("NA", ci)] <- NA
				ci
			},
			p = as.character(formatC(round(p, 2), format = "f", digits = 2)),
			`SE` = as.character(formatC(round(SE, 3), format = "f", digits = 3)),
			events,
			df = as.character(formatC(round(df, 2), format = "f", digits = 2))
		)]

		coef.tab[grepl("NA", p), p := NA]

		if (spline_year | spline_yin) {
			coef.tab[grepl("spline", Covariate), Covariate := paste0(
				Covariate, " ($df = ", df, "$)"
			)]
		}

		# Save coefficient table ####
		if (is.null(directory)) {
			directory <- gsub("//", "/", here::here(paste(
				"./reports/resources",
				paste0(
					ifelse(hwse2 | hwse3,
								 paste0(ifelse(hwse2, "hwse 2", "hwse 3"),
								 			 ifelse(is.finite(year.max), paste0("/FU through ", year.max), ""),
								 			 "/Lag ", ifelse(!(hwse2 | hwse3), exposure.lag, 1) + additional.lag,
								 			 ifelse(employment_status.lag != 0,
								 			 			 paste0("/Employment status lagged ", employment_status.lag, " years"), "")),
								 paste0("Lag ", exposure.lag))),
				ifelse(grepl("age", time_scale), "indexed by age", "indexed by calendar"),
				ifelse(hwse2, employment.which, ""),
				sep = "/"
			))
			)
		}

		if (mi > 0) {
			directory <- paste(directory, ifelse(mi > 0, "mi race", ""), sep = "/")
		}

		dir.create(directory, showWarnings = F, recursive = T)

		if (is.null(file.prefix)) {
			file.prefix <- paste0(
				code,
				ifelse(is.finite(messy_sol), # & !(hwse2 | hwse3),
							 paste0("_sol", messy_sol %/% .01), ""),
				ifelse(spline_year,
							 ifelse(spline_yin, paste0("_splined"), "_splinedyear"),
							 ifelse(spline_yin, paste0("_splinedyin"), "")),
				ifelse(mi > 0, paste0(".M", mi), "")
			)

			file.prefix <- gsub("^_", "", file.prefix)}

		saveRDS(tmp.coef,
						paste(directory,
									paste0(file.prefix, ".coef.rds"),
									sep = "/")
		)

		saveRDS(coef.tab,
						paste(
							directory,
							paste0(file.prefix, ".tab.rds"),
							sep = "/")
		)

		coef.tab[is.na(coef.tab)] <- ""
		coef.tab[, (names(coef.tab)) := lapply(coef.tab, function(x) {
			gsub("\\\\,|\\$|NA", " ", x)
		})]

		coef.tab[, `:=`(` ` = ifelse(is.finite(as.numeric(p)), ifelse(
			p < 0.01, "**", ifelse(p < 0.05, "*", ifelse(p < 0.1, ".", ""))
		), ""))]

		assign(gsub(" ", "_", unlist(description)),
					 coef.tab,
					 envir = .GlobalEnv)

		# Print to console ####
		cat(paste0("\n", paste0(rep("_", 100), collapse = ""), "\n"))
		cat(unlist(description))
		cat("\n")
		cat(paste0("Exposure lagged ", ifelse(!(hwse2 | hwse3), exposure.lag, 1) + additional.lag, " years\n"))
		# Time elapsed
		since_start <- time_length(difftime(Sys.time(), start), "minutes")
		if (mi > 0) {
			if (since_start <= 90) {
				cat(paste0(round(since_start, 2), " minutes since get.coef() was called."))
			} else {
				cat(paste0(round(since_start/60, 3), " hours since get.coef() was called."))
			}}
		cat(paste0("\n", paste0(rep("_", 100), collapse = ""), "\n"))
		print(coef.tab[, -c("events")])
	}))
}

# Get figures ####
clean.coef.tab <- function(x = data.table::copy(coef.tab),
													 table.engine = "xtable",
													 additional.lag = 0) {
	tmp.tab <- data.table::copy(x)
	if (table.engine == 'xtable') {
		tmp.tab[is.na(`(95% CI)`), `(95% CI)` := "\\multicolumn{1}{c}{--}"]
		for (y in c("n", "p", "df")) {
			names(tmp.tab)[grep(paste0("^", y, "$"), names(tmp.tab))] <- paste0("$", y, "$")
		}
		names(tmp.tab)[grep("95%", names(tmp.tab))] <- gsub("%", "\\\\%", grep("95%", names(tmp.tab), value = T))
	} else {
		tmp.tab[is.na(`(95% CI)`), `(95% CI)` := "$-$"]
	}

	tmp.tab[!is.na(as.numeric(`$p$`)),`:=`(
		` ` = ifelse(`$p$` < 0.05 & !is.na(as.numeric(`$p$`)), "$*$", "")
	)]

	if (table.engine == 'xtable') {
		tmp.tab[!grepl("0.00", `$p$`) & !is.na(`$p$`), `:=`(
			`$p$` = paste0("", `$p$`, "\\phantom{0}"))]}
	tmp.tab$`$p$` <- gsub("0.00", "< 0.005", tmp.tab$`$p$`)

	return(tmp.tab)
}

# Set up ggplot-friendly data ####
get.ggtab <- function(mwf = "Straight",
											outcomes = outcomes.which,
											time_scale = "age",
											spline_year = F,
											spline_yin = F,
											messy_sol = 0.05,
											additional.lag = 0,
											coef.directory = NULL,
											mi = 0) {
	lapply(outcomes, function(i = outcomes[1]) {
		# Exposure-incidence model
		code <- unlist(incidence.key[i, 1]); description <- unlist(incidence.key[i, 2])

		if (is.null(coef.directory)) {
			coef.directory <- here::here(paste(
				"./reports/resources",
				paste0("Lag ", exposure.lag),
				ifelse(grepl("age", time_scale), "indexed by age", "indexed by calendar"),
				sep = "/"
			))
		}

		if (mi > 0) {
			coef.directory <- paste0(coef.directory, "/mi race")
		}

		file.prefix <- paste0(
			code,
			ifelse(is.finite(messy_sol),
						 paste0("_sol", messy_sol %/% .01), ""),
			ifelse(spline_year,
						 ifelse(spline_yin, paste0("_splined"), "_splinedyear"),
						 ifelse(spline_yin, paste0("_splinedyin"), "")),
			ifelse(mi > 0, paste0(".M", mi), "")
		)

		file.prefix <- gsub("^_", "", file.prefix)

		coef.tab <- readRDS(paste(
			coef.directory,
			paste0(file.prefix, ".tab.rds"),
			sep = "/"))

		coef.tab <- clean.coef.tab(coef.tab)

		# Clean up covariate name
		coef.tab[Covariate == "Cumulative soluble 5",
						 Covariate := "Cumulative soluble"]

		coef.tab <- coef.tab[grepl(mwf, zoo::na.locf(Covariate), ignore.case = T)]
		coef.tab[grepl("--", `(95\\% CI)`), `(95\\% CI)` := NA]
		coef.tab[, `:=`(
			Outcome = paste0(description),#, " (", events, " cases)"),
			cases = paste0("(", events, " cases)"),
			Covariate = zoo::na.locf(Covariate),
			HR = gsub("\\$", "", HR),
			level = paste0(level, " (", `$n$`, " cases)"),
			ci.lower = substr(`(95\\% CI)`, 3, unlist(gregexpr(",\\\\,", `(95\\% CI)`)) - 1),
			ci.upper = substr(`(95\\% CI)`, unlist(gregexpr(",\\\\,", `(95\\% CI)`)) + 3, nchar(`(95\\% CI)`) - 2)
		)]

		coef.tab[,.(
			Outcome,
			Covariate,
			level,
			HR,
			ci.lower,
			ci.upper,
			cases
		)]
		# data.table(
		# 	Outcome = unique(coef.tab$Outcome),
		# 	Covariate = unique(coef.tab$Covariate),
		# 	HR = NA,
		# 	ci.lower = NA,
		# 	ci.upper = NA
		# )
	})
}

# HWSE Condition 2 figures
get.hwse.ggtab <- function(outcomes = outcomes.which,
													 time_scale = "age",
													 messy_sol = 0.05,
													 spline_year = F,
													 spline_yin = F,
													 additional.lag = 0,
													 employment_status.lag = 0,
													 coef.directory = NULL,
													 mi = 0,
													 year.max = 1994) {
	lapply(outcomes, function(
		i = outcomes.which[1]
	) {
		# Exposure-incidence model
		code <- unlist(incidence.key[i, 1]); description <- unlist(incidence.key[i, 2])

		if (is.null(coef.directory)) {
			coef.directory <- here::here(paste0("./reports/resources/hwse 2",
																					ifelse(is.finite(year.max), paste0("/FU through ", year.max), ""),
																					paste0("/lag ", 1 + additional.lag),
																					ifelse(employment_status.lag == 0, "",
																								 paste0("/Employment status lagged ",
																								 			 employment_status.lag, " years")),
																					ifelse(!grepl("age", time_scale),
																								 "/indexed by calendar", "/indexed by age"),
																					"/"))
		}

		dir.names <- c("Binary",
									 "Age 50", "Age 55",
									 "Age 60")

		if (length(coef.directory) <= 1) {
			coef.directory <- paste(coef.directory, dir.names, sep = "/")
		}

		if (mi > 0) {
			coef.directory <- paste0(coef.directory, "/mi race")
		}

		coef.directory <- gsub("//", "/", coef.directory)

		levels.new <- list(c("Still at work", "At any age"),
											 "50 or younger", "55 or younger",
											 "60 or younger")
		rows.which <- list(1:2,
											 2, 2,
											 2)
		coef.tab <- rbindlist(
			lapply(1:length(dir.names), function(i = 1) {
				coef.tab <- readRDS(
					paste(coef.directory[i],
								paste0(code,
											 ifelse(
											 	is.finite(messy_sol),
											 	paste0("_sol", messy_sol %/% .01),
											 	""),
											 ifelse(spline_year,
											 			 ifelse(spline_yin,
											 			 			 paste0("_splined"),
											 			 			 "_splinedyear"),
											 			 ifelse(spline_yin,
											 			 			 paste0("_splinedyin"),
											 			 			 "")),
											 ifelse(mi > 0, paste0(".M", mi), ""),
											 ".tab.rds"), sep = "/"))
				names(coef.tab)[grepl("95", names(coef.tab))] <- "(95% CI)"
				coef.tab <- clean.coef.tab(coef.tab[rows.which[[i]],])
				coef.tab[,level := unlist(levels.new[i])]
			}))

		coef.tab[, `:=`(
			Outcome = paste0(description),#, " (", events, " cases)"),
			cases = paste0("(", events, " cases)"),
			HR = gsub("\\$", "", HR),
			level = paste0(level, " (", `$n$`, " cases)"),
			ci.lower = substr(`(95\\% CI)`, 3, unlist(gregexpr(",\\\\,", `(95\\% CI)`)) - 1),
			ci.upper = substr(`(95\\% CI)`, unlist(gregexpr(",\\\\,", `(95\\% CI)`)) + 3, nchar(`(95\\% CI)`) - 2)
		)]

		coef.tab[,.(
			Outcome,
			level,
			HR,
			ci.lower,
			ci.upper,
			cases
		)]
	})
}

# Tikz and lualatex####
get.tikz <- function(
	ggtab.prefix = c("str", "sol", "syn"),
	file.prefix = NULL,
	directory = NULL,
	width = 6.5, height = 8.5) {

	prefix.which <- 1:length(ggtab.prefix)

	sapply(prefix.which, function(i = 3) {

		gg.tab <- get(paste0(ggtab.prefix[i], ".ggtab"))

		# legend.tab <- gg.tab[,.(
		# 	I = c(quantile(I, 0.25), quantile(I, 0.25) + 1.5),
		# 	description = c(cases[1], Outcome[1])),
		# 	by = .(Outcome)]
		# legend.tab[grepl('cases', description), Outcome := NA]
		legend.tab <- gg.tab[,.(
			I = c(quantile(I, 0.5)),
			description = paste0(
				"\\begin{tabular}{l}",
				Outcome[1],	"\\\\",
				cases[1],
				"\\end{tabular}")
		),
		by = .(Outcome)]

		dir.create(directory, showWarnings = F, recursive = T)
		tikz(file = paste(directory, paste0(ifelse(is.null(file.prefix[i]), ggtab.prefix[i], file.prefix[i]), ".tex"), sep = "/"),
				 standAlone = T, width = width, height = height)
		gridExtra::grid.arrange(
			ggplot(gg.tab, aes(
				x = I,
				y = as.numeric(HR),
				ymin = as.numeric(ci.lower),
				ymax = as.numeric(ci.upper),
				shape =  Outcome
			)) + geom_pointrange(size = 0.07) +
				scale_shape_manual(values = c(1:7, 1:7, 1:4)) +
				scale_x_continuous(breaks = gg.tab$I,
													 labels = gg.tab$level) +
				geom_hline(aes(yintercept = 1), color = 'gray') +
				coord_flip(
					ylim = {
						if (ggtab.prefix[i] != "hwse") {
							c(0.4, 2.1)} else {c(0, 3.6)}},
					xlim = c(
						min(gg.tab$I),
						max(gg.tab$I))) +
				mytheme + labs(y = "") +
				theme(plot.margin = unit(c(0.1, 0, 0.1, 0.1), "cm"),
							panel.grid = element_blank(),
							axis.title.y = element_text(color = "white"),
							legend.position = "none"),
			# Legend
			ggplot(legend.tab, aes(
				x = I,
				y = 1,
				ymin = 1 - 0.01,
				ymax = 1 + 0.01,
				shape = Outcome,
				label = gsub(
					" cancer| cancers", "",
					gsub("nervous system", "neural", description))
			)) +
				geom_pointrange(size = 0.07) +
				scale_shape_manual(values = c(1:7, 1:7, 1:4)) +
				geom_text(aes(
					x = I,
					y = 1 + 0.05), position = 'identity', hjust = 0, size = 2.75) +
				# facet_wrap(. ~ prefix, scales = 'free_x', ncol = 1) +
				theme_classic() + coord_flip(ylim = c(0.99 , 1.4),
																		 xlim = c(
																		 	min(gg.tab$I),
																		 	max(gg.tab$I))) +
				labs(x = "", y = '') +
				scale_x_continuous(breaks = gg.tab$I) +
				theme(plot.margin = unit(c(0.1, 0.1, 0.1, -0.2), "cm"),
							legend.position = 'none',
							axis.text = element_text(colour = 'white'),
							axis.ticks = element_blank(),
							axis.line = element_line(color = "white"),
							axis.title.y = element_text(color = "white"),
							strip.text = element_text(colour = 'white'),
							strip.background = element_rect(color = "white")),
			ncol = 2, widths = c(0.6, 0.4)
		)
		dev.off()
	})
}

# Figures in facet view ####
get.facet_tikz <- function(
	ggtab.prefix = c("str", "sol", "syn", "hwse"),
	file.prefix = NULL,
	directory = NULL,
	messy_sol = 0.05,
	year.max = 1994,
	width = 11, height = 8.5) {

	lapply(ggtab.prefix, function(prefix = "sol") {

		i <- which(ggtab.prefix == prefix)

		if (is.null(directory)) {
			directory <- here::here(paste0("./reports/resources",
																		 ifelse(grepl("hwse", prefix), "/hwse 2", ""),
																		 ifelse(grepl("hwse", prefix) & year.max != 2015, paste0("/FU through ", year.max), ""),
																		 "/Lag ", exposure.lag))
		}


		gg.tab <- data.table::copy(get(paste0("og_", prefix, ".ggtab")))
		gg.tab <- gg.tab[!grepl(" 50| 55| 60", level)]
		gg.tab[grepl("any", level), level := gsub("At any age", "Left work", level)]

		gg.tab[is.na(as.numeric(ci.lower)), `:=`(
			ci.lower = "1.00", ci.upper = "1.00"
		)]

		gg.tab[, `:=`(
			Outcome = factor(paste(Outcome, cases),
											 levels = unique(paste(Outcome, cases))),
			HR = as.numeric(HR),
			ci.lower = as.numeric(ci.lower),
			ci.upper = as.numeric(ci.upper))]

		gg.tab[,`:=`(
			level.factor = factor(.N:1)
		), by = .(Outcome)]

		tikz(file = paste(
			directory,
			paste0(ifelse(is.null(file.prefix[i]), paste0(prefix, paste0("_sol", messy_sol %/% .01)), file.prefix[i]), "_facet.tex"), sep = "/"),
			standAlone = T, width = width, height = height)
		print(ggplot(gg.tab, aes(
			x = level.factor,
			y = HR,
			ymin = ci.lower,
			ymax = ci.upper
		)) + geom_pointrange(size = 0.15) +
			geom_hline(aes(yintercept = 1), color = 'gray') +
			geom_text(aes(x = level.factor, y = ifelse(grepl("hwse", prefix), 0, 0.33), label = level), hjust = 1, size = 3) +
			coord_flip(
				ylim = {if (grepl("hwse", prefix)) {
					c(-1.25, 2)
				} else {c(-1.225, 2.5)}},
				xlim = c(max(gg.tab[,.N, by = .(Outcome)]$N) + 1, 0)
			) +
			scale_y_continuous(breaks = c(0.5, 1, 1.5, 2, if (!grepl("hwse", prefix)) {2.5})) +
			mytheme + labs(y = "") +
			facet_wrap(. ~ Outcome, ncol = 3) +
			theme(
				# plot.margin = unit(c(0.1, 0, 0.1, 0.1), "cm"),
				axis.text.y = element_blank(),
				axis.ticks.y = element_blank(),
				panel.grid = element_blank(),
				axis.title.y = element_text(color = "white"),
				legend.position = "none"))
		dev.off()
	})
}

