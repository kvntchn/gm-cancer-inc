# mi-for-splines.R
# Running models with splined MWF exposure
# With MI for missing race
# July 23, 2021
# Kevin Chen

# Packages to load or install, if needed
packages <- c("here", "knitr", "tidyverse", "data.table", "lubridate", "Hmisc", "survival", "boxr", "tikzDevice", "xtable")

# Load or install packages
invisible(sapply(packages, function(package) {
	if (package %in% rownames(installed.packages())) {
		library(package, character.only = T)
	} else {
		install.packages(package)
		library(package, character.only = T)
	}
}))

# Authenticate to get data from Box
box_auth()

# # Get data and helper functions for building analytic data
# source(here("get-data.R"))
# box_save(
# 	list = ls(),
# 	file_name = "cohort_analytic.rdata",
# 	dir_id = 139229852587)
if (!"cohort_analytic" %in% ls()) {box_load(828354451667)}
source(here("incidence.R"))

# Outcomes to run
outcomes.which <- sapply(c("co", "es", "pa", "re", "st"), function(x) {
	which(incidence.key$code == x)
})

get.splined_mwf <- function(
	outcome = 3,
	spline_which = "straight",
	get.dat = F,
	mi = 18,
	width = 4, height = 3
) {

	if (get.dat) {
		# Build analytic data for Cox PH analysis for NHL
		get.coxph(outcomes = outcome, run_model = F, mi = mi)

		dat.og <- copy(get(paste0(
			incidence.key[outcome,"code"], ".dat")))

		mi_race <- copy(get(paste0(
			incidence.key[outcome,"code"], ".mi_race.M", mi)))

		# Save race imputations to Box
		box_write(
			mi_race,
			paste0(incidence.key[outcome,"code"], ".mi_race.M", mi, ".rds"),
			139229852587)
	}

	# Output directory
	output.dir <- to_drive_D(here::here(
		paste("resources", paste0("Lag ", exposure.lag), "indexed by age", "mi race",
					incidence.key[outcome, "code"], paste0("splined ", spline_which),
					sep = "/")))

	dir.create(output.dir, showWarnings = F, recursive = T)

	mod_mi.coxph <- lapply(1:mi, function(i = 1) {

		message("Running models for ", tolower(incidence.key[outcome, "description"]),
						" with splined ", spline_which, " (M =", mi, ")\r")
		pb <- txtProgressBar(min = 0, max = mi, style = 3)

		# Get imputed race
		dat <- merge(dat.og, mi_race[[i]], by = c("studyno", "year"), all.x = T, all.y = F)
		dat[finrace %in% c(0, 9), Race := Race.impute]
		dat[,Race := factor(Race, levels = c("White", "Black"))]

		# Run Cox PH model
		mod.coxph <- coxph(as.formula(paste(
			"Surv(age.year1, age.year2, status) ~",
			ifelse(grepl("str", spline_which), "pspline(cum_straight, df = 3) +", "`Cumulative straight` +"),
			ifelse(grepl("sol", spline_which), "pspline(cum_soluble, df = 3) +", "`Cumulative soluble 5` +"),
			ifelse(grepl("syn", spline_which), "pspline(cum_synthetic, df = 3) +", "`Cumulative synthetic` +"),
			"Year + `Year of hire` + Race + Plant + Sex"
		)),
		data = dat[immortal != 1 & right.censored != 1],
		method = "efron")

		# Saving models
		saveRDS(mod.coxph,
						file = paste(
							output.dir,
							paste0(incidence.key[outcome,'code'] , "_m", i, ".coxph.rds"),
							sep = "/"))

		setTxtProgressBar(pb, i)

		return(mod.coxph)

	})

	# Give data some within-individual indexing
	setorder(dat.og, studyno, year)
	dat.og[immortal != 1 & right.censored != 1, `:=`(I = 1:(.N), N = .N), studyno]

	outcome_mi.ggtab <- lapply(1:mi, function(i) {
		# Save data used to generate the spline plots
		outcome.termplot <- termplot(
			mod_mi.coxph[[i]],
			term = which(c("straight", "soluble", "synthetic") == spline_which),
			se = T, plot = F)

		# Extract relevant quantities
		outcome.ggtab <- data.frame(
			exposure = outcome.termplot[[1]]$x,
			`Point estimate` = outcome.termplot[[1]]$y,
			se = outcome.termplot[[1]]$se,
			Reference = mean(with(
				outcome.termplot[[1]],
				y[signif(x, 2) <= c(0, 0.05)[
					which(c("straight", "soluble") == spline_which)]
				])),
			df = mod_mi.coxph[[i]]$df[
				which(c("straight", "soluble") == spline_which)
			],
			mwf_ref = c(0, 0.05)[
				which(c("straight", "soluble") == spline_which)],
			mwf_max = as.numeric(quantile(
				unlist(dat.og[I == N, paste0("cum_", spline_which), with = F]),
				0.99)),
			check.names = F
		)
		setDT(outcome.ggtab)
		setorder(outcome.ggtab, exposure)

		# Contrast to reference
		outcome.ggtab[,`:=`(
			`Point estimate` = `Point estimate` - Reference
		)]

		# Thin out grid, if desired
		target_resolution <- max(50 * (max(outcome.ggtab$mwf_max) - min(outcome.ggtab$mwf_ref)), 1500)
		if (nrow(outcome.ggtab) > target_resolution) {
			outcome.ggtab <- outcome.ggtab[
				seq(1, nrow(outcome.ggtab), nrow(outcome.ggtab)/(50 * (mwf_max[1] - mwf_ref[1])))]
		}

		return(outcome.ggtab)
	})

	assign(paste0(incidence.key[outcome, "code"], "_", spline_which, "_mi", mi, ".ggtab"),
				 outcome_mi.ggtab,
				 envir = .GlobalEnv)

	# outcome_mi.ggtab <- co_straight_mi2.ggtab

	# Get pooled estimates
	outcome_mi.ggtab <- rbindlist(outcome_mi.ggtab)[, .(
		`Point estimate` = mean(`Point estimate`),
		se.pooled = sqrt(mean(se^2)),
		var.excess = sum((`Point estimate` - mean(`Point estimate`))^2/(
			5 - 1)),
		df = mean(df),
		mwf_ref = mwf_ref[1],
		mwf_max = mwf_max[1]
	), by = .(exposure)]

	# Calculate the standard error of the MI estimate
	outcome_mi.ggtab[, se := sqrt(se.pooled^2 + (1 + 1/5) * var.excess)]

	# Rugplot data
	rugplot.tab <- melt(
		dat.og[status == 1,.(Straight = cum_straight, Soluble = cum_soluble, Synthetic = cum_synthetic)],
		measure.vars = c("Straight", "Soluble", "Synthetic"),
		value.name = "exposure",
		variable.name = "mwf"
	)
	rugplot.tab$`Point estimate` <- -Inf

	rugplot.tab[tolower(mwf) == spline_which, mwf_max := unique(outcome_mi.ggtab$mwf_max)]

	rugplot.tab <- rugplot.tab[exposure < mwf_max]

	# Plot splines
	outcome_mi.ggtab[exposure > mwf_ref & exposure <= mwf_max] %>% ggplot(
		aes(x = exposure, y = exp(`Point estimate`))) +
		stat_function(fun = function(x) {1}, geom = 'line',
									color = 'grey') +
		geom_ribbon(aes(
			ymin = exp(`Point estimate` - qnorm(0.975) * se),
			ymax = exp(`Point estimate` + qnorm(0.975) * se)),
			alpha = 0.1) +
		geom_line() +
		geom_rug(data = rugplot.tab, length = unit(0.015, "npc")) +
		geom_label(data = outcome_mi.ggtab[,.(df = df[1], mwf_max = mwf_max[1])],
							 aes(x = 0.04 * mwf_max,
							 		y = 4.45,
							 		label = paste0("df = ", signif(df, 3))),
							 size = 2.5) +
		coord_cartesian(ylim = c(0.5, 4.5)) +
		# facet_wrap(. ~ mwf, scales = "free_x") +
		labs(
			x = paste0('Cumulative exposure to ', spline_which,  'MWF (mg/m$^3\\cdot$years)'),
			y = paste0("HR of ", tolower(incidence.key[outcome, "description"]))) +
		theme_bw() + theme(strip.text = element_text(size = 11)) -> spline.ggplot

	print(spline.ggplot)

	plot.dir <- here::here(paste0(
		"./reports/resources",
		"/Lag ", exposure.lag,
		"/indexed by age",
		"/mi race",
		"/splined mwf"))

	dir.create(plot.dir, showWarnings = F, recursive = T)

	tikzDevice::tikz(file = paste(
		plot.dir,
		paste0(incidence.key[outcome, "code"], "_", spline_which,  "_M", mi, ".tex"),
		sep = "/"),
		standAlone = T, width = width, height = height)
	print(spline.ggplot)
	dev.off()

}

for (outcome in outcomes.which) {
	get.splined_mwf(outcome, "straight", T)
}

for (outcome in outcomes.which) {
	get.splined_mwf(outcome, "synthetic", F)
}