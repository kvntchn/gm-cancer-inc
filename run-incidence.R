# Run Cancer Incidence functions ####
# Kevin Chen
# July 15, 2020

library(here)

# rm(list = ls())
# rm(list = ls()[-grep("outcome.selected", ls())])
# rm(list = ls()[-grep('cohort', ls())])

# outcomes.which <- c(3, 4, 6:8, 10:12, 15, 18, 19, 21, 25, 31, 32)
# Specify YOUT ####
yout.which <- "yout"
source(here::here("get-data.R"))
source(here::here("incidence.R"))
employment_status.lag <- c(0)
additional.lag <- c(0)
outcomes.which1 <- grep("male|stomach|hodgkin|lung|rectal|colon|breas|prostat|colorec", incidence.key$description, ignore.case = T)
outcomes.which <- grep("colon|rectal|pancreatic|esophageal|stomach|laryn|lung|breast|prostate|kidney|bladder|melanom|^leuk|hodgkin|all cancers", incidence.key$description, ignore.case = T)
outcomes.which <- outcomes.which[!outcomes.which %in% outcomes.which1[-(10:11)]]

incidence.key[outcomes.which,]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # Get imputed race ####
# gm_dta <- cohort_analytic[
# 	nohist == 0 & wh == 1 & immortal == 0 &
# 		right.censored != 1 & year >= 1941 & (
# 			year(yin) < 1938 | year >= year(yin + 365.25 * 3))]
#
# gm_dta[,yin.date := gm.to.date(yin.gm)]
#
# gm_dta[,`:=`(employment.years = time_length(difftime(year2, yin.date), 'years')), by = .(studyno)]
#
# gm_dta[, `:=`(
# 	`Cumulative_time_off` = cut(cum_off, c(-Inf, 0, quantile(cum_off[cum_off > 0], seq(0.2, 1, 0.2)))),
# 	`Duration_of_employment` = cut(employment.years, c(-Inf, quantile(employment.years, seq(0.2, 1, 0.2)))),
# 	canc_corec = factor(canc_corec),
# 	canc_co = factor(canc_co),
# 	canc_re = factor(canc_re),
# 	canc_st = factor(canc_st),
# 	canc_lu = factor(canc_lu),
# 	canc_br = factor(canc_br),
# 	canc_fe = factor(canc_fe),
# 	canc_ma = factor(canc_ma),
# 	canc_pr = factor(canc_pr),
# 	canc_nhl = factor(canc_nhl),
# 	Age = cut(age.year2/365, quantile(age.year2/365, seq(0, 1, 1/6)), include.lowest = T),
# 	Year = cut(year,
# 						 unique(quantile(year, seq(0, 1, 1/6))),
# 						 include.lowest = T, dig.lab = 4),
# 	yin = cut(yin.gm,
# 						unique(quantile(yin.gm, seq(0, 1, 1/6))),
# 						include.lowest = T),
# 	Straight = cut(cum_straight,
# 								 c(-Inf, 0, quantile(cum_straight[cum_straight > 0], seq(1/3, 1, 1/3))),
# 								 include.lowest = T),
# 	Soluble5 = cut(cum_soluble,
# 								 c(-Inf, 0.05, quantile(cum_soluble[cum_soluble > 0.05], seq(1/3, 1, 1/3))),
# 								 include.lowest = T),
# 	Synthetic = cut(cum_synthetic,
# 									unique(c(-Inf, 0, quantile(cum_synthetic[cum_synthetic > 0], seq(1/3, 1, 1/3)))),
# 									include.lowest = T),
# 	Sex = sex,
# 	Plant = plant,
# 	All_causes = `All causes`,
# 	Race = factor(
# 		finrace,
# 		levels = c(1, 2, 0, 9),
# 		labels = c("White", "Black", NA, NA)),
# 	# An indicator of having observed race
# 	R = ifelse(finrace %in% c(0, 9), 0, 1)
# )]
#
# source(here::here("../gm-race-imputation/imputation.R"))
# race.imputed.melt <- get.race_posterior(
# 	gm_dta,
# 	dir.id = 132157022198,
# 	file.name = "race.mcmc.rds"
# )
# box_write(race.imputed.melt,
# 					file_name = "race.imputed.melt.rds",
# 					dir_id = 132157022198,
# 					description = "Predictive distribution for race.")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Run model ####
# MWF-Cancer with Messy soluble ####
get.coxph(
	outcomes = outcomes.which,
	run_model = T,
	time_scale = "age",
	mi = 15)
get.coef(
	outcomes = outcomes.which,
	# new_dat = F,
	time_scale = "age",
	mi = 50)

rm(list = ls()[grepl("dat$", ls())]); Sys.sleep(0)

# HWSE 2 with Messy soluble ####
get.hwse2.coxph(
	outcomes = outcomes.which,
	run_model = T,
	spline_year = T,
	spline_yin = T,
	time_scale = "age",
	additional.lag = additional.lag,
	employment_status.lag = employment_status.lag,
	year.max = 1994,
	mi = 50
)
sapply(c("Binary", paste("Age", seq(50, 60, 5))),
			 function(x = "Age 55") {get.coef(
			 	outcomes = outcomes.which,
			 	new_dat = F,
			 	time_scale = "age",
			 	employment.which = x,
			 	spline_year = T,
			 	spline_yin = T,
			 	hwse2 = T,
			 	additional.lag = additional.lag,
			 	employment_status.lag = employment_status.lag,
			 	year.max = 1994,
			 	mi = 50)})

# #  HWSE 3 with Messy soluble ####
# get.hwse3.coxph(
# 	run_model = T,
# 	additional.lag = additional.lag,
# 	employment_status.lag = employment_status.lag,
# 	year.max = 1994,
# 	mi = 50)
# get.coef(hwse3 = T,
# 				 additional.lag = additional.lag,
# 				 employment_status.lag = employment_status.lag,
# 				 year.max = 1994,
# 				 new_dat = F,
# 				 mi = 50)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# # MWF-outcome HRs ####
# str.ggtab <- rbindlist(get.ggtab())
# str.ggtab[,`:=`(I = .N:1)]
# og_str.ggtab <- as.data.table(as.data.frame(str.ggtab))
# str.ggtab <- str.ggtab[
# 	!(grepl("^0|^\\$0", level) | is.na(level) | level == ""),
# 	]
#
# sol.ggtab <- rbindlist(get.ggtab("Soluble"))
# sol.ggtab[,`:=`(I = .N:1)]
# og_sol.ggtab <- as.data.table(as.data.frame(sol.ggtab))
# sol.ggtab <- sol.ggtab[
# 	!(grepl("^0|^\\$0", level) | level == "0 to 0.05" | is.na(level) | level == ""),
# 	]
#
# syn.ggtab <- rbindlist(get.ggtab("Synthetic"))
# syn.ggtab[,`:=`(I = .N:1)]
# og_syn.ggtab <- as.data.table(as.data.frame(syn.ggtab))
# syn.ggtab <- syn.ggtab[
# 	!(grepl("^0|^\\$0", level) | is.na(level) | level == ""),
# 	]
#
# # Compile figure ####
# get.tikz(
# 	ggtab.prefix = c("str", "sol", "syn"),
# 	file.prefix = c("str_sol5", "sol_sol5", "syn_sol5"),
# 	directory = here::here(paste0("./reports/resources/Lag ", exposure.lag)))
# lualatex(pattern = "*sol5\\.tex",
# 				 directory = here::here(paste0("./reports/resources/Lag ", exposure.lag)))
#
# # HWSE 2 HRs ####
# for (j in 1:length(employment_status.lag)) {
# 	og_hwse.ggtab <- rbindlist(get.hwse.ggtab(
# 		outcomes = outcomes.which,
# 		# # Change messy_sol for clean soluble referent group
# 		# messy_sol = 0.05,
# 		time_scale = "age",
# 		additional.lag = additional.lag[j],
# 		employment_status.lag = employment_status.lag[j]))
# 	og_hwse.ggtab[,`:=`(I = .N:1)]
# 	hwse.ggtab <- as.data.table(as.data.frame(og_hwse.ggtab))
#
#
# 	hwse.ggtab <- hwse.ggtab[
# 		!(grepl("Still", level) | is.na(level) | level == ""),
# 	]
#
# 	hwse.ggtab <- hwse.ggtab[!grepl("Still| 50| 55", level)]
# 	hwse.ggtab <- hwse.ggtab[,.(
# 		level = c(NA, level),
# 		HR = c(NA, HR),
# 		ci.lower = c(NA, ci.lower),
# 		ci.upper = c(NA, ci.upper),
# 		cases = c(NA, cases),
# 		I = NA
# 	), by = .(Outcome)]
#
# 	hwse.ggtab[,`:=`(I = .N:1)]
# 	hwse.ggtab <- hwse.ggtab[!is.na(level)]
#
# 	# Compile for HWSE 2 ####
# 	get.tikz(ggtab.prefix = "hwse",
# 					 directory = here::here(paste0(
# 					 	"./reports/resources/hwse 2/Lag ", 1 + additional.lag[j],
# 					 	ifelse(employment_status.lag[j] != 0,
# 					 				 paste0("/Employment status lagged ", employment_status.lag[j], " years"),
# 					 				 "")
# 					 )))
# 	lualatex(pattern = "*\\.tex",
# 					 directory = here::here(paste0(
# 					 	"./reports/resources/hwse 2/Lag ", 1 + additional.lag[j],
# 					 	ifelse(employment_status.lag[j] != 0,
# 					 				 paste0("/Employment status lagged ", employment_status.lag[j], " years"),
# 					 				 "")
# 					 )))
# }
#
# # Facet view ###
# get.facet_tikz()
# # Compile
# lualatex(pattern = "*facet\\.tex",
# 				 directory = here::here(paste0("./reports/resources/Lag ", exposure.lag)))
# lualatex(pattern = "*facet\\.tex",
# 				 directory = here::here(paste0("./reports/resources/hwse 2/Lag ", exposure.lag)))