# Run Cancer Incidence functions ####
# Kevin Chen
# July 15, 2020

library(here)
set.seed(03032021)

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
outcomes.which <- grep("colon|^rectal|pancreatic|esophageal|stomach|laryn|lung|breast|prostate|kidney|bladder|melanom|^leuk|hodgkin|all cancers", incidence.key$description, ignore.case = T)

# # Clean up environment a little
# rm(cohort, dta, jobhist, jobhist_py, jobhist_py.cast, exposure)

incidence.key[outcomes.which,]

cohort_analytic <- cohort_analytic[wh == 0 & nohist == 1]
cohort2 <- cohort2[wh == 0 & nohist == 1]
mortality.cohort2 <- copy(cohort2)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # Get race model ####
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
# box_save(race.imputed.melt,
# 					file_name = "race.imputed.melt.rdata",
# 					dir_id = 132157022198,
# 					description = "Predictive distribution for race.")
# name        : race.imputed.melt.rdata
# file id     : 783387319906
# version     : V1
# size        : 2 GB
# modified at : 2021-03-04 23:02:10
# created at  : 2021-03-04 23:02:10
# uploaded by : kevchen@berkeley.edu
# owned by    : spa-ehsadmin@berkeley.edu
# shared link : None
#
# parent folder name :  MI Race
# parent folder id   :  132157022198
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Run models ####
# MWF-Cancer with Messy soluble ####
get.coxph(
	outcomes = outcomes.which,
	run_model = T,
	time_scale = "age",
	# mi = 50,
	age_under = 65,
	directory.name = to_drive_D(here::here(
		paste0("resources/under 65/Lag ", exposure.lag,
					 "/indexed by age",
					 ifelse(0, "/mi race", ""),
					 "/")))
)
get.coef(
	outcomes = outcomes.which,
	new_dat = F,
	time_scale = "age",
	# mi = 50,
	age_under = 65,
	mod.directory = to_drive_D(here::here(
		paste0("resources/under 65/Lag ", exposure.lag,
					 "/indexed by age",
					 "/"))),
	directory = here::here(paste(
		"./reports/resources/under 65",
		paste0("Lag ", exposure.lag),
		"indexed by age",
		sep = "/"
	))
)

# rm(list = ls()[grepl("dat$", ls())]); Sys.sleep(0)

# HWSE 2 with Messy soluble ####
get.hwse2.coxph(
	outcomes = outcomes.which,
	run_model = T,
	spline_year = F,
	spline_yin = F,
	time_scale = "age",
	additional.lag = additional.lag,
	employment_status.lag = employment_status.lag,
	year.max = 1994,
	# mi = 50,
	age_under = 65,
	directory.name = to_drive_D(gsub("//", "/", here::here(
		paste('./resources/under 65/hwse 2',
					paste0("FU through ", 1994),
					paste0("lag ", 1 + additional.lag),
					ifelse(employment_status.lag != 0, paste0(
						"Employment status lagged ", employment_status.lag, " years"), "" ),
					"indexed by age", sep = "/"))))
)
for (x in c("Binary", paste("Age", seq(50, 60, 5)))) {
	get.coef(
		outcomes = outcomes.which,
		new_dat = T,
		time_scale = "age",
		employment.which = x,
		spline_year = F,
		spline_yin = F,
		hwse2 = T,
		additional.lag = additional.lag,
		employment_status.lag = employment_status.lag,
		year.max = 1994,
		# mi = 50,
		age_under = 65,
		mod.directory = to_drive_D(here::here(
			paste0("./resources/under 65",
						 paste0("/hwse 2",
						 			 "/FU through 1994",
						 			 "/Lag ", 1 + additional.lag,
						 			 ifelse(employment_status.lag != 0, paste0(
						 			 	"/Employment status lagged ", employment_status.lag, " years"), "")),
						 paste("/indexed by age", x, sep = "/")))),
		directory = here::here(paste0(
			"./reports/resources/under 65",
			paste0(paste0("/hwse 2",
										paste0("/FU through 1994"),
										"/Lag ", 1 + additional.lag,
										ifelse(employment_status.lag != 0, paste0(
											"/Employment status lagged ", employment_status.lag, " years"), "")),
						 paste("/indexed by age", x, sep = "/")
			)))
	)
}

#  HWSE 3 with Messy soluble ####
get.hwse3.coxph(
	run_model = T,
	additional.lag = additional.lag,
	employment_status.lag = employment_status.lag,
	year.max = 1994,
	# mi = 50
	age_under = 65,
	directory.name = to_drive_D(here::here(paste0(
		'./resources/under 65/hwse 3',
		paste0("/FU through ", 1994),
		"/Lag ", 1 + additional.lag,
		ifelse(employment_status.lag != 0, paste0(
			"/Employment status lagged ", employment_status.lag, " years"), ""),
		"/indexed by age"
	)))
)
get.coef(hwse3 = T,
				 additional.lag = additional.lag,
				 employment_status.lag = employment_status.lag,
				 year.max = 1994,
				 new_dat = F,
				 # mi = 50,
				 age_under = 65,
				 mod.directory = to_drive_D(here::here(
				 	paste0("./resources/under 65",
				 				 paste0("/hwse 3",
				 				 			 "/FU through 1994",
				 				 			 "/Lag ", 1 + additional.lag,
				 				 			 ifelse(employment_status.lag != 0, paste0(
				 				 			 	"/Employment status lagged ", employment_status.lag, " years"), "")),
				 				 "/indexed by age"))),
				 directory = here::here(paste0(
				 	"./reports/resources/under 65",
				 	paste0(paste0("/hwse 3",
				 								paste0("/FU through 1994"),
				 								"/Lag ", 1 + additional.lag,
				 								ifelse(employment_status.lag != 0, paste0(
				 									"/Employment status lagged ", employment_status.lag, " years"), "")),
				 				 "/indexed by age"
				 	)))
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Figures ####
# MWF-outcome HRs ####
str.ggtab <- rbindlist(get.ggtab(
	# mi = 50,
	coef.directory = here::here(paste(
		"./reports/resources/under 65",
		paste0("Lag ", exposure.lag),
		"indexed by age",
		sep = "/"
	))
))
str.ggtab[,`:=`(I = .N:1)]
og_str.ggtab <- as.data.table(as.data.frame(str.ggtab))
str.ggtab <- str.ggtab[
	!(grepl("^0|^\\$0", level) | is.na(level) | level == ""),
]

sol.ggtab <- rbindlist(get.ggtab(
	"Soluble",# mi = 50,
	coef.directory = here::here(paste(
		"./reports/resources/under 65",
		paste0("Lag ", exposure.lag),
		"indexed by age",
		sep = "/"
	))))
sol.ggtab[,`:=`(I = .N:1)]
og_sol.ggtab <- as.data.table(as.data.frame(sol.ggtab))
sol.ggtab <- sol.ggtab[
	!(grepl("^0|^\\$0", level) | level == "0 to 0.05" | is.na(level) | level == ""),
]

syn.ggtab <- rbindlist(get.ggtab(
	"Synthetic",
	# mi = 50,
	coef.directory = here::here(paste(
		"./reports/resources/under 65",
		paste0("Lag ", exposure.lag),
		"indexed by age",
		sep = "/"
	))))
syn.ggtab[,`:=`(I = .N:1)]
og_syn.ggtab <- as.data.table(as.data.frame(syn.ggtab))
syn.ggtab <- syn.ggtab[
	!(grepl("^0|^\\$0", level) | is.na(level) | level == ""),
]

# Compile figure ####
get.tikz(
	ggtab.prefix = c("str", "sol", "syn"),
	file.prefix = paste0(c("str_sol5", "sol_sol5", "syn_sol5"),
											 # ".M50"
											 NULL),
	directory = here::here(paste0("./reports/resources/under 65/Lag ", exposure.lag)))
lualatex(pattern = paste0("*sol5",
													# "\\.M50",
													"\\.tex"),
				 directory = here::here(paste0("./reports/resources/under 65/Lag ", exposure.lag)))

# HWSE 2 HRs ####
for (j in 1:length(employment_status.lag)) {
	og_hwse.ggtab <- rbindlist(get.hwse.ggtab(
		outcomes = outcomes.which,
		# # Change messy_sol for clean soluble referent group
		# messy_sol = 0.05,
		time_scale = "age",
		additional.lag = additional.lag[j],
		employment_status.lag = employment_status.lag[j],
		# mi = 50,
		coef.directory = here::here(paste0(
			"./reports/resources/under 65",
			paste0(paste0("/hwse 2",
										paste0("/FU through 1994"),
										"/Lag ", 1 + additional.lag,
										ifelse(employment_status.lag != 0, paste0(
											"/Employment status lagged ", employment_status.lag, " years"), "")),
						 "/indexed by age"
			)))
		))
	og_hwse.ggtab[,`:=`(I = .N:1)]
	hwse.ggtab <- data.table::copy(og_hwse.ggtab)

	hwse.ggtab[grepl("Still", level), `:=`(ci.lower = 1, ci.upper = 1)]
	hwse.ggtab[grepl("any", level), level := gsub("At any age", "Left work", level)]

	hwse.ggtab <- hwse.ggtab[
		!(#grepl("Still", level) |
			is.na(level) | level == ""),
	]

	hwse.ggtab <- hwse.ggtab[!grepl(" 50| 55| 60", level)]
	hwse.ggtab <- hwse.ggtab[,.(
		level = c(NA, level),
		HR = c(NA, HR),
		ci.lower = c(NA, ci.lower),
		ci.upper = c(NA, ci.upper),
		cases = c(NA, cases),
		I = NA
	), by = .(Outcome)]

	hwse.ggtab[,`:=`(I = .N:1)]
	hwse.ggtab <- hwse.ggtab[!is.na(level)]

	# Compile for HWSE 2 ####
	get.tikz(ggtab.prefix = "hwse",
					 file.prefix = paste0("hwse_sol5",
					 										 # ".M50"
					 										 NULL),
					 directory = here::here(paste0(
					 	"./reports/resources/under 65/hwse 2/FU through 1994/Lag ", 1 + additional.lag[j],
					 	ifelse(employment_status.lag[j] != 0,
					 				 paste0("/Employment status lagged ", employment_status.lag[j], " years"),
					 				 ""),
					 	"/indexed by age"
					 )))
	lualatex(pattern = "*\\.tex",
					 directory = here::here(paste0(
					 	"./reports/resources/under 65/hwse 2/FU through 1994/Lag ", 1 + additional.lag[j],
					 	ifelse(employment_status.lag[j] != 0,
					 				 paste0("/Employment status lagged ", employment_status.lag[j], " years"),
					 				 ""),
					 	"/indexed by age")))
}

# Facet view ####
get.facet_tikz(
	ggtab.prefix = c("str", "sol", "syn"),
	file.prefix = paste0(c("str_sol5", "sol_sol5", "syn_sol5"),
											 # ".M50"
											 NULL),
	directory = here::here(paste0("./reports/resources/under 65/Lag ", exposure.lag)))

get.facet_tikz(
	ggtab.prefix = "hwse",
	file.prefix = paste0("hwse_sol5",
											 # ".M50"
											 NULL),
	directory = here::here(paste0(
		"./reports/resources/under 65/hwse 2/FU through 1994/Lag ", 1 + additional.lag[j],
		ifelse(employment_status.lag[j] != 0,
					 paste0("/Employment status lagged ", employment_status.lag[j], " years"),
					 ""),
		"/indexed by age"
	)),
	height = 6)

# Compile
lualatex(pattern = "*facet\\.tex", directory = here::here(paste0("./reports/resources/under 65/Lag ", exposure.lag)))
lualatex(pattern = "*facet\\.tex",
				 directory = here::here(paste0(
				 	"./reports/under 65/resources/hwse 2/FU through 1994/Lag ", 1 + additional.lag[j],
				 	ifelse(employment_status.lag[j] != 0,
				 				 paste0("/Employment status lagged ", employment_status.lag[j], " years"),
				 				 ""),
				 	"/indexed by age")))
