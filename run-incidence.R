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
race.imputed.melt <- F
fix_discrepancies <- T
hwse_age_under <- 75
source(here::here("get-data.R"))
source(here("breaks.R"))
source(here::here("incidence.R"))
employment_status.lag <- c(0)
additional.lag <- c(0)
lazy_indexing <- F

incidence.key <- rbindlist(list(
	incidence.key,
	data.table(
		code = "first_mort",
		description = "All cancer mortality",
		var.name = "All cancers",
		date.name = "yod")
))

outcomes.which <- grep("colon|^rectal|pancreatic|esophageal|stomach|laryn|lung|breast|prostate|kidney|bladder|melanom|^leuk|hodgkin|all cancers", incidence.key$description, ignore.case = T)
# outcomes.which <- grep("stomach", incidence.key$description, ignore.case = T)
length(outcomes.which)

# outcomes.which <- grep("mort", incidence.key$description, ignore.case = T)


# # Clean up environment a little
# rm(cohort, dta, jobhist, jobhist_py, jobhist_py.cast, exposure)

incidence.key[outcomes.which,]

# # Tweak yout (left consider left job if )
# cohort_analytic[yod < jobloss.date, jobloss.date := yod]
# cohort2[yod < jobloss.date, jobloss.date := yod]

# Unknown race as separate category? ####
cohort_analytic[finrace %in% c(0, 9), race := "Unknown"]
cohort2[finrace %in% c(0, 9), race := "Unknown"]

# Filter ####
cohort_analytic <- cohort_analytic[wh == 1 & nohist == 0]
cohort2 <- cohort2[wh == 1 & nohist == 0]
# mortality.cohort2 <- copy(cohort2)
# cohort2_85 <- cohort2[year >= 1985]
# cohort2_38 <- cohort2[yin.gm >= 1938]

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
# # MWF-Cancer with Messy soluble ####
# get.coxph(
# 	outcomes = outcomes.which,
# 	cohort_name = "cohort2_38",
# 	run_model = T,
# 	spline_year = F,
# 	spline_yin = F,
# 	time_scale = "age",
# 	mi = 50
# )
# get.coef(
# 	outcomes = outcomes.which,
# 	cohort_name = "cohort2_38",
# 	new_dat = F,
# 	time_scale = "age",
# 	spline_year = F,
# 	spline_yin = F,
# 	mi = 50
# )
#
# rm(list = ls()[grepl("dat$", ls())]); Sys.sleep(0)

# HWSE 2 with Messy soluble ####
get.hwse2.coxph(
	outcomes = outcomes.which,
	run_model = T,
	spline_year = T,
	year.df = 3,
	spline_yin = T,
	yin.df = 3,
	time_scale = "age",
	additional.lag = additional.lag,
	employment_status.lag = employment_status.lag,
	year.max = 1994,
	age_under = hwse_age_under,
	lazy_indexing = lazy_indexing,
	# mi = 50,
)
for (x in c("Binary", paste("Age", seq(50, 60, 5)))) {
	# x = "Binary"
	get.coef(
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
		age_under = hwse_age_under,
		lazy_indexing = lazy_indexing,
		# mi = 50,
	)
}

#  HWSE 3 with Messy soluble ####
get.hwse3.coxph(
	run_model = T,
	additional.lag = additional.lag,
	employment_status.lag = employment_status.lag,
	spline_year = T,
	include_yin = T,
	spline_yin = T,
	year.df = 0,
	year.max = 1994,
	age_under = hwse_age_under
	# mi = 50
)
get.coef(hwse3 = T,
				 additional.lag = additional.lag,
				 employment_status.lag = employment_status.lag,
				 new_dat = F,
				 spline_year = T,
				 spline_yin = T,
				 year.max = 1994,
				 age_under = hwse_age_under
				 # mi = 50
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# # Figures ####
# # MWF-outcome HRs ####
# str.ggtab <- rbindlist(get.ggtab(
# 	mi = 50
# ))
# str.ggtab[,`:=`(I = .N:1)]
# og_str.ggtab <- as.data.table(as.data.frame(str.ggtab))
# str.ggtab <- str.ggtab[
# 	!(grepl("^0|^\\$0", level) | is.na(level) | level == ""),
# ]
#
# sol.ggtab <- rbindlist(get.ggtab(
# 	"Soluble",
# 	mi = 50
# 	))
# sol.ggtab[,`:=`(I = .N:1)]
# og_sol.ggtab <- as.data.table(as.data.frame(sol.ggtab))
# sol.ggtab <- sol.ggtab[
# 	!(grepl("^0|^\\$0", level) | level == "0 to 0.05" | is.na(level) | level == ""),
# ]
#
# syn.ggtab <- rbindlist(get.ggtab(
# 	"Synthetic",
# 	mi = 50
# 	))
# syn.ggtab[,`:=`(I = .N:1)]
# og_syn.ggtab <- as.data.table(as.data.frame(syn.ggtab))
# syn.ggtab <- syn.ggtab[
# 	!(grepl("^0|^\\$0", level) | is.na(level) | level == ""),
# ]
#
# # Compile figure ####
# get.tikz(
# 	ggtab.prefix = c("str", "sol", "syn"),
# 	file.prefix = paste0(c("str_sol5", "sol_sol5", "syn_sol5"),
# 											 ".M50",
# 											 NULL),
# 	directory = here::here(paste0(
# 		"./reports/resources/",
# 		"Lag ", exposure.lag))
# 	)
# lualatex(pattern = paste0(
# 	"*sol5",
# 	"\\.M50",
# 	"\\.tex"),
# 	directory = here::here(paste0(
# 		"./reports/resources/",
# 		"Lag ", exposure.lag)))
#
# HWSE 2 HRs ####
og_hwse.ggtab <- rbindlist(get.hwse.ggtab(
	outcomes = outcomes.which,
	# Change messy_sol for clean soluble referent group
	messy_sol = 0.05,
	time_scale = "age",
	additional.lag = additional.lag,
	employment_status.lag = employment_status.lag,
	# mi = 50,
	spline_year = T,
	spline_yin = T,
	year.max = 1994,
	age_under = hwse_age_under,
	lazy_indexing = lazy_indexing
))
og_hwse.ggtab[,`:=`(I = .N:1)]

og_hwse.ggtab[grepl("Still", level), `:=`(ci.lower = 1, ci.upper = 1)]
og_hwse.ggtab[grepl("any", level), level := gsub("At any age", "Left work", level)]

og_hwse.ggtab <- og_hwse.ggtab[
	!(#grepl("Still", level) |
		is.na(level) | level == ""),
]

hwse.ggtab <- data.table::copy(og_hwse.ggtab)

hwse.ggtab <- hwse.ggtab[!grepl("^50|^55|^60", level)]


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
				 	"./reports/resources/",
				 	# "garcia yout/",
				 	ifelse(is.finite(hwse_age_under), paste0("under ", hwse_age_under, "/"), ""),
				 	"hwse 2/",
				 	"FU through 1994/",
				 	# "FU 1985 through 1994/",
				 	"Lag ", 1 + additional.lag,
				 	ifelse(employment_status.lag != 0,
				 				 paste0("/Employment status lagged ", employment_status.lag, " years"),
				 				 ""),
				 	"/indexed by age",
				 	if (lazy_indexing) {"/left work for entire year"} else  {NULL}
				 )),
				 hr.min = 0.3, hr.max = 7)
lualatex(pattern = "hwse_sol5\\.tex",
				 directory = here::here(paste0(
				 	"./reports/resources/",
				 	# "garcia yout/",
				 	ifelse(is.finite(hwse_age_under), paste0("under ", hwse_age_under, "/"), ""),
				 	"hwse 2/",
				 	"FU through 1994/",
				 	# "FU 1985 through 1994/",
				 	"Lag ", 1 + additional.lag,
				 	ifelse(employment_status.lag != 0,
				 				 paste0("/Employment status lagged ", employment_status.lag, " years"),
				 				 ""),
				 	"/indexed by age",
				 	if (lazy_indexing) {"/left work for entire year"} else  {NULL})))

# Selected outcomes ####
hwse.ggtab.selected <- og_hwse.ggtab[
	as.vector(sapply(c("All cancers", "Prostate", "Lung", "Colon", "Stomach"),
									 grep, Outcome))]
hwse.ggtab.selected[,`:=`(I = 1:.N)]
hwse.ggtab.selected[, Outcome := gsub(" cancer.*", "", Outcome)]
hwse.ggtab.selected[grepl("Lung", Outcome), Outcome := "Lung and\nbronchial"]
hwse.ggtab.selected[,`:=`(
	n = {
		paren1.where <- unlist(gregexpr("\\(", level))
		paren2.where <- unlist(gregexpr(" cases\\)", level))
		as.numeric(substr(level, paren1.where + 1, paren2.where - 1))
	},
	level = gsub(" \\(.*", "", level)
)]

label.ggtab <- hwse.ggtab.selected[,.(
	I.median = median(I[!grepl("Still at work", level)])
), Outcome]

ggplot(hwse.ggtab.selected[!grepl("Still at work", level)], aes(
	x = I, y = as.numeric(HR),
	ymin = as.numeric(ci.lower), ymax = as.numeric(ci.upper),
	shape = level)) +
	geom_pointrange(size = 0.5) +
	geom_hline(yintercept = 1, linetype = 2, size = 0.5, color = "grey") +
	scale_y_continuous(breaks = seq(0.5, 3, 0.5), trans = "log") +
	scale_x_continuous(breaks = label.ggtab$I.median,
										 labels = label.ggtab$Outcome) +
	labs(shape = "Age at leaving work ",
			 x = "Cancer type", y = "HR") +
	mytheme + theme(
		panel.grid = element_blank(),
		legend.position = "bottom"
	) + guides(shape = guide_legend(nrow = 2, byrow = TRUE)) -> hwse2.ggplot

hwse2.ggplot

# Compile for HWSE 2 (selected outcomes) ####
tikz(here::here(
	paste0("reports/resources/",
				 # "garcia yout/",
				 ifelse(is.finite(hwse_age_under), paste0("under ", hwse_age_under, "/"), ""),
				 "hwse 2/",
				 "FU through 1994/",
				 # "FU 1985 through 1994/",
				 "Lag ", 1 + additional.lag,
				 ifelse(employment_status.lag != 0,
				 			 paste0("/Employment status lagged ", employment_status.lag, " years"),
				 			 ""),
				 "/indexed by age",
				 if (lazy_indexing) {"/left work for entire year"} else  {NULL}),
	"hwse2_selected.tex"), standAlone = T, width = 4, height = 3.5)
print(hwse2.ggplot)
dev.off()

lualatex(pattern = "hwse2_selected\\.tex",
				 directory = here::here(
				 	paste0("reports/resources/",
				 				 # "garcia yout/",
				 				 ifelse(is.finite(hwse_age_under), paste0("under ", hwse_age_under, "/"), ""),
				 				 "hwse 2/",
				 				 "FU through 1994/",
				 				 # "FU 1985 through 1994/",
				 				 "Lag ", 1 + additional.lag,
				 				 ifelse(employment_status.lag != 0,
				 				 			 paste0("/Employment status lagged ", employment_status.lag, " years"),
				 				 			 ""),
				 				 "/indexed by age",
				 				 if (lazy_indexing) {"/left work for entire year"} else  {NULL})))

# # Facet view ####
# get.facet_tikz(
# 	ggtab.prefix = c("str", "sol", "syn"),
# 	file.prefix = paste0(c("str_sol5", "sol_sol5", "syn_sol5"),
# 											 # ".M50"
# 											 NULL),
# 	directory = here::here(paste0("./reports/resources/",
# 																"Lag ", exposure.lag)))
#
get.facet_tikz(
	ggtab.prefix = "hwse",
	file.prefix = paste0("hwse_sol5",
											 # ".M50"
											 NULL),
	directory = here::here(paste0(
		"./reports/resources/",
		# "garcia yout/",
		ifelse(is.finite(hwse_age_under), paste0("under ", hwse_age_under, "/"), ""),
		"hwse 2/",
		"FU through 1994/",
		# "FU 1985 through 1994/",
		"Lag ", 1 + additional.lag,
		ifelse(employment_status.lag != 0,
					 paste0("/Employment status lagged ", employment_status.lag, " years/"),
					 ""),
		"/indexed by age",
		if (lazy_indexing) {"/left work for entire year"} else  {NULL}
	)),
	hr.min = 0.01, hr.max = 5.5,
	width = 9, height = 6)

# # Compile
# lualatex(pattern = "*facet\\.tex", directory = here::here(paste0(
# 	"./reports/resources/",
# 	"Lag ", exposure.lag)))
lualatex(pattern = "*facet\\.tex",
				 directory = here::here(paste0(
				 	"./reports/resources/",
				 	# "garcia yout/",
				 	ifelse(is.finite(hwse_age_under), paste0("under ", hwse_age_under, "/"), ""),
				 	"hwse 2/",
				 	"FU through 1994/",
				 	# "FU 1985 through 1994/",
				 	"Lag ", 1 + additional.lag,
				 	ifelse(employment_status.lag != 0,
				 				 paste0("/Employment status lagged ", employment_status.lag, " years/"),
				 				 ""),
				 	"/indexed by age",
				 	if (lazy_indexing) {"/left work for entire year"} else  {NULL})))
