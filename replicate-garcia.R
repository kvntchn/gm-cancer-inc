# Run Cancer Incidence functions ####
# Kevin Chen
# July 15, 2020

library(boxr); box_auth()
library(lubridate)
library(tidyverse)
library(data.table)
# library(date)
library(xtable)
library(sas7bdat)
library(Hmisc)
library(tikzDevice)
source("~/headRs/00-my-theme.R")
library(here)

if (!('cohort' %in% ls(envir = .GlobalEnv))) {
	source(here::here('../gm-wrangling/wrangling', '00-hello.R'))
}
# get.exposure()
# no_exposure.who <- exposure[,.(no_straight = sum(!is.na(straight)) == 0), studyno][no_straight == T]$studyno

# rm(list = ls())
# rm(list = ls()[-grep("outcome.selected", ls())])
# rm(list = ls()[-grep('cohort', ls())])

# outcomes.which <- c(3, 4, 6:8, 10:12, 15, 18, 19, 21, 25, 31, 32)

source(here::here("incidence.R"))

incidence.key <- fread(here::here("resources", 'cancer-key-expanded.csv'))

outcomes.which <- grep("prostate|lung|colorectal|All cancers", incidence.key$description, ignore.case = T)

# # Clean up environment a little
# rm(cohort, dta, jobhist, jobhist_py, jobhist_py.cast, exposure)

incidence.key[outcomes.which,]

yout.which <- "yout"
use_seer <- F
race.imputed.melt <- NULL
employment_status.lag <- c(0)
additional.lag <- c(0)
deathage.max <- Inf

enforce_3_year <- F
fix_discrepancies <- F
restrict_by_exposure <- T

# # Get data ####
# rm(cohort_analytic); rm(cohort2)
# source(here::here("get-data.R"))
# box_write(cohort_analytic,
# 				 dir_id = 118903632077,
# 				 file_name = "cohort_full_replicate_garcia_hwse3.rds")
# #  name        : cohort_full_replicate_garcia_hwse3.rds
# #  file id     : 811171425676
# #  version     : V1
# #  size        : 30 MB
# #  modified at : 2021-05-15 18:25:29
# #  created at  : 2021-05-15 18:25:29
# #  uploaded by : kevchen@berkeley.edu
# #  owned by    : spa-ehsadmin@berkeley.edu
# #  shared link : None
# #
# #  parent folder name :  data
# #  parent folder id   :  118903632077
#
# box_write(cohort2,
# 				 dir_id = 118903632077,
# 				 file_name = "cohort_replicate_garcia_hwse3.rds")
# #  name        : cohort_replicate_garcia_hwse3.rds
# #  file id     : 811171083241
# #  version     : V1
# #  size        : 30 MB
# #  modified at : 2021-05-15 18:26:08
# #  created at  : 2021-05-15 18:26:08
# #  uploaded by : kevchen@berkeley.edu
# #  owned by    : spa-ehsadmin@berkeley.edu
# #  shared link : None
# #
# #  parent folder name :  data
# #  parent folder id   :  118903632077

# Download data
# box_auth()
cohort2_big <- box_read(811171083241)
# n_distinct(cohort[studyno %in% exposure$studyno])
# n_distinct(cohort2.og$studyno)
# cohort2_big <- copy(cohort2.og)

cohort2_big[, jobloss.date := gm.to.date(yout15[1]), studyno]

cohort2_big[,yoc := min(yoc[1], yob[1] + years(75), na.rm = T), studyno]
cohort2_big <- cohort2_big[year <= year(yoc) | is.na(yoc)]

cohort2_big <- cohort2_big[
	time_length(difftime(jobloss.date, yin), "year") > 3 &
		studyno %in% cohort[
			cancinccoh == 1 & cancinccoh09 == 1 &
				cancinccoh_new == 1 & cancinccoh09_new == 1 &
				cancinccoh15 & cancinccoh15_new == 1 & possdiscr == 0,
			studyno] &
		yin.gm >= 1938 &
		yin.gm <= 1985 &
		flag77 == 0 &
		possdiscr_new == 0 &
		oddend == 0 &
		nohist == 0 &
		wh == 1 &
		year < 1995 &
		# finrace %in% c(1, 2) &
		yob.gm > 1910 &
		(yod >= as.Date("1985-01-01") | is.na(yod)) &
		year <- year(jobloss.date)
]

n_distinct(cohort2_big$studyno)
# table(cohort2_big[year == year(jobloss.date), sex])
# sum(table(cohort2_big[year == year(jobloss.date), sex]))

cohort2_big_male <- cohort2_big[sex == "M"]
cohort2_big_female <- cohort2_big[sex == "F"]

cohort2 <- copy(cohort2_big)
cohort2 <- cohort2[year >= 1985]
# cohort2[,yoc := min(yoc[1], yob[1] + years(75), na.rm = T), studyno]
# cohort2 <- cohort2[year <= year(yoc) | is.na(yoc)]

cohort2_male <- cohort2[sex == "M"]
cohort2_female <- cohort2[sex == "F"]

cohort2_male[, .(
	`All cancers` = sum(canc_first == 1 & (ddiag_first <= yoc | is.na(yoc))),
	`Prostate cancer` = sum(canc_pr == 1 & (ddiag_pr <= yoc | is.na(yoc))),
	`Lung cancer` = sum(canc_lu == 1 & (ddiag_lu <= yoc | is.na(yoc))),
	`Colorectal cancer` = sum(canc_corec == 1 & (ddiag_corec <= yoc | is.na(yoc)))
),
studyno][,-1] %>% apply(2, table) %>% apply(2, function(x) {c(x, total = sum(x))})

cohort2_female[,.(
	`All cancers` = sum(canc_first == 1)
),
studyno][,-1] %>% apply(2, table)

# HWSE 2 ####
# Male
get.hwse2.coxph(
	outcomes = outcomes.which,
	cohort_name = "cohort2_male",
	run_model = T,
	spline_year = T,
	year.df = 3,
	spline_yin = T,
	yin.df = 3,
	time_scale = "age",
	# messy_sol = NA,
	additional.lag = additional.lag,
	employment_status.lag = employment_status.lag,
	year.max = 1994,
	# mi = 50,
	# age_under = 75,
	directory.name = to_drive_D(gsub("//", "/", here::here(
		paste('./resources/replicate garcia 2017/hwse 2',
					paste0("FU through ", 1994),
					paste0("lag ", 1 + additional.lag),
					ifelse(employment_status.lag != 0, paste0(
						"Employment status lagged ", employment_status.lag, " years"), "" ),
					"indexed by age", sep = "/"))))
)

for (x in c("Binary", paste("Age", seq(50, 60, 5)))) {
	get.coef(
		cohort_name = "cohort2_male",
		outcomes = outcomes.which,
		new_dat = F,
		# messy_sol = NA,
		time_scale = "age",
		employment.which = x,
		spline_year = T,
		spline_yin = T,
		hwse2 = T,
		additional.lag = additional.lag,
		employment_status.lag = employment_status.lag,
		year.max = 1994,
		# mi = 50,
		# age_under = 75,
		mod.directory = to_drive_D(here::here(
			paste0("./resources/replicate garcia 2017",
						 paste0("/hwse 2",
						 			 "/FU through 1994",
						 			 "/Lag ", 1 + additional.lag,
						 			 ifelse(employment_status.lag != 0, paste0(
						 			 	"/Employment status lagged ", employment_status.lag, " years"), "")),
						 paste("/indexed by age", x, sep = "/")))),
		directory = here::here(paste0(
			"./reports/resources/replicate garcia 2017",
			paste0(paste0("/hwse 2",
										paste0("/FU through 1994"),
										"/Lag ", 1 + additional.lag,
										ifelse(employment_status.lag != 0, paste0(
											"/Employment status lagged ", employment_status.lag, " years"), "")),
						 paste("/indexed by age", x, sep = "/")
			)))
	)
}

# Female
get.hwse2.coxph(
	outcomes = outcomes.which[length(outcomes.which)],
	cohort_name = "cohort2_female",
	run_model = T,
	spline_year = T,
	year.df = 3,
	spline_yin = T,
	yin.df = 3,
	time_scale = "age",
	# messy_sol = NA,
	additional.lag = additional.lag,
	employment_status.lag = employment_status.lag,
	year.max = 1994,
	# mi = 50,
	# age_under = 75,
	directory.name = to_drive_D(gsub("//", "/", here::here(
		paste('./resources/replicate garcia 2017/hwse 2',
					paste0("FU through ", 1994),
					paste0("lag ", 1 + additional.lag),
					ifelse(employment_status.lag != 0, paste0(
						"Employment status lagged ", employment_status.lag, " years"), "" ),
					"indexed by age",
					"female",
					sep = "/"))))
)

for (x in c("Binary", paste("Age", seq(50, 60, 5)))) {
	get.coef(
		cohort_name = "cohort2_female",
		outcomes = outcomes.which[length(outcomes.which)],
		new_dat = F,
		# messy_sol = NA,
		time_scale = "age",
		employment.which = x,
		spline_year = T,
		spline_yin = T,
		hwse2 = T,
		additional.lag = additional.lag,
		employment_status.lag = employment_status.lag,
		year.max = 1994,
		# mi = 50,
		# age_under = 75,
		mod.directory = to_drive_D(here::here(
			paste0("./resources/replicate garcia 2017",
						 paste0("/hwse 2",
						 			 "/FU through 1994",
						 			 "/Lag ", 1 + additional.lag,
						 			 ifelse(employment_status.lag != 0, paste0(
						 			 	"/Employment status lagged ", employment_status.lag, " years"), "")),
						 paste("/indexed by age/female", x, sep = "/")))),
		directory = here::here(paste0(
			"./reports/resources/replicate garcia 2017",
			paste0(paste0("/hwse 2",
										paste0("/FU through 1994"),
										"/Lag ", 1 + additional.lag,
										ifelse(employment_status.lag != 0, paste0(
											"/Employment status lagged ", employment_status.lag, " years"), "")),
						 paste("/indexed by age/female", x, sep = "/")
			)))
	)
}

# HWSE 2 Figure ####
og_hwse.ggtab <- rbindlist(get.hwse.ggtab(
	outcomes = outcomes.which[c(4, 3, 2, 1)],
	# Change messy_sol for clean soluble referent group
	messy_sol = 0.05,
	time_scale = "age",
	spline_year = T,
	spline_yin = T,
	additional.lag = additional.lag,
	employment_status.lag = employment_status.lag,
	# mi = 50,
	coef.directory = here::here(paste0(
		"./reports/resources/replicate garcia 2017",
		paste0(paste0("/hwse 2",
									paste0("/FU through 1994"),
									"/Lag ", 1 + additional.lag,
									ifelse(employment_status.lag != 0, paste0(
										"/Employment status lagged ", employment_status.lag, " years"), "")),
					 "/indexed by age"
		)))
))

og_hwse.ggtab[, Outcome := gsub(" cancer.*", "", Outcome)]
og_hwse.ggtab[grepl("All", Outcome), Outcome := "All cancers\namong men"]
og_hwse.ggtab[grepl("Lung", Outcome), Outcome := "Lung and\nbronchial"]

og_hwse_fe.ggtab <- rbindlist(get.hwse.ggtab(
	outcomes = outcomes.which[length(outcomes.which)],
	# Change messy_sol for clean soluble referent group
	messy_sol = 0.05,
	spline_year = T,
	spline_yin = T,
	time_scale = "age",
	additional.lag = additional.lag,
	employment_status.lag = employment_status.lag,
	# mi = 50,
	coef.directory = here::here(paste0(
		"./reports/resources/replicate garcia 2017",
		paste0(paste0("/hwse 2",
									paste0("/FU through 1994"),
									"/Lag ", 1 + additional.lag,
									ifelse(employment_status.lag != 0, paste0(
										"/Employment status lagged ", employment_status.lag, " years"), "")),
					 "/indexed by age/female"
		)))
))

og_hwse_fe.ggtab$Outcome <- "All cancers\namong women"

og_hwse.ggtab <- rbindlist(list(og_hwse.ggtab, og_hwse_fe.ggtab))

hwse.ggtab <- data.table::copy(og_hwse.ggtab)

hwse.ggtab[grepl("Still", level), `:=`(ci.lower = 1, ci.upper = 1)]
# hwse.ggtab[grepl("any", level), level := gsub("At any age", "Left work", level)]

hwse.ggtab <- hwse.ggtab[
	!(#grepl("Still", level) |
		is.na(level) | level == ""),
]

hwse.ggtab <- hwse.ggtab[,.(
	level = c(NA, gsub(" \\(.*\\)", "", level)),
	HR = as.numeric(c(NaN, HR)),
	ci.lower = as.numeric(c(NaN, ci.lower)),
	ci.upper = as.numeric(c(NaN, ci.upper)),
	cases = c(NA, cases),
	I = as.numeric(NaN)
), by = .(Outcome)]

hwse.ggtab[,`:=`(I = 1:.N)]
hwse.ggtab <- hwse.ggtab[!is.na(level)]

hwse.ggtab[,`:=`(
	I.median = median(I[!grepl("Still at work", level)])
), Outcome]

label.ggtab <- hwse.ggtab[,.(
	I.median = median(I[!grepl("Still at work", level)])
), Outcome]

ggplot(hwse.ggtab[!grepl("Still at work", level)], aes(
	x = I, y = HR, ymin = ci.lower, ymax = ci.upper,
	shape = level)) +
	geom_pointrange(size = 0.5) +
	geom_hline(yintercept = 1, linetype = 2, size = 0.5, color = "grey") +
	scale_y_continuous(breaks = seq(0.5, 3, 0.5), trans = "log") +
	scale_x_continuous(breaks = label.ggtab$I.median,
										 labels = label.ggtab$Outcome) +
	labs(shape = "Age at leaving work ",
			 x = "Cancer type") +
	mytheme + theme(
		panel.grid = element_blank(),
		legend.position = "bottom"
	) + guides(shape = guide_legend(nrow = 2, byrow = TRUE)) -> hwse2.ggplot

hwse2.ggplot

# Compile for HWSE 2 ####
tikz(here::here(
	"reports/resources/replicate garcia 2017/hwse 2/FU through 1994/Lag 1/indexed by age",
	"hwse2.tex"), standAlone = T, width = 4, height = 3.5)
print(hwse2.ggplot)
dev.off()

lualatex(pattern = "*\\.tex",
				 directory = here::here(
	"reports/resources/replicate garcia 2017/hwse 2/FU through 1994/Lag 1/indexed by age"))

#  HWSE 3 with Messy soluble ####
get.hwse3.coxph(
	cohort_name = "cohort2_big_male",
	run_model = T,
	additional.lag = additional.lag,
	employment_status.lag = employment_status.lag,
	year.max = 1994,
	# mi = 50
	spline_year = T,
	year.df = 3,
	start.year = -Inf,
	messy_sol = 1.2742,
	directory.name = to_drive_D(here::here(paste0(
		'./resources/replicate garcia 2017/hwse 3',
		paste0("/FU through ", 1994),
		"/Lag ", 1 + additional.lag,
		ifelse(employment_status.lag != 0, paste0(
			"/Employment status lagged ", employment_status.lag, " years"), ""),
		"/indexed by age"
	)))
)

get.coef(
	cohort_name = "cohort2_big_male",
	hwse3 = T,
	additional.lag = additional.lag,
	employment_status.lag = employment_status.lag,
	year.max = 1994,
	new_dat = F,
	messy_sol = 1.2742,
	spline_year = T,
	# mi = 50,
	mod.directory = to_drive_D(here::here(
		paste0("./resources/replicate garcia 2017",
					 paste0("/hwse 3",
					 			 "/FU through 1994",
					 			 "/Lag ", 1 + additional.lag,
					 			 ifelse(employment_status.lag != 0, paste0(
					 			 	"/Employment status lagged ", employment_status.lag, " years"), "")),
					 "/indexed by age"))),
	directory = here::here(paste0(
		"./reports/resources/replicate garcia 2017",
		paste0(paste0("/hwse 3",
									paste0("/FU through 1994"),
									"/Lag ", 1 + additional.lag,
									ifelse(employment_status.lag != 0, paste0(
										"/Employment status lagged ", employment_status.lag, " years"), "")),
					 "/indexed by age"
		)))
)

get.hwse3.coxph(
	cohort_name = "cohort2_big_female",
	run_model = T,
	additional.lag = additional.lag,
	employment_status.lag = employment_status.lag,
	year.max = 1994,
	start.year = -Inf,
	# mi = 50
	spline_year = T,
	year.df = 3,
	messy_sol = 0.11,
	directory.name = to_drive_D(here::here(paste0(
		'./resources/replicate garcia 2017/hwse 3',
		paste0("/FU through ", 1994),
		"/Lag ", 1 + additional.lag,
		ifelse(employment_status.lag != 0, paste0(
			"/Employment status lagged ", employment_status.lag, " years"), ""),
		"/indexed by age",
		"/female"
	)))
)

get.coef(
	cohort_name = "cohort2_big_female",
	hwse3 = T,
	additional.lag = additional.lag,
	employment_status.lag = employment_status.lag,
	year.max = 1994,
	new_dat = T,
	messy_sol = 0.11,
	spline_year = T,
	# mi = 50,
	mod.directory = to_drive_D(here::here(
		paste0("./resources/replicate garcia 2017",
					 paste0("/hwse 3",
					 			 "/FU through 1994",
					 			 "/Lag ", 1 + additional.lag,
					 			 ifelse(employment_status.lag != 0, paste0(
					 			 	"/Employment status lagged ", employment_status.lag, " years"), "")),
					 "/indexed by age",
					 "/female"))),
	directory = here::here(paste0(
		"./reports/resources/replicate garcia 2017",
		paste0(paste0("/hwse 3",
									paste0("/FU through 1994"),
									"/Lag ", 1 + additional.lag,
									ifelse(employment_status.lag != 0, paste0(
										"/Employment status lagged ", employment_status.lag, " years"), "")),
					 "/indexed by age",
					 "/female"
		)))
)
