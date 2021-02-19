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
outcomes.which <- grep("male|stomach|hodgkin|lung|rectal|colon|breas|prostat|colorec", incidence.key$description, ignore.case = T)

incidence.key[outcomes.which,]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Run model ####
# MWF-Cancer with Messy soluble ####
get.coxph(
	# outcomes = c(25, 31),
	run_model = T,
	time_scale = "age")
get.coef(
	# outcomes = c(25, 31),
	# new_dat = F,
	time_scale = "age")

rm(list = ls()[grepl("dat$", ls())]); Sys.sleep(0)

# HWSE 2 with Messy soluble ####
for (j in 1:length(employment_status.lag)) {
get.hwse2.coxph(
	# outcomes = grep("lung|rectal|colon|breast|all", incidence.key$description, ignore.case = T),
	run_model = T,
	spline_year = T,
	spline_yin = T,
	time_scale = "age",
	additional.lag = additional.lag[j],
	employment_status.lag = employment_status.lag[j],
	year.max = 1994
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
			 	additional.lag = additional.lag[j],
			 	employment_status.lag = employment_status.lag[j],
			 	year.max = 1994)})
}

#  HWSE 3 with Messy soluble ####
get.hwse3.coxph(
	run_model = T,
	additional.lag = additional.lag[j],
	employment_status.lag = employment_status.lag[j],
	year.max = 1994)
get.coef(time_scale = "age", hwse3 = T,
				 additional.lag = additional.lag[j],
				 employment_status.lag = employment_status.lag[j],
				 year.max = 1994)

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