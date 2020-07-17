# Run Cancer Incidence functions (race stratified) ####
# Kevin Chen
# July 15, 2020

library(here)

# rm(list = ls())
# rm(list = ls()[-grep("outcome.selected", ls())])
# rm(list = ls()[-grep('cohort', ls())])

outcomes.which <- c(3, 4, 6:8, 10:12, 15, 18, 19, 21, 25, 31, 32)
# Specify YOUT ####
yout.which <- "yout"
source(here::here("incidence.R"))
cohort_analytic_black <- cohort_analytic[race == "Not white"]
cohort2_black <- cohort2[race == "Not white"]
mortality.cohort2_black <- mortality.cohort2[race == "Not white"]

nevent.black <- get.nevent(cohort_name = "cohort_analytic_black")

nevent.black[outcomes.which,]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Run model ####
# MWF-Cancer with Messy soluble ####
get.coxph(
	cohort_name = "cohort_analytic_black",
	# outcomes = c(25, 31),
	run_model = T,
	time_scale = "age")
get.coef(
	cohort_name = "cohort_analytic_black",
	# outcomes = c(25, 31),
	# new_dat = F,
	time_scale = "age")

rm(list = ls()[grepl("dat$", ls())]); Sys.sleep(0)

# HWSE 2 with Messy soluble ####
get.hwse2.coxph(
	cohort_name = "cohort2_black",
	# outcomes = c(25, 31),
	run_model = T,
	time_scale = "age",
	additional.lag = 0,
	employment_status.lag = employment_status.lag
)
sapply(c("Binary", paste("Age", seq(50, 60, 5))),
			 function(x) {get.coef(
			 	cohort_name = "cohort2_black",
			 	# outcomes = c(25, 31),
			 	new_dat = F,
			 	time_scale = "age",
			 	employment.which = x,
			 	spline_year = T,
			 	spline_yin = T,
			 	hwse2 = T,
			 	additional.lag = 0,
			 	employment_status.lag = employment_status.lag)})

#  HWSE 3 with Messy soluble ####
get.hwse3.coxph(
	cohort_name = "cohort2_black",
	run_model = T,
	additional.lag = 0,
	employment_status.lag = employment_status.lag)
get.coef(
	cohort_name = "cohort2_black",
	time_scale = "age", hwse3 = T,
	additional.lag = 0,
	employment_status.lag = employment_status.lag)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# MWF-outcome HRs ####
str.ggtab <- rbindlist(get.ggtab())
str.ggtab[,`:=`(I = .N:1)]
og_str.ggtab <- as.data.table(as.data.frame(str.ggtab))
str.ggtab <- str.ggtab[
	!(grepl("^0|^\\$0", level) | is.na(level) | level == ""),
	]

sol.ggtab <- rbindlist(get.ggtab("Soluble"))
sol.ggtab[,`:=`(I = .N:1)]
og_sol.ggtab <- as.data.table(as.data.frame(sol.ggtab))
sol.ggtab <- sol.ggtab[
	!(grepl("^0|^\\$0", level) | level == "0 to 0.05" | is.na(level) | level == ""),
	]

syn.ggtab <- rbindlist(get.ggtab("Synthetic"))
syn.ggtab[,`:=`(I = .N:1)]
og_syn.ggtab <- as.data.table(as.data.frame(syn.ggtab))
syn.ggtab <- syn.ggtab[
	!(grepl("^0|^\\$0", level) | is.na(level) | level == ""),
	]

# Compile figure ####
get.tikz(
	ggtab.prefix = c("str", "sol", "syn"),
	file.prefix = c("str_sol5", "sol_sol5", "syn_sol5"),
	directory = here::here(paste0("./reports/resources/Lag ", exposure.lag)))
lualatex(pattern = "*sol5\\.tex",
				 directory = here::here(paste0("./reports/resources/Lag ", exposure.lag)))

# HWSE 2 HRs ####
og_hwse.ggtab <- rbindlist(get.hwse.ggtab(
	additional.lag = 0,
	# # Change messy_sol for clean soluble referent group
	# messy_sol = 0.05,
	time_scale = "age",
	employment_status.lag = employment_status.lag))
og_hwse.ggtab[,`:=`(I = .N:1)]
hwse.ggtab <- as.data.table(as.data.frame(og_hwse.ggtab))


hwse.ggtab <- hwse.ggtab[
	!(grepl("Still", level) | is.na(level) | level == ""),
	]

hwse.ggtab <- hwse.ggtab[!grepl("Still| 50| 55", level)]
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
				 directory = here::here(paste0(
				 	"./reports/resources/hwse 2/Lag ", exposure.lag,
				 	ifelse(employment_status.lag != 0,
				 				 paste0("/Employment status lagged ", employment_status.lag, " years"),
				 				 "")
				 )))
lualatex(pattern = "*\\.tex",
				 directory = here::here(paste0(
				 	"./reports/resources/hwse 2/Lag ", exposure.lag,
				 	ifelse(employment_status.lag != 0,
				 				 paste0("/Employment status lagged ", employment_status.lag, " years"),
				 				 "")
				 )))


# Facet view ###
get.facet_tikz()
# Compile
lualatex(pattern = "*facet\\.tex",
				 directory = here::here(paste0("./reports/resources/Lag ", exposure.lag)))
lualatex(pattern = "*facet\\.tex",
				 directory = here::here(paste0("./reports/resources/hwse 2/Lag ", exposure.lag)))