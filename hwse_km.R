# Survival curves -- small cells? ####
invisible(sapply(nevent[, which(cases >= 100)], function(
	i = which(grepl("Lung", incidence.key$description, ignore.case = T))
) {
	# New yout
	dat <- get(paste0(incidence.key[i, 1], ".dat2"))
	names(dat)[names(dat) == "Employment status"] <- "Employment"
	km.mod <- survfit(Surv(age.year1, age.year2, status) ~ Employment,
										data = dat)
	km.dat <- data.table(
		# Time in days to years
		Age = km.mod$time / 365,
		Survival = km.mod$surv,
		Employment = gsub("Employment=", "", rep(names(km.mod$strata),
																	 km.mod$strata)),
		type = "Kaplan-Meier survival function",
		yout.which = "New YOUT"
	)

	# OG yout
	dat <- get(paste0(incidence.key[i, 1], ".dat2og"))
	names(dat)[names(dat) == "Employment status"] <- "Employment"
	km.mod <- survfit(Surv(age.year1, age.year2, status) ~ Employment,
										data = dat)
	km.dat <- rbindlist(list(
		km.dat,
		data.table(
		# Time in days to years
		Age = km.mod$time / 365,
		Survival = km.mod$surv,
		Employment = gsub("Employment=", "", rep(names(km.mod$strata),
																						 km.mod$strata)),
		type = "Kaplan-Meier survival function",
		yout.which = "Old YOUT"
	)))

	km.dat[, Employment := factor(Employment,
																levels = 0:1,
																labels = c("At work", "Left work"))]

	tikz(
		file = here::here(
			"cancer incidence/resources/survival plots/by employment status",
			paste0("survplot_", incidence.key[i, 1],
						 '.tex')
		),
		standAlone = T,
		width = 6.5,
		height = 3
	)
	print(
		ggplot(km.dat, aes(
			x = Age,
			y = Survival,
			color = Employment
		)) +
			geom_line() + xlab('Age (years)') +
			ylab(paste0(
				'Survival until ',
				tolower(incidence.key[i, 2])
			)) + labs(color = 'Left work?') +
			facet_wrap(. ~ yout.which, ncol = 2) +
			mytheme + theme(legend.position = "bottom")
	)
	dev.off()
}))

