# Comparing under 65 vs no age restriction ####
# Kevin Chen
# April 8, 2021

# Figures ####
# MWF-outcome HRs ####
str.ggtab <- rbindlist(list(
	rbindlist(get.ggtab(mi = 50)),
	rbindlist(get.ggtab(
	# mi = 50,
	coef.directory = here::here(paste(
		"./reports/resources/under 65",
		paste0("Lag ", exposure.lag),
		"indexed by age",
		sep = "/"
	))))
), idcol = c("Age restriction"))
str.ggtab[,`:=`(`Age restriction` = factor(`Age restriction`, 1:2, c("No restriction", "Under 65 at start of FU")))]
str.ggtab[,`:=`(I = .N:1), by = `Age restriction`]
og_str.ggtab <- as.data.table(as.data.frame(str.ggtab))
str.ggtab <- str.ggtab[
	!(grepl("^0|^\\$0", level) | is.na(level) | level == ""),
]

sol.ggtab <- rbindlist(list(
	rbindlist(get.ggtab("Soluble", mi = 50)),
	rbindlist(get.ggtab(
	"Soluble",# mi = 50,
	coef.directory = here::here(paste(
		"./reports/resources/under 65",
		paste0("Lag ", exposure.lag),
		"indexed by age",
		sep = "/"
	))))
	), idcol = "Age restriction")
sol.ggtab[,`:=`(`Age restriction` = factor(`Age restriction`, 1:2, c("No restriction", "Under 65 at start of FU")))]
sol.ggtab[,`:=`(I = .N:1), by = `Age restriction`]
og_sol.ggtab <- as.data.table(as.data.frame(sol.ggtab))
sol.ggtab <- sol.ggtab[
	!(grepl("^0|^\\$0", level) | level == "0 to 0.05" | is.na(level) | level == ""),
]

syn.ggtab <- rbindlist(list(
	rbindlist(get.ggtab("Synthetic", mi = 50)),
	rbindlist(get.ggtab(
	"Synthetic",# mi = 50,
	coef.directory = here::here(paste(
		"./reports/resources/under 65",
		paste0("Lag ", exposure.lag),
		"indexed by age",
		sep = "/"
	))))
	), idcol = "Age restriction")
syn.ggtab[,`:=`(`Age restriction` = factor(`Age restriction`, 1:2, c("No restriction", "Under 65 at start of FU")))]
syn.ggtab[,`:=`(I = .N:1), by = `Age restriction`]
og_syn.ggtab <- as.data.table(as.data.frame(syn.ggtab))
syn.ggtab <- syn.ggtab[
	!(grepl("^0|^\\$0", level) | is.na(level) | level == ""),
]

# Make long figure ####
ggtab.prefix <- c("str", "sol", "syn")
directory <- here::here(paste0("./reports/resources/under 65/Lag ", exposure.lag))
sapply(1:3, function(i = 1) {

		gg.tab <- get(paste0(ggtab.prefix[i], ".ggtab"))

		# legend.tab <- gg.tab[,.(
		# 	I = c(quantile(I, 0.25), quantile(I, 0.25) + 1.5),
		# 	description = c(cases[1], Outcome[1])),
		# 	by = .(Outcome)]
		# legend.tab[grepl('cases', description), Outcome := NA]
		legend.tab <- gg.tab[,.(
			I = c(quantile(I, 0.6)),
			description = paste0(
				"\\begin{tabular}{l}",
				Outcome[1],	"\\\\",
				cases[1],
				"\\end{tabular}")
		),
		by = .(Outcome)]

		dir.create(directory, showWarnings = F, recursive = T)
		tikz(file = paste(directory, paste0(ggtab.prefix[i], "_age-restriction.tex"), sep = "/"),
				 standAlone = T, width = 8.5, height = 8.5)
		gridExtra::grid.arrange(
			ggplot(gg.tab[grepl("No", `Age restriction`)], aes(
				x = I,
				y = as.numeric(HR),
				ymin = as.numeric(ci.lower),
				ymax = as.numeric(ci.upper),
				shape =  Outcome
			)) + geom_pointrange(size = 0.07) +
				scale_shape_manual(values = c(1:7, 1:7, 1:4)) +
				scale_x_continuous(breaks = gg.tab[grepl("No", `Age restriction`)]$I,
													 labels = gg.tab[grepl("No", `Age restriction`)]$level) +
				geom_hline(aes(yintercept = 1), color = 'gray') +
				coord_flip(ylim = {
					if (ggtab.prefix[i] != "hwse") {
						c(0.4, 2.1)} else {c(0, 3.6)}},
					xlim = c(
						min(gg.tab[grepl("No", `Age restriction`)]$I) + 1,
						max(gg.tab[grepl("No", `Age restriction`)]$I) - 1)) +
				facet_wrap(. ~ `Age restriction`, scales = "free_y") +
				mytheme + labs(y = "") +
				theme(plot.margin = unit(c(0.1, 0, 0.1, 0.1), "cm"),
							panel.grid = element_blank(),
							axis.title.y = element_text(color = "white"),
							legend.position = "none"),
			ggplot(gg.tab[!grepl("No", `Age restriction`)], aes(
				x = I,
				y = as.numeric(HR),
				ymin = as.numeric(ci.lower),
				ymax = as.numeric(ci.upper),
				shape =  Outcome
			)) + geom_pointrange(size = 0.07) +
				scale_shape_manual(values = c(1:7, 1:7, 1:4)) +
				scale_x_continuous(breaks = gg.tab[!grepl("No", `Age restriction`)]$I,
													 labels = gg.tab[!grepl("No", `Age restriction`)]$level) +
				geom_hline(aes(yintercept = 1), color = 'gray') +
				coord_flip(ylim = {
					if (ggtab.prefix[i] != "hwse") {
						c(0.4, 2.1)} else {c(0, 3.6)}},
					xlim = c(
						min(gg.tab[!grepl("No", `Age restriction`)]$I) + 1,
						max(gg.tab[!grepl("No", `Age restriction`)]$I) - 1)) +
				facet_wrap(. ~ `Age restriction`, scales = "free_y") +
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
																		 	min(gg.tab$I) + 1,
																		 	max(gg.tab$I) - 1)) +
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
			ncol = 3, widths = c(0.35, 0.35, 0.3)
		)
		dev.off()
	})

lualatex(pattern = "*age-restriction\\.tex", directory = directory)

# # Facet view ####
# messy_sol <- 0.05
# year.max <- 1994
# width <- 11
# height <- 8.5
#
# lapply(ggtab.prefix, function(prefix = "sol") {
#
# 	i <- which(ggtab.prefix == prefix)
#
# 	if (is.null(directory)) {
# 		directory <- here::here(paste0("./reports/resources",
# 																	 ifelse(grepl("hwse", prefix), "/hwse 2", ""),
# 																	 ifelse(grepl("hwse", prefix) & year.max != 2015, paste0("/FU through ", year.max), ""),
# 																	 "/Lag ", exposure.lag))
# 	}
#
#
# 	gg.tab <- data.table::copy(get(paste0("og_", prefix, ".ggtab")))
# 	gg.tab <- gg.tab[!grepl(" 50| 55| 60", level)]
# 	gg.tab[grepl("any", level), level := gsub("At any age", "Left work", level)]
#
# 	gg.tab[is.na(as.numeric(ci.lower)), `:=`(
# 		ci.lower = "1.00", ci.upper = "1.00"
# 	)]
#
# 	gg.tab[, `:=`(
# 		Outcome = factor(paste(Outcome, cases),
# 										 levels = unique(paste(Outcome, cases))),
# 		HR = as.numeric(HR),
# 		ci.lower = as.numeric(ci.lower),
# 		ci.upper = as.numeric(ci.upper))]
#
# 	gg.tab[,`:=`(
# 		level.factor = factor(.N:1)
# 	), by = .(Outcome)]
#
# 	tikz(file = paste(
# 		directory,
# 		paste0(ifelse(is.null(file.prefix[i]), paste0(prefix, paste0("_sol", messy_sol %/% .01)), file.prefix[i]), "_facet.tex"), sep = "/"),
# 		standAlone = T, width = width, height = height)
# 	print(ggplot(gg.tab, aes(
# 		x = level.factor,
# 		y = HR,
# 		ymin = ci.lower,
# 		ymax = ci.upper
# 	)) + geom_pointrange(size = 0.15) +
# 		geom_hline(aes(yintercept = 1), color = 'gray') +
# 		geom_text(aes(x = level.factor, y = ifelse(grepl("hwse", prefix), 0, 0.33), label = level), hjust = 1, size = 3) +
# 		coord_flip(
# 			ylim = {if (grepl("hwse", prefix)) {
# 				c(-1.25, 2)
# 			} else {c(-1.225, 2.5)}},
# 			xlim = c(max(gg.tab[,.N, by = .(Outcome)]$N) + 1, 0)
# 		) +
# 		scale_y_continuous(breaks = c(0.5, 1, 1.5, 2, if (!grepl("hwse", prefix)) {2.5})) +
# 		mytheme + labs(y = "") +
# 		facet_wrap(Outcome ~ `Age restriction`, ncol = 2) +
# 		theme(
# 			# plot.margin = unit(c(0.1, 0, 0.1, 0.1), "cm"),
# 			axis.text.y = element_blank(),
# 			axis.ticks.y = element_blank(),
# 			panel.grid = element_blank(),
# 			axis.title.y = element_text(color = "white"),
# 			legend.position = "none"))
# 	dev.off()
# })
#
# # Compile
# lualatex(pattern = "*facet\\.tex", directory = here::here(paste0("./reports/resources/under 65/Lag ", exposure.lag)))
# lualatex(pattern = "*facet\\.tex",
# 				 directory = here::here(paste0(
# 				 	"./reports/under 65/resources/hwse 2/FU through 1994/Lag ", 1 + additional.lag[j],
# 				 	ifelse(employment_status.lag[j] != 0,
# 				 				 paste0("/Employment status lagged ", employment_status.lag[j], " years"),
# 				 				 ""),
# 				 	"/indexed by age")))
