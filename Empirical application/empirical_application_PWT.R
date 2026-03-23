setwd(this.path::here())

library(sparseVARboot)
library(readxl)

data <- read_xlsx("pwt110.xlsx", sheet = 3)

countries <- unique(data$countrycode)
years <- unique(data$year)
real_gdp <- matrix(data$rgdpna, ncol = length(countries), 
                   dimnames = list(year = years, countries = countries))

library(bootUR)

bootUR::plot_missing_values(real_gdp)

real_gdp_clean <- real_gdp[-(1:20), ]
nT <- nrow(real_gdp_clean)
real_gdp_clean <- t(na.omit(t(real_gdp_clean)))

country_names <- unique(data$country)
country_names <- country_names[countries %in% colnames(real_gdp_clean)]
country_key <- cbind(colnames(real_gdp_clean), country_names)

gdp_growth <- 100 * (log(real_gdp_clean[-1, ]) - log(real_gdp_clean[-nT,]))
nN <- ncol(gdp_growth)

plot(colMeans(gdp_growth))

save(gdp_growth, country_key, file = "PWT_gdpgrowth.RData")

load("PWT_gdpgrowth.RData")

source("StepM.R")
names_units <- country_key[, 1]
boot_vec <- c("VAR-MB", "VAR-EB", "BWB", "MBB", "MB", "EB")

set.seed(13682)
test_mu1 <- StepM(x = gdp_growth, names_units = names_units, mu0 = 1, boot = boot_vec)
test_mu2 <- StepM(x = gdp_growth, names_units = names_units, mu0 = 2, boot = boot_vec)

save(test_mu1, test_mu2, file = "PWT_gdp_tests.RData")

load("PWT_gdp_tests.RData")

pval1 <- reshape2::melt(test_mu1$p_val[, -2], id.vars = "unit")
pval2 <- reshape2::melt(test_mu2$p_val[, -2], id.vars = "unit")
colnames(pval2) <- c("Country", "Bootstrap", "Pvalue")
pval2$Country <- factor(rep(1:nN, length(boot_vec)), labels = test_mu2$p_val[, 1])
pval2$index <- as.numeric(pval2$Country)
test_mu2$details$`VAR-MB`$coef_pre


library(ggplot2)

pvalsub <- pval2[pval2$Country %in% pval2$Country[1:20], ]

palette <- c("#097e50", "#66a61e", "#a54043", "#e5694a", "#464b92", "#698fb2")

ggplot(pvalsub, 
       aes(x = Country, y = Pvalue, fill = Bootstrap)) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 0.05, ymax = 0.1), 
            fill = "grey90", alpha = 1) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 0.01, ymax = 0.05), 
            fill = "grey95", alpha = 1) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 0, ymax = 0.01), 
            fill = "grey90", alpha = 1) +
  geom_hline(yintercept = 0.1, colour = "grey70") +
  geom_hline(yintercept = 0.05, colour = "grey70") +
  geom_hline(yintercept = 0.01, colour = "grey70") +
  geom_hline(yintercept = 0, colour = "grey70") +
  geom_col(position = "dodge") +
  scale_fill_manual(values = palette) +
  coord_cartesian(ylim = c(0, 0.12)) +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  scale_y_continuous(name = "p-value", breaks = c(0.01, 0.05, 0.1),
                     expand = expansion(mult = c(0, 0))) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE)) + 
  theme(legend.position = "bottom",
        panel.background = element_rect(fill = "grey95"),
        plot.background =  element_rect(fill = "transparent"),
        panel.grid.major.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_text(size = rel(1.2)),
        axis.title.y = element_text(size = rel(1.2)),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = rel(1.2)),
        legend.key = element_blank())

ggsave("pval_gdp.pdf", width = 8, height = 4)
ggsave("../../Presentations/Figures/pval_gdp.pdf", width = 8, height = 4)

r