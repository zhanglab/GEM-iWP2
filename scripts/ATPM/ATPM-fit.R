# Set to working directory
wd <- commandArgs(trailingOnly = TRUE)
setwd(wd)

df_4 <- read.csv('./04C-ATPM-Robustness.tsv', header=F, sep='\t')
df_15 <- read.csv('./15C-ATPM-Robustness.tsv', header=F, sep='\t')
df_20 <- read.csv('./20C-ATPM-Robustness.tsv', header=F, sep='\t')

# Plot out values of ATP maintenance flux (x) vs biomass flux (y)
# Horizontal lines indicate the derived experimental biomass flux
# Linear region around desired biomass is from ~0 to ~10.
# For 4C the linear fitting will be done on the points from ATP maintenance flux of 1 to 9.
fit_4 <- lm(df_4$V2[df_4$V1 > 1 & df_4$V1 < 9] ~ df_4$V1[df_4$V1 > 1 & df_4$V1 < 9])
# Intercept 0.15350
# Slope -0.01535
# Value of X when y = 0.1108573
atpm_4 = (0.1108573 - fit_4$coefficients[1]) / fit_4$coefficients[2]
# 2.778026

# Linear region around desired biomass is from ~1 to ~4
# For 15C the linear fitting will be done on the points from 1 to 3
fit_15 <- lm(df_15$V2[df_15$V1 > 1 & df_15$V1 < 3] ~ df_15$V1[df_15$V1 > 1 & df_15$V1 < 3])
# Intercept 0.27913
# -0.01833
# Value of X when y = 0.254725
atpm_15 = (0.254725 - fit_15$coefficients[1]) / fit_15$coefficients[2]
# 1.331424

# Linear region around the desired biomass is from ~4 to ~9.5
# For 20C the linear fitting will be done on the points from 5 to 9.5
fit_20 <- lm(df_20$V2[df_20$V1 > 5 & df_20$V1 < 9.5] ~ df_20$V1[df_20$V1 > 5 & df_20$V1 < 9.5])
# Intercept 0.28131
# Slope -0.02001
# Value of X when y = 0.09309851
atpm_20 = (0.09309851 - fit_20$coefficients[1]) / fit_20$coefficients[2]
# 9.406421

df <- rbind(c('4C', atpm_4), c('15C', atpm_15), c('20C', atpm_20))
df


pdf('./S5 Fig.pdf', height=12, width=4)
par(mfrow=c(3,1))
plot(df_4, pch='o', cex=0.25, main='4 Degree ATP Maintenance vs Biomass', xlab='ATP Maintenance Flux', ylab='Biomass Flux', ylim=c(0,0.30))
abline(h = 0.11085728, col='lightblue')
abline(v = 2.78, col='lightblue')
text(x = 3.4, y = 0, '2.78', col='black')
text(x = 0.2, y = 0.119, '0.111', col='black')
clip(1,9, -100, 100)
abline(a=0.15350, b=-0.01535, col='red')

plot(df_15, pch='o', cex=0.25, main='15 Degree ATP Maintenance vs Biomass', xlab='ATP Maintenance Flux', ylab='Biomass Flux', ylim=c(0,0.30))
abline(h = 0.25472496, col='lightblue')
abline(v = 1.33, col='lightblue')
text(x = 1.9, y = 0, '1.33', col='black')
text(x = 0.19, y = 0.245, '0.255', col='black')
clip(1,3, -100, 100)
abline(a= 0.27913, b=-0.01833, col='red')

plot(df_20, pch='o', cex=0.25, main='20 Degree ATP Maintenance vs Biomass', xlab='ATP Maintenance Flux', ylab='Biomass Flux', ylim=c(0,0.30))
abline(h = 0.09309851, col='lightblue')
abline(v = 9.41, col='lightblue')
text(x = 0.19, y = 0.1, '0.093', col='black')
text(x = 8.52, y = 0, '9.41', col='black')
clip(5,9.5, -100, 100)
abline(a=0.28131, b=-0.02001, col='red')
dev.off()
