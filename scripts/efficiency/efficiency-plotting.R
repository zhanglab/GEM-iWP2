wd <- commandArgs(trailingOnly = TRUE)
setwd(wd)

df <- read.csv('../TMFA-results/TMFA-Condensed_Flux_Variability.tsv', header=T, sep='\t')
car_ex <- read.csv('./carbon-exchange', header=F, sep='\t')
colnames(car_ex) <- c('ID', 'count')

ex_df <- df[startsWith(as.character(df$ID), 'EX_'),]
sub_ex_df <- ex_df[,c(2,3,4,6,7,9,10)]

sub_ex_carbnorm <- merge(car_ex, sub_ex_df, by='ID', all.y=F, all.x=T)
sub_ex_carbnorm$lower_4 <- sub_ex_carbnorm$lower_4 * sub_ex_carbnorm$count
sub_ex_carbnorm$lower_15 <- sub_ex_carbnorm$lower_15 * sub_ex_carbnorm$count
sub_ex_carbnorm$lower_20 <- sub_ex_carbnorm$lower_20 * sub_ex_carbnorm$count


four_uptake <- abs(sum(sub_ex_carbnorm$lower_4[sub_ex_carbnorm$lower_4 < 0]))
four_production <- abs(sum(sub_ex_carbnorm$lower_4[sub_ex_carbnorm$lower_4 > 0]))

fifteen_uptake <- abs(sum(sub_ex_carbnorm$lower_15[sub_ex_carbnorm$lower_15 < 0]))
fifteen_production <- abs(sum(sub_ex_carbnorm$lower_15[sub_ex_carbnorm$lower_15 > 0]))

twenty_uptake <- abs(sum(sub_ex_carbnorm$lower_20[sub_ex_carbnorm$lower_20 < 0]))
twenty_production <- abs(sum(sub_ex_carbnorm$lower_20[sub_ex_carbnorm$lower_20 > 0]))

four_cue <- (four_uptake - four_production) / four_uptake
fifteen_cue <- (fifteen_uptake - fifteen_production) / fifteen_uptake
twenty_cue <- (twenty_uptake - twenty_production) / twenty_uptake

pdf('./Fig 2.pdf', height=5, width=16)
par(mfrow=c(1,3))
barplot(height = c(four_cue, fifteen_cue, twenty_cue),
        names.arg = c('4', '15', '20'),
        ylab = 'Carbon use efficiency (CUE)',
        xlab='Temperature (°C)', col=c('gray', 'gray', 'gray'), ylim = c(0,0.3), cex.lab=1.5, cex.axis=1.5, cex.names=1.5)

df <- read.csv('../TMFA-results/TMFA-Condensed_absMin_Variability.tsv', header=T, sep='\t')

df <- df[,c(2,3,4,6,7,9,10)]
atpm <- df[df$ID == 'ATPM',]
atpm <- c(atpm$upper_4, atpm$upper_15, atpm$upper_20)
atp <- read.csv('./atp-rxns', header=F, sep='\t')

atp_fluxes <- merge(atp, df, by.x = 1, by.y =1)
atp_fluxes$lower_4[is.na(atp_fluxes$lower_4)] <- 0
atp_fluxes$upper_4[is.na(atp_fluxes$upper_4)] <- 0
atp_fluxes$lower_15[is.na(atp_fluxes$lower_15)] <- 0
atp_fluxes$upper_15[is.na(atp_fluxes$upper_15)] <- 0
atp_fluxes$lower_20[is.na(atp_fluxes$lower_20)] <- 0
atp_fluxes$upper_20[is.na(atp_fluxes$upper_20)] <- 0

atp_fluxes$lower_4 <- atp_fluxes$lower_4  * atp_fluxes$V3
atp_fluxes$upper_4 <- atp_fluxes$upper_4  * atp_fluxes$V3

atp_fluxes$lower_15 <- atp_fluxes$lower_15  * atp_fluxes$V3
atp_fluxes$upper_15 <- atp_fluxes$upper_15  * atp_fluxes$V3

atp_fluxes$lower_20 <- atp_fluxes$lower_20  * atp_fluxes$V3
atp_fluxes$upper_20 <- atp_fluxes$upper_20  * atp_fluxes$V3

four_min_prod <- sum(atp_fluxes$lower_4[atp_fluxes$lower_4 > 0])
four_min_con <- abs(sum(atp_fluxes$lower_4[atp_fluxes$lower_4 < 0]))
four_max_prod <- sum(atp_fluxes$upper_4[atp_fluxes$upper_4 > 0])
four_max_con <- abs(sum(atp_fluxes$upper_4[atp_fluxes$upper_4 < 0]))

fifteen_min_prod <- sum(atp_fluxes$lower_15[atp_fluxes$lower_15 > 0])
fifteen_min_con <- abs(sum(atp_fluxes$lower_15[atp_fluxes$lower_15 < 0]))
fifteen_max_prod <- sum(atp_fluxes$upper_15[atp_fluxes$upper_15 > 0])
fifteen_max_con <- abs(sum(atp_fluxes$upper_15[atp_fluxes$upper_15 < 0]))

twenty_min_prod <- sum(atp_fluxes$lower_20[atp_fluxes$lower_20 > 0])
twenty_min_con <- abs(sum(atp_fluxes$lower_20[atp_fluxes$lower_20 < 0]))
twenty_max_prod <- sum(atp_fluxes$upper_20[atp_fluxes$upper_20 > 0])
twenty_max_con <- abs(sum(atp_fluxes$upper_20[atp_fluxes$upper_20 < 0]))

x <- rbind(c(four_max_prod, fifteen_max_prod, twenty_max_prod), c(four_min_prod, fifteen_min_prod, twenty_min_prod))

atp_per_carbon <- x


atp_vals <- t(atp_per_carbon)
atp_vals <- as.data.frame(cbind(c('4', '15', '20'), atp_vals))
colnames(atp_vals) <- c('temp', 'lower', 'upper')
atp_vals$lower <- as.numeric(atp_vals$lower)
atp_vals$upper <- as.numeric(atp_vals$upper)
atp_vals$temp <- factor(atp_vals$temp, levels = c('4', '15', '20'), ordered = T)

atp_vals$lower
carb <- df[df$ID == 'EX_cpd_acgam[e]',]
carb <- c(carb$upper_4, carb$upper_15, carb$upper_20)*-1

atp_vals$production <- (atp_vals$lower) / carb


(carb) * c(0.11071657, 0.25472549, 0.09337582)

barplot(atp_vals$production, names.arg = atp_vals$temp, ylim = c(0, 3.5),
        xlab='Temperature (°C)', ylab='ATP Produced / GlcNac', col=c('gray', 'gray', 'gray'), cex.lab=1.5, cex.axis=1.5, cex.names=1.5)

atp_vals$effic <- atp_vals$lower - atpm


barplot(atp_vals$effic, names.arg = atp_vals$temp, ylim = c(0, 100),
        xlab='Temperature (°C)', ylab='ATP Consumed / gDW Biomass', col=c('gray', 'gray', 'gray'), cex.lab=1.5, cex.axis=1.5, cex.names=1.5)
dev.off()
