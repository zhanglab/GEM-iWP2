wd <- commandArgs(trailingOnly = TRUE)
setwd(wd)

library(ggplot2)
library(ggpubr)

df <- read.csv('../TMFA-results/TMFA-Condensed_Flux_Variability.tsv', sep='\t', header=T)
df <- df[df$range_4 != 0 | df$range_15 != 0 | df$range_20 != 0,]

four <- as.data.frame(df[c('ID', 'range_4')])
fifteen <- as.data.frame(df[c('ID', 'range_15')])
twenty <- as.data.frame(df[c('ID', 'range_20')])
colnames(four) <- c('id','flux')
colnames(fifteen) <- c('id','flux')
colnames(twenty) <- c('id','flux')
four$flux <- four$flux / 0.11082563
fifteen$flux <- fifteen$flux / 0.25474175
twenty$flux <- twenty$flux / 0.0930269
four$temp <- '4'
fifteen$temp <- '15'
twenty$temp <- '20'
comb <- rbind(four, fifteen, twenty)
comb$temp <- factor(comb$temp, levels = c('4', '15', '20'))
comb$flux <- as.numeric(comb$flux)
temp_cols <- c('Blue', 'Yellow', 'Red')

comb$id <- NULL
g <- ggplot(comb, aes(x=temp, y=flux)) + geom_boxplot(varwidth=F, fill=temp_cols) + labs(title=NULL)+
  scale_y_continuous(limits = c(0,32))+theme_bw(base_size = 16)+
  xlab('Temperature (°C)')+ylab('Flux Range \n (Maximum Normalized Flux - Minimum Normalized Flux)')+
  stat_compare_means(comparisons = list(c('4', '15'), c('15', '20'), c('4', '20')), label.y = c(22, 24, 30), method ='wilcox.test')

ggsave(filename = './S3 Fig.pdf')
