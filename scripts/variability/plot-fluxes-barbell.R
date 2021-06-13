library(ggplot2)
library(gtable)
library(gridExtra)
library(grid)
library(ggalt)

wd <- commandArgs(trailingOnly = TRUE)
setwd(wd)
######################


gibbs <- read.csv('../TMFA-results/TMFA-Condensed_deltaG_Variability.tsv', sep='\t')

four <- gibbs[,c('ID', 'lower_4', 'upper_4')]
colnames(four) <- c('ID', 'lower', 'upper')
four$temp <- '4°C'
fifteen <-gibbs[,c('ID', 'lower_15', 'upper_15')]
colnames(fifteen) <- c('ID', 'lower', 'upper')
fifteen$temp <- '15°C'
twenty <- gibbs[,c('ID', 'lower_20', 'upper_20')]
colnames(twenty) <- c('ID', 'lower', 'upper')
twenty$temp <- '20°C'


gibbs_combined <- rbind(four, fifteen, twenty)


x <- c('ENO', 'TPI', 'GAPD', 'PGK', 'PGM', 'PYK', 'PDH', 'G6PDHy', 'PGDH', 'EDA', 'DRPA', 'PPM2', 'ICDHy', 'ICL', 'FUM', 'ME2', 'PTAr')

gibbs_combined <- gibbs_combined[gibbs_combined$ID %in% x,]


gibbs_combined$color <- '#FF0000'
gibbs_combined$color[gibbs_combined$temp == '15°C'] <- '#bdbf11'
gibbs_combined$color[gibbs_combined$temp == '4°C'] <- '#002EFF'

gibbs_combined$temp <- factor(gibbs_combined$temp, levels = c('4°C', '15°C', '20°C'))
gibbs_combined$ID <- factor(gibbs_combined$ID, levels=c('ENO', 'TPI', 'GAPD', 'PGK', 'PGM', 'PYK', 'PDH', 'G6PDHy', 'PGDH', 'EDA', 'DRPA', 'PPM2', 'ICDHy', 'ICL', 'FUM', 'ME2', 'PTAr'))


g <- ggplot(gibbs_combined, aes(x=lower, xend=upper, y=temp))+
  geom_dumbbell(color=gibbs_combined$color, size_x=3, size_xend = 3, size=2)+
  facet_wrap(ID~., scales='free_y', strip.position='top', ncol = c(3), dir = 'v')+
  labs(y=NULL, x='Gibbs Free Energy (kJ/mol)')+theme_bw()+
  theme(strip.text.x = element_text(size = 14, angle = 0), axis.text.y = element_text(size=14), axis.text.x = element_text(size=12), axis.title = element_text(size=18))+geom_vline(xintercept = 0, color='black')

ggsave(filename = 'S2 Fig.pdf', height=8, width=16)




####################

metabolites <- read.csv('../TMFA-results/TMFA-Condensed_Compund_Variability.tsv', sep='\t')

four <- metabolites[,c('compound', 'lower_4', 'upper_4')]
colnames(four) <- c('ID', 'lower', 'upper')
four$temp <- '4°C'
fifteen <-metabolites[,c('compound', 'lower_15', 'upper_15')]
colnames(fifteen) <- c('ID', 'lower', 'upper')
fifteen$temp <- '15°C'
twenty <- metabolites[,c('compound', 'lower_20', 'upper_20')]
colnames(twenty) <- c('ID', 'lower', 'upper')
twenty$temp <- '20°C'


metab_combined <- rbind(four, fifteen, twenty)




metab_combined <- metab_combined[metab_combined$ID %in% c('cpd_atp[c]', 'cpd_adp[c]', 'cpd_pi[c]', 'cpd_nadh[c]', 'cpd_nad[c]', 'cpd_h2o2[c]', 'cpd_g3p[c]', 'cpd_13dpg[c]', 'cpd_2pg[c]', 'cpd_3pg[c]', 'cpd_pep[c]', 'cpd_r1p[c]', 'cpd_xu5p-D[c]', 'cpd_fum[c]', 'cpd_icit[c]', 'cpd_mal-L[c]', 'cpd_akg[c]', 'cpd_methf[c]', 'cpd_thf[c]'),]
metab_combined$color <- '#FF0000'
metab_combined$color[metab_combined$temp == '15°C'] <- '#bdbf11'
metab_combined$color[metab_combined$temp == '4°C'] <- '#002EFF'

metab_combined$temp <- factor(metab_combined$temp, levels = c('4°C', '15°C', '20°C'))

metab_combined$ID <- factor(metab_combined$ID, levels=c('cpd_g3p[c]', 'cpd_13dpg[c]', 'cpd_2pg[c]', 'cpd_3pg[c]', 'cpd_pep[c]', 'cpd_r1p[c]', 'cpd_xu5p-D[c]', 'cpd_fum[c]', 'cpd_icit[c]', 'cpd_mal-L[c]', 'cpd_akg[c]', 'cpd_methf[c]', 'cpd_thf[c]', 'cpd_atp[c]', 'cpd_adp[c]', 'cpd_pi[c]', 'cpd_nadh[c]', 'cpd_nad[c]', 'cpd_h2o2[c]'),
                            labels=c('Glyceraldehyde 3-phosphate', '3-Phospho-D-glyceroyl phosphate', 'D-Glycerate 2-phosphate', '3-Phospho-D-glycerate', 'Phosphoenolpyruvate', 'alpha-D-Ribose 1-phosphate', 'D-Xylulose 5-phosphate', 'Fumarate', 'Isocitrate', 'Malate', '2-Oxoglutarate', '5,10-Methenyltetrahydrofolate', '5,6,7,8-Tetrahydrofolate', 'ATP', 'ADP', 'Phosphate', 'NADH', 'NAD+', 'Hydrogen Peroxide'))

g <- ggplot(metab_combined, aes(x=lower, xend=upper, y=temp))+
  geom_dumbbell(color=metab_combined$color, size_x=3, size_xend = 3, size=2)+
  facet_wrap(ID~., scales='free_y', ncol = c(3), dir='v')+
  labs(y=NULL, x='Concentration (M)')+theme_bw()+
  theme(strip.text.x = element_text(size = 14, angle = 0), axis.text.y = element_text(size=14), axis.text.x = element_text(size=14), axis.title = element_text(size=18))

ggsave(filename = 'S1 Fig.pdf', height=8, width = 16)
