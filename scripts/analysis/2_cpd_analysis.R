library("DESeq2")
library("rcompanion")
library("stringr")
library("ggplot2")
library("ggdendro")
library("ggExtra")
library("grid")
library("gridExtra")
library("tidyr")
library("tidyverse")
library("plyr")
library("scales")


ccrand4<-'../TMFA-results/cpd_random_04C.tsv'
ccrand15<-'../TMFA-results/cpd_random_15C.tsv'
ccrand20<-'../TMFA-results/cpd_random_20C.tsv'

############ Read compound concentrations from random simulation #############
#-read cc
ccrand4<-read.table(ccrand4,header=FALSE,sep="\t",quote="\"",fill=FALSE,
                    row.names='V1')
ccrand15<-read.table(ccrand15,header=FALSE,sep="\t",quote="\"",fill=FALSE, 
                     row.names='V1')
ccrand20<-read.table(ccrand20,header=FALSE,sep="\t",quote="\"",fill=FALSE,
                     row.names='V1')

#-round numbers to 0.001
ccrand4<-round(ccrand4,digits=3)
ccrand15<-round(ccrand15,digits=3)
ccrand20<-round(ccrand20,digits=3)

cpdlist<-intersect(rownames(ccrand4),
                   intersect(rownames(ccrand15),rownames(ccrand20)))

############ initiate summary stats from random simulations #############
#-compound concentrations
cpdall<-data.frame(matrix(nrow=length(cpdlist),ncol=0))
rownames(cpdall)<-cpdall$cpdid<-cpdlist
#--4C
cpdall$t4median<-apply(ccrand4,1,median)
cpdall$t4q25<-apply(ccrand4,1,quantile,probs=0.25)
cpdall$t4q75<-apply(ccrand4,1,quantile,probs=0.75)
cpdall$t4min<-apply(ccrand4,1,min)
cpdall$t4max<-apply(ccrand4,1,max)
#--15C
cpdall$t15median<-apply(ccrand15,1,median)
cpdall$t15q25<-apply(ccrand15,1,quantile,probs=0.25)
cpdall$t15q75<-apply(ccrand15,1,quantile,probs=0.75)
cpdall$t15min<-apply(ccrand15,1,min)
cpdall$t15max<-apply(ccrand15,1,max)
#--20C
cpdall$t20median<-apply(ccrand20,1,median)
cpdall$t20q25<-apply(ccrand20,1,quantile,probs=0.25)
cpdall$t20q75<-apply(ccrand20,1,quantile,probs=0.75)
cpdall$t20min<-apply(ccrand20,1,min)
cpdall$t20max<-apply(ccrand20,1,max)
#--write to output
write.table(cpdall,file=paste("cpdall_summary.tsv",sep=''),quote=FALSE,
            row.names=FALSE,col.names=TRUE,sep="\t")

######## define cpd.sig, significantly differ across temperature ########
#-cpd.sig record the significance of the Kruskal-Wallis test on concentration
cpd.sig<-data.frame(matrix(nrow=0,ncol=6))
colnames(cpd.sig)<-c("cpdid","p_kw","es","p_4v15","p_4v20","p_20v15")
#--perform kwt and posthoc tests for each compound
for(cid in cpdlist){
  #-temporary data frame for storing the melt data of concentrations in 3 temp
  cf<-data.frame(concentration=c(as.numeric(ccrand4[cid,]),
                                 as.numeric(ccrand15[cid,]),
                                 as.numeric(ccrand20[cid,])),
                 temperature=c(rep("T4",ncol(ccrand4)),
                               rep("T15",ncol(ccrand15)),
                               rep("T20",ncol(ccrand20))))
  
  #-Kruskal-Wallis test
  kw<-kruskal.test(concentration~temperature, data=cf)
  if(is.na(kw$p.value)){ kw$p.value<-1000 }
  # if(kw$p.value<=0.05){
  rowid=nrow(cpd.sig)+1
  cpd.sig[rowid,]$cpdid<-cid
  cpd.sig[rowid,]$p_kw<-kw$p.value
  #-calculate effect size, epsilon-squared > 0.64 suggest very strong effect
  #-https://peterstatistics.com/CrashCourse/3-TwoVarUnpair/NomOrd/NomOrd3c.html
  cpd.sig[rowid,]$es<-epsilonSquared(x = cf$concentration, g = cf$temperature)
  #-posthoc test
  pwt<-pairwise.wilcox.test(cf$concentration,g=cf$temperature,
                            p.adjust.method = "bonferroni")
  cpd.sig[rowid,]$p_4v15<-pwt$p.value["T4","T15"]
  cpd.sig[rowid,]$p_4v20<-pwt$p.value["T4","T20"]
  cpd.sig[rowid,]$p_20v15<-pwt$p.value["T20","T15"]
  # }
}
cpdoverall <- merge(cpdall, cpd.sig, by=c(1,1), all.x=T, all.y=T)
write.table(cpdoverall, file='./cpd_all_with_sig.tsv', sep='\t', row.names=F, col.names=T)
#-rxns with significant concentration variance based on Kruskal-Wallis test: 425
write.table(cpd.sig,file=paste("cpd_kwt_sig_r3.tsv",sep=''),quote=FALSE,
            row.names=FALSE,col.names=TRUE,sep="\t")

######### filter cpd.sig by posthoc analysis #########################
#-The threshold for effect size, epsilon-squared was based on below:
#-https://peterstatistics.com/CrashCourse/3-TwoVarUnpair/NomOrd/NomOrd3c.html
#--0.00 < 0.01 - Negligible
#--0.01 < 0.04 - Weak
#--0.04 < 0.16 - Moderate
#--0.16 < 0.36 - Relatively strong
#--0.36 < 0.64 - Strong
#--0.64 < 1.00 - Very strong
#-remove entries in cpd.sig if failed effect size (>=0.36) or posthoc test
tmplist<-c()
for(r in rownames(cpd.sig)){
  if(cpd.sig[r,]$es < 0.36){ tmplist[length(tmplist)+1]=r }
  if(is.na(cpd.sig[r,]$p_4v15) | is.na(cpd.sig[r,]$p_4v20) | 
     is.na(cpd.sig[r,]$p_20v15)){ tmplist[length(tmplist)+1]=r }
  else if(cpd.sig[r,]$p_4v15 > 0.05 & cpd.sig[r,]$p_4v20 > 0.05 & 
          cpd.sig[r,]$p_20v15 > 0.05){ tmplist[length(tmplist)+1]=r }
}
#-rxns with significant flux variance with ES and posthoc test: 302
cpd.sig.ph<-cpd.sig[!(rownames(cpd.sig) %in% tmplist),]$cpdid
write.table(cpd.sig[cpd.sig$cpdid %in% cpd.sig.ph,],
            file=paste("cpd_kwt_sig_r3_posthoc.tsv",sep=''),quote=FALSE,
            row.names=FALSE,col.names=TRUE,sep="\t")

######### Data prep for visualization #########################
#-initiate column map
#--"median": median of 1000 random simulation
#--"q25": 25th quantile of 1000 random simulation
#--"q75": 75th quantile of 1000 random simulation
#--"min": minimum of 1000 random simulation
#--"max": maximum of 1000 random simulation
#--"label": label for each temperature to be used in plotting
#--"color": color for each temperature to be used in plotting
colmap<-data.frame(median=c("t4median","t15median","t20median"),
                   q25=c("t4q25","t15q25","t20q25"),
                   q75=c("t4q75","t15q75","t20q75"),
                   min=c("t4min","t15min","t20min"),
                   max=c("t4max","t15max","t20max"),
                   label=c("4C","15C","20C"),
                   color=c("blue","gold","firebrick1"))
rownames(colmap)<-c("T4","T15","T20")

#-prepare data for scaling
cpd.scale<-cpdall[cpd.sig$cpdid,c("t4median", "t15median", "t20median")]
# for(i in rownames(cpd.scale)){
#   cpd.scale[i,]<-scale(t(cpd.scale[i,colmap$median]),center=FALSE,scale=TRUE)
# }
cpd.scale<-exp(cpd.scale)
cpd.scale$cpdid<-cpd.sig$cpdid
cpd.scale<-cpd.scale[!rownames(cpd.scale) %in% c("cpd_cl[c]"),]
colnames(cpd.scale)<-c("T4","T15","T20")

#-cc is the core data for visualization of gene expression
#--make cc based on scaled values
cc<-pivot_longer(data=cpd.scale,
                 cols=c(T4, T15, T20),
                 names_to = "temperature",
                 values_to = "concentration")

#-focus on compounds in [c]
cc<-cc[!str_detect(cc$cpdid,"\\[e\\]|\\[p\\]"),]

######### Visualization: histogram #####################
#-Define a function for plotting pairs of columns in fvgv 
#--data: the fvgv data frame
#--t: definition of a group value, should be retrieved from the column 'grp'
#--grp: the column name for matching the variable t
#--x: the column name used for the x-axis
#--y: the column name used for the y-axis
#--xlab: the label used for the x-axis
#--ylab: the label used for the y-axis
#--colmap: a dataframe that stores the mappign of grp to colors,
#          colmap should at least contain columns named "color" and "label"
pfunc<-function(data,x,y,xlabel,ylabel){
  ggplot(data[!is.na(data[,x] & !is.na(data[,y])),],
         aes_string(x=x, y=y)) +
    geom_point(size = 1,position=position_jitter(w=0.05,h=0)) +
    # xlim(-0.1,1.55) +
    # ylim(-0.1,1.55) +
    # scale_colour_manual(values=t(colmap)["color",]) +
    labs(title="", x=xlabel, y=ylabel) +
    theme(legend.position="none", title = element_text(colour = "#57585a"),
          legend.text = element_text(colour="#57585a", size = 9))
}

#-Define labels corresponding to column names
lablist<-data.frame(name=c("T4","T15","T20"),
                    label=c("4C","15C","20C"))
rownames(lablist)<-lablist$name

#-main figure: Plot pairwise relationship between deltaG, flux, expression, 
#-with all growth phases combined in the expression
#--ax: a list of all columns for pairwise comparison
#---main figure, plot 15v4 and 15v20
ax<-c("T15","T4","T20")
outpdf<-"cpd_concentration_all.pdf"
#---only use first index as x-axis, others as y-axis
pcnt<-length(ax)-2
#--pi: the index of the plot serial
pi<-0 
plots<-c()
for(i in 1:2){
  xs<-ax[i]
  for(j in (i+1):length(ax)){
    ys<-ax[j]
    pm<-pfunc(cpd.scale,xs,ys,lablist[xs,]$label,
                lablist[ys,]$label)
    pi<-pi+1
    plots[[pi]]<-ggMarginal(pm, type="boxplot")
  }
}
pdf(outpdf,height=4*pcnt, width=8)
grid.arrange(grobs=plots,nrow=pcnt)
dev.off()

######### visualization 1: Plotting individual copds or ratios ###################
#-Define cpds of interest
#--full list of compounds and their names in original Fig S1
coi<-data.frame(cpdid=c("cpd_gam6p[c]","cpd_gam1p[c]",
                        "cpd_g3p[c]","cpd_pep[c]",
                        "cpd_e4p[c]","cpd_s7p[c]",
                        "cpd_dxyl5p[c]","cpd_etoh[c]",
                        "cpd_mal-L[c]","cpd_succ[c]","cpd_amet[c]"),
                cpdname=c("D-Glucosamine 6-phosphate","D-Glucosamine 1-phosphate",
                          "Glyceraldehyde 3−phosphate","Phosphoenolpyruvate",
                          "Erythrose 4-phosphate","Sedoheptulose 7-phosphate",
                          "1-deoxy-D-xylulose 5-phosphate","Ethanol",
                          "L-Malate","Succinate","S-Adenosyl-L-methionine"))
#-making plots per compound
outpdf<-"cpd_concentration_examples.pdf"
#---set the total number of rows
pcnt<-3
#--initiate the plot list
plots=c()
for(i in 1:length(coi$cpdid)){
  c<-coi[i,]$cpdid
  #-temporary data fram for storing the melt data of CC in 3 temperature
  rf<-data.frame(concentration=c(exp(as.numeric(ccrand4[c,])),
                        exp(as.numeric(ccrand15[c,])),
                        exp(as.numeric(ccrand20[c,]))),
                 temperature=c(rep("T4",ncol(ccrand4)),
                               rep("T15",ncol(ccrand15)),
                               rep("T20",ncol(ccrand20))),
                 cpdid=rep(c,ncol(ccrand4)+ncol(ccrand15)+ncol(ccrand20)))
  
  #-plot compound concentration
  plots[[i]]<-ggplot(rf,aes(x=factor(temperature,level = c('T4','T15','T20')), 
                    y=concentration, fill=temperature)) +
    geom_boxplot() +
    scale_fill_manual(values=t(colmap)["color",]) +
    labs(title="", x="", y="Concentration (M/L)", fill="Temperature") +
    theme(legend.position="none") +
    facet_wrap( ~ cpdid,nrow=1)
}
pdf(outpdf,height=4*pcnt, width=20)
grid.arrange(grobs=plots,nrow=pcnt)
dev.off()

######### visualization 2: Plotting compound ratios ###################
#-Define cpds of interest
#--full list of compounds to be used in the calculation of CC ratios
coi_ratio<-data.frame(divident=c("cpd_nadh[c]","cpd_nadph[c]","cpd_atp[c]"),
                      divisor=c("cpd_nad[c]","cpd_nadp[c]","cpd_adp[c]"),
                      name=c("NADH/NAD+","NADPH/NADP+","ATP/ADP"))
#-making plots per ratio
outpdf<-"Figure4.pdf"
#---set the total number of rows
pcnt<-1
#--initiate the plot list
plots=c()
for(i in 1:length(coi_ratio$divident)){
  dt<-coi_ratio[i,]$divident
  dr<-coi_ratio[i,]$divisor
  #-temporary data fram for storing the melt data of CC in 3 temperature
  rf<-data.frame(ratio=c(t(exp(ccrand4[dt,])/exp(ccrand4[dr,])),
                            t(exp(ccrand15[dt,])/exp(ccrand15[dr,])),
                            t(exp(ccrand20[dt,])/exp(ccrand20[dr,]))),
                 temperature=c(rep("T4",ncol(ccrand4)),
                               rep("T15",ncol(ccrand15)),
                               rep("T20",ncol(ccrand20))),
                 name=rep(coi_ratio[i,]$name,
                          ncol(ccrand4)+ncol(ccrand15)+ncol(ccrand20)))
  
  #-plot ratio of metabolite concentrations
  plots[[i]]<-ggplot(rf,aes(x=factor(temperature,level = c('T4','T15','T20'),
                                     labels = c('4°C', '15°C', '20°C')), 
                            y=ratio, fill=temperature)) +
    geom_boxplot() +
    scale_y_log10(labels = comma) +
    scale_fill_manual(values=t(colmap)["color",]) +
    labs(title="", x="", y="Concentration Ratio", fill="Temperature") +
    theme(legend.position="none") +
    facet_wrap( ~ name,nrow=1)
}
pdf(outpdf,height=4, width=10)
grid.arrange(grobs=plots,nrow=pcnt)
dev.off()


