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

dgrand4<-'../TMFA-results/dg_random_04C.tsv'
dgrand15<-'../TMFA-results/dg_random_15C.tsv'
dgrand20<-'../TMFA-results/dg_random_20C.tsv'
fluxrand4<-'../TMFA-results/flux_random_04C.tsv'
fluxrand15<-'../TMFA-results/flux_random_15C.tsv'
fluxrand20<-'../TMFA-results/flux_random_20C.tsv'
rtgf<-'../TMFA-results/WP2-gene-associations.tsv'
ctsf<-'../TMFA-results/full_count_table_summary.csv'
ctsannof<-'../TMFA-results/full_count_table_annotation.tsv'
def_early<-"../TMFA-results/de_pairwise_early.tsv"
def_late<-"../TMFA-results/de_pairwise_late.tsv"       
def_stat<-"../TMFA-results/de_pairwise_stationary.tsv"

#-read rxn to gene mapping
rtg<-read.table(rtgf,header=TRUE,sep="\t",quote="\"",fill=FALSE)
rownames(rtg)<-rtg$id
for(i in rownames(rtg)){
  rtg[i,]$genes<-str_replace_all(rtg[i,]$genes,"\\(|\\)","")
  rtg[i,]$genes<-strsplit(as.character(rtg[i,]$genes)," and | or ")
}

#-create gene to reaction mapping
rtg.guniq<-unique(unlist(rtg$genes))
gtr<-data.frame(matrix(nrow=length(rtg.guniq),ncol=2))
colnames(gtr)<-c('gid','rxns')
rownames(gtr)<-gtr$gid<-rtg.guniq
for(r in rtg$id){
  for(g in rtg[r,]$genes){
    gtr[g,]$rxns<-list(unique(c(unlist(gtr[g,]$rxns),r)))
  }
}
for(g in rownames(gtr)){
  tmp=unlist(gtr[g,]$rxns)
  gtr[g,]$rxns<-list(tmp[!is.na(tmp)])
}

############ Initiating gene expression values #############
#-read raw count data
cts<-as.matrix(read.table(ctsf,header=TRUE,sep=",",quote="\"",fill=FALSE,
                          row.names="gene_id"))
cts<-cts[,colnames(cts)!="X"]
coldata<-read.table(ctsannof,header=TRUE,sep="\t",quote="\"",fill=FALSE,
                    row.names="sample")
coldata$condition <- factor(coldata$condition)
coldata$growthphase <- factor(coldata$growthphase)
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition + growthphase)
#--median of ratios, 
#--https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html
#--https://scienceparkstudygroup.github.io/research-data-management-lesson/median_of_ratios_manual_normalization/index.html
dds<-estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)
cts.mrn<-counts(dds, normalized = TRUE)

#-read pairwise DE from exponential phase (Analysis for Data S2)
de_late<-read.table(def_late,header=TRUE,sep="\t",quote="\"",fill=FALSE)
colnames(de_late)<-c("gene", "lfc4v15","padj4v15","lfc15v20","padj15v20",
                      "lfc4v20","padj4v20","gname","pathway")
rownames(de_late)<-de_late$gene
deglist<-de_late[!de_late$padj4v15>0.05 | !de_late$padj4v20>0.05 | 
           !de_late$padj15v20>0.05,]$gene

############ Read flux and deltaG from random simulation #############
#-read deltaG
dgrand4<-read.table(dgrand4,header=FALSE,sep="\t",quote="\"",fill=FALSE,
                    row.names='V1')
dgrand15<-read.table(dgrand15,header=FALSE,sep="\t",quote="\"",fill=FALSE, 
                     row.names='V1')
dgrand20<-read.table(dgrand20,header=FALSE,sep="\t",quote="\"",fill=FALSE,
                     row.names='V1')
#-read flux
frand4<-read.table(fluxrand4,header=FALSE,sep="\t",quote="\"",fill=FALSE,
                   row.names='V1')
frand15<-read.table(fluxrand15,header=FALSE,sep="\t",quote="\"",fill=FALSE, 
                    row.names='V1')
frand20<-read.table(fluxrand20,header=FALSE,sep="\t",quote="\"",fill=FALSE,
                    row.names='V1')

#-round numbers to 0.001
dgrand4<-round(dgrand4,digits=3)
dgrand15<-round(dgrand15,digits=3)
dgrand20<-round(dgrand20,digits=3)
frand4<-round(frand4,digits=3)
frand15<-round(frand15,digits=3)
frand20<-round(frand20,digits=3)

rxnlist<-intersect(intersect(rownames(frand4),rownames(frand15)),
                   rownames(frand20))

############ initiate summary stats from random simulations #############
#-flux
fluxall<-data.frame(matrix(nrow=length(rxnlist),ncol=0))
rownames(fluxall)<-fluxall$rxnid<-rxnlist
#--4C
fluxall$t4median<-apply(frand4,1,median)
fluxall$t4q25<-apply(frand4,1,quantile,probs=0.25)
fluxall$t4q75<-apply(frand4,1,quantile,probs=0.75)
fluxall$t4min<-apply(frand4,1,min)
fluxall$t4max<-apply(frand4,1,max)
#--15C
fluxall$t15median<-apply(frand15,1,median)
fluxall$t15q25<-apply(frand15,1,quantile,probs=0.25)
fluxall$t15q75<-apply(frand15,1,quantile,probs=0.75)
fluxall$t15min<-apply(frand15,1,min)
fluxall$t15max<-apply(frand15,1,max)
#--20C
fluxall$t20median<-apply(frand20,1,median)
fluxall$t20q25<-apply(frand20,1,quantile,probs=0.25)
fluxall$t20q75<-apply(frand20,1,quantile,probs=0.75)
fluxall$t20min<-apply(frand20,1,min)
fluxall$t20max<-apply(frand20,1,max)
#--write to output
write.table(fluxall,file=paste("fluxall_summary.tsv",sep=''),quote=FALSE,
            row.names=FALSE,col.names=TRUE,sep="\t")

#-deltaG
dgall<-data.frame(matrix(nrow=length(rxnlist),ncol=0))
rownames(dgall)<-dgall$rxnid<-rxnlist
#--4C
dgall$t4median<-apply(dgrand4,1,median)
dgall$t4q25<-apply(dgrand4,1,quantile,probs=0.25)
dgall$t4q75<-apply(dgrand4,1,quantile,probs=0.75)
dgall$t4min<-apply(dgrand4,1,min)
dgall$t4max<-apply(dgrand4,1,max)
#--15C
dgall$t15median<-apply(dgrand15,1,median)
dgall$t15q25<-apply(dgrand15,1,quantile,probs=0.25)
dgall$t15q75<-apply(dgrand15,1,quantile,probs=0.75)
dgall$t15min<-apply(dgrand15,1,min)
dgall$t15max<-apply(dgrand15,1,max)
#--20C
dgall$t20median<-apply(dgrand20,1,median)
dgall$t20q25<-apply(dgrand20,1,quantile,probs=0.25)
dgall$t20q75<-apply(dgrand20,1,quantile,probs=0.75)
dgall$t20min<-apply(dgrand20,1,min)
dgall$t20max<-apply(dgrand20,1,max)
#--write to output
write.table(dgall,file=paste("dgall_summary.tsv",sep=''),quote=FALSE,
            row.names=FALSE,col.names=TRUE,sep="\t")

######## overview of all flux and deltaG ################
#-Fluxes
#--get rows with non-zero values
#--rxns showing non-zero in at least one temperature: 371
fnozero<-fluxall[fluxall$t4median!=0 | fluxall$t15median!=0 | 
                   fluxall$t20median!=0,]$rxnid
#--get rows where flux differ between at least two temperatures
#---based on the data, fv completely overlaps with fnozero, 
#---so we will use fv for futher analysis/plotting
#--rxns showing variability between at least two temperatures: 366
fv<-fluxall[fluxall$t4median!=fluxall$t15median | 
              fluxall$t20median!=fluxall$t15median | 
              fluxall$t4median!=fluxall$t20median,]$rxnid

#-deltaG
#--Initially trying to do a parallel analysis with deltaG, 
#--but meaning of zeros in deltaG is not clear, so skip nonzero analysis for dg
#--deltaG of non-zero in at least one temperature: 954
# dgnozero<-dgall[dgall$t4median!=0 | dgall$t15median!=0 | 
#                    dgall$t20median!=0,]$rxnid
#--gv is a much bigger list, showing more variable deltaG but not in flux
#--deltaG of showing variability between at least two temperatures: 700
#--overlap between fv and gv is only 230 rxns
gv<-dgall[dgall$t4median!=dgall$t15median | 
            dgall$t20median!=dgall$t15median | 
            dgall$t4median!=dgall$t20median,]$rxnid


######## define fv.sig, fluxes significantly differ across temperature ########
#-fv.sig record the significance of the Kruskal-Wallis test on flux
fv.sig<-data.frame(matrix(nrow=0,ncol=6))
colnames(fv.sig)<-c("rxnid","p_kw","es","p_4v15","p_4v20","p_20v15")
#--perform kwt and posthoc tests for each reaction
for(r in rxnlist){
  #-temporary data fram for storing the melt data of fluxes in 3 temperature
  rf<-data.frame(flux=c(as.numeric(frand4[r,]),
                        as.numeric(frand15[r,]),
                        as.numeric(frand20[r,])),
                 temperature=c(rep("T4",ncol(frand4)),
                               rep("T15",ncol(frand15)),
                               rep("T20",ncol(frand20))))
  
  #-Kruskal-Wallis test
  kw<-kruskal.test(flux~temperature, data=rf)
  if(is.na(kw$p.value)){ kw$p.value<-1000 }
  if(kw$p.value<=0.05){
    rowid=nrow(fv.sig)+1
    fv.sig[rowid,]$rxnid<-r
    fv.sig[rowid,]$p_kw<-kw$p.value
    #-calculate effect size, epsilon-squared > 0.64 suggest very strong effect
    #-https://peterstatistics.com/CrashCourse/3-TwoVarUnpair/NomOrd/NomOrd3c.html
    fv.sig[rowid,]$es<-epsilonSquared(x = rf$flux, g = rf$temperature)
    #-posthoc test
    pwt<-pairwise.wilcox.test(rf$flux,g=rf$temperature,
                         p.adjust.method = "bonferroni")
    fv.sig[rowid,]$p_4v15<-pwt$p.value["T4","T15"]
    fv.sig[rowid,]$p_4v20<-pwt$p.value["T4","T20"]
    fv.sig[rowid,]$p_20v15<-pwt$p.value["T20","T15"]
  }
}
#-rxns with significant flux variance based on Kruskal-Wallis test: 425
write.table(fv.sig,file=paste("fv_kwt_sig_r3.tsv",sep=''),quote=FALSE,
            row.names=FALSE,col.names=TRUE,sep="\t")

######## define gv.sig, dges significantly differ across temperature ########
#-gv.sig record the significance of the Kruskal-Wallis test on deltaG
gv.sig<-data.frame(matrix(nrow=0,ncol=6))
colnames(gv.sig)<-c("rxnid","p_kw","es","p_4v15","p_4v20","p_20v15")
#-gv2.sig record the deltaG that flipped signs among different temperature 
gv2.sig<-data.frame(matrix(nrow=0,ncol=1))
colnames(gv2.sig)<-c("rxnid")
#--perform kwt, posthoc test, and flipped signs test for each reaction
for(r in rxnlist){
  #-temporary data fram for storing the melt data of dges in 3 temperature
  gf<-data.frame(dg=c(as.numeric(dgrand4[r,]),
                        as.numeric(dgrand15[r,]),
                        as.numeric(dgrand20[r,])),
                 temperature=c(rep("T4",ncol(dgrand4)),
                               rep("T15",ncol(dgrand15)),
                               rep("T20",ncol(dgrand20))))
  
  #-Kruskal-Wallis test
  kw<-kruskal.test(dg~temperature, data=gf)
  if(is.na(kw$p.value)){ kw$p.value<-1000 }
  if(kw$p.value<=0.05){
    rowid<-nrow(gv.sig)+1
    gv.sig[rowid,]$rxnid<-r
    gv.sig[rowid,]$p_kw<-kw$p.value
    #-calculate effect size, epsilon-squared > 0.64 suggest very strong effect
    #-https://peterstatistics.com/CrashCourse/3-TwoVarUnpair/NomOrd/NomOrd3c.html
    gv.sig[rowid,]$es<-epsilonSquared(x = gf$dg, g = gf$temperature)
    #-posthoc test
    pwt<-pairwise.wilcox.test(gf$dg,g=gf$temperature,
                              p.adjust.method = "bonferroni")
    gv.sig[rowid,]$p_4v15<-pwt$p.value["T4","T15"]
    gv.sig[rowid,]$p_4v20<-pwt$p.value["T4","T20"]
    gv.sig[rowid,]$p_20v15<-pwt$p.value["T20","T15"]
  }
  
  #-check flipped signs
  x<-data.frame(median=c(median(as.numeric(dgrand4[r,])),
       median(as.numeric(dgrand15[r,])),
       median(as.numeric(dgrand20[r,]))))#
  if(product(x$median[c(1,2)])<0 | product(x$median[c(1,3)])<0 | 
     product(x$median[c(2,3)])<0){
    gv2.sig[nrow(gv2.sig)+1,"rxnid"]=r
  }
}
#-write gv.sig to output
#-rxns with significant deltaG variance based on Kruskal-Wallis test: 753
write.table(gv.sig,file=paste("gv_kwt_sig_r3.tsv",sep=''),quote=FALSE,
            row.names=FALSE,col.names=TRUE,sep="\t")
#-write gv2.sig to output
#-rxns with flipped deltaG: 32
write.table(gv2.sig,file=paste("gv_flipped_r3.tsv",sep=''),quote=FALSE,
            row.names=FALSE,col.names=TRUE,sep="\t")

######### filter fv.sig and gv.sig by posthoc analysis #########################
#-The threshold for effect size, epsilon-squared was based on below:
#-https://peterstatistics.com/CrashCourse/3-TwoVarUnpair/NomOrd/NomOrd3c.html
#--0.00 < 0.01 - Negligible
#--0.01 < 0.04 - Weak
#--0.04 < 0.16 - Moderate
#--0.16 < 0.36 - Relatively strong
#--0.36 < 0.64 - Strong
#--0.64 < 1.00 - Very strong
#-remove entries in fv.sig if failed effect size (>=0.36) or posthoc test
tmplist<-c()
for(r in rownames(fv.sig)){
  if(fv.sig[r,]$es < 0.36){ tmplist[length(tmplist)+1]=r }
  if(is.na(fv.sig[r,]$p_4v15) | is.na(fv.sig[r,]$p_4v20) | 
     is.na(fv.sig[r,]$p_20v15)){ tmplist[length(tmplist)+1]=r }
  else if(fv.sig[r,]$p_4v15 > 0.05 & fv.sig[r,]$p_4v20 > 0.05 & 
     fv.sig[r,]$p_20v15 > 0.05){ tmplist[length(tmplist)+1]=r }
}
#-rxns with significant flux variance with ES and posthoc test: 302
fv.sig.ph<-fv.sig[!(rownames(fv.sig) %in% tmplist),]$rxnid
write.table(fv.sig[fv.sig$rxnid %in% fv.sig.ph,],
            file=paste("fv_kwt_sig_r3_posthoc.tsv",sep=''),quote=FALSE,
            row.names=FALSE,col.names=TRUE,sep="\t")
#-remove entries in gv.sig if failed effect size or posthoc test
tmplist<-c()
for(r in rownames(gv.sig)){
  if(gv.sig[r,]$es < 0.36){ tmplist[length(tmplist)+1]=r }
  if(is.na(gv.sig[r,]$p_4v15) | is.na(gv.sig[r,]$p_4v20) | 
     is.na(gv.sig[r,]$p_20v15)){ tmplist[length(tmplist)+1]=r }
  else if(gv.sig[r,]$p_4v15 > 0.05 & gv.sig[r,]$p_4v20 > 0.05 & 
          gv.sig[r,]$p_20v15 > 0.05){ tmplist[length(tmplist)+1]=r }
}
#-rxns with significant deltaG variance with ES and posthoc test: 265
gv.sig.ph<-gv.sig[!(rownames(gv.sig) %in% tmplist),]$rxnid
write.table(gv.sig[gv.sig$rxnid %in% gv.sig.ph,],
            file=paste("gv_kwt_sig_r3_posthoc.tsv",sep=''),quote=FALSE,
            row.names=FALSE,col.names=TRUE,sep="\t")

######### Subsetting to gene associated reactions #########################
#-update rxn lists to keep only gene associated reactions
#-updates were applied to: rxnlist, fvlist, gvlist
gtr<-gtr[!(rownames(gtr) %in% c("Diffusion","Gap","Sink","spontaneous")),]
#--total # of gene associated rxns (GAR): 890
garlist<-unique(unlist(gtr$rxns))

#------------write significant differences only for GAR------------------
write.table(fluxall[garlist,],file=paste("fluxall_gar_summary.tsv",sep=''),
            quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")
write.table(dgall[garlist,],file=paste("dgall_gar_summary.tsv",sep=''),
            quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")
write.table(fv.sig[fv.sig$rxnid %in% fv.sig.ph & fv.sig$rxnid %in% garlist,],
            file=paste("fv_kwt_sig_r3_posthoc_gar.tsv",sep=''),quote=FALSE,
            row.names=FALSE,col.names=TRUE,sep="\t")
write.table(gv.sig[gv.sig$rxnid %in% gv.sig.ph & gv.sig$rxnid %in% garlist,],
            file=paste("gv_kwt_sig_r3_posthoc_gar.tsv",sep=''),quote=FALSE,
            row.names=FALSE,col.names=TRUE,sep="\t")
#------------------------------------------------------------------------

#--total # of GAR with different flux between at least 2 temperatures: 329
garfv<-intersect(garlist,fv)
#--total # of GAR that showed significant flux variability: 274
garfvsig<-intersect(garlist,fv.sig.ph)
#--total # of GAR that showed significant deltaG variability: 254
gargvsig<-intersect(garlist,gv.sig.ph)
#-take intersection from fvlist & gvlist
#--total # of GAR with sig variability in both flux & deltaG: 78
garfvgvsig<-intersect(garfvsig,gargvsig)
# #-remove OMP associated reactions
# #--total # of GAR (non OMP) with sig variability in both flux & deltaG: 77
# garfvgvsig<-garfvgvsig[!(str_detect(garfvgvsig,"^OMP_"))]

######### Subsetting to GAR and significant DE #########################
#-DEGAR: differentially expressed gene (DEG)-associated reactions
#--total # of DEGAR: 517
degar<-unique(unlist(gtr[intersect(deglist,gtr$gid),]$rxns))
#--total # of DEGAR with different flux between at least 2 temperatures: 208
degarfv<-intersect(degar,fv)

######### Prepare data for visualization #########################
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
                   label=c("4°C","15°C","20°C"),
                   color=c("blue","gold","firebrick1"))
rownames(colmap)<-c("T4","T15","T20")

#-Define labels corresponding to column names
lablist<-data.frame(name=c("gvmedian","fvmedian","gemedian_all",
                           "gemedian_early","gemedian_late","gemedian_stat"),
                    label=c("delta G","Metabolic Flux",
                            "Gene Expression","Gene Expression - early",
                            "Gene Expression - late",
                            "Gene Expression - stationary"))
rownames(lablist)<-lablist$name

#-Prepare scaled data, based on median, with abs() or not
#-data maybe raw or scaled; scaling better shows the changes
#-data can be based on subsetting with garfv or garfvgvsig:
#--garfv: GAR & variable flux regardless of deltaG
rselect<-garfv
#rselect<-gargvsig
# #--degarfv: GAR & variable flux regardless of deltaG & Gene is significant DE
# rselect<-degarfv
#--perform scaling
fv.scale<-fluxall[rselect,colmap$median]
gv.scale<-dgall[rselect,colmap$median]
#---absolute values may be taken to facilitate the comparison
#----taking absolute values allow the comparison of differences in magnitude
#----only two reactions in degarfv ("GAPD" and "PTAr") had flipped deltaG,
#----so it will have a minimal effect on the interpretation of the deltaG
fv.scale<-abs(fv.scale)
gv.scale<-abs(gv.scale)
#---set "center" as FALSE in scale() when abs() was taken
for(i in rselect){
  fv.scale[i,]<-scale(t(fv.scale[i,colmap$median]),center=FALSE,scale=TRUE)
  gv.scale[i,]<-scale(t(gv.scale[i,colmap$median]),center=FALSE,scale=TRUE)
}

#-fvgv is the core data for visualization of flux vs delta G vs gene expression
#--make fvgv based on scaled absolute values
fvgv<-data.frame(matrix(nrow=length(rselect)*3,ncol=0))
fvgv$rxnid<-rep(rselect,3)
fvgv$temperature<-c(rep("T4",length(rselect)),rep("T15",length(rselect)),
                 rep("T20",length(rselect)))
for(i in colnames(colmap)[1:1]){
  fvgv[1:nrow(fvgv),paste("fv",i,sep='')]<-c(fv.scale[rselect,colmap["T4",i]],
                                             fv.scale[rselect,colmap["T15",i]],
                                             fv.scale[rselect,colmap["T20",i]])
  fvgv[1:nrow(fvgv),paste("gv",i,sep='')]<-c(gv.scale[rselect,colmap["T4",i]],
                                             gv.scale[rselect,colmap["T15",i]],
                                             gv.scale[rselect,colmap["T20",i]])
}

#-Add gene expression data to fvgv
#--mrnmap provides the column numbers mapping for the cts.mrn
mrnmap<-data.frame(matrix(nrow=3,ncol=0))
rownames(mrnmap)<-rownames(colmap)
mrnmap$gemedian_all<-c(list(c(1,2,3,4,5,6)),list(c(7,8,9,10,11,12)),
              list(c(13,14,15,16,17,18)))
mrnmap$gemedian_early<-c(list(c(1,2)),list(c(7,8)),list(c(13,14)))
mrnmap$gemedian_late<-c(list(c(3,4)),list(c(9,10)),list(c(15,16)))
mrnmap$gemedian_stat<-c(list(c(5,6)),list(c(11,12)),list(c(17,18)))
#-populate gene expression data for all growth phases or each growth phase
#--Median across all samples was used as the gene expression MRN
#--if multiple genes associate with one rxn, MRN sum across all genes was used
#---so for multiple gene associated reactions, sum MRN was calculated per sample
#---before taking the median
fvgv$gemedian_all<-rep(NA,length(rselect)*3)
fvgv$gemedian_early<-rep(NA,length(rselect)*3)
fvgv$gemedian_late<-rep(NA,length(rselect)*3)
fvgv$gemedian_stat<-rep(NA,length(rselect)*3)
for(i in rownames(fvgv)){
  r<-fvgv[i,]$rxnid
  for(j in colnames(mrnmap)){
    x<-cts.mrn[rownames(cts.mrn) %in% unlist(rtg[r,]$genes),
               unlist(mrnmap[fvgv[i,]$temperature,j])]
    fvgv[i,j]<-median(apply(as.matrix(x),2,sum))
  }
}
#-scaling the gene expression columns across all temperatures for each reaction
for(r in rselect){
  for(j in colnames(mrnmap)){
    rs<-rownames(fvgv[fvgv$rxnid==r,]) #rows involved in the scaling
    fvgv[rs,j]=scale(fvgv[rs,j],center=FALSE,scale=TRUE)
  }
}

######### Visualization: deltaG vs Flux vs Gene Expression #####################
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
pfunc<-function(data,t,grp,x,y,xlabel,ylabel,colmap){
  ggplot(data[data[,grp]==t & !is.na(data[,x]) & !is.na(data[,y]),],
         aes_string(x=x, y=y)) +
    geom_point(aes_string(colour=grp),size = 1,
               position=position_jitter(w=0.05,h=0)) +
    xlim(-0.1,1.55) +
    ylim(-0.1,1.55) +
    scale_colour_manual(values=t(colmap)["color",]) +
    labs(title=colmap[t,]$label, x=xlabel, y=ylabel) +
    theme(legend.position="none", title = element_text(colour = "#57585a"),
          legend.text = element_text(colour="#57585a", size = 9))
}

#-main figure: Plot pairwise relationship between deltaG, flux, expression, 
#-with all growth phases combined in the expression
#--ax: a list of all columns for pairwise comparison
#---main figure
ax<-c("gvmedian","fvmedian","gemedian_all")
outpdf<-"deltaG_Flux_GE_all_comparisons.pdf"
# #---Supplemental figure
# ax<-c("gvmedian",colnames(mrnmap)[2:ncol(mrnmap)])
# outpdf<-"FigSX_deltaG_GE_phases.pdf"
#--pcnt: number of pairwise combinations, equal to # or rows in final plot
# #---all combinations, lower index as x-axis; higher index as y-axis
# pcnt<-length(ax)*(length(ax)-1)/2
#---only use first index as x-axis, others as y-axis
pcnt<-length(ax)-1
#--pi: the index of the plot serial
pi<-0 
plots<-c()
# #---all combinations, lower index as x-axis; higher index as y-axis
# for(i in 1:(length(ax)-1)){
#---only use first index as x-axis, others as y-axis
for(i in 1:1){
  xs<-ax[i]
  for(j in (i+1):length(ax)){
    ys<-ax[j]
    for(k in 1:nrow(colmap)){
      t<-rownames(colmap)[k] #temperature label
      pm<-pfunc(fvgv,t,"temperature",xs,ys,lablist[xs,]$label,
                lablist[ys,]$label,colmap)
      pi<-pi+1
      plots[[pi]]<-ggMarginal(pm, type="boxplot")
    }  
  }
}
pdf(outpdf,height=4*pcnt, width=12)
grid.arrange(grobs=plots,nrow=pcnt)
dev.off()

#-Alternative FigX: Plot boxplot of the scaled flux, expression, and deltaG
outpdf<-"deltaG_Flux_GE_all_boxplot.pdf"
fvgv.long<-pivot_longer(data=fvgv[!is.na(fvgv$gvmedian),],
                        cols=c("fvmedian","gemedian_all","gvmedian"),
                        names_to = "type",
                        values_to = "value") %>% 
  mutate(type = recode(type, 
                       "fvmedian" = lablist["fvmedian",]$label, 
                       "gemedian_all" = lablist["gemedian_all",]$label, 
                       "gvmedian" = lablist["gvmedian",]$label))

ggplot(fvgv.long,
       aes(x=factor(temperature,level = c('T4','T15','T20')),y=value,
           fill=temperature)) +
  geom_boxplot() +
  scale_fill_manual(values=t(colmap)["color",]) +
  labs(title='', x="", y="") +
  scale_x_discrete(breaks=rownames(colmap), labels=colmap$label) +
  theme(legend.position="none", title = element_text(colour = "#57585a"),
        legend.text = element_text(colour="#57585a", size = 9)) +
  facet_wrap( ~ factor(type,level=c(lablist["fvmedian",]$label,
                                    lablist["gemedian_all",]$label,
                                    lablist["gvmedian",]$label)),nrow=1,)
ggsave(outpdf,height=5,width=12)

########## Visualization: temperature comparison with dG, Flux, GE #############
#-Define a function for plotting pairs of columns in fvgv2 
#--data: the fvgv2 data frame
#--t: definition of a group value, should be retrieved from the column 'grp'
#--grp: the column name for matching the variable t
#--x: the column name used for the x-axis
#--y: the column name used for the y-axis
#--xlab: the label used for the x-axis
#--ylab: the label used for the y-axis
#--colmap: a dataframe that stores the mappign of grp to colors,
#          colmap should at least contain columns named "color" and "label"
lm_eqn <- function(df){
  m <- lm(y ~ x, df)
  eq <- substitute(italic(y) == b %.% italic(x) + a*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq))
}

pfunc_2<-function(data,t,grp,x,y,xlabel,ylabel,colmap){
  df=data[data[,grp]==t & !is.na(data[,x]) & !is.na(data[,y]),c(x,y)]
  colnames(df)=c('x','y')
  title<-lablist[t,]$label
  if(title == 'delta G'){ title<-expression(paste(Delta[r],"G'",sep='')) }
  ggplot(data[data[,grp]==t & !is.na(data[,x]) & !is.na(data[,y]),],
         aes_string(x=x, y=y)) +
    geom_point(aes_string(colour=grp),size = 1,
               position=position_jitter(w=0.05,h=0)) +
    geom_text(x = 0.82, y = 1.2, label = lm_eqn(df), parse = TRUE) +
    # xlim(0,1.5) +
    # ylim(0,1.55) +
    stat_smooth(method=lm) +
    #    stat_smooth(formula = y ~ s(x, bs = "cs"))+
    scale_colour_manual(values='grey') +
    labs(title=title, x=xlabel, y=ylabel) +
    theme(legend.position="none", title = element_text(colour = "#57585a"),
          legend.text = element_text(colour="#57585a", size = 9))
}

#-fvgv2 is the core data for visualization of flux, dG, GE between temperatures
#--make fvgv2 based on scaled absolute values from fvgv
fvgv2<-data.frame(matrix(nrow=length(rselect)*6,ncol=0))
fvgv2$rxnid<-rep(rselect,6)
fvgv2$metric<-c(rep("fvmedian",length(rselect)),
                rep("gvmedian",length(rselect)),
                rep("gemedian_all",length(rselect)),
                rep("gemedian_early",length(rselect)),
                rep("gemedian_late",length(rselect)),
                rep("gemedian_stat",length(rselect))
)
for(i in rownames(colmap)){
  fvgv2[1:nrow(fvgv2),i]<-c(fvgv[fvgv$temperature==i,]$fvmedian,
                            fvgv[fvgv$temperature==i,]$gvmedian,
                            fvgv[fvgv$temperature==i,]$gemedian_all,
                            fvgv[fvgv$temperature==i,]$gemedian_early,
                            fvgv[fvgv$temperature==i,]$gemedian_late,
                            fvgv[fvgv$temperature==i,]$gemedian_stat
  )
}

#-plotting and save to file
ax<-c("T4","T15","T20")
# mt<-c("gvmedian","fvmedian","gemedian_all",
#       "gemedian_early","gemedian_late","gemedian_stat")
mt<-c("gvmedian","fvmedian","gemedian_all")
pcnt<-length(mt)
#--pi: the index of the plot serial
pi<-0 
plots<-c()
outpdf<-'FigureS2.pdf'
for(m in mt){ #metric to plot
  for(i in 1:(length(ax)-1)){
    #  for(i in 1:1){
    #    xs<-colmap[ax[i],]$median
    xs<-ax[i]
    for(j in (i+1):length(ax)){
      #      ys<-colmap[ax[j],]$median
      ys<-ax[j]
      pm<-pfunc_2(fvgv2,m,"metric",ys,xs,colmap[ys,]$label,
                  colmap[xs,]$label,colmap)
      pi<-pi+1
      plots[[pi]]<-ggMarginal(pm, type="violin",
                              xparams=list(fill=colmap[ys,]$color),
                              yparams=list(fill=colmap[xs,]$color))
    }  
  }
}
pdf(outpdf,height=4*pcnt, width=12)
grid.arrange(grobs=plots,nrow=pcnt)
dev.off()

######### visualization: Plotting individual genes of interest #################
#-Define genes of interest
# #--full list of reactions in Fig 1
# goi<-c("ACGAMK","AGDC","G6PDA","FBP","FBA","PGI","G6PDHy","PGL","PGDH","PGDHY",
#        "EDA","TPI","TKT2","TAL","TKT1","RPI","RPE","XPK","GAPD","PGK","PGM",
#        "ENO","PPS","PPC","PYK","PFL","PDH","PTAr","ACKr","CS","ACONT","ICDHy",
#        "ICDHxi","ICL","SUCOAS","SUCD7","FRD10","FRD11","FUM","ME2")
#-list of reactions in Fig 1 and highly interesting, all has DE in late phase
#--TCA/glyoxylate shunt: "CS","ICL"
#--Glycolysis: "PFL","FBA"
#--Substrate-level phosphorylation: "PTAr"
#--ED: "PGDHY"
#--Oxidative PPP: "G6PDHy"
#--Non-oxidative PPP: "TKT2"
#--energy metabolism: "ATPS4r",CYTBO3_4
#--e-transport: "NADH4", "NADH13" (same genes as NADH11/16, only included one)
#--NDH: NADH dehydrogenase
goi<-data.frame(rxnid=c("PFL","FBA","G6PDHy","PGDHY","ICL","PTAr","ATPS4r",
                        "SUCD7","CYTBO3_4","","",""),#"TKT2","PGDHY","CS","CYTBD","ATPS4r","CYTBO3_4","SUCD7","FRD10","FRD11","SUCOAS")
                rname=c("PFL","FBA","G6PDHy","PGDHY","ICL","PTAr","ATPS4r",
                        "SUCD7", "CYTBO3_4","FRD10/11",
                        "NADH11/13/16",
                        "NADH4"),
                main=c("Pyruvate formate-lyase",
                       "Fructose-bisphosphate aldolase",
                       "G6P dehydrogenase","phosphogluconate dehydratase",
                       "Isocitrate lyase","Phosphate Acetyltransferase",
                       "ATP Synthase","Succinate dehydrogenase", 
                       "cytochrome bb3 ubiquinol oxidase","Fumarate reductase",
                       "NADH:quinone oxidoreductase",
                       "NADH dehydrogenase"), #"Transketolase"
                deg_late=rep("",12) #deg will be populated with list of DE genes
                )
goi[nrow(goi)-2,]$rxnid=list(c("FRD10","FRD11"))
goi[nrow(goi)-1,]$rxnid=list(c("NADH11","NADH13","NADH16"))
#goi[nrow(goi),]$rxnid=list(c("NADH4","NADH12","NADH14"))
goi[nrow(goi),]$rxnid=list(c("NADH4"))
# ed<-c("FBA","FBP","PGI","G6PDHy","PGL","PGDH","PGDHY","EDA","TPI")
#-making plots per reaction
plots=c()
for(i in 1:nrow(goi)){
  r<-unlist(goi$rxnid[i])
  goi[i,]$deg_late<-list(intersect(deglist,unlist(rtg[r,]$genes)))
  #-temporary data fram for storing the melt data of fluxes in 3 temperature
  #-flux, abs() was applied to look at the flux magnitude
  #-abs() is no problem because NONE of the "goi" has flipped deltaG
  rf<-data.frame(flux=c(apply(abs(frand4[r,]),2,sum),
                        apply(abs(frand15[r,]),2,sum),
                        apply(abs(frand20[r,]),2,sum)),
                 temperature=c(rep("T4",ncol(frand4)),
                               rep("T15",ncol(frand15)),
                               rep("T20",ncol(frand20))),
                 rxnid=rep(goi$rname[i],ncol(frand4)+ncol(frand15)+ncol(frand20)))
  #-deltaG
  rdg<-data.frame(dg=c(apply(dgrand4[r,],2,sum),
                       apply(dgrand4[r,],2,sum),
                       apply(dgrand4[r,],2,sum)),
                  temperature=c(rep("T4",ncol(dgrand4)),
                                rep("T15",ncol(dgrand15)),
                                rep("T20",ncol(dgrand20))),
                  rxnid=rep(goi$rname[i],ncol(dgrand4)+ncol(dgrand15)+ncol(dgrand20)))
  #-gene expression
  cts.genelist<-data.frame(matrix(nrow=0,ncol=4))
  for(g in unique(unlist(rtg[r,]$genes))){
    d<-data.frame(mrn=cts.mrn[g,],
                  temperature=c(rep("T4",6),rep("T15",6),rep("T20",6)),
                  growthphase=rep(c("early","early","late","late",
                                    "stat","stat"),3),
                  gene=rep(g,18))
    cts.genelist<-rbind(cts.genelist,d)
  }
  
  #-plot reaction fluxes
  g1<-ggplot(rf,aes(x=factor(temperature,level = c('T4','T15','T20'),
                             labels=c('4°C', '15°C', '20°C')), 
                    y=flux, fill=temperature)) +
    geom_boxplot() +
    scale_y_continuous(expand = expansion(c(0,0.2))) +
    #geom_jitter(width=0.35,height=0) +
    scale_fill_manual(values=t(colmap)["color",]) +
    labs(title="", x="", y="Metabolic Flux", fill="Temperature") +
    theme(legend.position="none", plot.margin=margin(1,0,1,1, unit='cm')) +
    facet_wrap( ~ rxnid,nrow=1)
  
  #-plot deltaG
  g2<-ggplot(rdg,aes(x=factor(temperature,level = c('T4','T15','T20'),
                              labels=c('4°C', '15°C', '20°C')), 
                     y=dg, fill=temperature)) +
    geom_boxplot() +
    scale_fill_manual(values=t(colmap)["color",]) +
    labs(title="", x="", y="delta G", fill="Temperature") +
    theme(legend.position="none") +
    facet_wrap( ~ rxnid,nrow=1)
  
  #-plot gene expression, only plot the late exponential phase
  # g3<-ggplot(cts.genelist[cts.genelist$growthphase=="late",],
  g3<-ggplot(cts.genelist[cts.genelist$mrn!=0,],
                        aes(x=factor(temperature,level = c('T4','T15','T20'),
                                     labels=c('4°C', '15°C', '20°C')), 
                              y=mrn, fill=temperature)) + 
    geom_boxplot() +
    scale_fill_manual(values=t(colmap)["color",]) +
    scale_y_log10(expand=expansion(c(0,0.2))) +
    labs(title="", x="", y="Gene Expression", fill="Temperature") +
    theme(legend.position="none", axis.text.x = element_text(angle = 90, vjust = 1, hjust=1), plot.margin=margin(1,1,1,0, unit='cm')) +
    facet_wrap( ~ growthphase,nrow=1)
  
  #-save combined plot to plots[[]]
  title=textGrob(goi$main[i], just="top",gp=gpar(fontsize=12,fontface="bold")) 
  plots[[i]]<-grid.arrange(g1,g3,nrow=1,top=title)
}
pdf(paste("FigureS1.pdf",sep=""),height=20, width=20)
grid.arrange(grobs=plots,nrow=4)
dev.off()


######################### ######################### #########################
######################### ######################### #########################
######################### ######################### #########################


