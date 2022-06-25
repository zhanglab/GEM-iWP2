library(stringr)

colors<-data.frame(T4='blue', T15='gold', T20='firebrick1')

################# Process Experimental Data ################################
#-read experimental growth data
gc4 <- read.csv('./4C-gcurve-replicates', sep='\t', header=T)
gc15 <- read.csv('./15C-gcurve-replicates', sep='\t', header=T)
gc20 <- read.csv('./20C-gcurve-replicates', sep='\t', header=T)

#-ignore initial inoculation (time 0) in the plot
gc4 <- gc4[-1,]
gc15 <- gc15[-1,]
gc20 <- gc20[-1,]

#-convert growth to gDW, time to Days, calculate mean and sd
#--4C
gc4_gdw <- gc4[,2:4]*0.69
gc4_gdw$time<-gc4[,1]/24
gc4_gdw$mean <- apply(gc4_gdw[,c("X1","X2","X3")], 1, mean)
gc4_gdw$sd <- apply(gc4_gdw[,c("X1","X2","X3")], 1, sd)
gc4_gdw$temperature<-rep("T4",nrow(gc4_gdw))
#--15C
gc15_gdw <-  gc15[,2:4]*0.69
gc15_gdw$time<-gc15[,1]/24
gc15_gdw$mean <- apply(gc15_gdw[,c("X1","X2","X3")], 1, mean)
gc15_gdw$sd <- apply(gc15_gdw[,c("X1","X2","X3")], 1, sd)
gc15_gdw$temperature<-rep("T15",nrow(gc15_gdw))
#--20C
gc20_gdw <-  gc20[,2:4]*0.69
gc20_gdw$time<-gc20[,1]/24
gc20_gdw$mean <- apply(gc20_gdw[,c("X1","X2","X3")], 1, mean)
gc20_gdw$sd <- apply(gc20_gdw[,c("X1","X2","X3")], 1, sd)
gc20_gdw$temperature<-rep("T20",nrow(gc20_gdw))

#-calculate growth rate
#--range of data used for the growth rate calculation
gr_range<-data.frame(T4=c(5,21,12,21,22), T15=c(1,7,4,6,8), T20=c(1,8,5,8,11))
rownames(gr_range)<-c("start","end","tr1","tr2","tr3") #tr, transcriptome
#--combine experimental data from all three temperature
gcall<-rbind(gc4_gdw,gc15_gdw,gc20_gdw)
#--define function that calculates growth rates from two time points
gr_calculate<-function(x,y,s,e){ 
  #---function to calculate growth rate
  #---x: vector contain biomass data from one experimental replicate
  #---y: vector contain the corresponding times of measurement
  #---s: start position in x for growth rate calculation
  #---e: end position in x for growth rate calculation
  (x[e]-x[s])/(y[e]-y[s])
}
#--populate the growth rate matrix, save to output
gr<-data.frame(matrix(nrow=3,ncol=7))
rownames(gr)<-c("T4","T15","T20")
colnames(gr)<-c("X1","X2","X3","min","max","mean","sd")
for(i in rownames(gr)){
  for(j in c("X1","X2","X3")){
    gr[i,j]<-gr_calculate(gcall[gcall$temperature==i,j],
                          gcall[gcall$temperature==i,]$time,
                          gr_range["start",i],
                          gr_range["end",i])
  }
  gr[i,"min"]<-min(gr[i,c("X1","X2","X3")])
  gr[i,"max"]<-max(gr[i,c("X1","X2","X3")])
  gr[i,"mean"]<-mean(as.numeric(gr[i,c("X1","X2","X3")]))
  gr[i,"sd"]<-sd(as.numeric(gr[i,c("X1","X2","X3")]))
}
write.table(gr,"growth_rates.tsv",quote=FALSE,sep="\t",row.names=TRUE)

################## Fit ATPM ###################################################
#-read results from robustness simulations
simu4f<-'04C-ATPM-Robustness.tsv'
simu15f<-'15C-ATPM-Robustness.tsv'
simu20f<-'20C-ATPM-Robustness.tsv'
#--4C
simu4<-read.table(simu4f,header=FALSE,sep="\t")
simu4$temperature<-rep("T4",nrow(simu4))
#--15C
simu15<-read.table(simu15f,header=FALSE,sep="\t")
simu15$temperature<-rep("T15",nrow(simu15))
#--20C
simu20<-read.table(simu20f,header=FALSE,sep="\t")
simu20$temperature<-rep("T20",nrow(simu20))

#-combine data from all temperatures
simuall<-rbind(simu4,simu15,simu20)
colnames(simuall)=c("atpm","biomass","temperature")

#-range of data used for the ATPM fitting
# For 4C the linear fitting was on the points from ATPM of 1 to 9.
# For 15C the linear fitting was on the points from ATPM of 1 to 3
# For 20C the linear fitting was on the points from ATPM of 5 to 9.5
simu_range<-data.frame(T4=c(1,9), T15=c(1,3), T20=c(5,9.5))
rownames(simu_range)<-c("lb","ub") 

#-populate the ATPM settings based on fitting results
atpm<-data.frame(matrix(nrow=3,ncol=3))
rownames(atpm)<-c("T4","T15","T20")
colnames(atpm)<-c("X1","X2","X3")
#--record the fitting parameters from each simulation
fit_all<-data.frame(matrix(nrow=0,ncol=4))
colnames(fit_all)<-c("temperature","replicate","slope","intercept")
for(i in rownames(atpm)){
  s<-simuall[simuall$temperature==i,]
  lb<-gr[i,]$min-gr[i,]$min*0.05
  ub<-gr[i,]$max+gr[i,]$max*0.05
  fit<-lm(s[s$biomass>=lb & s$biomass<=ub,]$biomass ~ 
            s[s$biomass>=lb & s$biomass<=ub,]$atpm)
  for(j in colnames(atpm)){
    b<-gr[i,j]
    atpm[i,j]<-(b-fit$coefficients[1])/fit$coefficients[2]
    fit_all[nrow(fit_all)+1,c("temperature","replicate")]<-c(i,j)
    fit_all[nrow(fit_all),c("slope","intercept")]<-c(fit$coefficients[2],
                                                    fit$coefficients[1])
  }
}
atpm$min<-apply(atpm[,c("X1","X2","X3")],1,min)
atpm$max<-apply(atpm[,c("X1","X2","X3")],1,max)
atpm$mean<-apply(atpm[,c("X1","X2","X3")],1,mean)
atpm$sd<-apply(atpm[,c("X1","X2","X3")],1,sd)
write.table(atpm,"atpm.tsv",quote=FALSE,sep="\t",row.names=TRUE)
write.table(fit_all,"atpm_fit.tsv",quote=FALSE,sep="\t",row.names=FALSE)

############ Plot Fig.1 ######################################################
pdf('./Figure1.pdf', height=8, width=8)
par(mfrow=c(2,2))
#-panel A, experimental data
plot(NULL,ylim=c(0.005,0.7),xlim=c(0,10),
     xlab='Time (Days)', ylab='Biomass, log(gDW)',
     main="Experimental Data",log='y')
for(i in rownames(gr)){
  gci<-gcall[gcall$temperature==i,]
  points(x=gci$time,y=gci$mean,
         pch='o', cex=0.75, col=as.character(colors[i]))
  arrows(gci$time, gci$mean-gci$sd, 
         gci$time, gci$mean+gci$sd, 
         length=0.025, angle=90, code=3, col=as.character(colors[i]))
  points(gci[gr_range[c("tr1","tr2","tr3"),i],]$time, 
         gci[gr_range[c("tr1","tr2","tr3"),i],]$mean*1.1,
         pch='*', cex=1, col=as.character(colors[i]))
  lines(x=gci[gr_range[c("start","end"),i],]$time,
        y=gci[gr_range[c("start","end"),i],]$mean, col=as.character(colors[i]))
}
text(gc4_gdw[18,]$time+0.6,gc4_gdw[17,]$mean-0.02,
     labels=paste(round(gr["T4",]$mean,2),"gDW/Day",sep=' '),
     cex=0.8, col=as.character(colors["T4"]))
text(gc15_gdw[6,]$time-0.9,gc15_gdw[5,]$mean+0.1,
     labels=paste(round(gr["T15",]$mean,2),"gDW/Day",sep=' '),
     cex=0.8, col=as.character(colors["T15"]))
text(gc20_gdw[10,]$time+0.6,gc20_gdw[9,]$mean-0.03,
     labels=paste(round(gr["T20",]$mean,2),"gDW/Day",sep=' '),
     cex=0.8, col=as.character(colors["T20"]))
legend('bottomright', c('4°C', '15°C', '20°C'), 
       fill=c('blue', 'gold', 'firebrick1'), bty='n', cex=1)

#-Panel B, C, D, ATPM fitting for each temperature
for(i in rownames(atpm)){
  s<-simuall[simuall$temperature==i,]
  lb<-gr[i,]$min
  ub<-gr[i,]$max
  plot(s$atpm, s$biomass, pch='o', cex=0.25, 
       main=paste(str_replace(i,"T",''),'C ATPM Calibration',sep=''), 
       xlab='ATP Maintenance Flux', ylab='Biomass Flux')
  clip(0, 20, lb, ub)
  abline(a=mean(fit_all[fit_all$temperature==i,]$intercept),
         b=mean(fit_all[fit_all$temperature==i,]$slope),col="red",lwd=4)
  clip(-10, atpm[i,]$min, lb-0.001, ub+0.001)
  abline(h = gr[i,]$max, col='lightblue',lwd=2)
  clip(-10, atpm[i,]$max, lb-0.001, ub+0.001)
  abline(h = gr[i,]$min, col='lightblue',lwd=2)
  clip(-10, 20, -10, gr[i,]$max)
  abline(v = atpm[i,]$min, col='lightblue',lwd=2)
  clip(-10, 20, -10, gr[i,]$min)
  abline(v = atpm[i,]$max, col='lightblue',lwd=2)
  clip(-10, 20, -10, 10)
  text(atpm[i,]$mean+2,gr[i,]$mean+0.1*gr[i,]$mean,"ATPM flux:")
  text(atpm[i,]$mean+2.5,gr[i,]$mean-0.005,
       paste("[",round(atpm[i,]$min,2),",",round(atpm[i,]$max,2),"]",sep=''))
}
dev.off()

