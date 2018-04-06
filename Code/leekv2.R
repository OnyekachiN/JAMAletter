.libPaths(c("~/Documents/Rpackagesv2",.libPaths()))
library(readr)
library(ggridges)
library(viridis)
library("gridExtra")
library(pROC)
###################################################################################################
# IMPORTING DATA 
###################################################################################################
dat <- read_csv("/Volumes/Lab_Gerke/pval/data/abstractPvals/leekDatav2.csv")

dat$journalName <- gsub("\\..*","",dat$journal)
dat$journalName <- ifelse(dat$journalName=="JAMA1","JAMA",dat$journalName)

dat$CIwidth <- dat$U95 - dat$L95
dat$CIwidth2 <- log(dat$U95) - log(dat$L95)

dat$pGroup1 <- ifelse(dat$pvalue<=0.001,"<=0.001",
                      ifelse(dat$pvalue<=0.01,"(0.001,0.01]",
                             ifelse(dat$pvalue<=0.05,"(0.01,0.05]",">0.05")))
dat$pGroup1 <- factor(dat$pGroup1, levels = c("<=0.001","(0.001,0.01]","(0.01,0.05]",">0.05"))
dat$pGroup2 <- ifelse(dat$pvalue<=0.005,"<=0.005",
                      ifelse(dat$pvalue<=0.05,"(0.005,0.05]",">0.05"))
dat$pGroup2 <- factor(dat$pGroup2, levels = c("<=0.005","(0.005,0.05]",">0.05"))

dat$effectEst2 <- ifelse(dat$effectEst<1,1/dat$effectEst,dat$effectEst)
dat$effectEst3 <- log(dat$effectEst2)

dat$L95v2 <- ifelse(dat$effectEst<1,1/dat$U95,dat$L95)

dat$sampleSize2 <- log(dat$sampleSize)

dat$missingp <- ifelse(!is.na(dat$pvalue),NA,ifelse(dat$L95<1 & dat$U95>1,1,0))

dat_major <- dat[dat$journalName %in% c("American Journal of Epidemiology","Lancet","JAMA"),]

dat_p <- dat_major[!is.na(dat_major$pvalue),]
p_bonferroni <- 0.05/nrow(dat_p)
dat_bon <- dat_p[dat_p$pvalue<=p_bonferroni,]
###################################################################################################
# BASIC INFO
###################################################################################################

# ALL DATA 

table(dat$journalName)

tapply(dat$pvalue,dat$journalName,summary)
ggplot(dat, aes(x = pvalue, y = journalName)) +
  geom_density_ridges() 

tapply(dat$effectEst,dat$journalName,summary)
ggplot(dat, aes(x = effectEst, y = journalName)) +
  geom_density_ridges()

tapply(dat$sampleSize,dat$journalName,summary)
ggplot(dat, aes(x = sampleSize, y = journalName)) +
  geom_density_ridges()

###################################################################################################
# AJE, LANCET, AND JAMA ONLY

# p-value by journal
tapply(dat_major$pvalue,dat_major$journalName,summary)
kruskal.test(dat_major$pvalue~as.factor(dat_major$journalName))

# effect estimate(s) by journal
tapply(dat_major$effectEst,dat_major$journalName,summary)
kruskal.test(dat_major$effectEst~as.factor(dat_major$journalName))
tapply(dat_major$effectEst2,dat_major$journalName,summary)
kruskal.test(dat_major$effectEst2~as.factor(dat_major$journalName))
tapply(dat_major$effectEst3,dat_major$journalName,summary)
kruskal.test(dat_major$effectEst3~as.factor(dat_major$journalName))

# sample size by journal
tapply(dat_major$sampleSize,dat_major$journalName,summary)
kruskal.test(dat_major$sampleSize~as.factor(dat_major$journalName))
tapply(dat_major$sampleSize2,dat_major$journalName,summary)
kruskal.test(dat_major$sampleSize2~as.factor(dat_major$journalName))

# plotting effect estimates by binned p-values
ggplot(dat_major,aes(x=as.factor(pGroup1),y=effectEst)) +
  geom_boxplot()

# plotting effect estimates (<10) by binned p-values
ggplot(dat_major[dat_major$effectEst2<=10,],aes(x=as.factor(pGroup1),y=effectEst2,fill=as.factor(pGroup1))) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=rep("white",5)) +
  geom_jitter(position=position_jitterdodge(jitter.width=.5, jitter.height=0, dodge.width=.75),
              aes(fill=as.factor(pGroup1),col=as.factor(pGroup1)),alpha=0.25) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  guides(fill=FALSE, colour=guide_legend(title="P-value")) +
  labs(x="P-value",y="Effect Estimate") +
  geom_hline(yintercept=1)

# plotting log of effect estimates by binned p-values
ggplot(dat_major,aes(x=as.factor(pGroup1),y=effectEst3,fill=as.factor(pGroup1))) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=rep("white",5)) +
  geom_jitter(position=position_jitterdodge(jitter.width=.5, jitter.height=0, dodge.width=.75),
              aes(fill=as.factor(pGroup1),col=as.factor(pGroup1)),alpha=0.25) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  guides(fill=FALSE, colour=guide_legend(title="P-value")) +
  labs(x="P-value",y="Log(Effect Estimate)")

# effect estimates by binned p-values 
tapply(dat_major$effectEst2,dat_major$pGroup1,summary)
kruskal.test(dat_major$effectEst2~as.factor(dat_major$pGroup1))
# tapply(dat_major$effectEst3,dat_major$pGroup1,summary)
# kruskal.test(dat_major$effectEst3~as.factor(dat_major$pGroup1))
tapply(dat_major$effectEst3,ifelse(is.na(dat_major$pGroup1),1,0),summary)
wilcox.test(dat_major$effectEst3~ifelse(is.na(dat_major$pGroup1),1,0))

# plotting pvalues by effect estimates (<10) by log of the CI width 
ggplot(dat_major[dat_major$effectEst2<=10 & !is.na(dat_major$effectEst2) & !is.na(dat_major$pvalue),],
       aes(x=pvalue,y=effectEst2)) +
  geom_point(aes(colour=CIwidth2)) +
  scale_colour_gradientn(colours = topo.colors(20)) +
  geom_smooth(color="grey",na.rm = TRUE,se=FALSE) +
  guides(fill=FALSE, colour=guide_colorbar(title="Log CI Width"),size=FALSE) +
  labs(x="P-value",y="Effect Estimates") +
  theme_classic()
# plotting pvalues (<=0.05) by effect estimates (<10) by log of the CI width 
ggplot(dat_major[dat_major$effectEst2<=10 & !is.na(dat_major$effectEst2) & !is.na(dat_major$pvalue) & dat_major$pvalue<=0.05,],
       aes(x=pvalue,y=effectEst2)) +
  geom_point(aes(colour=CIwidth2)) +
  scale_colour_gradientn(colours = topo.colors(20)) +
  geom_smooth(color="grey",na.rm = TRUE,se=FALSE) +
  guides(fill=FALSE, colour=guide_colorbar(title="Log CI Width"),size=FALSE) +
  labs(x="P-value",y="Effect Estimates") +
  theme_classic()

# plotting pvalues by effect estimates (<10) by log of the sample size
ggplot(dat_major[dat_major$effectEst2<=10 & !is.na(dat_major$effectEst2) & !is.na(dat_major$pvalue),],
       aes(x=pvalue,y=effectEst2)) +
  geom_point(aes(colour=sampleSize2)) +
  scale_colour_gradientn(colours = topo.colors(20)) +
  geom_smooth(color="grey",na.rm = TRUE,se=FALSE) +
  guides(fill=FALSE, colour=guide_colorbar(title="Log Sample Size")) +
  labs(x="P-value",y="Effect Estimates") +
  theme_classic()
# plotting pvalues (<=0.05) by effect estimates (<10) by log of the sample size
ggplot(dat_major[dat_major$effectEst2<=10 & !is.na(dat_major$effectEst2) & !is.na(dat_major$pvalue) & dat_major$pvalue<=0.05,],
       aes(x=pvalue,y=effectEst2)) +
  geom_point(aes(colour=sampleSize2)) +
  scale_colour_gradientn(colours = topo.colors(20)) +
  geom_smooth(color="grey",na.rm = TRUE,se=FALSE) +
  guides(fill=FALSE, colour=guide_colorbar(title="Log Sample Size")) +
  labs(x="P-value",y="Effect Estimates") +
  theme_classic()

###################################################################################################
#
# "traditional" p-value bins
#
# log of the effect estimate by the log of the CI width 
ggplot(dat_major,
       aes(x=effectEst3,y=CIwidth2)) +
  geom_point(aes(colour=pGroup1,alpha=0.01)) +
  guides(fill=FALSE, colour=guide_legend(title="P-value"), alpha=FALSE) +
  labs(x="log(Effect Estimate)",y="log(CI width)") +
  theme_classic() + scale_colour_manual(values=viridis(12)[c(2,5,9,12)])

# log of the effect estimate by the log of the CI width by binned p-value group 
ggplot(dat_major,
       aes(x=effectEst3,y=CIwidth2)) +
  geom_point(aes(colour=CIwidth2)) +
  scale_colour_gradientn(colours = topo.colors(20)) +
  guides(fill=FALSE, colour=FALSE) +
  labs(x="log(Effect Estimate)",y="log(CI width)") +
  theme_classic() +
  facet_wrap( ~ pGroup1, nrow = 5)

# log of the effect estimates by p value bins 
ggplot(dat_major,aes(x=as.factor(pGroup1),y=effectEst3,fill=as.factor(pGroup1))) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=rep("white",5)) +
  geom_jitter(position=position_jitterdodge(jitter.width=.5, jitter.height=0, dodge.width=.75),
              aes(fill=as.factor(pGroup1),col=as.factor(pGroup1)),alpha=0.25) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  guides(fill=FALSE, colour=guide_legend(title="P-value")) +
  labs(x="P-value",y="Log(Effect Estimate)") 

# log of the CI width by p value bins 
ggplot(dat_major,aes(x=as.factor(pGroup1),y=CIwidth2,fill=as.factor(pGroup1))) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=rep("white",5)) +
  geom_jitter(position=position_jitterdodge(jitter.width=.5, jitter.height=0, dodge.width=.75),
              aes(fill=as.factor(pGroup1),col=as.factor(pGroup1)),alpha=0.25) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  guides(fill=FALSE, colour=guide_legend(title="P-value")) +
  labs(x="P-value",y="Log(CI width)") 

###################################################################################################
#
# "ionnidis" p-value bins
#
# log of the effect estimate by the log of the CI width 
p1 <- ggplot(dat_major,
       aes(x=effectEst3,y=CIwidth2)) +
  geom_point(aes(colour=pGroup2,alpha=0.01)) +
  guides(fill=FALSE, colour=guide_legend(title="P-value"), alpha=FALSE) +
  labs(x="log(Effect Estimate)",y="log(U95) - log(L95)") +
  theme_classic() +
  scale_colour_manual(values=viridis(10)[c(2,5,9)], na.value="grey")

# log of the effect estimates by p value bins 
p2 <- ggplot(dat_major,aes(x=as.factor(pGroup2),y=effectEst3,fill=as.factor(pGroup2))) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=rep("white",4)) +
  geom_jitter(position=position_jitterdodge(jitter.width=.5, jitter.height=0, dodge.width=.75),
              aes(fill=as.factor(pGroup2),col=as.factor(pGroup2)),alpha=0.25) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  guides(fill=FALSE, colour=FALSE) +
  labs(x="P-value",y="log(Effect Estimate)") +
  scale_colour_manual(values=viridis(10)[c(2,5,9)], na.value="grey")

# log of the CI width by p value bins 
p3 <- ggplot(dat_major,aes(x=as.factor(pGroup2),y=CIwidth2,fill=as.factor(pGroup2))) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=rep("white",5)) +
  geom_jitter(position=position_jitterdodge(jitter.width=.5, jitter.height=0, dodge.width=.75),
              aes(fill=as.factor(pGroup2),col=as.factor(pGroup2)),alpha=0.25) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  guides(fill=FALSE, colour=FALSE) +
  labs(x="P-value",y="log(U95) - log(L95)") +
  scale_colour_manual(values=viridis(10)[c(2,5,9)], na.value="grey")

grid.arrange(arrangeGrob(p2, p3, ncol = 2),p1,nrow = 2) 

kruskal.test(dat_major$effectEst3~dat_major$pGroup2)
kruskal.test(dat_major$CIwidth2~dat_major$pGroup2)

###################################################################################################
#
# "ionnidis" p-value bins
#
# effect estimate by the log of the CI width 
p1 <- ggplot(dat_major[dat_major$effectEst2<=10 & dat_major$L95v2<=10,],
             aes(x=effectEst2,y=L95v2)) +
  geom_point(aes(colour=pGroup2,alpha=0.01)) +
  guides(fill=FALSE, colour=guide_legend(title="P-value"), alpha=FALSE) +
  labs(x="Effect Estimate",y="L95") +
  theme_classic() +
  scale_colour_manual(values=viridis(10)[c(2,5,9)], na.value="grey")

# effect estimates by p value bins 
p2 <- ggplot(dat_major[dat_major$effectEst2<=10,],aes(x=as.factor(pGroup2),
                                                      y=effectEst2,fill=as.factor(pGroup2))) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=rep("white",4)) +
  geom_jitter(position=position_jitterdodge(jitter.width=.5, jitter.height=0, dodge.width=.75),
              aes(fill=as.factor(pGroup2),col=as.factor(pGroup2)),alpha=0.25) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  guides(fill=FALSE, colour=FALSE) +
  labs(x="P-value",y="Effect Estimate") +
  scale_colour_manual(values=viridis(10)[c(2,5,9)], na.value="grey") 

# L95 by p value bins 
p3 <- ggplot(dat_major[dat_major$L95v2<=8,],
             aes(x=as.factor(pGroup2),y=L95v2,fill=as.factor(pGroup2))) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=rep("white",5)) +
  geom_jitter(position=position_jitterdodge(jitter.width=.5, jitter.height=0, dodge.width=.75),
              aes(fill=as.factor(pGroup2),col=as.factor(pGroup2)),alpha=0.25) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  guides(fill=FALSE, colour=FALSE) +
  labs(x="P-value",y="L95") +
  scale_colour_manual(values=viridis(10)[c(2,5,9)], na.value="darkgrey") +
  geom_hline(yintercept = 1, col="grey",linetype = 2)

grid.arrange(arrangeGrob(p2, p3, ncol = 2),p1,nrow = 2) 

kruskal.test(dat_major$effectEst3~dat_major$pGroup2)
kruskal.test(dat_major$CIwidth2~dat_major$pGroup2)

###################################################################################################
# dat_major$pGroup3 <- ifelse(dat_major$pvalue<=0.05,1,0)
# dat_major$pGroup3 <- ifelse(dat_major$pvalue<=0.005,2,dat_major$pGroup3)
# 
# roc.multi <- multiclass.roc(dat_major$pGroup2,dat_major$effectEst2)
# rs <- roc.multi[['rocs']]
# plot.roc(rs[[1]])
# sapply(2:length(rs),function(i) lines.roc(rs[[i]],col=i))
# legend("bottomright", legend=c("<=0.005 vs (0.005,0.05]",
#                                "<=0.005 vs >0.05","(0.005,0.05] vs >0.05"),
#        col=c("black","red","green"), lwd=2)
# 
# 
# roc.multi <- multiclass.roc(dat_major$pGroup2,dat_major$L95v2)
# rs <- roc.multi[['rocs']]
# plot.roc(rs[[1]])
# sapply(2:length(rs),function(i) lines.roc(rs[[i]],col=i))
# legend("bottomright", legend=c("<=0.005 vs (0.005,0.05]",
#                                "<=0.005 vs >0.05","(0.005,0.05] vs >0.05"),
#        col=c("black","red","green"), lwd=2)

thresholds <- c(1.2,1.5,2,2.5)
roc1 <- roc(ifelse(dat_major$pvalue<=0.005,1,0),dat_major$effectEst2)
roc2 <- roc(ifelse(dat_major$pvalue<=0.05,1,0),dat_major$effectEst2)
roc1$thresholds <- round(roc1$thresholds,2)
# roc.multi <- multiclass.roc(dat_major$pGroup3,dat_major$effectEst2)
# rs <- roc.multi[['rocs']]
plot.roc(roc1, print.thres = thresholds, print.thres.cex=par("cex"))
plot.roc(roc2,col="blue", print.thres = thresholds, add=TRUE, print.thres.cex=par("cex")/1.5, print.thres.col = "blue")
# lines.roc(rs[[3]],col="green")
legend("bottomright", legend=c("p<=0.005","p<=0.05"),
       col=c("black","blue"), lwd=2)
roc.test(roc1,roc2)
wilcox.test(dat_major$effectEst2~ifelse(dat_major$pvalue<=0.005,1,0))
wilcox.test(dat_major$effectEst2~ifelse(dat_major$pvalue<=0.05,1,0))

rocD <- as.data.frame(cbind(sens=roc1$sensitivities,
                            spec=roc1$specificities,
                            thresh=round(roc1$thresholds,1)))
rocD <- rocD[!duplicated(rocD$thresh),]
rocD <- rocD[rocD$thresh %in% thresholds,]

plot(roc1$specificities,roc1$sensitivities,type="line", xlim = c(1,0),lwd=2, xlab="Specificity",
     ylab="Sensitivity")
points(x=c(0.35,0.62,0.81,0.88),y=c(0.93,0.74,0.50,0.35),pch=20)
text(x=c(0.35,0.62,0.81,0.88),y=c(0.93,0.74,0.50,0.35),labels=c("1.2","1.5","2.0","2.5"), pos=4,
     cex = 0.75)
lines(roc2$specificities,roc2$sensitivities,type="line", lwd=2, col="blue")
points(x=c(0.65,0.87,0.97,0.98),y=c(0.91,0.68,0.43,0.29),pch=20,col="blue")
text(x=c(0.65,0.87,0.97,0.98),y=c(0.91,0.68,0.43,0.29),labels=c("1.2","1.5","2.0","2.5"),
     pos=2, col="blue",cex = 0.75)
lines(x = c(1,0), y = c(0,1), col="grey")
legend("bottomright", legend=c("p<=0.005","p<=0.05"),
       col=c("black","blue"), lwd=2, bty = "n")

roc1.1 <- roc(ifelse(dat_major$pvalue<=0.005,1,0),dat_major$L95v2)
roc2.1 <- roc(ifelse(dat_major$pvalue<=0.05,1,0),dat_major$L95v2)
roc.multi <- multiclass.roc(dat_major$pGroup3,dat_major$L95v2)
rs <- roc.multi[['rocs']]
plot.roc(roc1.1)
lines.roc(roc2.1,col="blue")
lines.roc(rs[[3]],col="green")
legend("bottomright", legend=c("p<=0.005","p<=0.05","<0.05 vs < 0.005"),
       col=c("black","blue","green"), lwd=2)
wilcox.test(dat_major$L95v2~ifelse(dat_major$pvalue<=0.005,1,0))
wilcox.test(dat_major$L95v2~ifelse(dat_major$pvalue<=0.05,1,0))

roc1.2 <- roc(ifelse(dat_major$pvalue<=0.005,1,0),dat_major$CIwidth2)
roc2.2 <- roc(ifelse(dat_major$pvalue<=0.05,1,0),dat_major$CIwidth2)
plot.roc(roc1.2)
lines.roc(roc2.2,col="blue")
legend("bottomright", legend=c("p<=0.005","p<=0.05"),
       col=c("black","blue"), lwd=2)
wilcox.test(dat_major$CIwidth2~ifelse(dat_major$pvalue<=0.005,1,0))
wilcox.test(dat_major$CIwidth2~ifelse(dat_major$pvalue<=0.05,1,0))

summary(dat_major[dat_major$pvalue<=0.05,]$effectEst2)
summary(dat_major[dat_major$pvalue<=0.005,]$effectEst2)
summary(dat_major[dat_major$pGroup2=="(0.005,0.05]",]$effectEst2)

summary(dat_major[dat_major$pvalue<=0.05,]$L95v2)
summary(dat_major[dat_major$pvalue<=0.005,]$L95v2)
summary(dat_major[dat_major$pGroup2=="(0.005,0.05]",]$L95v2)

summary(dat_major[dat_major$pvalue<=0.05,]$CIwidth2)
summary(dat_major[dat_major$pvalue<=0.005,]$CIwidth2)
wilcox.test(dat_major[dat_major$pvalue<=0.05,]$CIwidth2,dat_major[dat_major$pvalue<=0.005,]$CIwidth2)
