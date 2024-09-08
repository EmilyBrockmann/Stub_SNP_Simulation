library(updog)
library(dplyr)
library(seqinr)
library(tibble)
library(ldsep)
library(purrr)
library(latex2exp)
library(truncnorm)
library(ggplot2)
library(cowplot)
library(xtable)
#Import path
setwd("c:/Users/mimib/Documents/Uni/Masterarbeit/Stub_NGS_simulation/Simulation_true_positives/src/")
pathDat <- "../../../Data/simulatedData/vcfSimulation/Parameters/"

#PLOT PARAMETERS
width <- 15
height <- 10
dpi <- 500



#####_P1_#####
load(file = paste0(pathDat, "p1/Analysis.RData"))
changeIndvProbBinomVec <- seq(0.0125,0.1,0.0125)

corrPlots <- list()
for(i in seq_along(corrList)){
  titlePlot <- paste0("$p_{shift} = ", changeIndvProbBinomVec[i], "$")
  p <- ggplot(data=corrList[[i]],aes(corr))+
    geom_histogram(alpha=0.8,bins=floor(sqrt(1000)),aes(y = after_stat(count / sum(count))), fill="cadetblue")+
    labs(title = TeX(titlePlot))+xlim(c(0,1))+theme_classic()+
    xlab(TeX("$r^2$"))+ylab("density")
  print(p)
  corrPlots[[length(corrPlots)+1]] <- p
}
corrPlot <- plot_grid(plotlist = corrPlots, ncol=2)
ggsave(plot=corrPlot, file=paste0("../Results/r1_linkage.jpeg"), width=width, height = height, unit="cm", dpi=dpi)


#densityPlot
plotCorrDensities <- function(dataList, adjust){
  # :param dataList: named list or 1D vectors
  densityTable <- tibble(densi=numeric(), label=character())
  for(i in names(dataList)){
    # tmp0 <- unlist(dataList[[i]])
    # print(tmp0)
    tmp <- tibble(densi=unlist(dataList[[i]]$corr))
    tmp <- tmp %>% mutate(label=i)
    densityTable <- densityTable %>% add_row(tmp)
  }
  densityTable$label <- factor(densityTable$label)
  
  #colors in plot 
  blues <- colorRampPalette(c("#E0FFFF","cadetblue","#2F4F4F"))(8)
  
  #plot densities overlapping
  p <- ggplot(data=densityTable, aes(x=densi, color=label, fill=label))+
    geom_density(alpha=0.1,adjust=adjust, size=.9)+theme_classic()+
    xlab(TeX("$r^2$"))+ylab("density")+scale_colour_manual(values=blues)+
    scale_fill_manual(values=blues)+labs(fill=TeX("$p_{shift}$"),color=TeX("$p_{shift}$"))
  return(p)
}

names(corrList) <- changeIndvProbBinomVec
p <- plotCorrDensities(corrList,6)
ggsave(plot=p, file=paste0("../Results/r1_linkage.jpeg"), width=width, height = height, unit="cm", dpi=dpi)



names(datComp) <- c("$p_{shift}$", "mean linkage", "median linkage", "# Haplotypes")

sink("../Results/determine_p1.tex")
texTab <- xtable(datComp, digits = c(0,4,rep(3,2), 0))
print(texTab, include.rownames = FALSE,  # Exclude row names
      label = "tab:sample")
sink()


##linear regression
lm(datComp$mean~datComp$p1)
cor(datComp$mean,datComp$p1)
pLinreg <- ggplot(datComp, aes(x=p1, y=mean))+geom_point(fill="cadetblue", color="cadetblue", size=1.5)+ 
  ylab(TeX("mean $r^2$"))+xlab(TeX("$p_{shift}$"))+theme_classic()+
  geom_smooth(method='lm', formula= y~x, color="darkslategrey", alpha=0.1)
pLinreg
ggsave(plot=pLinreg, file=paste0("../Results/p1_linreg.jpeg"), width=width, height = height, unit="cm", dpi=dpi)

  
# write.csv2(datComp, file = "../Results/determine_p1.csv", row.names = F)
#####_P4_#####
load(file = paste0(pathDat, "p4/Analysis.RData"))
pUncorrelatedVec <- seq(0,0.02,0.001)
corrPlots <- list()
for(i in seq_along(corrList)){
  titlePlot <- paste0("$p_{uncorr} = ", pUncorrelatedVec[i], "$")
  p <- ggplot(data=corrList[[i]],aes(corr))+
    geom_histogram(alpha=0.8,bins=floor(sqrt(1000)),aes(y = after_stat(count / sum(count))), fill="cadetblue")+
    labs(title = TeX(titlePlot))+xlim(c(0,1))+theme_classic()+
    xlab(TeX("$r^2$"))+ylab("density")
  print(p)
  corrPlots[[length(corrPlots)+1]] <- p
}

corrPlot <- plot_grid(plotlist = corrPlots[seq(1,length(corrPlots), 2)], ncol=2)
ggsave(plot=corrPlot, file=paste0("../Results/r2_linkage.jpeg"), width=width, height = height, unit="cm", dpi=dpi)

#p4 is renamed p2 in thesis due to removal of inversion concept in simulation
sink("../Results/determine_p2.tex")
# names(datComp) <- c("$\displaystyle p_{uncorr}$", "mean linkage", "median linkage", "# Haplotypes")
names(datComp) <- c("$p_{uncorr}$", "mean linkage", "median linkage", "# Haplotypes")


texTab <- xtable(datComp, digits = c(rep(3,4), 0))
print(texTab, include.rownames = FALSE,  # Exclude row names
      label = "tab:sample")
sink()
# write.csv2(datComp, file = "../Results/determine_p4.csv", row.names = F)