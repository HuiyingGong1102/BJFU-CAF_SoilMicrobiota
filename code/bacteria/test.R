setwd("E:/CAF/2022-10-20/Bacteria")
library("openxlsx")
library(mvtnorm)
properties<-read.xlsx("E:/CAF/CAF/1 soil_properties.xlsx")[-(181:210),]
bacteria_otu <- read.xlsx("E:/CAF/CAF/3 bacteria_out.xlsx")[-(181:210),]
Treatment <- bacteria_otu$Treatment
#row sum of OTU
#apply(bacteria_otu[,-(1:6)], 1, sum)

relative_abu <- apply(bacteria_otu[,-(1:6)], 2, function(x) sum(x))
#20535*180*0.0005

#Screening samples based on relative abundance
#data normalization
f <- function(x){
  if(length(which(x!=0))>0){
    return((x - min(x))/(max(x) - min(x)))
  }else{
    return(x)
  }
}

yi <- list()
for (i in 1:length(table(Treatment))) {
  index <- which(Treatment==names(table(Treatment))[i])
  yi[[i]] <- apply(bacteria_otu[index,-(1:6)][,which(relative_abu>=1848)], 2, f)
  
}


exp_index <- c()
xy <- list()
order_properties <- list()
for (i in 1:length(table(Treatment))) {
  index <- which(Treatment==names(table(Treatment))[i])
  total<-apply(yi[[i]],1,sum)
  ei_index<-total[order(total)]
  X<-log10(ei_index)
  yin <- cbind(bacteria_otu[index,1:6],yi[[i]])
  data_gene_exp_order<-yin[order(total),]
  order_properties[[i]] <- properties[index,][order(total),]
  xy[[i]]<-cbind(X,data_gene_exp_order)
  exp_index <- c(exp_index,ei_index)
}


l <- dim(xy[[1]])[2]
x <- c(xy[[1]]$X,xy[[2]]$X,xy[[3]]$X,xy[[4]]$X,xy[[5]]$X,xy[[6]]$X)
all.data <- rbind(xy[[1]][,8:l],xy[[2]][,8:l],xy[[3]][,8:l],
                  xy[[4]][,8:l],xy[[5]][,8:l],xy[[6]][,8:l])
all.data <- t(all.data)
dim(all.data)






