setwd("F:/CAF/2022-11-10/fungi")
library("openxlsx")
library(mvtnorm)
set.seed(1324)

properties<-read.xlsx("F:/CAF/CAF/1 soil_properties.xlsx")[-(181:210),]
fungi_otu <- read.xlsx("F:/CAF/CAF/2 fungi_otu.xlsx")[-(181:210),]
Treatment <- fungi_otu$Treatment
#row sum of OTU
#apply(fungi_otu[,-(1:6)], 1, sum)

relative_abu <- apply(fungi_otu[,-(1:6)], 2, function(x) sum(x))
#33951*180*0.0005

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
  yi[[i]] <- apply(fungi_otu[index,-(1:6)][,which(relative_abu>=3055)], 2, f)
  
}

exp_index <- c()
xy <- list()
order_properties <- list()
for (i in 1:length(table(Treatment))) {
  index <- which(Treatment==names(table(Treatment))[i])
  total<-apply(yi[[i]],1,sum)
  ei_index<-total[order(total)]
  X<-log10(ei_index)
  yin <- cbind(fungi_otu[index,1:6],yi[[i]])
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



