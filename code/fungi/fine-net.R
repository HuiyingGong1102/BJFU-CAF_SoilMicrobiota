library(orthopolynom)
library(pbapply)
library(glmnet)
library(ggplot2)
source("F:/CAF/2022-11-10/All.function3.R")
source("F:/CAF/2022-11-10/lop.fit.R")

step_fct<-function(x,step)
{
  t<-length(x)
  x_new<-x[1]
  for (i in 1:(t-1))
  {
    x_int<-rev(seq(x[i+1],x[i],by=-step))
    x_new<-c(x_new,x_int)
  }
  x_new
}

get_after <- function(i){
  temp <- matrix(NA,nrow = length(i[[2]]),ncol=3)
  temp[,1] <- i[[2]]
  temp[,2] <- i[[1]]
  temp[,3] <- i[[5]][2:(length(i[[2]])+1)]
  
  colnames(temp) <- c('source','target','weight')
  temp <- data.frame(temp)
  temp[,3] <- as.numeric(as.character(temp[,3]))
  return(temp)
}

get_extra <- function(i){
  temp <- i[[5]][1]
  return(temp)
}

k <- 11
load("F:/CAF/2022-10-20/fungi/cluster-result/k11.RData")
module<-return_object$clustered_data[,181]

#Calculate the sum of OTUs for each category
data_clustered<-rep()
for (i in 1:k)
{
  new_i<-apply(all.data[which(module==i),],2,sum)
  data_clustered<-rbind(data_clustered,new_i)
}
data_clustered<-t(data_clustered)

load("F:/CAF/2022-11-10/fungi/lop_par_all.RData")

category <- 6
x_int <- xy[[category]]$X
#x_smooth<-seq(x_int[1],x_int[30],length.out = 50)
step<-0.001
x_smooth <- step_fct(x_int,step)
x_num <- c()
for (i in 1:length(x_int)) {
  x_num <- c(x_num,which(x_smooth==x_int[i]))
}

data_whole <- data_clustered[(1+30*(category-1)):(30*category),]
colnames(data_whole) <- paste("M",1:k,sep = "")

#(least square method (The initial value is positive number))
data_smooth <- rep()
data_orig <- rep()
for (i in 1:k) {
  yy <- f(data_whole[,i])
  plot(x_int,yy,ylim = c(-0.1,1))
  
  l_par_i <- lop_par_all[[category]][i,]
  
  data_i <- legendre_fit(l_par_i,x_int)
  lines(x_int,data_i)
  data_i2 <- legendre_fit(l_par_i,x_smooth)
  lines(x_smooth,data_i2)
  
  data_smooth <- cbind(data_smooth,data_i2)
  data_orig <- cbind(data_orig,data_i)
}

dim(data_smooth)
dim(data_orig)

colnames(data_smooth) <- paste("M",1:k,sep = "")
colnames(data_orig) <- paste("M",1:k,sep = "")

relationship_all <- pblapply(1:ncol(data_smooth),function(c)
  get_interaction(data_smooth, c, reduction = FALSE))


par_all_obj <- list(lop_par = lop_par,
                    relationship = relationship_all,
                    dataset = data_smooth,
                    times = x_smooth,
                    order = 3)

module_net1 <- get_all_net(par_all_obj)
get_decomposition_plot(par_all_obj, module_net1, 2)


extra <- sapply(module_net1,get_extra)
links_gene <- do.call(rbind,lapply(module_net1, get_after))
links_gene<-as.data.frame(links_gene)
links_gene$edge_type<-ifelse(links_gene$weight>0,1,2)
links_gene$weight <- abs(as.numeric(links_gene$weight))
rownames(links_gene)<-rep(1:nrow(links_gene))

snode <- c()
tnode <- c()
for (i in 1:dim(links_gene)[1]) {
  snode <- c(snode,extra[which(colnames(data_whole)==links_gene$source[i])])
  tnode <- c(tnode,extra[which(colnames(data_whole)==links_gene$target[i])])
}
links_gene$snode <- snode
links_gene$tnode <- tnode
write.csv(links_gene, file = "E:/CAF/2022-11-10/fungi/fig6/T45R/links_gene.csv",row.names = F,quote=F)
table(c(links_gene$source,links_gene$target))
dim(links_gene)





######OTU network#####

k <- 11
load("F:/CAF/2022-10-20/fungi/cluster-result/k11.RData")
module<-return_object$clustered_data[,181]
table(module)

index <- which(module==6)
rownames(all.data)[index]
data_clustered <- t(all.data[index,])

category <- 6
x_int <- xy[[category]]$X
#x_smooth<-seq(x_int[1],x_int[30],length.out = 50)
step<-0.001
x_smooth <- step_fct(x_int,step)
x_num <- c()
for (i in 1:length(x_int)) {
  x_num <- c(x_num,which(x_smooth==x_int[i]))
}

data_whole <- data_clustered[(1+30*(category-1)):(30*category),]

#(least square method (The initial value is positive number))
data_smooth <- rep()
data_orig <- rep()
for (i in 1:ncol(data_clustered)) {
  yy <- f(data_whole[,i])
  plot(x_int,yy,ylim = c(-0.1,1))
  
  l_par_i <- get_legendre_par2(yy,3,x_int)
  
  data_i <- legendre_fit(l_par_i,x_int)
  lines(x_int,data_i)
  data_i2 <- legendre_fit(l_par_i,x_smooth)
  lines(x_smooth,data_i2)
  
  data_smooth <- cbind(data_smooth,data_i2)
  data_orig <- cbind(data_orig,data_i)
}

dim(data_smooth)
dim(data_orig)

otu_names <- c()
for (i in 1:ncol(data_clustered)) {
  otu_names <- c(otu_names,substr(colnames(data_whole)[i],4,nchar(colnames(data_whole)[i])))
}
colnames(data_smooth) <- as.numeric(otu_names)
colnames(data_orig) <- as.numeric(otu_names)

relationship_all <- pblapply(1:ncol(data_smooth),function(c)
  get_interaction(data_smooth, c, reduction = FALSE))

#T45R M6
relationship_all[[1]] <- get_interaction(data_smooth[,-6], 1, reduction = FALSE)
relationship_all[[2]] <- get_interaction(data_smooth[,-6], 2, reduction = FALSE)
relationship_all[[3]] <- get_interaction(data_smooth[,-6], 3, reduction = FALSE)
relationship_all[[4]] <- get_interaction(data_smooth[,-6], 4, reduction = FALSE)
relationship_all[[5]] <- get_interaction(data_smooth[,-6], 5, reduction = FALSE)
relationship_all[[6]] <- get_interaction(data_smooth, 6, reduction = FALSE)
relationship_all[[7]] <- get_interaction(data_smooth[,-6], 6, reduction = FALSE)
save(relationship_all,file = "F:/CAF/2022-11-10/fungi/fig6/T45R/M6-relationship.RData")


#T45R M2
# relationship_all[[1]] <- get_interaction(data_smooth[,-c(8,11,13)], 1, reduction = FALSE)
# relationship_all[[2]] <- get_interaction(data_smooth[,-c(8,11,13)], 2, reduction = FALSE)
# relationship_all[[3]] <- get_interaction(data_smooth[,-c(8,11,13)], 3, reduction = FALSE)
# relationship_all[[4]] <- get_interaction(data_smooth[,-c(8,11,13)], 4, reduction = FALSE)
# relationship_all[[5]] <- get_interaction(data_smooth[,-c(8,11,13)], 5, reduction = FALSE)
# relationship_all[[6]] <- get_interaction(data_smooth[,-c(8,11,13)], 6, reduction = FALSE)
# relationship_all[[7]] <- get_interaction(data_smooth[,-c(8,11,13)], 7, reduction = FALSE)
# relationship_all[[8]] <- get_interaction(data_smooth[,-c(11,13)], 8, reduction = FALSE)
# relationship_all[[9]] <- get_interaction(data_smooth[,-c(8,11,13)], 8, reduction = FALSE)
# relationship_all[[10]] <- get_interaction(data_smooth[,-c(8,11,13)], 9, reduction = FALSE)
# relationship_all[[11]] <- get_interaction(data_smooth[,-c(8,13)], 10, reduction = FALSE)
# relationship_all[[12]] <- get_interaction(data_smooth[,-c(8,11,13)], 10, reduction = FALSE)
# relationship_all[[13]] <- get_interaction(data_smooth[,-c(8,11)], 11, reduction = FALSE)
# relationship_all[[14]] <- get_interaction(data_smooth[,-c(8,11,13)], 11, reduction = FALSE)


lop_par <- pblapply(1:ncol(data_smooth),function(c)
  get_value(data_smooth,relationship_all[[c]],x_smooth,3))
save(lop_par,file = "F:/CAF/2022-11-10/fungi/fig6/T45R/M6-0.00001-par6.RData")
#load("E:/CAF/2022-10-1/fungi/10-11/fig6/M10/a=0.00001-lop_par5.RData")

order <- 3
library(parallel)
core.number <- detectCores()
cl <- makeCluster(getOption("cl.cores", 8))
clusterEvalQ(cl, {require(orthopolynom)})
clusterExport(cl, c("data_smooth","x_smooth","order","relationship_all",
                    "get_value","ode_optimize","get_effect"),envir=environment())
lop_par <- pblapply(1:ncol(data_smooth),function(c)
  get_value(data_smooth,relationship_all[[c]],x_smooth,order),cl=cl)
stopCluster(cl)


par_all_obj <- list(lop_par = lop_par,
                    relationship = relationship_all,
                    dataset = data_smooth,
                    times = x_smooth,
                    order = 3)

module_net1 <- get_all_net(par_all_obj)
get_decomposition_plot(par_all_obj, module_net1, 2)

extra <- sapply(module_net1,get_extra)
links_gene <- do.call(rbind,lapply(module_net1, get_after))
links_gene<-as.data.frame(links_gene)
links_gene$edge_type<-ifelse(links_gene$weight>0,1,2)
links_gene$weight <- abs(as.numeric(links_gene$weight))
rownames(links_gene)<-rep(1:nrow(links_gene))

snode <- c()
tnode <- c()
for (i in 1:dim(links_gene)[1]) {
  snode <- c(snode,extra[which(otu_names==links_gene$source[i])])
  tnode <- c(tnode,extra[which(otu_names==links_gene$target[i])])
}
links_gene$snode <- snode
links_gene$tnode <- tnode
write.csv(links_gene, file = "F:/CAF/2022-11-10/fungi/fig6/T45R/M6-links_gene.csv", row.names = F,quote=F)
table(c(links_gene$source,links_gene$target))
dim(links_gene)
