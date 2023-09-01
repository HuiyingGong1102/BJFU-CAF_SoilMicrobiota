library(orthopolynom)
library(pbapply)
library(glmnet)
library(ggplot2)
source("E:/CAF/2022-10-1/fungi/005-figure/network/All.function3.R")
source("E:/CAF/2022-10-1/fungi/005-figure/network/lop.fit.R")
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


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

k <- 11
load("E:/CAF/2022-10-20/fungi/cluster-result/k11.RData")
module<-return_object$clustered_data[,181]
table(module)

#Calculate the sum of OTUs for each category
data_clustered<-rep()
for (i in 1:k)
{
  new_i<-apply(all.data[which(module==i),],2,sum)
  data_clustered<-rbind(data_clustered,new_i)
}
data_clustered<-t(data_clustered)

category <- 3
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

colnames(data_smooth) <- paste("M",1:k,sep = "")
colnames(data_orig) <- paste("M",1:k,sep = "")

relationship_all <- pblapply(1:ncol(data_smooth),function(c)
  get_interaction(data_smooth, c, reduction = FALSE))

lop_par <- pblapply(1:ncol(data_smooth),function(c)
  get_value(data_smooth,relationship_all[[c]],x_smooth,3))

par_all_obj <- list(lop_par = lop_par,
                    relationship = relationship_all,
                    dataset = data_smooth,
                    times = x_smooth,
                    order = 3)

module_net1 <- get_all_net(par_all_obj)
get_decomposition_plot(par_all_obj, module_net1, 2)

col <- c('#ff7171','#0dceda','#9ede73')
gene_whole<-rep()
self_size<-rep()
plist <- list()
for (i in 1:ncol(data_whole)) {
  # cluster_mean <- par_all_obj$dataset[x_num,]
  # plot_data <- data.frame(module_net1[[i]][[6]][x_num,],check.names = F)
  # colnames(plot_data)[-1] <- module_net1[[i]][[2]]
  # plot_data$sum <- rowSums(plot_data)
  # plot_data$origin <- as.numeric(cluster_mean[,i])
  # plot_data$x <- x_int
  
  plot_data <- get_curve_data(par_all_obj, module_net1, i)
  ind_name <- module_net1[[i]][[1]]
  
  fitted_gene_matrix <- data.frame(module_net1[[i]][[6]][x_num,2:length(module_net1[[i]][[5]])],check.names = F)
  gene_cluster<-cbind(module_net1[[i]][[2]],module_net1[[i]][[1]],t(fitted_gene_matrix))
  gene_whole<-rbind(gene_whole,gene_cluster)
  self_size<-rbind(self_size,module_net1[[i]][[6]][x_num,1])
  
  d <- plot_data
  p1 <- ggplot(d,aes(x=x))+geom_line(d,mapping=aes(x=x,y=origin),color=col[2],size=1.5)+
    geom_line(d,mapping=aes(x=x,y=ind_effect),color=col[1],size=1.4)
  p1 <- p1+geom_line(d,mapping=aes(x=x,y=sum),color="black",size=1.4)
  for (j in 1:(ncol(d)-4)) {
    p1 <- p1+geom_line(aes_string(y=d[,1+j]),color=col[3],size=1.2,alpha=1)+
      annotate('text',x=max(as.numeric(d$x))+0.03,y=d[nrow(d),(j+1)],
               label=colnames(d)[j+1],size=3)
    
  }
  p1 <- p1+ theme(panel.grid = element_blank(),
                  panel.background = element_rect(color = 'black', fill = 'transparent'),
                  legend.title = element_blank()) +
    xlab("Habitat Index") + ylab('Effect')+
    theme(plot.margin = unit(c(0,0,0,0),"lines"))+
    geom_hline(yintercept = 0, size = 0.6) +
    annotate('text',x=mean(as.numeric(d$x)),y=round(max(d[,-ncol(d)])),label=ind_name,size=5)
  plist[[i]] <- p1
}


# pdf("E:/Doctoral/CAF/2022-10-1/Bacteria/fig3/T45R/T45R-a=0.001.pdf",height=2,width=15)
# multiplot(plotlist=plist,cols=5)
# dev.off()

dim(gene_whole)
dim(self_size)

affect_value_April<-apply(gene_whole[,which(xy[[category]][,7]=="April")+2],2,as.numeric)
value_April<-apply(affect_value_April,1,mean)

affect_value_May<-apply(gene_whole[,which(xy[[category]][,7]=="May")+2],2,as.numeric)
value_May<-apply(affect_value_May,1,mean)

affect_value_June<-apply(gene_whole[,which(xy[[category]][,7]=="June")+2],2,as.numeric)
value_June<-apply(affect_value_June,1,mean)

affect_value_July<-apply(gene_whole[,which(xy[[category]][,7]=="July")+2],2,as.numeric)
value_July<-apply(affect_value_July,1,mean)

affect_value_August<-apply(gene_whole[,which(xy[[category]][,7]=="August")+2],2,as.numeric)
value_August<-apply(affect_value_August,1,mean)

affect_value_September<-apply(gene_whole[,which(xy[[category]][,7]=="September")+2],2,as.numeric)
value_September<-apply(affect_value_September,1,mean)



ave_self_April<-apply(self_size[,which(xy[[category]][,7]=="April")],1,mean)
ave_self_May<-apply(self_size[,which(xy[[category]][,7]=="May")],1,mean)
ave_self_June<-apply(self_size[,which(xy[[category]][,7]=="June")],1,mean)
ave_self_July<-apply(self_size[,which(xy[[category]][,7]=="July")],1,mean)
ave_self_August<-apply(self_size[,which(xy[[category]][,7]=="August")],1,mean)
ave_self_September<-apply(self_size[,which(xy[[category]][,7]=="September")],1,mean)

####################################################
####construct the data we need to plot for  April->May of CKN   
####################################################

value_diff<-value_May-value_April
gene_whole_threhold<-gene_whole
range(abs(value_diff))
gene_whole_threhold<-gene_whole_threhold[abs(value_diff)>0.01,]
dim(gene_whole_threhold)

links_gene<-cbind(gene_whole_threhold[,1:2],value_diff[abs(value_diff)>0.01])
colnames(links_gene)<-c("source","target","weight")
# links_gene<-c(gene_whole_threhold[1:2],value_diff[abs(value_diff)>0.01])
# links_gene <- as.matrix(links_gene)
# links_gene <- t(links_gene)
colnames(links_gene)<-c("source","target","weight")
links_gene<-as.data.frame(links_gene)
links_gene$edge_type<-ifelse(links_gene$weight>0,1,2)
links_gene$weight <- abs(as.numeric(links_gene$weight))
rownames(links_gene)<-rep(1:nrow(links_gene))

# self_index<-unique(c(gene_whole_threhold[,1],gene_whole_threhold[,2]))
self_index<-unique(c(gene_whole_threhold[1],gene_whole_threhold[2]))
snode <- c()
tnode <- c()
for (i in 1:dim(links_gene)[1]) {
  snode <- c(snode,ave_self_May[which(colnames(data_whole)==links_gene$source[i])])
  tnode <- c(tnode,ave_self_May[which(colnames(data_whole)==links_gene$target[i])])
}
links_gene$snode <- snode
links_gene$tnode <- tnode
# links_gene <- paste("M",1:9,sep = "")
# tnode <- c()
# for (i in 1:length(links_gene)) {
#   tnode <- c(tnode,ave_self_May[which(colnames(data_whole)==links_gene[i])])
# }
# links_gene <- cbind(links_gene,tnode)
# colnames(links_gene)<-c("source","tnode")
# links_gene<-as.data.frame(links_gene)
write.csv(links_gene, file = "E:/CAF/2022-10-20/fungi/fig3/T45R/links_gene_May.csv", row.names = F,quote=F)
table(gene_whole_threhold[,1:2])
dim(gene_whole_threhold)
####################################################
####construct the data we need to plot for  April->June of CKN   
####################################################


value_diff<-value_June-value_April
gene_whole_threhold<-gene_whole
range(abs(value_diff))
gene_whole_threhold<-gene_whole_threhold[abs(value_diff)>0.01,]
dim(gene_whole_threhold)

links_gene<-cbind(gene_whole_threhold[,1:2],value_diff[abs(value_diff)>0.01])
colnames(links_gene)<-c("source","target","weight")
links_gene<-as.data.frame(links_gene)
links_gene$edge_type<-ifelse(links_gene$weight>0,1,2)
links_gene$weight <- abs(as.numeric(links_gene$weight))
rownames(links_gene)<-rep(1:nrow(links_gene))

self_index<-unique(c(gene_whole_threhold[,1],gene_whole_threhold[,2]))
snode <- c()
tnode <- c()
for (i in 1:dim(links_gene)[1]) {
  snode <- c(snode,ave_self_June[which(colnames(data_whole)==links_gene$source[i])])
  tnode <- c(tnode,ave_self_June[which(colnames(data_whole)==links_gene$target[i])])
}
links_gene$snode <- snode
links_gene$tnode <- tnode
write.csv(links_gene, file = "E:/CAF/2022-10-20/fungi/fig3/T45R/links_gene_June.csv", row.names = F,quote=F)
table(gene_whole_threhold[,1:2])
dim(gene_whole_threhold)
####################################################
####construct the data we need to plot for  April->July of CKN   
####################################################

value_diff<-value_July-value_April
gene_whole_threhold<-gene_whole
range(abs(value_diff))
gene_whole_threhold<-gene_whole_threhold[abs(value_diff)>0.01,]
dim(gene_whole_threhold)

links_gene<-cbind(gene_whole_threhold[,1:2],value_diff[abs(value_diff)>0.01])
colnames(links_gene)<-c("source","target","weight")
links_gene<-as.data.frame(links_gene)
links_gene$edge_type<-ifelse(links_gene$weight>0,1,2)
links_gene$weight <- abs(as.numeric(links_gene$weight))
rownames(links_gene)<-rep(1:nrow(links_gene))

self_index<-unique(c(gene_whole_threhold[,1],gene_whole_threhold[,2]))
snode <- c()
tnode <- c()
for (i in 1:dim(links_gene)[1]) {
  snode <- c(snode,ave_self_July[which(colnames(data_whole)==links_gene$source[i])])
  tnode <- c(tnode,ave_self_July[which(colnames(data_whole)==links_gene$target[i])])
}
links_gene$snode <- snode
links_gene$tnode <- tnode
# links_gene <- paste("M",1:11,sep = "")
# tnode <- c()
# for (i in 1:length(links_gene)) {
#   tnode <- c(tnode,ave_self_July[which(colnames(data_whole)==links_gene[i])])
# }
# links_gene <- cbind(links_gene,tnode)
# colnames(links_gene)<-c("source","tnode")
# links_gene<-as.data.frame(links_gene)
write.csv(links_gene, file = "E:/CAF/2022-10-20/fungi/fig3/T45R/links_gene_July.csv", row.names = F,quote=F)
table(gene_whole_threhold[,1:2])
dim(gene_whole_threhold)
####################################################
####construct the data we need to plot for  April->August of CKN   
####################################################

value_diff<-value_August-value_April
gene_whole_threhold<-gene_whole
range(abs(value_diff))
gene_whole_threhold<-gene_whole_threhold[abs(value_diff)>0.01,]
dim(gene_whole_threhold)

links_gene<-cbind(gene_whole_threhold[,1:2],value_diff[abs(value_diff)>0.01])
#links_gene<-c(gene_whole_threhold[1:2],value_diff[abs(value_diff)>0.01])
colnames(links_gene)<-c("source","target","weight")
#names(links_gene)<-c("source","target","weight")
links_gene<-as.data.frame(links_gene)
#links_gene<-as.data.frame(t(links_gene))
links_gene$edge_type<-ifelse(links_gene$weight>0,1,2)
links_gene$weight <- abs(as.numeric(links_gene$weight))
rownames(links_gene)<-rep(1:nrow(links_gene))

self_index<-unique(c(gene_whole_threhold[,1],gene_whole_threhold[,2]))
#self_index<-unique(c(gene_whole_threhold[1],gene_whole_threhold[2]))
snode <- c()
tnode <- c()
for (i in 1:dim(links_gene)[1]) {
  snode <- c(snode,ave_self_August[which(colnames(data_whole)==links_gene$source[i])])
  tnode <- c(tnode,ave_self_August[which(colnames(data_whole)==links_gene$target[i])])
}
links_gene$snode <- snode
links_gene$tnode <- tnode
# links_gene <- paste("M",1:9,sep = "")
# tnode <- c()
# for (i in 1:length(links_gene)) {
#   tnode <- c(tnode,ave_self_August[which(colnames(data_whole)==links_gene[i])])
# }
# links_gene <- cbind(links_gene,tnode)
# colnames(links_gene)<-c("source","tnode")
# links_gene<-as.data.frame(links_gene)
write.csv(links_gene, file = "E:/CAF/2022-10-20/fungi/fig3/T45R/links_gene_August.csv", row.names = F,quote=F)
table(gene_whole_threhold[,1:2])
dim(gene_whole_threhold)
####################################################
####construct the data we need to plot for  April->September of CKN   
####################################################

value_diff<-value_September-value_April
gene_whole_threhold<-gene_whole
range(abs(value_diff))
gene_whole_threhold<-gene_whole_threhold[abs(value_diff)>0.01,]
dim(gene_whole_threhold)

links_gene<-cbind(gene_whole_threhold[,1:2],value_diff[abs(value_diff)>0.01])
colnames(links_gene)<-c("source","target","weight")
links_gene<-as.data.frame(links_gene)
links_gene$edge_type<-ifelse(links_gene$weight>0,1,2)
links_gene$weight <- abs(as.numeric(links_gene$weight))
rownames(links_gene)<-rep(1:nrow(links_gene))

self_index<-unique(c(gene_whole_threhold[,1],gene_whole_threhold[,2]))
snode <- c()
tnode <- c()
for (i in 1:dim(links_gene)[1]) {
  snode <- c(snode,ave_self_September[which(colnames(data_whole)==links_gene$source[i])])
  tnode <- c(tnode,ave_self_September[which(colnames(data_whole)==links_gene$target[i])])
}
links_gene$snode <- snode
links_gene$tnode <- tnode
write.csv(links_gene, file = "E:/CAF/2022-10-20/fungi/fig3/T45R/links_gene_September.csv", row.names = F,quote=F)
table(gene_whole_threhold[,1:2])
dim(gene_whole_threhold)









