for (j in 1:8) {
  filename <- paste("E:/CAF/2022-10-1/Bacteria/all-fit-curves/fitcurve",j,".tiff",sep = "")
  tiff(filename,width=32,height=20,units = "cm",res=300)
  
  par(mfrow=c(5,8))
  for (example_index in (1+40*(j-1)):(40*j)) {
    power_par <- c()
    initial_power_par <- c(0.5,0.5)
    for (i in 1:6) {
      power_par <- rbind(power_par,get_power_par(initial_power_par,xy[[i]][,8:l][,example_index],xy[[i]]$X))
    }
    
    par(mar=c(1,1,1,1))
    plot(xy[[1]]$X,xy[[1]][,8:l][,example_index],type = "p",pch="+",col="lightcoral",
         xlim = c(range(x)[1]-0.01,range(x)[2]+0.01),ylim = c(0,1.1))
    points(xy[[2]]$X,xy[[2]][,8:l][,example_index],pch="+",col="steelblue1")
    points(xy[[3]]$X,xy[[3]][,8:l][,example_index],pch=19,col="lightcoral")
    points(xy[[4]]$X,xy[[4]][,8:l][,example_index],pch=19,col="steelblue1")
    points(xy[[5]]$X,xy[[5]][,8:l][,example_index],pch=17,col="lightcoral")
    points(xy[[6]]$X,xy[[6]][,8:l][,example_index],pch=17,col="steelblue1")
    
    lines(xy[[1]]$X,power.equation(power_par[1,],xy[[1]]$X),col="red",lwd=2)
    lines(xy[[2]]$X,power.equation(power_par[2,],xy[[2]]$X),col="blue",lwd=2)
    lines(xy[[3]]$X,power.equation(power_par[3,],xy[[3]]$X),col="red",lwd=2)
    lines(xy[[4]]$X,power.equation(power_par[4,],xy[[4]]$X),col="blue",lwd=2)
    lines(xy[[5]]$X,power.equation(power_par[5,],xy[[5]]$X),col="red",lwd=2)
    lines(xy[[6]]$X,power.equation(power_par[6,],xy[[6]]$X),col="blue",lwd=2)
    
    text(range(x)[1]-0.01+(range(x)[2]-range(x)[1]+0.02)/2,1.08,colnames(xy[[2]][,8:l])[example_index])
    
  }
  
  dev.off()
}


filename <- paste("E:/CAF/2022-10-1/Bacteria/all-fit-curves/fitcurve",9,".tiff",sep = "")
tiff(filename,width=32,height=20,units = "cm",res=300)

par(mfrow=c(5,8))
for (example_index in 321:356 ) {
  power_par <- c()
  initial_power_par <- c(0.5,0.5)
  for (i in 1:6) {
    power_par <- rbind(power_par,get_power_par(initial_power_par,xy[[i]][,8:l][,example_index],xy[[i]]$X))
  }
  
  par(mar=c(1,1,1,1))
  plot(xy[[1]]$X,xy[[1]][,8:l][,example_index],type = "p",pch="+",col="lightcoral",
       xlim = c(range(x)[1]-0.01,range(x)[2]+0.01),ylim = c(0,1.1))
  points(xy[[2]]$X,xy[[2]][,8:l][,example_index],pch="+",col="steelblue1")
  points(xy[[3]]$X,xy[[3]][,8:l][,example_index],pch=19,col="lightcoral")
  points(xy[[4]]$X,xy[[4]][,8:l][,example_index],pch=19,col="steelblue1")
  points(xy[[5]]$X,xy[[5]][,8:l][,example_index],pch=17,col="lightcoral")
  points(xy[[6]]$X,xy[[6]][,8:l][,example_index],pch=17,col="steelblue1")
  
  lines(xy[[1]]$X,power.equation(power_par[1,],xy[[1]]$X),col="red",lwd=2)
  lines(xy[[2]]$X,power.equation(power_par[2,],xy[[2]]$X),col="blue",lwd=2)
  lines(xy[[3]]$X,power.equation(power_par[3,],xy[[3]]$X),col="red",lwd=2)
  lines(xy[[4]]$X,power.equation(power_par[4,],xy[[4]]$X),col="blue",lwd=2)
  lines(xy[[5]]$X,power.equation(power_par[5,],xy[[5]]$X),col="red",lwd=2)
  lines(xy[[6]]$X,power.equation(power_par[6,],xy[[6]]$X),col="blue",lwd=2)
  
  text(range(x)[1]-0.01+(range(x)[2]-range(x)[1]+0.02)/2,1.08,colnames(xy[[2]][,8:l])[example_index])
  
}

dev.off()


#Fig1 Power fitting######################

tiff("E:/CAF/2022-10-20/Bacteria/fig1/fig1.tiff",width=20,height=13,units = "cm",res=300)
par(mfrow=c(5,6))
par(mar = c(0, 0, 0,0),oma = c(1, 1, 1, 1),mgp=c(1,0.5,0))
col1 <- c("lightcoral","steelblue1","lightcoral","steelblue1","lightcoral","steelblue1")
col2 <- c("red","blue","red","blue","red","blue")
#col2 <- c("red","blue","red3","blue3","red4","blue4")


examples <- c(6,45,72,226,352)
index <- 1
for (i in 1:5) {
  example_index <- examples[index]
  initial_power_par <- c(0.5,0.5)
  for (t in 1:6) {
    power_par <- get_power_par(initial_power_par,xy[[t]][,8:l][,example_index],xy[[t]]$X)
    plot(xy[[t]]$X,xy[[t]][,8:l][,example_index],pch=16,cex=1,col=col1[t],
         type = "p",ylab = "",xlab = "",axes=FALSE,ylim = c(-0.05,1.05),
         xlim = c(range(xy[[t]]$X)[1]-0.005,range(xy[[t]]$X)[2]+0.005))
    lines(xy[[t]]$X,power.equation(power_par,xy[[t]]$X),col=col2[t],lwd=2)
    box(lwd=1)
    if(t==1&&i!=5){
      axis(2,at = seq(0,1,0.2),labels =FALSE,cex.axis=1)
    }else if(t==1&&i==5){
      axis(2,at = seq(0,1,0.2),labels =FALSE,cex.axis=1)
      cat("################",round(seq(range(xy[[t]]$X)[1],range(xy[[t]]$X)[2],length=4),2),"\n")
      axis(1,at=seq(range(xy[[t]]$X)[1],range(xy[[t]]$X)[2],length=4),labels =FALSE,cex.axis=1)
    }else if(t!=1&&i==5){
      cat("################",round(seq(range(xy[[t]]$X)[1],range(xy[[t]]$X)[2],length=4),2),"\n")
      axis(1,at=seq(range(xy[[t]]$X)[1],range(xy[[t]]$X)[2],length=4),labels =FALSE,cex.axis=1)
    }
  }
  index <- index+1
}

dev.off()




tiff("E:/CAF/2022-10-20/Bacteria/fig1/Fig1-Residuals.tiff",width=20,height=15,units = "cm",res=300)
par(mfrow=c(5,6))
par(mar = c(0, 0, 0,0),oma = c(2.5, 2.5, 0.5, 0.5),mgp=c(1,0.5,0),mai = c(0.2,0.2,0,0))

examples <- c(6,45,72,226,352)
index <- 1
for (i in 1:5) {
  example_index <- examples[index]
  initial_power_par <- c(0.5,0.5)
  for (t in 1:6) {
    power_par <- get_power_par(initial_power_par,xy[[t]][,8:l][,example_index],xy[[t]]$X)
    fit <- power.equation(power_par,xy[[t]]$X)
    raw <- xy[[t]][,8:l][,example_index]
    plot(fit,fit-raw,pch=1,cex=1.5,col=col1[t],xlab="",ylab="",cex.lab=1.5,
         axes=FALSE,xlim=c(range(fit)[1]*0.95,range(fit)[2]*1.05),
         ylim=c(-max(abs(range(fit)))-0.05,max(abs(range(fit)))+0.05))
    abline(h=0,lwd=2,lty=2)
    axis(1,at=round(seq(range(fit)[1],range(fit)[2],length=4),2),
         labels =round(seq(range(fit)[1],range(fit)[2],length=4),2),cex.axis=0.8,gap.axis = 0.2)
    axis(2,at = seq(-max(abs(range(fit))),max(abs(range(fit))),length=4),
         labels =round(seq(-max(abs(range(fit))),max(abs(range(fit))),length=4),2),cex.axis=0.8,gap.axis = 0.2)
    box(lwd=1)
  }
  index <- index+1
}
dev.off()

#Fig2 FunClu######################

k <- 9
# filename<-paste("E:/CAF/2022-10-1/Bacteria/cluster-result/k",k,".RData",sep = "")
# load(filename)
load("E:/CAF/2022-10-20/Bacteria/data2-clu/k9.RData")
module<-return_object$clustered_data[,181]
table(module)
par<-return_object$par[-c(1:33)]

tiff("E:/Doctoral/CAF/2022-10-20/Bacteria/fig2/Fig2-k9.tiff",width=20,height=20,units = "cm",res=300)
par(mfrow=c(k,6))
par(mar = c(0, 0, 0,0),oma = c(1, 1, 1, 1),mgp=c(1,0.5,0))
col1 <- c("pink1","skyblue1","pink1","skyblue1","pink1","skyblue1")
col2 <- c("red","blue","red","blue","red","blue")
#col2 <- c("red","blue","red3","blue3","red4","blue4")


for (i in 1:k) {
  index<-which(module==i)
  
  for (t in 1:6) {
    initial_power_par <- c(0.5,0.5)
    a1 <- get_power_par(initial_power_par,rowMeans(xy[[t]][,8:l][,index]),xy[[t]]$X)
    curve <- power.equation(a1,xy[[t]]$X)
    plot(xy[[t]]$X,curve,type = "n",ylab = "",xlab = "",axes=FALSE,ylim = c(0.05,0.8),
         xlim = c(range(xy[[t]]$X)[1]-0.005,range(xy[[t]]$X)[2]+0.005))
    
    for (ii in index){
      parr <- get_power_par(initial_power_par,xy[[t]][,8:l][,ii],xy[[t]]$X)
      lines(xy[[t]]$X,power.equation(parr,xy[[t]]$X),col=col1[t])
      #points(xy[[t]]$X,xy[[t]][,8:l][,ii],col=col1[t],pch=20)
    }
    lines(xy[[t]]$X,curve,col=col2[t],lwd=2)
    box(lwd=1)
    if(t==1&&i!=k){
      axis(2,at = seq(0,1,0.2),labels =FALSE,cex.axis=1)
    }else if(t==1&&i==k){
      axis(2,at = seq(0,1,0.2),labels =FALSE,cex.axis=1)
      cat("################",round(seq(range(xy[[t]]$X)[1],range(xy[[t]]$X)[2],length=4),2),"\n")
      axis(1,at=seq(range(xy[[t]]$X)[1],range(xy[[t]]$X)[2],length=4),labels =FALSE,cex.axis=1)
    }else if(t!=1&&i==k){
      cat("################",round(seq(range(xy[[t]]$X)[1],range(xy[[t]]$X)[2],length=4),2),"\n")
      axis(1,at=seq(range(xy[[t]]$X)[1],range(xy[[t]]$X)[2],length=4),labels =FALSE,cex.axis=1)
    }
  }
}

dev.off()

#Fig3 Bar statistical graph######################
library(RColorBrewer)
#"CKN","CKR","T30N","T30R","T45N","T45R"
k=9

i="T45R"
plotfilename <- paste("E:/CAF/2022-10-20/Bacteria/fig4/",i,".tiff",sep = "")
df <- c()
for (j in c("May","June","July","August","September")) {
  filename <- paste("./fig3/",i,"/links_gene_",j,".csv",sep = "")
  links <- read.csv(filename)[,1:2]
  if(i=="T30R"&&j=="May"){
    df_j <- data.frame(val.Var1=paste("M",1:k,sep = ""),val.Freq=0.00001,direction=rep("incoimg"),months=rep(j))
  }else if(i=="T30R"&&j=="August"){
    df_j <- data.frame(val.Var1=paste("M",1:k,sep = ""),val.Freq=0.00001,direction=rep("incoimg"),months=rep(j))
  }else{
    df_j <- rbind(data.frame(val=table(links$source),direction=rep("outgoing"),months=rep(j)),
                  data.frame(val=-table(links$target),direction=rep("incoimg"),months=rep(j)))
    if(length(setdiff(paste("M",1:k,sep = ""),names(table(links$target))))>0){
      df_j <- rbind(df_j,data.frame(val.Var1=setdiff(paste("M",1:k,sep = ""),names(table(links$target))),val.Freq=0.00001,direction=rep("incoimg"),months=rep(j)))
    }
  }
  df <- rbind(df,df_j)
}

df$class <- paste(df$direction,"-",df$months,sep = "")
df$class <- factor(df$class,levels = c(paste("incoimg","-",c("May","June","July","August","September"),sep = ""),
                                       paste("outgoing","-",c("May","June","July","August","September"),sep = "")))
# df$direction <- factor(df$direction)
df$months <- factor(df$months,levels = rev(c("May","June","July","August","September")))
df$val.Var1 <- factor(df$val.Var1,levels = paste("M",k:1,sep = ""))

p <- ggplot(df)+
  geom_bar(aes(x=val.Var1,y=val.Freq,group=months,fill=class),
           stat = 'identity',width = 0.8,position = position_dodge(width = 0.8))+
  scale_fill_manual(values = c(brewer.pal(9, 'Blues')[4:8],brewer.pal(9, 'Reds')[3:7]))+
  coord_flip()+
  scale_y_continuous(breaks = seq(-10, 10, 2), labels = as.character(abs(seq(-10, 10, 2))), limits = range(df[,2])) +
  guides(fill="none",alpha="none")+
  theme(panel.grid = element_blank(),
        panel.background = element_rect(color = 'black', fill = 'transparent'),
        legend.title = element_blank(),
        plot.margin=unit(c(0.2,0.2,0.2,0.5), 'lines'),
        axis.title = element_text(color='black',size=15),
        axis.title.x = element_text(vjust=1),
        axis.title.y = element_text(vjust=2.5),
        axis.text.x=element_text(vjust=0.5,size=10,color='black'),
        axis.text.y=element_text(hjust=0.5,size=10,color='black')) +
  xlab("") + ylab('')+
  geom_hline(yintercept = 0, size = 0.6,col="gray20")
tiff(plotfilename,width=12,height=8,units = "cm",res=300)
p
dev.off()



#Fig4 Effect curves######################

library(RColorBrewer)
#ncolor <- brewer.pal(7,"Accent")
ncolor <- brewer.pal(9, 'Blues')[2:7]

k <- 9
load("E:/CAF/2022-10-20/Bacteria/data2-clu/k9.RData")
module<-return_object$clustered_data[,181]

#Calculate the sum of OTUs for each category
data_clustered<-rep()
for (i in 1:k)
{
  new_i<-apply(all.data[which(module==i),],2,sum)
  data_clustered<-rbind(data_clustered,new_i)
}
data_clustered<-t(data_clustered)

tiff("E:/Doctoral/CAF/2022-10-20/Bacteria/fig5/Fig5.tiff",width=20,height=20,units = "cm",res=300)
par(mfcol=c(k,6))
par(mar = c(0, 0, 0,0),oma = c(1, 1, 1, 1),mgp=c(1,0.5,0))

for (j in 1:6) {
  category <- j
  x_int <- xy[[category]]$X
  step<-0.001
  x_smooth <- step_fct(x_int,step)
  nt <- length(x_smooth)
  Months <- xy[[j]]$Month
  name_mon <- names(table(Months))
  
  data_whole <- data_clustered[(1+30*(category-1)):(30*category),]
  colnames(data_whole) <- paste("M",1:k,sep = "")
  
  data_smooth <- rep()
  for (ii in 1:k) {
    yy <- f(data_whole[,ii])
    l_par_i <- lop_par_all[[j]][ii,]
    data_i2 <- legendre_fit(l_par_i,x_smooth)
    data_smooth <- cbind(data_smooth,data_i2)
  }
  colnames(data_smooth) <- paste("M",1:k,sep = "")
  
  relationship_all <- pblapply(1:ncol(data_smooth),function(c)
    get_interaction(data_smooth, c, reduction = FALSE))
  par_dat <- paste("E:/Doctoral/CAF/2022-10-20/Bacteria/fig3/a=0.001-lop_par",j,".RData",sep = "")
  load(par_dat)
  par_all_obj <- list(lop_par = lop_par,
                      relationship = relationship_all,
                      dataset = data_smooth,
                      times = x_smooth,
                      order = 3)
  module_net1 <- get_all_net(par_all_obj)
  
  for (i in 1:k) {
    plot_data <- get_curve_data(par_all_obj, module_net1, i)
    plot(plot_data$x,plot_data$sum,col="#0dceda",type = "n",lwd=2,ylab = "",xlab = "",axes=FALSE,
         ylim = c(-0.3,1),xlim = c(range(plot_data$x)[1]-0.005,range(plot_data$x)[2]+0.03))
    
    for (m in 1:length(name_mon)) {
      points(x=(x_int)[which(Months==name_mon[m])],
             y=f(data_whole[,i])[which(Months==name_mon[m])],
             pch=19,cex=1.2,col=ncolor[m])
    }
    lines(plot_data$x,plot_data$sum,col="#08306B",lwd=2)
    lines(plot_data$x,plot_data$ind_effect,col="#ff7171",lwd=2)
    for (d in 1:(ncol(plot_data)-4)) {
      lines(plot_data$x,plot_data[,1+d],col="#9ede73",lwd=2)
      #text(max(as.numeric(plot_data$x))+0.2+1*(j-1),plot_data[nrow(plot_data),(d+1)]+1.25*(i-1),colnames(plot_data)[d+1])
    }
    dep_name <- colnames(plot_data)[2:(ncol(plot_data)-3)]
    if(length(dep_name)==1){
      text(x=max(as.numeric(plot_data$x))+0.02, 
           y=plot_data[nt,2:(ncol(plot_data)-3)]+0.1, 
           labels=dep_name,cex=1.3)
    }else{
      text(x=rep(max(as.numeric(plot_data$x))+0.02, length(dep_name)), 
           y=seq(min(plot_data[,2:(ncol(plot_data)-3)])+0.05,max(0.2,max(plot_data[,2:(ncol(plot_data)-3)])),length=length(dep_name)), 
           labels=dep_name,cex=1.3)
    }
    
    box(lwd=1)
    abline(h=0,lty=2)
    
    if(j==1&&i==k){
      cat(round(seq(range(plot_data$x)[1],range(plot_data$x)[2],length=4),2),"\n")
      axis(2,at = seq(-0.2,1,0.3),labels =FALSE,cex.axis=1)
      axis(1,at=seq(range(plot_data$x)[1],range(plot_data$x)[2],length=4),labels =FALSE,cex.axis=1)
    }else if(j==1&&i!=k){
      axis(2,at = seq(-0.2,1,0.3),labels =FALSE,cex.axis=1)
    }else if(j!=1&&i==k){
      cat(round(seq(range(plot_data$x)[1],range(plot_data$x)[2],length=4),2),"\n")
      axis(1,at=seq(range(plot_data$x)[1],range(plot_data$x)[2],length=4),labels =FALSE,cex.axis=1)
    }
  }
}

dev.off()


