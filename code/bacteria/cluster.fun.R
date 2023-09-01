
get_cluster <- function(data,k){
  #set.seed(1324)
  init_cluster <- kmeans(data,k,iter.max = 100)
  init_curve_par<-c()
  init_sd_par<-c(1.01198516,-0.01671540,1.00107770,-0.02469911,0.97814544,-0.02846965,0.99866566,-0.02355896,
                 1.01502251,-0.01813548,1.00581784,-0.01475097,0.43450061,0.44796524,0.42083692,0.41300408,
                 0.28627527,0.63312904,0.54376242,0.41906569,0.41783621,0.65774494,0.41466895,0.37359392,
                 0.41544035,0.30733342,0.28137631,0.50000000,0.50000000,0.50000000,0.50000000,0.50000000,
                 0.50000000)
  initial_power_par <- c(0.5,0.5)
  init_curve_par <- c(sapply(1:k, function(c)categorys_get_power_par(initial_power_par,data[which(init_cluster$cluster==c),],x)))
  init_pro <- table(init_cluster$cluster)/nrow(data)
  
  Delta <- 1000; iter <- 1; itermax <- 1000;
  pro <- init_pro
  init_para <- c(init_sd_par,init_curve_par)
  LL_mem <- mle(par=init_para,data=data,prob=pro)
  
  while ( Delta > 0.1 && iter <= itermax ) {
    
    #E step, calculate the posterior probability
    old_par <- init_para
    
    par1<-init_para[1:33]
    par2<-init_para[-c(1:33)]
    covM<-get_SAD1_matrix(par1,n=30,all.time=x,traits=6)
    mu <-  sapply(1:k, function(c) categorys_fit_values(par2[(1+12*(c-1)):(12*c)],x))
    mvn.c <- sapply(1:k, function(c) dmvnorm(data,
                                             mu[,c],
                                             covM)*pro[c] )
    omega <- mvn.c/rowSums(mvn.c)
    #M step, calculate parameters
    pro <- colSums(omega)/sum(omega)
    
    mle.covar <- function(npar){
      
      nnpar <- c(npar,init_para[-c(1:33)])
      AA <- mle(nnpar,data,pro)
      AA
    }
    r1.covar <- optim(init_para[1:33],mle.covar,method = "BFGS",control=list(maxit=32000))
    init_para[1:33] <- r1.covar$par
    
    par1<-init_para[1:33]
    par2<-init_para[-c(1:33)]
    covM<-get_SAD1_matrix(par1,n=30,all.time=x,traits=6)
    mu <-  sapply(1:k, function(c) categorys_fit_values(par2[(1+12*(c-1)):(12*c)],x))
    mvn.c <- sapply(1:k, function(c) dmvnorm(data,
                                             mu[,c],
                                             covM)*pro[c] )
    omega <- mvn.c/rowSums(mvn.c)
    pro <- colSums(omega)/sum(omega)
    
    mle.mean <- function(npar){
      
      nnpar <- c(init_para[1:33],npar)
      AA <- mle(nnpar,data,pro)
      AA
    }
    r1.curve <- optim(init_para[-c(1:33)],mle.mean,method = "BFGS",control=list(maxit=32000))
    
    L_Value <- r1.curve$value
    
    init_para[-c(1:33)] <- r1.curve$par
    Delta <- abs(L_Value-LL_mem)
    #Delta <- max(abs( old_par - new_par$par) );
    cat('\n',"iter=",iter,"LL=",L_Value,'\n')
    iter <- iter+1; LL_mem <- L_Value
  } 
  BIC <- 2*L_Value+log(nrow(data))*length(init_para)
  cluster <- apply(omega,1,which.max)
  clustered_df <- cbind(data,cluster)
  init_para <- init_para
  
  return_object <- list(init_para,pro,LL_mem,BIC,clustered_df)
  names(return_object)<-c("par", "pro", "LL", "BIC","clustered_data")
  
  return(return_object)
}

initial_power_par <- c(0.5,0.5)
all_par <- sapply(1:356, function(c)categorys_get_power_par(initial_power_par,all.data[c,],x))
data2 <-  t(sapply(1:356, function(c) categorys_fit_values(all_par[,c],x)))

for (k in c(2:15)) {
  return_object <- get_cluster(data=data2,k=k)
  table(return_object$clustered_data[,181])
  
  filename<-paste("E:/CAF/2022-10-20/Bacteria/data2-clu/k",k,".RData",sep = "")
  save(return_object,file=filename)
}

norder<-c(2:15)
bic<-c()
for (i in norder) {
  filename<-paste("E:/CAF/2022-10-20/Bacteria/data2-clu/k",i,".RData",sep = "")
  load(filename)
  BIC <- return_object$BIC
  L_Value<-return_object$LL
  BIC <- 2*L_Value+log(nrow(data2)*30)*6*(length(return_object$par)-6)
  
  bic<-c(bic,BIC)
}
norder[which(bic==min(bic))]

plot(norder,bic,type = "l",xlab="modules",ylab = "BIC",yaxt = "n",
     lwd=2,col="red",lty=1,mgp=c(2.5,1,0),cex.lab=1.5)
