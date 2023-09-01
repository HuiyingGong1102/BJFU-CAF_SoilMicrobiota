get_SAD1_varmatrix <- function(par,times){
  phi <- par[1]; gamma <- par[2];
  n <- length(times)
  sigma <- array(dim=c(n,n))
  #formula 1, diag element
  diag(sigma) <- sapply(1:n, function(c)(1-phi^(2*c))/(1-phi^2) )
  #formula 2, non-diag element
  sigma[lower.tri(sigma)] <- do.call(c,lapply(1:(n-1),function(c)phi^(((c+1):n)-c)*diag(sigma)[c]))
  sigma[upper.tri(sigma)] <- t(sigma)[upper.tri(t(sigma))]
  sigma <- gamma^2*sigma
  # sigma[1,] <- 1e-10
  # sigma[,1] <- 1e-10
  return(sigma)
}

get_SAD1_covmatrix <- function(par,times){
  phi1 <- par[1]; gamma1 <- par[2]; phi2 <- par[3]; gamma2 <- par[4]; rho <- par[5]
  n <- length(times)
  sigma <- array(dim=c(n,n))
  #formula 1, diag element
  diag(sigma) <- sapply(1:n, function(c)(1-(phi1*phi2)^c)/(1-phi1*phi2) )
  #formula 2, non-diag element
  sigma[lower.tri(sigma)] <- do.call(c,lapply(1:(n-1),function(c)
    (phi1^(((1+c):n)-c)-(phi2^c)*phi1^((1+c):n))/(1-phi1*phi2)))
  sigma[upper.tri(sigma)] <- do.call(c,lapply(1:(n-1),function(c)
    (phi2^((c+1)-(1:c))-(phi2^(c+1))*phi1^(1:c))/(1-phi1*phi2)))
  sigma <- gamma1*gamma2*rho*sigma
  # sigma[1,] <- 1e-10
  # sigma[,1] <- 1e-10
  return(sigma)
}


get_SAD1_matrix <- function(par,n=35,all.time,traits=6){
  rho.par <- par[-(1:(2*traits))]
  sigma <- array(dim=c(traits*n,traits*n))
  r <- 1
  for (i in 1:traits) {
    times <- all.time[(1+(i-1)*n):(i*n)]
    sigma[(1+(i-1)*n):(i*n),(1+(i-1)*n):(i*n)] <- get_SAD1_varmatrix(par[(1+(i-1)*2):(2*i)],times)
    if(i<6){
      for (j in (i+1):traits) {
        times <- all.time[(1+(j-1)*n):(j*n)]
        sigma[(1+(i-1)*n):(i*n),(1+(j-1)*n):(j*n)] <- get_SAD1_covmatrix(c(par[(1+(i-1)*2):(2*i)],par[(1+(j-1)*2):(2*j)],rho.par[r]),times)
        r <- r+1
      }
    }
  }
  sigma[lower.tri(sigma)] <- t(sigma)[lower.tri(t(sigma))]
  return(sigma)
}


power.equation <- function(par,x){
  par[1]*x^par[2]
}

derivative.power.equation <- function(par,x){
  par[1]*par[2]*x^(par[2]-1)
}

s.mle<-function(par,data,x){
  y <- power.equation(par,x)
  yi <- data
  res <- sum((yi-y)^2)
  return(res)
}

get_power_par<-function(initial_power_par,data,x){
  
  a <- optim(initial_power_par,s.mle,data=data,x=x,method = "Nelder-Mead")
  cat("value=",a$value, "\n")
  return(a$par)
}

categorys_get_power_par<-function(initial_power_par,tmp,x){
  if(is.null(dim(tmp))){
    pp <- sapply(1:6, function(c)get_power_par(initial_power_par,tmp[(1+30*(c-1)):(30*c)],x[(1+30*(c-1)):(30*c)])) 
  }else{
    pp <- sapply(1:6, function(c)get_power_par(initial_power_par,colMeans(tmp)[(1+30*(c-1)):(30*c)],x[(1+30*(c-1)):(30*c)])) 
  }
  return(as.vector(pp))
}

categorys_fit_values <- function(par,x){
  c(sapply(1:6, function(c)power.equation(par[(1+2*(c-1)):(2*c)],x[(1+30*(c-1)):(30*c)])))
}

mle <- function(par,data,prob){
  par1<-par[1:27]
  par2<-par[-c(1:27)]
  covM<-get_SAD1_matrix(par1,n=30,all.time=x,traits=6)
  mu <-  sapply(1:k, function(c) categorys_fit_values(par2[(1+12*(c-1)):(12*c)],x))
  temp_S <- sapply(1:k,function(c) dmvnorm(data,
                                           mu[,c],
                                           covM)*prob[c] )
  LL <- sum(-log(rowSums(temp_S)))
  return(LL)
}


get_cluster <- function(data,k){
  #set.seed(1324)
  init_cluster <- kmeans(data,k,iter.max = 100)
  init_curve_par<-c()
  init_sd_par<-c(0.99698404,0.04259229,1.02503600,0.05190239,0.99909588,0.05613654,1.03011437,0.05428003,
                 0.97937515,0.05533830,0.99749095,0.05061320,0.46967808,0.55160289,0.46488999,0.52887021,
                 0.51731941,0.62517417,0.64365041,0.61037882,0.60995562,0.64869961,0.72432177,0.73948490,
                 0.68079467,0.67100283,0.76733228)
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
    
    par1<-init_para[1:27]
    par2<-init_para[-c(1:27)]
    covM<-get_SAD1_matrix(par1,n=30,all.time=x,traits=6)
    mu <-  sapply(1:k, function(c) categorys_fit_values(par2[(1+12*(c-1)):(12*c)],x))
    mvn.c <- sapply(1:k, function(c) dmvnorm(data,
                                             mu[,c],
                                             covM)*pro[c] )
    omega <- mvn.c/rowSums(mvn.c)
    #M step, calculate parameters
    pro <- colSums(omega)/sum(omega)
    
    mle.covar <- function(npar){
      
      nnpar <- c(npar,init_para[-c(1:27)])
      AA <- mle(nnpar,data,pro)
      AA
    }
    r1.covar <- optim(init_para[1:27],mle.covar,method = "BFGS",control=list(maxit=32000))
    init_para[1:27] <- r1.covar$par
    
    par1<-init_para[1:27]
    par2<-init_para[-c(1:27)]
    covM<-get_SAD1_matrix(par1,n=30,all.time=x,traits=6)
    mu <-  sapply(1:k, function(c) categorys_fit_values(par2[(1+12*(c-1)):(12*c)],x))
    mvn.c <- sapply(1:k, function(c) dmvnorm(data,
                                             mu[,c],
                                             covM)*pro[c] )
    omega <- mvn.c/rowSums(mvn.c)
    pro <- colSums(omega)/sum(omega)
    
    mle.mean <- function(npar){
      
      nnpar <- c(init_para[1:27],npar)
      AA <- mle(nnpar,data,pro)
      AA
    }
    r1.curve <- optim(init_para[-c(1:27)],mle.mean,method = "BFGS",control=list(maxit=32000))
    
    L_Value <- r1.curve$value
    
    init_para[-c(1:27)] <- r1.curve$par
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
all_par <- sapply(1:207, function(c)categorys_get_power_par(initial_power_par,all.data[c,],x))
data2 <-  t(sapply(1:207, function(c) categorys_fit_values(all_par[,c],x)))

for (k in c(2:15)) {
  return_object <- get_cluster(data=data2,k=k)
  table(return_object$clustered_data[,181])
  
  filename<-paste("E:/CAF/2022-10-20/fungi/cluster-result/k",k,".RData",sep = "")
  save(return_object,file=filename)
}

norder<-c(2:15)
bic<-c()
for (i in norder) {
  filename<-paste("E:/CAF/2022-10-20/fungi/cluster-result/k",i,".RData",sep = "")
  load(filename)
  BIC <- return_object$BIC
  L_Value<-return_object$LL
  BIC <- 2*L_Value+log(nrow(data2)*30)*6*length(return_object$par)
  
  bic<-c(bic,BIC)
}
norder[which(bic==min(bic))]

plot(norder,bic,type = "l",xlab="modules",ylab = "BIC",yaxt = "n",
     lwd=2,col="red",lty=1,mgp=c(2.5,1,0),cex.lab=1.5)




