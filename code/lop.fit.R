# get_legendre_par <- function(y,legendre_order,x){
#   #lm_method
#   legendre_par <- as.numeric(coef(lm(y~get_legendre_matrix(x,legendre_order))))
#   return(legendre_par)
# }

get_legendre_par2 <- function(y,legendre_order,x){
  #lm_method
  initpar <- rep(0.1,legendre_order)
  smle <- function(y,legendre_order,x){
    if(min(legendre_fit(legendre_order,x))<0){
      return(100)
    }else{
      sum((y-legendre_fit(legendre_order,x))^2)
    }
  }
  legendre_par <- optim(initpar,smle,y=y,x=x,method = "Nelder-Mead")
  return(legendre_par$par)
}

get_legendre_matrix <- function(x,legendre_order){
  legendre_coef <- legendre.polynomials(n = legendre_order, normalized=F)
  legendre_matrix <- as.matrix(as.data.frame(polynomial.values(
    polynomials = legendre_coef, x = scaleX(x, u = -1, v = 1))))
  colnames(legendre_matrix) <- paste0("legendre_",0:legendre_order)
  return(legendre_matrix[,2:(legendre_order+1)])
}

legendre_fit <- function(par,x){
  legendre_order = length(par)
  fit <- sapply(1:length(par), function(c)
    par[c]*legendre.polynomials(n=legendre_order, normalized=F)[[c]])
  legendre_fit <- as.matrix(as.data.frame(polynomial.values(
    polynomials = fit, x = scaleX(x, u = -1, v = 1))))
  x_interpolation <- rowSums(legendre_fit)
  return(x_interpolation)
}

