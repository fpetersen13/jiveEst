#' JIVE Estimator
#'
#'
#' This function calculates estimands by using the JIVE estimator
#' as proposed by Angrist, Imbens, Krueger 1999. The output includes the coefficents,
#' variance matrix, standard errors and test statistics.
#' The type of estimator (i.e. JIVE1 vs. JIVE2) and standard errors
#' (i.e. homoskedastic: HC0 and heteroskedastic robust: HCE) can be given as parameters.
#'
#' Version: 0.1
#'
#' @param y: a n x 1 matrix of the dependent variable,
#' @param X: a n x k matrix of independent variables,
#' @param Z: a n x l matrix of instrumental variables with l >= k,
#' @param type: JIVE1 or JIVE2 for the type of estimator that should be used,
#' @param vcov: HC0 for homoskedastic or HCE for heteroskedastic robust standard errors.
#' @return A list of regression results
#'
#' URL: https://github.com/fpetersen13/
#' BugReports: https://github.com/fpetersen13/jiveEst/issues
#'
#'

#' @export
jive_est <- function(y,X,Z,type=JIVE1, vcov = HCE){
  n  <- length(y);
  k  <- dim(X)[2];
  ret <- list();
  if(type==JIVE1){
    jive <- estimate_jive(y,X,Z,type,vcov)
    ret$coe <- jive$coe
    ret$var <- jive$var
    ret$se  <- sqrt(diag(jive$var))
    ret$t   <- ret$coe / ret$se
    ret$p   <- 2*pt(ret$t,n-1)
  }

  return(ret)
}#jive_est


##################################################
### JIVE
estimate_jive <- function(y,X,Z,type,vcov){
  n <- length(y);
  k <- dim(X)[2];
  ret <- list();
  X.hat <- matrix(0,n,k);
  h <- vector("numeric",n);

  if(type==JIVE1){
    for(i in 1:n){
      h   <- solve(t(Z[-i,])%*%Z[-i,])%*%(t(Z[-i,])%*%X[-i,])
      X.hat[i,] <- Z[i,]%*%h
    }#n
    b_jive <- solve(t(X.hat)%*%X)%*%t(X.hat)%*%y
    ret$coe <- b_jive
  }

  if(type==JIVE2){
    for(i in 1:n){
      h   <- solve(t(Z)%*%Z)%*%(t(Z[-i,])%*%X[-i,])*(n/(n-1))
      X.hat[i,] <- Z[i,]%*%h
    }#n
    b_jive <- solve(t(X.hat)%*%X)%*%t(X.hat)%*%y
    ret$coe <- b_jive
  }




  ### Heteroskedasticity Consistent Standard Erros
  if(vcov==HCE){
    e.hat2 <- matrix(0,n,1);
    for(i in 1:n){
      e.hat2[i,] <- (y[i] - X[i,] %*% b_jive)^2
    }#n

    cov.ma <- matrix(0,k,k);
    for(i in 1:n){
      cov.ma <- cov.ma + e.hat2[i,]*X.hat[i,]%*%t(X.hat[i,])
    }#n

    var.hat <- solve(t(X.hat)%*%X)%*%cov.ma%*%solve(t(X)%*%X.hat)
    ret$var  <- var.hat
  }

  ### Homoskedastic Standard Erros
  if(vcov==HC0){
    sigma.hat2 <- sum((y - X%*%b_jive)^2)/n
    var.hat <- sigma.hat2 * solve(t(X.hat)%*%X) %*% t(X.hat)%*%X.hat %*%solve(t(X)%*%X.hat)
    ret$var  <- var.hat
  }
}#JIVE1
