"OLS for R:
*** make sure to have all functions included here when running OLS_with_robust_standard_errors_f
INPUTS:
y takes a vector of dependent variable
X takes a matrix (or data.frame) of regressors
print_flag takes either 1 (to print estimates and a number of measures) or 0 (not to print them)

OUTPUTS:
 b carries the estimates for the regressors' coefficients
 se carries the standard error
 t carries the t-values
 v carries the Asymptotic Variance of the estimates
 nobs carries the number of observations
 mean_y carries the mean of the dependent variable
 R_sq carries the R Squared measure for the regression
 see carries the standard error of the equation
 ssr  carries the sum of squared residuals"

OLS_with_robust_standard_errors_f <- function(y,X,flag_print) {
     
     # Do OLS algebra
     nobs <- nrow(X)  # number of observations
     df <- nobs-ncol(X)  # degrees of freedom
     S_xx <- as.matrix(t(X)%*%X/nobs)  # sample cross moments
     s_xy <- t(X)%*%y/nobs  # x"y/n
     s_yy <- t(y)%*%y/nobs  # y"y/n
     mean_y <- mean(y)     # mean of y
     b <- solve(S_xx)%*%s_xy   # the OLS formula in sample averages
     ssr <- nobs*(s_yy-t(s_xy)%*%b)   # sum of squared residuals
     var_y <- s_yy-mean_y^2  # variance of y
     R_sq <- 1-(ssr/nobs)/var_y   # R-squared
     see <- sqrt(ssr/df)     # standard error of the equation
     
     # Calculate standard errors
     e <- y-X%*%b
     e <- c(e)   # convert e to a vector
     temp <- e*X
     S_hat <- t(temp)%*%temp/nobs
     S_xx_inv <- solve(S_xx)
     v <- S_xx_inv%*%S_hat%*%S_xx_inv/nobs
     se <- sqrt(diag(v))  # standard errors
     t <- b/se   # t-value
     
     # Print results if desired
     if (flag_print==1){
          cat("************* OLS with Robust Standard Errors ***************\n")
          cat(" Number of Observations:",nobs);cat("\n")
          cat(" Mean of the Dependent Variable:",mean_y);cat("\n")
          cat(" Variance of the Dependent Variable:",var_y);cat("\n")
          cat(" Centered R-squared:",R_sq);cat("\n")
          cat(" Standard Error of the Equation:",see);cat("\n")
          cat(" Sum of Squared Residuals:",ssr);cat("\n")
          cat("--------------------------------------------------------------\n")
          cat("variable            estimate       s.e.        t\n")
          cat("--------------------------------------------------------------\n")
          for (i in 1:length(b)){
               if (is.null(colnames(X))) {
                    cat(sprintf('%6.0f  %20.6f %12.6f %12.6f\n', i, b[i], se[i], t[i]))
               }
               else {
                    cat(sprintf('%6s  %20.6f %12.6f %12.6f\n', colnames(X)[i], b[i], se[i], t[i]))
               }
          }
          cat("--------------------------------------------------------------\n")
     }
     return(list
            (b=b,se=se,t=t,v=v,nobs=nobs,mean_y=mean_y,R_sq=R_sq,see=see,ssr=ssr))
}

sample_j_order_autocovariance <- function(sample,j,Mean=NULL) {
     z <- sample
     l_z <- length(z)
     if (is.null(Mean)) {
          Mean<-mean(z)
     }
     return(sum((z[(j+1):l_z]-Mean)*((r[1:(l_z-j)]-Mean)))/l_z)
}

sample_j_order_autocovariance_coefficient <- function(sample,j,Mean=NULL) {
     return(sample_j_order_autocovariance(sample=sample,j=j,Mean)/sample_j_order_autocovariance(sample=sample,j=0,Mean))
}

annual_inflation_rate_over_t <- function(X) {
     return(c(0, ((X[2:length(X)]/X[1:(length(X)-1)])^12 -1)*100 )
     )
     
}

ljung_box_q_statistic <- function(n,rho_vector) {
     return(n*(n+2)*sum((rho_vector^2)/(n-(1:length(rho_vector)))))
}

standard_error_estimator <- function(x,epsilon) {
     S <- 0
     for (i in 1:length(epsilon)) {
          S <- S + (epsilon[i]^2)*(x[i]^2)
     }
     return(S/i)
}
