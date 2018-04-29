"2 Stage Least Square function for R:
INPUTS:
 y takes a vector of dependent variable
 X takes a matrix (or data.frame) of instruments
 Z takes a matrix (or data.frame) of regressors
 print_flag takes either 1 (to print estimates and a number of measures) or 0 (not to print them)

OUTPUTS:
 b carries the estimates for the regressors' coefficients
 se carries the standard error
 t carries the t-values
 Avar carries the Asymptotic Variance of the estimates
 df carries the degrees of freedom of the model
 nobs carries the number of observations
 mean_y carries the mean of the dependent variable
 R_sq carries the R Squared measure for the regression
 see carries the standard error of the equation
 ssr  carries the sum of squared residuals
 S_hat carries the estimate for the variance of the error
 Sargan carries the value of the Sargan Statistic for the regression"

twostageLS_f <- function(y, X, Z, flag_print) {
     
     cat(" By invoking this function you are assuming conditional homoskedasticity\n")
     
     # Checking for the value types
     if (!is.matrix(X)) {
          X <- as.matrix(X)
     }
     if (! is.matrix(Z)) {
          Z <- as.matrix(Z)
     }
     
     # Do GMM algebra
     nobs <- nrow(Z)  # number of observations
     K <- ncol(X)
     L <- ncol(Z)
     S_xx <- as.matrix((t(X)%*%X)/nobs)  # sample cross moments
     df <- K-L  # degrees of freedom
     s_xy <- t(X)%*%y/nobs  # x"y/n
     s_yy <- t(y)%*%y/nobs  # y"y/n
     S_xz <- t(X)%*%Z/nobs  # x"z/n
     S_zz <- t(Z)%*%Z/nobs  # z"z/n
     mean_y <- mean(y)     # mean of y
     W <- inv(S_xx)  # setting W to be the S_xx
     
     # 2SLS Consistent Estimation
     delta_hat <- solve(t(S_xz)%*%W%*%S_xz)%*%(t(S_xz)%*%W%*%s_xy)   # the 2SLS estimator
     
     # Calculate Control Statistics
     e <- c(y-Z%*%delta_hat)
     ssr <- (t(e)%*%(e))/nobs  # sum of squared residuals
     var_y <- s_yy-mean_y^2  # variance of y
     R_sq <- 1-(ssr/nobs)/var_y   # R-squared
     see <- sqrt(ssr/df)     # standard error of the equation
     
     # Calculate robust standard errors
     sigma_hat_sqr <- sum(e^2)/nobs
     g <- e*X
     S_hat <- t(g)%*%g/nobs
     W_hat <- solve(S_hat)
     
     Avar_delta_hat <- sigma_hat_sqr * solve(t(S_xz) %*% solve(S_xx) %*% S_xz)
     se <- diag(sqrt(Avar_delta_hat/nobs))  # standard errors
     g_bar <- s_xy - S_xz%*%delta_hat
     temp_inv <- solve(t(S_xz)%*%W_hat%*%S_xz)
     t <- delta_hat/se   # t-value
     # Sargan's Statistic
     Sargan <- nobs*t(s_xy-S_xz%*%delta_hat)%*%inv(S_xx)%*%(s_xy-S_xz%*%delta_hat)/(ssr)
     
     # Print results if desired
     if (flag_print==1){
          cat("******************** GMM **********************\n")
          cat(" Number of Observations:",nobs);cat("\n")
          cat(" Mean of the Dependent Variable:",mean_y);cat("\n")
          cat(" Variance of the Dependent Variable:",var_y);cat("\n")
          cat(" Centered R-squared:",R_sq);cat("\n")
          cat(" Standard Error of the Equation:",see);cat("\n")
          cat(" Sum of Squared Residuals:",ssr);cat("\n")
          cat(" The degrees of freedom:",df);cat("\n")
          cat(" Sargan's Statistic/Chi-sqrt:",Sargan,"/",qchisq(.95,df));cat("\n")
          cat("--------------------------------------------------------------\n")
          cat("variable            estimate       s.e.        t\n")
          cat("--------------------------------------------------------------\n")
          for (i in 1:length(delta_hat)){
               if (is.null(colnames(Z))) {
                    cat(sprintf('%6.0f  %20.6f %12.6f %12.6f\n', i, delta_hat[i], se[i], t[i]))
               }
               else {
                    cat(sprintf('%6s  %20.6f %12.6f %12.6f\n', colnames(Z)[i], delta_hat[i], se[i], t[i]))
               }
          }
          cat("--------------------------------------------------------------\n")
          
          
     }
     return(list(b=delta_hat,se=se,t=t,Avar_delta_hat=Avar_delta_hat, df= df,
                 nobs=nobs,mean_y=mean_y,R_sq=R_sq,see=see,ssr=ssr,S_hat=S_hat, Sargan = Sargan))
     
     
}
