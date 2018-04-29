gmm_f<- function(y, X ,Z, W=NULL, flag_print) {
     
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
     
     control_W_was_null <- FALSE
     # Check for W
     if (is.null(W)){
          cat(" No W is provided\n");cat(" W = S_xx is used\n")
          W <- inv(S_xx)
          control_W_was_null <- TRUE
     }
     else {
          W <- as.matrix(W)
     }
     
     # GMM Consistent Estimation
     delta_hat <- solve(t(S_xz)%*%W%*%S_xz)%*%(t(S_xz)%*%W%*%s_xy)   # the GMM estimator
     
     
     # Calculate Control Statistics
     e <- c(y-Z%*%delta_hat)
     ssr <- (t(e)%*%(e))/nobs  # sum of squared residuals
     var_y <- s_yy-mean_y^2  # variance of y
     R_sq <- 1-(ssr/nobs)/var_y   # R-squared
     see <- sqrt(ssr/df)     # standard error of the equation
     
     # Calculate robust standard errors
     cat("dim(x) ",dim(X)); cat("\n")
     cat("dim(e) ", dim(e))
     g <- e*X
     S_hat <- t(g)%*%g/nobs
     
     W_hat <- W
     
     temp_inv <- solve(t(S_xz)%*%W_hat%*%S_xz)
     Avar_delta_hat <- temp_inv%*%t(t(S_xz)%*%W_hat%*%S_hat%*%W_hat%*%S_xz)%*%temp_inv
     se <- diag(sqrt(Avar_delta_hat/nobs))  # standard errors
     g_bar <- s_xy - S_xz%*%delta_hat
     t <- delta_hat/se   # t-value
     
     # Hansen's J Statistic
     Hansen_J <- nobs*t(s_xy-S_xz%*%delta_hat)%*%S_hat%*%(s_xy-S_xz%*%delta_hat)
     
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
          cat(" Hansen's J Statistic/Chi-sqrt:",Hansen_J,"/",qchisq(.95,K-L));cat("\n")
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
                 nobs=nobs,mean_y=mean_y,R_sq=R_sq,see=see,ssr=ssr,S_hat=S_hat, Hansen_J = Hansen_J))
     
     
}
