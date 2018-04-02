
efficient2step_gmm_f <- function(y, X, Z, flag_print) {
     
     # Do Initial GMM Algebra
     nobs <- nrow(Z)  # number of observations
     K <- ncol(X)
     L <- ncol(Z)
     S_xx <- as.matrix((t(X)%*%X)/nobs)  # sample cross moments
     df <- nobs-K  # degrees of freedom
     s_xy <- t(X)%*%y/nobs  # x"y/n
     s_yy <- t(y)%*%y/nobs  # y"y/n
     S_xz <- t(Z)%*%X/nobs  # x"z/n
     S_zz <- t(Z)%*%Z/nobs  # z"z/n
     mean_y <- mean(y)     # mean of y
     
     #========================= Step I ==========================#
     # Initializing W
     W <- S_xx
     # Use GMM estimator to find S_hat based on this W
     delta_hat <- solve(t(S_xz)%*%W%*%S_xz)%*%(t(S_xz)%*%W%*%s_xy)   # the GMM estimator
     e <- y-Z%*%delta_hat
     e <- c(e)   # convert e to a vector
     g <- e*X
     S_hat <- t(g)%*%g/nobs
     
     #========================= Step II ==========================#
     # Set W to be the S_hat found and do GMM
     return(gmm_f(y=y,X=X,Z=Z,W=S_hat,flag_print=flag_print))
     
}