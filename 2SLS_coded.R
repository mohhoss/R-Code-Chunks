
twostageLS_f <- function(y, X, Z, flag_print) {
     
     cat(" By invoking this function you are assuming conditional homoskedasticity\n")
     # I can simply call my coded GMM as in the case W is not provided, my function simply W <- S_xx
     # This is what's required for 2SLS
     return(gmm_f(y=y,X=X,Z=Z,W=NULL,flag_print=flag_print))
}