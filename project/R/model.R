dlm.model <- function(y, alpha0=2, beta0=20, order=4, discount.factor=0.9) {
    mod <- dlmModPoly(order, dV = 1)
    modFilt <- dlmFilterDF(y, mod, DF=discount.factor)

    ## Filtering estimates
    out <- residuals(modFilt)
    beta <- beta0 + cumsum(out$res^2) / 2
    alpha <- alpha0 + (1 : length(y)) / 2
    Ctilde <- unlist(dlmSvd2var(modFilt$U.C, modFilt$D.C))[-1]
    prob <- 0.95
    tt <- qt(prob, df = 2 * alpha)
#     lower <-  dropFirst(modFilt$m) - tt * sqrt(Ctilde * beta / alpha)
#     upper <-  dropFirst(modFilt$m) + tt * sqrt(Ctilde * beta / alpha)
#     plot(y, ylab = "Filtering level estimates", type = "o", 
#             ylim = c(18, 40), col = "darkgray")
#     lines(dropFirst(modFilt$m), type = "o")
#     lines(lower, lty = 2, lwd = 2)
#     lines(upper, lty = 2, lwd = 2)
    ## One-step-ahead forecasts
    sigma2 <- c(beta0 / (alpha0-1), beta / (alpha-1))
    Qt <- out$sd^2 * sigma2[-length(sigma2)]
    alpha0T = c(alpha0,alpha)
    tt <- qt(prob, df = 2 * alpha0T[-length(alpha0T)])
    parf <- c(beta0 / alpha0, beta / alpha)
    parf <- parf[-length(parf)] * out$sd^2
    lower <- dropFirst(modFilt$f) - tt * sqrt(parf)
    upper <- dropFirst(modFilt$f) + tt * sqrt(parf)
#     plot(y, ylab = "One-step-ahead forecasts", type = "o",
#          ylim = c(20, 40), col = "darkgray")
#     lines(window(modFilt$f, start=1902),  type="o") 
#     lines(lower, lty = 2, lwd = 2)
#     lines(upper, lty = 2, lwd = 2)
    ## Smoothing estimates
#     modFilt$mod$JW <- matrix(1)
#     X <- unlist(dlmSvd2var(modFilt$U.W, modFilt$D.W))[-1]
#     modFilt$mod$X <- matrix(X)
#     modSmooth <- dlmSmooth(modFilt)
#     Stildelist <- dlmSvd2var(modSmooth$U.S, modSmooth$D.S)
#     TT <- length(y)
#     pars <- unlist(Stildelist) * (beta[TT] / alpha[TT])
#     tt <- qt(prob, df = 2 * alpha[TT])
#     plot(y, ylab = "Smoothing level estimates", 
#          type = "o", ylim = c(20, 40), col = "darkgray")
#     lines(dropFirst(modSmooth$s),  type = "o")
#     lines(dropFirst(modSmooth$s - tt * sqrt(pars)), lty = 3, lwd = 2)
#     lines(dropFirst(modSmooth$s + tt * sqrt(pars)), lty = 3, lwd = 2)

    ## Choosing the discount factor
#     vecDF <- c(1,0.9,0.8,0.3)
#     mat <- matrix(0, nrow=length(vecDF), ncol=5,dimnames = 
#         list(c("", "","",""), c("DF", "MAPE", "MAD","MSE","SD")))
#     par(mar = c(2, 4, 1, 1) + 0.1, cex = 0.6)
#     par(mfrow=c(2,2)) 
#     for (i in 1:length(vecDF)){
#        mod <- dlmModPoly(1,dV=1)
#        modFilt <- dlmFilterDF(y, mod, DF=vecDF[i])
#        out <- residuals(modFilt)
#        beta0 <- 20
#        alpha0 <-2
#        beta <- beta0 + 1/2 * cumsum(out$res^2)
#        alpha <- alpha0 + (1: length(y))*(1/2)
#        sigmaqhat <- beta/(alpha-1) 
#        sigmaqhat <- beta/(alpha-1) 
#        # Plot the estimates of sigma^2
#        # ts.plot(sigmaqhat, main=paste("DF=", as.character(vecDF[i])), ylab="")  
#        # Plot the one-step-ahead forecasts   
#        plot(window(y, start=1960), ylab="", ylim=c(20,45), lty=1, 
#             main=paste("DF=", as.character(vecDF[i])), col="darkgray") 
#        lines(window(modFilt$f, start=1960), lty=1)  
#        mat[i,1] <- vecDF[i]
#        mat[i,2] <- mean(abs(modFilt$f - y)/y)
#        mat[i,3] <- mean(abs(modFilt$f - y))
#        mat[i,4] <- mean((modFilt$f - y)^2)
#        mat[i,5] <- sigmaqhat[length(y)]
#     }
# 
#     round(mat,4)

    return(cbind(obs=y, ll=lower, ul=upper))
}
