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

    ## One-step-ahead forecasts
    sigma2 <- c(beta0 / (alpha0-1), beta / (alpha-1))
    Qt <- out$sd^2 * sigma2[-length(sigma2)]
    alpha0T = c(alpha0,alpha)
    tt <- qt(prob, df = 2 * alpha0T[-length(alpha0T)])
    parf <- c(beta0 / alpha0, beta / alpha)
    parf <- parf[-length(parf)] * out$sd^2
    lower <- dropFirst(modFilt$f) - tt * sqrt(parf)
    upper <- dropFirst(modFilt$f) + tt * sqrt(parf)

    return(cbind(obs=y, ll=lower, ul=upper))
}
