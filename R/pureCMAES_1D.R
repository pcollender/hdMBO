pureCMAES1D = 
  function (par, fun, lower = NULL, upper = NULL, sigma = 0.5, 
            stopfitness = -Inf, stopeval = 1000 * length(par)^2, ...) 
  {
    stopifnot(is.numeric(par))
    if (sigma < 0.1 || sigma > 0.9) 
      sigma <- max(min(sigma, 0.9), 0.1)
    N <- length(par)
    if (length(lower) == 1) 
      lower <- rep(lower, N)
    if (length(upper) == 1) 
      upper <- rep(upper, N)
    fct <- match.fun(fun)
    fun <- function(x) fct(x, ...)
    xmean <- par
    sigma <- sigma * (upper - lower)
    lambda <- 4 + floor(3 * log(N))
    mu <- lambda/2
    weights <- log(mu + 1/2) - log(1:mu)
    mu <- floor(mu)
    weights <- weights/sum(weights)
    mueff <- sum(weights)^2/sum(weights^2)
    cc <- (4 + mueff/N)/(N + 4 + 2 * mueff/N)
    cs <- (mueff + 2)/(N + mueff + 5)
    c1 <- 2/((N + 1.3)^2 + mueff)
    cmu <- min(1 - c1, 2 * (mueff - 2 + 1/mueff)/((N + 2)^2 + 
                                                    mueff))
    damps <- 1 + 2 * max(0, sqrt((mueff - 1)/(N + 1)) - 1) + 
      cs
    pc <- ps <- numeric(N)
    B <- diag(N, nrow = length(N))
    D <- rep(1, N)
    C <- as.matrix(B) %*% diag(D^2, nrow = length(D)) %*% B
    invsqrtC <- as.matrix(B) %*% diag(D^-1, nrow = length(D)) %*% B
    eigeneval <- 0
    chiN <- N^0.5 * (1 - 1/(4 * N) + 1/(21 * N^2))
    ml.triu <- function(M, k = 0) {
      if (k == 0) 
        M[lower.tri(M, diag = FALSE)] <- 0
      else M[col(M) <= row(M) + k - 1] <- 0
      return(M)
    }
    counteval <- 0
    while (counteval < stopeval) {
      arx <- matrix(0, nrow = N, ncol = lambda)
      arfitness <- numeric(lambda)
      for (k in 1:lambda) {
        arxk <- xmean + sigma * B %*% (D * rnorm(N))
        arxk <- ifelse(arxk > lower, ifelse(arxk < upper, 
                                            arxk, upper), lower)
        arx[, k] <- arxk
        arfitness[k] <- fun(arx[, k])
        counteval <- counteval + 1
      }
      arindex <- order(arfitness, decreasing = FALSE)
      arfitness <- arfitness[arindex]
      xold <- xmean
      xmean <- arx[, arindex[1:mu]] %*% weights
      ps <- (1 - cs) * ps + sqrt(cs * (2 - cs) * mueff) * 
        invsqrtC %*% (xmean - xold)/sigma
      hsig <- norm(ps, "F")/sqrt(1 - (1 - cs)^(2 * counteval/lambda))/chiN < 
        1.4 + 2/(N + 1)
      pc <- (1 - cc) * pc + hsig * sqrt(cc * (2 - cc) * mueff) * 
        (xmean - xold)/sigma
      artmp <- (1/sigma) * (arx[, arindex[1:mu]] - matrix(1, 
                                                          1, mu) %x% xold)
      C <- (1 - c1 - cmu) * C + c1 * (pc %*% t(pc) + (1 - 
                                                        hsig) * cc * (2 - cc) * C) + cmu * artmp %*% diag(weights,
                                                                                                          nrow=length(weights)) %*% 
        t(artmp)
      sigma <- sigma * exp((cs/damps) * (norm(ps, "F")/chiN - 
                                           1))
      if (counteval - eigeneval > lambda/(c1 + cmu)/N/10) {
        eigeneval <- counteval
        C <- ml.triu(C) + t(ml.triu(C, 1))
        if (any(is.nan(C))) 
          break
        EigVal <- eigen(C, symmetric = TRUE)
        B <- EigVal$vectors
        D <- sqrt(EigVal$values)
        invsqrtC <- B %*% diag(D^-1,nrow=length(D)) %*% t(B)
      }
      if (arfitness[1] <= stopfitness || max(D) > 1e+07 * 
          min(D)) 
        break
    }
    xmin <- arx[, arindex[1]]
    return(list(xmin = xmin, fmin = fun(xmin)))
  }
