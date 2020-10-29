#' Compute a design's augmented expected improvement (AEI)
#'
#' @description
#' `getAEI` returns the AEI for a proposed design points based
#' on a provided Gaussian Process model, the previously tested
#' design points, and the effective best solution (EBS) thus far.
#'
#' @details
#' TODO: PROVIDE DETAILS ON AEI COMPUTATION
#'
#' @param point Numeric vector of same dimension as `model` and `designs`
#' @param model Object of class `km` with attr(,"package") = "DiceKriging"
#' @param designs Matrix with columns corresponding to design parameters
#' @param ebs Output of internal function `getEBS`
#' @param se.threshold Setting AEI to zero for nearby tested values
#' @export
getAEI = function(point, model, designs, ebs, se.threshold = 1e-6) {
  # Predict surrogate function at test point
  p = DiceKriging::predict.km(model, newdata =
                                matrix(point,nrow=1,
                                       dimnames = list(
                                         NULL, names(designs))),
                              type = "SK", checkNames = F)
  p.mu = p$mean
  p.se = p$sd

  # compute AEI input values
  d = ebs$mu - p.mu
  xcr = d / p.se
  xcr.prob = pnorm(xcr)
  xcr.dens = dnorm(xcr)

  # get nugget effect
  pure.noise.var = model@covariance@nugget

  tau = sqrt(pure.noise.var)

  # final AEI computation (automatically set to zero if enough
  #   confidence re nearby measurements)
  res =
    ifelse(p.se < se.threshold, 0,
           (d * xcr.prob + p.se * xcr.dens) * (1 - tau / sqrt(tau^2 + p.se^2)))

  return(res)

}




getEBS = function(designs, model,
                  is.noisy, jitter) {

  preds = DiceKriging::predict.km(model, newdata = designs +
                                    ifelse(is.noisy,
                                           sample(c(-1,1),1)*jitter, 0),
                                  type = "SK",checkNames = F)
  mu = preds$mean
  se = preds$sd

  u.x = -1*mu-1*se
  j = which.max(u.x)
  return(list(index = j, des = designs[j, ],
              mu = mu[[j]], se = se[[j]]))
}
