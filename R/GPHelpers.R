### i. Fit GP model via DiceKriging package -----------------------------------
FitGPModel <- function(obj.evals, designs,
                       bfgs.init.pts.per.d, kernel, is.noisy){
  d = ncol(designs)
  # run Gaussian process model for each objective function
  apply(obj.evals, 2,
        function(resp) {
          DiceKriging::km(design       = designs,
                          covtype      = kernel,
                          response     = resp,
                          coef.trend   = 0,
                          nugget.estim = is.noisy,
                          control = list(trace = 0,
                                         pop.size =
                                           max(20, d*bfgs.init.pts.per.d)))
        })
}


### ii. Sample approx. function paths from the GP posterior distribution ------
SampleGPwRandomFeatures <- function(designs, obj.evals,
                                    model, nFeatures = 1000) {
  ##
  d      = model@d
  N_data = model@n

  sigma2 = model@covariance@sd2
  nu2    = model@covariance@nugget
  range.pars  = model@covariance@range.val

  ## Only supporting the squared exponential and Matern52 kernel thus far
  if(model@covariance@name == "matern5_2") {
    m = 5.0/2.0
    # unscaled random normal matrix
    W = matrix(rnorm(d*nFeatures), ncol=d)
    # scale values
    W = t(t(W) / range.pars) / sqrt(rgamma(nFeatures, shape=m, scale=1/m))

  } else if(model@covariance@name == "gauss") {
    # unscaled random normal matrix
    W = matrix(rnorm(d*nFeatures), ncol=d)
    # scale values
    W = t(t(W) / range.pars)
  } else {stop("Provide one of the supported kernels")}
  ## Random draws from 0 to 2pi
  b = runif(n=nFeatures, min=0,max=2*pi)

  ## Standard Normal draws
  rand.norm.nfeat = rnorm(nFeatures)

  # create design matrix
  phi.x = sqrt(2*sigma2/nFeatures) * cos(W %*% t(designs) + replicate(N_data,b))

  chol_Sigma_inv = chol(phi.x %*% t(phi.x) + nu2*diag(1,nFeatures))

  Sigma = chol2inv(chol_Sigma_inv)


  ### cholesky inverse faster + stable
  ## sample from posterior of approx linear model of GP
  m = solve( t(chol_Sigma_inv) %*% chol_Sigma_inv, phi.x %*% obj.evals)
  theta = m + t(rand.norm.nfeat %*% chol(Sigma*nu2))

  ### Create functional output
  approx.func = function(x, inputs) {
    ### inputs is a list containing:  sigma2, nFeatures, W, b
    result = t(theta) %*%
      (sqrt(2.0 * inputs$sigma2 / inputs$nFeatures) *
         cos(inputs$W %*% t(x) + replicate(nrow(x), inputs$b)))
    return(t(result))
  }
  return(list(approx.func=approx.func, inputs=list(sigma2=sigma2,nFeatures=nFeatures,W=W,b=b)))


}
