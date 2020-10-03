### i. AEI, SINGLE OBJECTIVE optimization
#' @export
MaxInfillNoisySingleObj <- function(model, designs, search.space,
                                    is.noisy = TRUE, jitter) {
  # Using CMAES algorithm

  # Get effective best design so far
  ebs = getEBS(designs, model, is.noisy, jitter)

  # AEI function using only point vector as arg. (minimize negative AEI)
  fn = function(x){getAEI(x,
                          model   = model,
                          designs = designs,
                          ebs     = ebs) * (-1)}

  opt.result =
    if(length(search.space)==1){
      pureCMAES1D(par   = as.numeric(designs[sample(nrow(designs),1),]),
                  fun   = fn,
                  lower = sapply(search.space, function(x)x[1]),
                  upper = sapply(search.space, function(x)x[2]))
    } else{
      adagio::pureCMAES(par   =
                          as.numeric(designs[sample(nrow(designs),1),]),
                        fun   = fn,
                        lower = sapply(search.space, function(x)x[1]),
                        upper = sapply(search.space, function(x)x[2]))}


  return(opt.result)
}

### ii.  surrogate-fct-based MULTI OBJECTIVE optimization
#' @export
MaxInfillNoisyMultiObj <- function(model_list, designs, search.space,
                                   points.per.iter, nsga2.hyper.pars, alpha,
                                   is.noisy = TRUE, jitter) {
  o = length(model_list)
  d = ncol(designs)
  # get infill criteria function for each obj dimension
  ## Note: this is MSPOT acquisition function..
  fn = function(x){sapply(seq_len(o),
                          function(i) {
                            pred=predict.km(object  = model_list[[i]],
                                            newdata =
                                              matrix(x,nrow=1,
                                                     dimnames =
                                                       list(NULL,
                                                            names(designs))),
                                            type = "SK")
                            return(pred$mean -
                                     ifelse(!is.null(alpha), alpha,
                                            -qnorm(0.5 * (0.5^(1/o))))
                                   * pred$sd)}
  )}

  # Program function for NSGAII is SMS-EGO criterion

  # run NSGA-II multi-criteria algorithm on above defined MO fn.
  opt.result =
    mco::nsga2(fn,
               idim = d,
               odim = o,
               lower.bounds = sapply(search.space, function(x)x[1]),
               upper.bounds = sapply(search.space, function(x)x[2]),
               vectorized = FALSE,
               popsize =
                 pmax(100, nsga2.hyper.pars$pop.per.dim *
                        d * o),
               generations = nsga2.hyper.pars$generations,
               cprob       = nsga2.hyper.pars$cprob,
               cdist       = nsga2.hyper.pars$cdist,
               mprob       = nsga2.hyper.pars$mprob,
               mdist       = nsga2.hyper.pars$mdist)

  # Choose m points from last generation via hypervolume criteria
  # uses nsgaii pareto front and designs
  hv.df = rbind(t(apply(designs, 1, fn)))
  props = matrix(nrow=0,ncol=ncol(designs))
  # proposal loop
  for (i in seq_len(points.per.iter)) {
    # Total Old HV
    hv.old = emoa::dominated_hypervolume(t(hv.df))
    # New HV for addition of each NSGAii candidate
    hvs.new = apply(opt.result$value, 1,
                    function(x) emoa::dominated_hypervolume(t(rbind(x,hv.df))))
    # max HV Contribution
    best.ind = which.max(hvs.new - hv.old)

    # add best to proposed points, remove from candidate list.
    props = rbind(props, opt.result$par[best.ind, ])
    hv.df = rbind(opt.result$value[best.ind, ])
    opt.result$par   = opt.result$par[-best.ind,, drop = FALSE]
    opt.result$value = opt.result$value[-best.ind,, drop = FALSE]
  }
  return(props)
}

