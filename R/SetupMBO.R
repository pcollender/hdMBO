# Initial design, evaluations, search space
#' @export
SetupMBO <- function(d.pars, bb.fn, target.fns, hyper.pars) {
  ### Purposes:
  ###   i)   Setup search space according to d.pars and hyper.pars
  ###   ii)  Generate initial designs and evaluate
  ###   iii) Initialize results list

  # Initialize: Recast BB function as only input of design parameters
  bb.fn.mbo      = function(x) {
    return(bb.fn(x, parent.frame$target.fns))
  }

  # check if parameter box boundaries are manually initiated
  if(length(setdiff(unique(sapply(d.pars, length)), 2)) != 0) {
    # assign bounds to search space
    design.pars.mbo =
      lapply(d.pars, function(x) x=c(-1*hyper.pars$designParamBound,
                                     hyper.pars$designParamBound))
    names(design.pars.mbo) = names(d.pars)
  } else{
    design.pars.mbo = d.pars
  }



  ## Generate Initial Designs
  init.des = GenBoxDesign(design.pars.mbo,
                          hyper.pars$initDesignsPerParam *
                            length(design.pars.mbo))

  ## Generate Initial Objective Evaluations

  if(hyper.pars$parallelize == TRUE){
    new.objf = parallelEval(bb.fn, designs = init.des,
                            nSampleAvg = hyper.pars$nSampleAvg,
                            nCores = hyper.pars$nCores)
  } else{new.objf = as.matrix(apply(init.des, 1, bb.fn.mbo),
                                nrow = ncol(results.mbo$outcomes$obj.evals))}

  ## initialize outcomes object
  outcomes.mbo = list(
    index     = 1:nrow(init.des),
    iteration = rep(0,nrow(init.des)),
    designs   = init.des,
    obj.evals = new.objf)

  ## outputs
  results.mbo <-
    list(search.space = design.pars.mbo,
         outcomes     = outcomes.mbo)

  return(results.mbo)


}
