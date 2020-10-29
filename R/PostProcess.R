# Single Objective Best Predicted
#' @export
finalEvalBestPred <-
  function(bb.fn, results.mbo, num.evals,
           parallelize, no.export,
           jitter, is.noisy = TRUE) {
    # get predictions at evaluated points
    preds=DiceKriging::predict.km(results.mbo$gp.models$obj1,
                                  newdata = results.mbo$outcomes$designs +
                                    ifelse(is.noisy, sample(c(-1,1),1)*jitter, 0),
                                  type = "SK",
                                  checkNames = F)$mean
    ind = which.min(preds)
    des = results.mbo$outcomes$designs[ind,]

    if(parallelize == TRUE) {
      eval = parallelEval(bb.fn, designs = des, nSampleAvg=num.evals)
    } else{
      eval = mean(sapply(seq_len(num.evals),
                         function(i) bb.fn(des)))
    }

    return(list(index     = results.mbo$outcomes$index[ind],
                iteration = results.mbo$outcomes$iteration[ind],
                designs   = results.mbo$outcomes$designs[ind,],
                obj.evals = eval,
                obj.pred  = preds[ind]))
  }


### Computes PF approximations using posterior simulations, summarizes results.
## SEE: "Quantifying uncertainty on Pareto fronts with Gaussian processconditional simulations" (Binois et al)
### TODO: implement CDF summarization for >2 objectives
#' @export
ComputeCondPF <- function(bb.fn.mbo, results.mbo,
                          n.resample.pf, final.pf.random.searches.per.d,
                          parallelize,
                          search.function, maxSeed = 1e9) {
  ## of points per simulation
  d = ncol(results.mbo$outcomes$designs)
  o = ncol(results.mbo$outcomes$obj.evals)
  npointssim = final.pf.random.searches.per.d * d

  ## store realized PF
  PF = t(emoa::nondominated_points(t(results.mbo$outcomes$obj.evals)))

  ## initialize simulation arrays
  sims.array = array(0, dim = c(n.resample.pf, npointssim, o))
  design.sim = array(0, dim = c(npointssim, d, n.resample.pf))
  ## start progress bar
  pb = txtProgressBar(0, n.resample.pf, style = 3)
  prog = function(i) setTxtProgressBar(pb, i)
  opts = list(progress = prog)

  ## run simulations on random designs from posterior distribution
  if(parallelize == TRUE) {
    ## need to establish seeds for reproducibility
    parallel_seeds <- sample(maxSeed, size = n.resample.pf)

    ## run posterior simulation loop
    sim.list =
      foreach::foreach(i = seq_len(n.resample.pf),
                       .packages = c("stats", "data.table", "DiceKriging"),
                       .options.snow = opts) %dopar% {
                         set.seed(parallel_seeds[i])
                         design.sim[,,i] = as.matrix(search.function())
                         for(j in seq_len(o)) {
                           sims.array[i,,j] =
                             DiceKriging::simulate(
                               results.mbo$gp.models[[j]],
                               nsim = 1, newdata = design.sim[,,i], cond = TRUE,
                               checkNames = FALSE, nugget.sim = 10^-8)
                           }
                return(sims.array[i,,])
              }

    ## assign the sim list to the simulation array
    for(i in seq_len(length(sim.list))){
      sims.array[i,,]=sim.list[[i]]
    }
  }
  if(parallelize == FALSE){ for(i in seq_len(n.resample.pf)){
    ## print progress bar update
    prog(i)
    ## sample simulation loop as above
    design.sim[,,i] = as.matrix(search.function())
    for(j in seq_len(o)) {
      sims.array[i,,j] =
        simulate(results.mbo$gp.models[[j]],
                 nsim = 1, newdata = design.sim[,,i], cond = TRUE,
                 checkNames = FALSE, nugget.sim = 10^-8) }
  }}
  ## end progress bar
  close(pb)


  ### TODO: extend CPF function to >2 objectives....
  if(o==2) {
    CPF.res =
      GPareto::CPF(sims.array[,,1], sims.array[,,2],
                   results.mbo$outcomes$obj.evals, paretoFront = PF)

  } else{CPF.res = NULL
  }
  return(list(sims.array = sims.array, CPF.res = CPF.res))
}


# Multi Objective Pareto Front (for noisy obj functions)
#' @export
finalEvalParetoFront <-
  function(bb.fn.mbo, results.mbo,
           hyper.pars) {

    ## Intialize results list
    solution = list()

    ## Creates random designs over space
    d = ncol(results.mbo$outcomes$designs)
    random_search_data = function(){
      dat=data.table::data.table(
        do.call(cbind,
                lapply(seq_len(ncol(results.mbo$outcomes$designs)),
                       function(i) runif(d * final.pf.random.searches.per.d,
                                         min=results.mbo$search.space[[i]][1],
                                         max=results.mbo$search.space[[i]][2]))
                ))
      names(dat)=names(results.mbo$outcomes$designs)
      return(dat)
    }

    if(hyper.pars$doPostProcessing == TRUE) {
      ## call needed hypers
      n.resample.pf                  = hyper.pars$n.resample.pf
      final.pf.method                = hyper.pars$final.pf.method
      final.pf.random.searches.per.d = hyper.pars$final.pf.random.searches.per.d
      parallelize                    = hyper.pars$parallelize
      nsga2.hyper.pars               = hyper.pars$nsga2.hyper.pars

      if(final.pf.method == "post.fct.approx") {
        ### Actual functional approximations from Posterior
        ### TODO: this functionality is not finished yet.
        for(i in seq_len(n.resample.pf)) {
          cat("resampling ", i)
          ## Sample from posterior at all designs
          post.sample.fns = lapply(
            as.list(seq_len(length(results.mbo$gp.models))),
            function(i) SampleGPwRandomFeatures(
              designs   = results.mbo$outcomes$designs,
              obj.evals = results.mbo$outcomes$obj.evals[,i],
              model     = results.mbo$gp.models[[i]],
              nFeatures = 1000))
          ## making this a two-step function creation process such that SampleGPwRandomFeatures is only called once
          get.post.sample.fns.output = function(x) {
            sapply(as.list(seq_len(length(results.mbo$gp.models))),
                   function(i)
                     post.sample.fns[[i]]$approx.func(matrix(x, nrow=1),
                                                      post.sample.fns[[i]]$inputs))
          }


          ## Find pareto front of sampled posterior functions using NSGAii
          opt.result =
            mco::nsga2(get.post.sample.fns.output,
                       idim = ncol(results.mbo$outcomes$designs),
                       odim = length(post.sample.fns),
                       lower.bounds = sapply(results.mbo$search.space, function(x)x[1]),
                       upper.bounds = sapply(results.mbo$search.space, function(x)x[2]),
                       vectorized = FALSE,
                       popsize =
                         pmax(100, nsga2.hyper.pars$pop.per.dim *
                                ncol(results.mbo$outcomes$designs) * length(post.sample.fns)),
                       generations = nsga2.hyper.pars$generations,
                       cprob       = nsga2.hyper.pars$cprob,
                       cdist       = nsga2.hyper.pars$cdist,
                       mprob       = nsga2.hyper.pars$mprob,
                       mdist       = nsga2.hyper.pars$mdist)

          ## Reduce to Pareto Optimal Solutions if necessary
          pf.df =  rbind(pf.df, opt.result$par[opt.result$pareto.optimal==TRUE, ])

          ### TODO: (low-dimensional) KDE of Histogram density estimation of Pareto Front
          ### TODO: (high-dimensional) KDE or Histogram density estimation of Pareto Front
        }
      }
      else if(final.pf.method == "post.simulation") {
        solution$post.sim =
          ComputeCondPF(bb.fn.mbo, results.mbo,
                        n.resample.pf, final.pf.random.searches.per.d,
                        parallelize = parallelize,
                        search.function = random_search_data)

      }
      else{stop("Provide a supported method for finalizing the multi-obj optimization")}
    }


    ## compute actual Pareto Front
    is_dom = emoa::is_dominated(t(results.mbo$outcomes$obj.evals))
    solution$obs.pf = !is_dom

    return(solution)
  }
