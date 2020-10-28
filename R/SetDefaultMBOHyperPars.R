## List of hyperparameters and settings for Bayesian Optimization
## note:hyperparameters stored as integer only if required to be an integer.
#' @export
SetDefaultMBOHyperPars <- function() {
  mbo.hyperparams = list(
    ### Extra Functionality
    progress.upd.settings = list(
      iters.per.progress.save    = 3L, ## how many iterations to run before saving progress?
      doSaveProgPlot             = TRUE, ## save a plot of optimization path or PF?
      filename.tag               = "my_opt_run", ## name of intmdt. save file (system time auto-added to the end)
      save.filedir               = "", ## directory for intmdt. saves
      save.time                  = TRUE, ## save the time?
      doSaveFinal                = TRUE, ## saves the completed run, deletes intermediate results
      max.cuts                   = 4L, ## how many pareto front layers to display on plots?
      plot.dims = list(obj.val.plot.dims  = list(width = 10, height = 5),
                       des.hist.plot.dims = list(width = 10, height = 8),
                       pf.plot.dims       = list(width = 10, height = 10))

    ),
    no.export = "", #any objs from your environment (ex. Rcpp functions) to not send to the cores
    
    ### Hyperparameters determining # of obj function evaluations
    initDesignsPerParam = 4L,
    iterations          = 8L,
    pointsPerIter       = 8L, #should each proposal involve multiple points?
    nSampleAvg          = 10L, #specify taking mean of size nSampleavg samples of obj fct at each design
    ### Parameters for Finalizing Results
    finalEvals          = 10L, #eval final selection multiple times to reduce noise?
    doPostProcessing    = TRUE, #run posterior simulations to estimate quality of observed Pareto Front?
    n.resample.pf       = 100L, #for MO noisy opt: how many samples of posterior functions to estimate prob of pareto front?
    final.pf.method     = "post.simulation", #options: "post.simulation", "post.fct.approx"
    final.pf.random.searches.per.d = 500L, #how many points per dimension for random search on posterior dist?
    ### Numeric parameter bounds
    designParamBound    = 10,
    ### Surrogate function hyperparameters
    kernel              = "matern5_2",
    bfgs.init.pts.per.d = 3, #Initial points for BFGS km optimization for high-dim search space
    bb.fn.noisy         = TRUE,
    jitter              = 1e-8, #very small jitter to designs during kriging predictions
    #(since predict.km will return exact values of training points)
    ### parallelization:
    parallelize         = TRUE,
    nCores              = NULL, # if NULL, the default is # of cores on machine -1
    ### High-Dimensional Optimization parameters
    max.d               = 10, #at what design space dimension does high-dimensional method deploy?
    high.dim.method     = "dropout", #options: "dropout", "elastic_gauss", "rand_embed"
    high.dim.mo.pts.per.opt = 2, #ex. Take top "i" HV contributors per NSGA2 run.
    # if set to 1, then NSGA2 will run pointsPerIter times per fct eval iteration.
    dropout.pars = list(
      n.dims.opt = 10, #anything btw 5 and 15 likely reasonable.  Higher dimension => more des vars have interacting impact
      #Higher dimension => more difficult optimization problem for the infill criterion.
      rand.prob   = 0.1, #A maximum of 0.2 is recommended and even 0 (i.e. full dropout copy strategy) can be OK fo SO.
      hv.rank.dc.factor = 1e-3 #weighted sample from pareto front by HV contribution,
      # weight_i = gamma^((rank_i - 1)/size_pf)
      # Dc.factor = 1 => all PF has equal weights, Dc.factor = 0 => always choose max HV pt.
    ),
    ### NSGA2 Parameters
    nsga2.hyper.pars = list(
      ## these are essentially default NSGA2 parameters provided by other MBO packages.
      pop.per.dim     = 4, ## population size per dim(obj)*dim(decision)
      generations = 100,
      cprob = 0.7,
      cdist = 5,
      mprob = 0.2,
      mdist = 10
    ),
    ### AF hyperparameters
    alpha = NULL #the gain provided to model standard error for MO AF. if NULL, uses default method sourced from GPareto package
  )
  return(mbo.hyperparams)
}
