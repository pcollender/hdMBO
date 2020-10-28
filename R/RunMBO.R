# Main BO Architecture
#' @export
RunMBO <- function(d.pars, bb.fn, hyper.pars,
                   results.mbo = NULL) {
  ### Purpose: Run Full Optimization Framework

  ## Initialize parallelization, if needed
  if(hyper.pars$parallelize == TRUE) {
    cat("Spinning up parallelization cores...")
    if(is.null(hyper.pars$nCores)){
      cl <- snow::makeCluster(detectCores() - 1)
      doSNOW::registerDoSNOW(cl)
    } else{
      cl <- snow::makeCluster(hyper.pars$nCores)
      doSNOW::registerDoSNOW(cl)
    }
  }

  ## Create Initial Designs
  if(is.null(results.mbo)) {
    cat("Generating Initial Designs...")
    results.mbo =
      SetupMBO(d.pars, bb.fn, hyper.pars)
    ## save the initial designs
    filepath = paste0(hyper.pars$progress.upd.settings$save.filedir,
                      hyper.pars$progress.upd.settings$filename.tag,
                      "_INITDES",
                      ifelse(hyper.pars$progress.upd.settings$save.time,
                             paste0("_", time), ""), ".RData")
    save(results.mbo, file=filepath)

  }
  d = ncol(results.mbo$outcomes$designs)
  o = ncol(results.mbo$outcomes$obj.evals)
  ## Run BO Loop
  iter=1
  while(iter <= hyper.pars$iterations) {
    warm.start.iter = max(results.mbo$outcomes$iteration) + 1

    if(iter > 1) {
      tdiff = round(as.double(difftime(tnew, told, units ="mins")), 2)
      told = tnew
    } else{told=Sys.time()
    tdiff = "NA"}

    cat("Starting Iteration", warm.start.iter, "...", "Iter time was",
        tdiff, "mins", "...")

    # INTERMEDIATE PROGRESS SAVES
    if((max(warm.start.iter,2)-1) %%
       hyper.pars$progress.upd.settings$iters.per.progress.save == 0) {
      time = gsub(" ", "_", gsub("\\:|\\-", "\\.", Sys.time()))
      filepath = paste0(hyper.pars$progress.upd.settings$save.filedir,
                        hyper.pars$progress.upd.settings$filename.tag,
                        "_ITER", warm.start.iter,
                        ifelse(hyper.pars$progress.upd.settings$save.time,
                               paste0("_", time), ""), ".RData")
      save(results.mbo, file=filepath)
      ## save a plot of progress?
      if(hyper.pars$progress.upd.settings$doSaveProgPlot) {
        if(o == 1) {
          PlotSOProgress(results.mbo,
                         plot.settings = hyper.pars$progress.upd.settings,
                         time = ifelse(hyper.pars$progress.upd.settings$save.time,
                                       paste0("_", time), ""), isFinal = FALSE)
        } else{## multi objective
          PlotMOProgress(results.mbo,
                         plot.settings = hyper.pars$progress.upd.settings,
                         iters.per.pf  =
                           hyper.pars$progress.upd.settings$iters.per.progress.save,
                         time = ifelse(hyper.pars$progress.upd.settings$save.time,
                                       paste0("_", time), ""), isFinal = FALSE)
        }
      }
      cat("Saved Intermediate Results")
    }

    # SINGLE OBJECTIVE INFILL CRITERION OPTIMIZATION:
    if(ncol(results.mbo$outcomes$obj.evals) == 1){

      for(j in 1:hyper.pars$pointsPerIter) {

        ## project design space to lower dimensional object, if needed
        if(d > hyper.pars$max.d) {
          results.for.ic.opt =
            projLowDimDes(results.mbo, max.d = hyper.pars$max.d,
                          high.dim.method    = hyper.pars$high.dim.method,
                          dropout.pars       = hyper.pars$dropout.pars)} else{
                            results.for.ic.opt = results.mbo}

        # run Gaussian process model for each objective function
        results.for.ic.opt$gp.models =
          FitGPModel(results.for.ic.opt$outcomes$obj.evals,
                     results.for.ic.opt$outcomes$designs,
                     hyper.pars$bfgs.init.pts.per.d,
                     hyper.pars$kernel, hyper.pars$bb.fn.noisy)
        names(results.for.ic.opt$gp.models) =
          paste0("obj", 1:length(results.for.ic.opt$gp.models))

        ### AEI optimization for each objective
        results.for.ic.opt$infill.opt =
          lapply(results.for.ic.opt$gp.models,
                 function(x) {
                   MaxInfillNoisySingleObj(
                     model = x,
                     designs = results.for.ic.opt$outcomes$designs,
                     search.space = results.for.ic.opt$search.space,
                     jitter = hyper.pars$jitter)})

        ## populate higher-dimensional suggestion, if necessary
        if(d > hyper.pars$max.d) {
          new.des =
            repopHighDimDesSO(results.mbo, results.for.ic.opt,
                              max.d           = hyper.pars$max.d,
                              high.dim.method = hyper.pars$high.dim.method,
                              dropout.pars    = hyper.pars$dropout.pars)
        }else{
          new.des = as.list(results.for.ic.opt$infill.opt$obj1$xmin)
        }

        ### GENERATING THE "LIE":
        ## append design selection to outcomes object
        ind = tail(results.mbo$outcomes$index, 1) + 1

        results.mbo$outcomes$index[ind]     = ind
        results.mbo$outcomes$iteration[ind] = warm.start.iter
        results.mbo$outcomes$designs =
          rbind(results.mbo$outcomes$designs,
                new.des)

        ## populate the low-dimensional optimization model for prediction
        results.mbo$gp.models$obj1 = results.for.ic.opt$gp.models$obj1
        ## for prediction, need to limit to optimized cols
        opt.cols = names(results.for.ic.opt$outcomes$designs)

        results.mbo$outcomes$obj.evals =
          ## the "lie" is set to be the GP prediction
          ## (i.e. variant of "constant liar" multi-point prop strategy)
          rbind(
            results.mbo$outcomes$obj.evals,
            DiceKriging::predict.km(
              results.mbo$gp.models$obj1,
              newdata = results.mbo$outcomes$designs[ind,..opt.cols],
              type = "SK")$mean
          )
      }
    } else{# MULTI OBJECTIVE INFILL CRITERION OPTIMIZATION:
      cum.props = 0
      while(cum.props < hyper.pars$pointsPerIter) {

        ## project design space to lower dimensional object, if needed
        if(d > hyper.pars$max.d) {
          results.for.ic.opt =
            projLowDimDes(results.mbo, max.d = hyper.pars$max.d,
                          high.dim.method    = hyper.pars$high.dim.method,
                          dropout.pars       = hyper.pars$dropout.pars)} else{
                            results.for.ic.opt = results.mbo}

        # run Gaussian process model for each objective function
        results.for.ic.opt$gp.models =
          FitGPModel(results.for.ic.opt$outcomes$obj.evals,
                     results.for.ic.opt$outcomes$designs,
                     hyper.pars$bfgs.init.pts.per.d,
                     hyper.pars$kernel, hyper.pars$bb.fn.noisy)
        names(results.for.ic.opt$gp.models) =
          paste0("obj", 1:length(results.for.ic.opt$gp.models))

        props =
          data.table(
            MaxInfillNoisyMultiObj(model_list = results.for.ic.opt$gp.models,
                                   designs = results.for.ic.opt$outcomes$designs,
                                   search.space = results.for.ic.opt$search.space,
                                   points.per.iter =
                                     ## differs if solving HighD problem.
                                     ifelse(d > hyper.pars$max.d,
                                            min(hyper.pars$high.dim.mo.pts.per.opt,
                                                hyper.pars$pointsPerIter - cum.props),
                                            hyper.pars$pointsPerIter),
                                   nsga2.hyper.pars = hyper.pars$nsga2.hyper.pars,
                                   alpha  = hyper.pars$alpha,
                                   jitter = hyper.pars$jitter)
          )
        names(props) = names(results.for.ic.opt$outcomes$designs)

        ## populate higher-dimensional suggestions, if necessary
        if(d > hyper.pars$max.d) {
          props =
            repopHighDimDesMO(results.mbo, results.for.ic.opt,
                              proposals       = props,
                              max.d           = hyper.pars$max.d,
                              high.dim.method = hyper.pars$high.dim.method,
                              dropout.pars    = hyper.pars$dropout.pars)
        }


        ## append design selections to outcomes object, generate "lies"
        ind = tail(results.mbo$outcomes$index, 1) + 1:(nrow(props))
        results.mbo$outcomes$index[ind]     = ind
        results.mbo$outcomes$iteration[ind] = warm.start.iter
        results.mbo$outcomes$designs =
          rbind(results.mbo$outcomes$designs,
                props)

        ## populate the low-dimensional optimization model for prediction
        results.mbo$gp.models = results.for.ic.opt$gp.models

        ## for prediction, need to limit to optimized cols
        opt.cols = names(results.for.ic.opt$outcomes$designs)

        preds = sapply(results.mbo$gp.models,
                       function(x)
                         DiceKriging::predict.km(
                           x, newdata = results.mbo$outcomes$designs[ind,..opt.cols],
                           type = "SK")$mean)
        colnames(preds)=NULL


        results.mbo$outcomes$obj.evals =
          ## the "lie" is set to be the GP prediction
          ## (i.e. variant of "constant liar" multi-point prop strategy)
          rbind(
            results.mbo$outcomes$obj.evals,
            preds
          )

        cum.props = cum.props + nrow(preds)
      }
    }

    ## Evaluate BB Objective Fct at Proposed Points
    ind = which(results.mbo$outcomes$iteration == max(results.mbo$outcomes$iteration))
    # fill in "lies" with actual objective function values
    # PARALLELIZABLE OBJ FCT EVALUATIONS
    new.x = as.matrix(results.mbo$outcomes$designs[ind])
    if(hyper.pars$parallelize == TRUE){
      new.objf = parallelEval(bb.fn, designs = new.x,
                              nSampleAvg = hyper.pars$nSampleAvg,
                              no.export  = hyper.pars$no.export)
    } else{new.objf = t(as.matrix(apply(new.x, 1, bb.fn),
                                  nrow = ncol(results.mbo$outcomes$obj.evals)))
    }
    results.mbo$outcomes$obj.evals[ind,] = new.objf


    iter=iter+1
    tnew=Sys.time()
  }

  ### Finalize Solution
  ## if a warm-started object, empty out previous solution module
  if(!is.null(results.mbo$solution)){results.mbo$solution = NULL}
  ## Run postprocessing scripts based on # of objectives
  if(ncol(results.mbo$outcomes$obj.evals) == 1){
    cat("Finalizing single-obj sol")
    # run final GP model
    results.mbo$gp.models =
      FitGPModel(results.mbo$outcome$obj.evals,
                 results.mbo$outcome$designs,
                 hyper.pars$bfgs.init.pts.per.d,
                 hyper.pars$kernel,
                 hyper.pars$bb.fn.noisy)
    names(results.mbo$gp.models) =
      paste0("obj", 1:length(results.mbo$gp.models))
    # Find best predicted design.  Evaluate finalEvals times.
    results.mbo$solution =
      finalEvalBestPred(bb.fn, results.mbo,
                        num.evals   = hyper.pars$finalEvals,
                        parallelize = hyper.pars$parallelize,
                        jitter      = hyper.pars$jitter)
  } else{## MULTIOBJECTIVE
    cat("Finalizing multi-obj sol")
    # Run GP model on final set of evaluated points.
    results.mbo$gp.models =
      FitGPModel(results.mbo$outcomes$obj.evals,
                 results.mbo$outcomes$designs,
                 hyper.pars$bfgs.init.pts.per.d,
                 hyper.pars$kernel, hyper.pars$bb.fn.noisy)
    names(results.mbo$gp.models) =
      paste0("obj", 1:length(results.mbo$gp.models))
    # Post-processing of the MO Bayesian Optimization Loop
    results.mbo$solution =
      finalEvalParetoFront(bb.fn.mbo        = bb.fn,
                           results.mbo      = results.mbo,
                           hyper.pars       = hyper.pars)
  }

  ## Save results
  if(hyper.pars$progress.upd.settings$doSaveFinal){## save results
    ## remove the intermediate saved results
    intmdt.files =
      list.files(path = hyper.pars$progress.upd.settings$save.filedir,
                 pattern = hyper.pars$progress.upd.settings$filename.tag,
                 recursive = TRUE)
    if(length(intmdt.files)>0){
      file.remove(paste0(hyper.pars$progress.upd.settings$save.filedir,
                         intmdt.files))}
    ## save final results
    time = gsub(" ", "_", gsub("\\:|\\-", "\\.", Sys.time()))
    filepath = paste0(hyper.pars$progress.upd.settings$save.filedir,
                      hyper.pars$progress.upd.settings$filename.tag,
                      "_FINAL",
                      time = ifelse(hyper.pars$progress.upd.settings$save.time,
                                    paste0("_", time), ""), ".RData")
    # save results
    save(results.mbo, file=filepath)
    # save plots of results
    if(o==1) {
      PlotSOProgress(results.mbo,
                     plot.settings = hyper.pars$progress.upd.settings,
                     time = ifelse(hyper.pars$progress.upd.settings$save.time,
                                   paste0("_", time), ""), isFinal = TRUE)
    } else{
      PlotMOProgress(results.mbo,
                     plot.settings = hyper.pars$progress.upd.settings,
                     iters.per.pf = hyper.pars$progress.upd.settings$iters.per.progress.save,
                     time = ifelse(hyper.pars$progress.upd.settings$save.time,
                                   paste0("_", time), ""), isFinal = TRUE)
    }

  }

  ## end parallelization, if needed
  if(hyper.pars$parallelize == TRUE) {
    cat("Spinning down parallelization cores...")
    snow::stopCluster(cl)
  }

  return(results.mbo)
}

