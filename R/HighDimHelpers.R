## i. Projecting Design Space to lower-dimension space for good IC optimization
projLowDimDes <-
  function(results.mbo, max.d, high.dim.method,
           dropout.pars) {
    d = ncol(results.mbo$outcomes$designs)
    ## a.  High-Dimensional Dropout
    if(high.dim.method == "dropout") {
      ## randomly select n.dims.opt dimensions to optimize AF on
      d.choice = sort(sample(d, dropout.pars$n.dims.opt, replace = FALSE))
      cols.chosen = names(results.mbo$outcomes$designs)[d.choice]

      low.dim.results = results.mbo
      low.dim.results$search.space     = low.dim.results$search.space[d.choice]
      low.dim.results$outcomes$designs = low.dim.results$outcomes$designs[,..cols.chosen]

    } else if(high.dim.method == "elastic_gauss") {### b. High-Dimensional Elastic Gaussian
      stop("Elastic Gaussian method not supported yet...")
    } else if(high.dim.method == "rand_embed") {### c. High-Dimensional Random Embeddings
      stop("Random Embeddings method not supported yet...")
    } else(stop("Please provide supported high-dimensional method"))
    return(low.dim.results)
  }

## ii. Repopulating full design post-Infill optimization
## Single Objective Case:
repopHighDimDesSO <-
  function(results.mbo, results.for.ic.opt, max.d, high.dim.method,
           dropout.pars) {
    low.dim.vars = names(results.for.ic.opt$search.space)
    dropped.vars =
      setdiff(names(results.mbo$search.space), low.dim.vars)
    dropped.search.space = results.mbo$search.space[dropped.vars]
    ## a. High-Dimensional Dropout
    if(high.dim.method == "dropout") {
      stopifnot(dropout.pars$rand.prob >=0 & dropout.pars$rand.prob <= 1)
      ## design of d-dimensions is a result of infill optimization
      low.dim.des =
        data.table::data.table(matrix(results.for.ic.opt$infill.opt$obj1$xmin,
                                      nrow=1))
      names(low.dim.des) = low.dim.vars

      ## create fill-in:
      ##   - with prob rand.prob, use a random selection of D-d design dimensions
      ##   - with prob 1-rand.prob, use best search so far.
      choice = runif(1)
      if(choice < dropout.pars$rand.prob) {## perform random fill-in
        dropped.des =
          data.table::data.table(do.call(
            cbind,
            lapply(seq_len(length(dropped.search.space)),
                   function(i) runif(1,
                                     min=dropped.search.space[[i]][1],
                                     max=dropped.search.space[[i]][2]))))
        names(dropped.des) = dropped.vars
      } else{## copy best search so far
        best.search.ind = which.min(results.mbo$outcomes$obj.evals)
        best.des    = results.mbo$outcomes$designs[best.search.ind,]
        dropped.des = best.des[, ..dropped.vars]
      }
      ## concatenate, order the optimized and dropped designs
      new.design = cbind(low.dim.des, dropped.des)
      data.table::setcolorder(new.design, names(results.mbo$search.space))

    } else if(high.dim.method == "elastic_gauss") {## b. High-Dimensional Elastic Gaussian
      stop("Elastic Gaussian method not supported yet...")
    } else if(high.dim.method == "rand_embed") {## c. High-Dimensional Random Embeddings
      stop("Random Embeddings method not supported yet...")
    } else(stop("Please provide supported high-dimensional method"))
    return(new.design)
  }

## Multi Objective Case:
repopHighDimDesMO <-
  function(results.mbo, results.for.ic.opt, proposals,
           max.d, high.dim.method,
           dropout.pars) {
    ## note that this function directly receives proposal DT as opposed to SO case.
    low.dim.vars = names(proposals)
    num.props = nrow(proposals)
    dropped.vars =
      setdiff(names(results.mbo$search.space), low.dim.vars)
    dropped.search.space = results.mbo$search.space[dropped.vars]
    ## a. High-Dimensional Dropout
    if(high.dim.method == "dropout") {
      stopifnot(dropout.pars$rand.prob >=0 & dropout.pars$rand.prob <= 1)
      ## design of d-dimensions is a result of infill optimization
      low.dim.des = proposals

      ## get current PF indices
      is_dom       = emoa::is_dominated(t(results.mbo$outcomes$obj.evals))
      is_par_front = !is_dom
      size_pf      = sum(is_par_front)

      ## Compute HV contributions of PF
      test_pts = results.mbo$outcomes$obj.evals[is_par_front,]
      full_hv  = emoa::dominated_hypervolume(results.mbo$outcomes$obj.evals)
      hv_contr = rep(NA, nrow(test_pts))
      # PF Sampling Weight Loop via HV contribution
      for (i in seq_len(size_pf)) {
        ind = which(is_par_front)[i]
        # Total Old HV
        hv_loo_pts = results.mbo$outcomes$obj.evals[-ind,]
        # New HV for addition of each NSGA2 candidate
        hv_loo    = emoa::dominated_hypervolume(hv_loo_pts)
        # max HV Contribution
        hv_contr[i]  = full_hv - hv_loo
      }
      hv_contr_rank     = rank(hv_contr)
      pf_sample_probs = dropout.pars$hv.rank.dc.factor^((hv_contr_rank-1)/size_pf) /
        sum(dropout.pars$hv.rank.dc.factor^((hv_contr_rank-1)/size_pf))
      ## create fill-in:
      ##   - with prob rand.prob, use a random selection of D-d design dimensions
      ##   - with prob 1-rand.prob, sample from PF using weights pf_sample_probs

      dropped.des =
        ## Note: similar process to SO but select from PF
        do.call(rbind,
                lapply(as.list(seq_len(num.props)),
                       function(row) {
                         choice = runif(1)
                         if(choice < dropout.pars$rand.prob) {## perform random fill-in
                           des.pt =
                             data.table::data.table(
                               do.call(cbind,
                                       lapply(seq_len(length(dropped.search.space)),
                                              function(i)
                                                runif(1,
                                                      min=dropped.search.space[[i]][1],
                                                      max=dropped.search.space[[i]][2]))))
                           names(des.pt) = dropped.vars
                         } else{## sample from PF according to pf_sample_probs
                           pf.pt.index = sample(seq_len(size_pf), size = 1,
                                                prob = pf_sample_probs, replace = FALSE)
                           des.pt      = results.mbo$outcomes$designs[pf.pt.index,]
                           des.pt      = des.pt[, ..dropped.vars]
                         }
                         return(des.pt)}
                ))


      ## concatenate, order the optimized and dropped designs
      new.design = cbind(low.dim.des, dropped.des)
      setcolorder(new.design, names(results.mbo$search.space))

    } else if(high.dim.method == "elastic_gauss") {## b. High-Dimensional Elastic Gaussian
      stop("Elastic Gaussian method not supported yet...")
    } else if(high.dim.method == "rand_embed") {## c. High-Dimensional Random Embeddings
      stop("Random Embeddings method not supported yet...")
    } else(stop("Please provide supported high-dimensional method"))
    return(new.design)
  }
