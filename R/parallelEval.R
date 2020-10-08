#' @export
parallelEval <- function(bb.fn, designs, nSampleAvg, nCores,
                         maxSeed = 1e9) {
  # This function simultaneously handles:
  #   i)  multipoint proposals
  #   ii) parallelized sample averages for noisy objective functions

  # helper function for means
  bind_and_sum = function(d1, d2) colSums(rbind(d1,d2))

  ## need to establish seeds for reproducibility
  parallel_seeds <- sample(maxSeed, size = nrow(designs)*nSampleAvg)

  # initialize backend
  if(is.null(nCores)){
    cl <- snow::makeCluster(detectCores() - 1)
    doSNOW::registerDoSNOW(cl)
  } else{
    cl <- snow::makeCluster(nCores)
    doSNOW::registerDoSNOW(cl)
  }

  if(nSampleAvg == 1) {
    new.objf = foreach::foreach(i = seq_len(nrow(designs)),
                                .packages = .packages(),
                                .combine = rbind,
                                .export   = ls(.GlobalEnv)) %dopar% {
                                  set.seed(parallel_seeds[i])
                                  bb.fn(designs[i,], parent.frame$target.fns)}
    stopCluster(cl)
  } else{
    new.objf =
      foreach::foreach(i = seq_len(nrow(designs)),
                       .packages = .packages(),
                       .combine  = rbind,
                       .export   = ls(.GlobalEnv)) %:%
      foreach::foreach(j = seq_len(nSampleAvg),
                       .combine='bind_and_sum',
                       .packages = c("stats"),
                       .export   = ls(.GlobalEnv)) %dopar%{
        set.seed(parallel_seeds[(i-1)*nSampleAvg+j])
        bb.fn(designs[i,], parent.frame$target.fns)}

    snow::stopCluster(cl)
    new.objf = new.objf / nSampleAvg
  }

  row.names(new.objf) = NULL
  return(new.objf)
}
