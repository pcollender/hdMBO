#' @export
parallelEval <- function(bb.fn, designs, nSampleAvg, no.export,
                         maxSeed = 1e9) {
  # This function simultaneously handles:
  #   i)  multipoint proposals
  #   ii) parallelized sample averages for noisy objective functions

  # helper function for means
  bind_and_sum = function(d1, d2) colSums(rbind(d1,d2))

  ## need to establish seeds for reproducibility
  parallel_seeds <- sample(maxSeed, size = nrow(designs)*nSampleAvg)


  if(nSampleAvg == 1) {
    new.objf = foreach::foreach(i = seq_len(nrow(designs)),
                                .packages = .packages(),
                                .combine  = rbind,
                                .export   = setdiff(ls(.GlobalEnv),
                                                    no.export),
                                .noexport = no.export) %dopar% {
                                  set.seed(parallel_seeds[i])
                                  bb.fn(designs[i,])}
    stopCluster(cl)
  } else{
    new.objf =
      foreach::foreach(i = seq_len(nrow(designs)),
                       .packages = .packages(),
                       .combine  = rbind,
                       .export   = setdiff(ls(.GlobalEnv),
                                           no.export),
                       .noexport = no.export) %:%
      foreach::foreach(j = seq_len(nSampleAvg),
                       .combine='bind_and_sum',
                       .packages = c("stats", "Rcpp"),
                       .export   = setdiff(ls(.GlobalEnv),
                                           no.export),
                       .noexport = no.export) %dopar%{
        set.seed(parallel_seeds[(i-1)*nSampleAvg+j])
        bb.fn(designs[i,])}

    new.objf = new.objf / nSampleAvg
  }

  row.names(new.objf) = NULL
  return(new.objf)
}
