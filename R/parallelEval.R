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
  
  #adding in mclapply option for Linux machines, since this can be faster (and plays much nicer with data.table)
  OS = Sys.info()['sysname']
  
  
  if(OS != 'Linux'){
    if(nSampleAvg == 1) {
    
    new.objf = foreach::foreach(i = seq_len(nrow(designs)),
                                .packages = .packages(),
                                .combine  = rbind,
                                .export   = setdiff(ls(.GlobalEnv),
                                                    no.export),
                                .noexport = no.export) %dopar% {
                                  set.seed(parallel_seeds[i])
                                  bb.fn(designs[i,])}
    
    }else{
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
      }}else{
    
      des.indices = rep(seq_len(nrow(designs)), nSampleAvg)
      
      new.objf = mclapply(X = seq_len(length(des.indices)),
                          FUN = function(i){set.seed(parallel_seeds[i]) 
                                            res = bb.fn(designs[des.indices[i],])
                                            pb = txtProgressBar(0, length(des.indices),initial = i,style = 3)
                                            return(res)},
                          mc.cores = getOption('mc.cores'))
      nulls = sapply(new.objf,is.null)
      while(any(nulls)){
        redoinds = which(nulls)
        new.objf[redoinds] = mclapply(X = redoinds,
                          FUN = function(i){set.seed(parallel_seeds[i]) 
                                            res = bb.fn(designs[des.indices[i],])
                                            pb = txtProgressBar(0, length(redoinds),initial = which(redoinds==i), style = 3)
                                            return(res)},
                          mc.cores = getOption('mc.cores'))
    
        nulls = sapply(new.objf,is.null)
        }
    new.objf = do.call(rbind, new.objf)
      new.objf = apply(new.objf, 2, tapply, des.indices, mean)
      }
                          
  row.names(new.objf) = NULL
  return(new.objf)
}
