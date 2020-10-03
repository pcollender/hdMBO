#' @export
PlotMOProgress <- function(results.mbo, plot.settings, iters.per.pf,
                           time, isFinal = FALSE){
  # clean data
  obj.df = data.table::data.table(results.mbo$outcomes$obj.evals)

  ## Pareto fronts over time
  obj.df$iteration = results.mbo$outcomes$iteration

  cuts = seq(iters.per.pf,
             max(obj.df$iteration), by = iters.per.pf)
  ## only plot max.cuts layers
  cuts = rev(cuts)[1:min(plot.settings$max.cuts, length(cuts))]
  cuts = rev(cuts)
  pf_dfs <- list()
  i=1
  for(c in cuts) {
    trim_df = obj.df[obj.df$iteration <= c]
    trim_df$iteration = NULL
    is_dom  = emoa::is_dominated(t(trim_df))
    is_pf = !is_dom
    pf_dfs[[i]] = trim_df[is_pf,]
    pf_dfs[[i]]$cut = as.character(c)
    i = i+1
  }
  pf_df = do.call(rbind, pf_dfs)
  ## paired objectives for PF plotting
  obj.names = setdiff(names(pf_df), "cut")
  obj.pairs = combn(obj.names,2)
  ## reverse order to full grid of paired objs
  obj.pairs = cbind(obj.pairs, apply(obj.pairs, 2, rev))

  pf_df_paired = apply(obj.pairs, 2, function(x)pf_df[,c(..x, "cut")])

  ## generate pair names as column
  pf_df_paired =
    lapply(pf_df_paired, function(x) {x[,pairname := paste(names(x)[1:2],
                                                           collapse="SPL")]
      names(x)[1:2] = c("val1", "val2")
      return(x)})
  ## collapse
  pf_df_paired = do.call(rbind, pf_df_paired)
  pf_df_paired[,':='(x_axis = gsub("SPL.*", "", pairname),
                     y_axis = gsub(".*SPL", "", pairname))]
  ## separate DF: designs for histogram of param vals
  des.plot.df = data.table::data.table(results.mbo$outcomes$designs)

  # plot the PFs
  # colors
  manual_colors = c("red", "purple", "blue",
                    "green", "yellow", "orange", "pink", "black")
  cut_inds = seq_len(length(unique(pf_df_paired$cut)))
  cut_inds = manual_colors[cut_inds %% (length(manual_colors)+1)]
  p <-
    ggplot2::ggplot() +
    geom_point(data = data.frame(pf_df_paired),
               aes_string(x=names(pf_df_paired)[1],
                          y=names(pf_df_paired)[2], col = "cut")) +
    scale_color_manual(values=cut_inds) +
    facet_grid(rows = y_axis ~ x_axis, scales = "free")
  # save info
  filepath = paste0(plot.settings$save.filedir, "plots/",
                    plot.settings$filename.tag,
                    ifelse(isFinal, "_FINAL_",
                           paste0("_ITER",
                                  parent.frame()$warm.start.iter, "_")),
                    "PFPlot",
                    time, ".png")
  print(p)
  ggplot2::ggsave(filepath,
                  width  = plot.settings$plot.dims$pf.plot.dims[[1]],
                  height = plot.settings$plot.dims$pf.plot.dims[[2]])



  # histograms of des par values
  p <-
    ggplot2::ggplot(data = suppressWarnings(data.table::melt(des.plot.df))) +
    geom_histogram(aes(x=value), bins=20, fill = "gold",
                   col="brown") +
    facet_wrap(~variable) +
    theme_bw()
  # save info
  filepath = paste0(plot.settings$save.filedir, "plots/",
                    plot.settings$filename.tag,
                    ifelse(isFinal, "_FINAL_",
                           paste0("_ITER",
                                  parent.frame()$warm.start.iter, "_")),
                    "DesParsHist",
                    time, ".png")
  print(p)
  ggplot2::ggsave(filepath,
                  width  = plot.settings$plot.dims$des.hist.plot.dims[[1]],
                  height = plot.settings$plot.dims$des.hist.plot.dims[[2]])

}

#' @export
PlotSOProgress <- function(results.mbo, plot.settings,
                           time, isFinal = FALSE){
  obj.plot.df = data.table::data.table(results.mbo$outcomes$obj.evals)
  names(obj.plot.df) = "ObjValue"

  ## minimum achieved value by iteration and thus far
  obj.plot.df$iteration = results.mbo$outcomes$iteration
  obj.plot.df =
    obj.plot.df[, .(min.value.iter = min(ObjValue)), by =  .(iteration)]
  obj.plot.df$min.thus.far = NA
  for(i in seq(0, max(obj.plot.df$iteration))) {
    obj.plot.df$min.thus.far[i+1] =
      min(obj.plot.df$min.value.iter[obj.plot.df$iteration <= i])
  }
  ## histogram of evaluations
  des.plot.df = data.table::data.table(results.mbo$outcomes$designs)

  ## make plots
  # line plot of values
  p <-
    ggplot2::ggplot(obj.plot.df, aes(x=iteration)) +
    geom_point(aes(y=min.value.iter, col="Iteration Minimum")) +
    geom_line(aes(y=min.thus.far, col="Cumulative Minimum")) +
    theme_bw() +
    labs(x="Iteration Num", y="Obj Value", col = "") +
    scale_color_manual(values =c("black","red"))

  # save info
  filepath = paste0(plot.settings$save.filedir, "plots/",
                    plot.settings$filename.tag,
                    ifelse(isFinal, "_FINAL_",
                           paste0("_ITER",
                                  parent.frame()$warm.start.iter, "_")),
                    "MinObjVals",
                    time, ".png")
  print(p)
  ggplot2::ggsave(filepath,
                  width  = plot.settings$plot.dims$obj.val.plot.dims[[1]],
                  height = plot.settings$plot.dims$obj.val.plot.dims[[2]])


  # histograms of des par values
  p <-
    ggplot2::ggplot(data = suppressWarnings(data.table::melt(des.plot.df))) +
    geom_histogram(aes(x=value), bins=20, fill = "gold",
                   col="brown") +
    facet_wrap(~variable) +
    theme_bw()
  # save info
  filepath = paste0(plot.settings$save.filedir, "plots/",
                    plot.settings$filename.tag,
                    ifelse(isFinal, "_FINAL_",
                           paste0("_ITER",
                                  parent.frame()$warm.start.iter, "_")),
                    "DesParsHist",
                    time, ".png")
  print(p)
  ggplot2::ggsave(filepath,
                  width  = plot.settings$plot.dims$des.hist.plot.dims[[1]],
                  height = plot.settings$plot.dims$des.hist.plot.dims[[2]])

}
