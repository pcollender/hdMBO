#' @export
GenBoxDesign <- function(params, n) {
  ### Purpose:
  ###   Creates random designs subject to box constraints
  ###     (i.e. uninteracting).
  ### Arguments:
  ###   params: Named list of numeric parameters with s
  ###             of form c(lower bound, upper bound)
  ###   n: Number of designs to generate

  des.qnt = lhs::randomLHS(n, length(params))

  range = sapply(params, function(t) (t[2]-t[1]))
  lb    = sapply(params, function(t) t[1])

  des.box = data.table::data.table()
  for(i in 1:length(params)) {
    des.box[,
            names(params)[i]:=des.qnt[,i]*range[i]+lb[i]]
  }

  return(des.box)

}


