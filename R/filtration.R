#' Construct a filtration
#' @description Creates a sequence of simplicial complexes connected by simplicial maps.
#' @param x the file name of the input simplicial complex, or a simplex tree. See details.
#' @param sm the file name of the the simplicial maps, or a vector of elementary operations. See details.
#' @details This function creates a dynamic representation of a sequence of simplicial complexes connected by simplicial maps.
#' Although the representation is not a filtration in the strict inclusion sense, for simplicity, this will be referred to as
#' a filtration. That is, in this context, a filtration \eqn{K'} is a sequence:
#' \deqn{K' = K -> K1 -> K2 ... -> KP}
#' where each \eqn{K*} is a simplicial complex and each \eqn{->} represents a simplicial map.
#' The simplicial maps are indexed by an arbitrary timestamp, i.e. a numeric indicating when the basis elements are born / die. \cr
#' \cr
#' Computationally, \code{filtration} parses the input simplicial maps into an ordered list of expressions representing
#' the appropriate elementary inclusions and/or vertex collapses. These expression are stored in a closure, which is returned to the user.
#' The resulting closure takes as input a timestamp and returns the simplical complex representing the
#' state of complex at that 'time' in sequence. Note that the inverse operations are not stored, so each time the closure is
#' called with a timestamp \eqn{a}, if the filtration is at state \eqn{b}, and \eqn{a < b}, then to get the complex to state \eqn{a}
#' it must be rebuilt by loading the initial complex and then applying the maps \eqn{1, 2, ..., a}.
#' @return A closure \eqn{f(t)}, which takes as input a timestamp \eqn{t} and returns a simplex tree object representing the complex at state \eqn{Kt}.
#' @export
filtration <- function(x, sm){
  stopifnot("Rcpp_SimplexTree" %in% class(x))
  if (is.vector(sm)){ stopifnot(all(substr(sm, 0, 1) %in% c("i", "c", "#"))) } else { stopifnot(file.exists(sm)) }
  si_maps <- if (is.vector(sm)){ sm } else { readLines(sm) }

  ## Read in the (elementary) simplicial operations as expressions for the simplex tree
  si_ops <- lapply(si_maps, function(si_map){
    op <- strsplit(si_map, split = " ")[[1]]
    if (op[1] == "i"){
      return(sprintf("st$insert(%s)", paste0("c(", paste0(as.integer(op[-1]), collapse = ","), ")")))
    } else if (op[1] == "c"){
      from <- as.integer(op[2:(which(op == "t")-1)])
      to <- as.integer(tail(op, 1))
      return(sapply(from, function(v) { sprintf("st$collapse(%d, %d, %d)", v, to, to) })) ## perform the collapse
    }
    else if (op[1] == "#") { return(as.numeric(op[2])) }
  })

  ## Reduce birth/death codes to a list
  bd_idx <- which(sapply(si_ops, is.numeric))
  si_filt <- structure(vector(mode = "list", length(bd_idx)))
  for (i in 1L:(length(bd_idx)-1)){
    if (i == 1L){ si_filt[[1L]] <- unlist(si_ops[1L:bd_idx[1]-1L]) }
    else { si_filt[[i]] <- unlist(si_ops[(bd_idx[i-1]+1):(bd_idx[i]-1L)]) }
  }
  si_filt <- lapply(si_filt, function(si_ops){ unlist(lapply(si_ops, function(op){ parse(text=op) }) ) })
  si_filt <- lapply(si_filt, unlist)

  ## Use closure to create the filtration function
  build_filtration <- function(base_st, filt, bd){
    # base_st_fn <- paste0(tempfile(), ".rds")
    st <- simplextree::simplex_tree()
    st$deserialize(base_st)
    c_idx <- 0L ## index representing what indices have already been computed (inclusive)
    function(eps){
      n_idx <- ifelse(eps < min(bd), 0, max(which(eps >= bd)))
      if (length(n_idx) == 0 || eps > max(bd)){ n_idx <- length(filt) }
      if (c_idx == n_idx){ return(st) }
      if (c_idx > n_idx){
        st <<- simplextree::simplex_tree()
        st$deserialize(base_st)
        c_idx <<- 0L
      }
      while(c_idx < n_idx){
        c_idx <<- c_idx + 1L
        for (expr in filt[[c_idx]]){ eval(expr, envir = parent.env(environment())) }
      }
      return(st)
    }
  }
  return(build_filtration(x$serialize(), si_filt, unlist(si_ops[bd_idx])))
}
