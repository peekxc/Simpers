#' SimPers
#' @description Simple R-wrapper for the SimPers software for computing topological persistence under simplicial maps (see 1).
#' @param x the file name for the input simplicial complex, or a simplex tree. See details.
#' @param elementary_mode whether to run in elementary mode or general mode.
#' @param validate_maps whether to check to ensure the simplicial maps are 'valid'. Defaults to true. See details.
#' @param dim maximum dimension of the barcodes.
#' @param sm input simplicial maps (either filename or character vector).
#' @param n_maps number of maps (general-mode only).
#' @param thresh barcode length threshold
#' @param ... Other simpers parameters.
#' @details Other parameters are converted to their <option>=<value> string equivalents and passed as arguments
#' to the main options parser. \cr
#' \cr
#' If \code{validate_maps} and \code{elementary_mode} are TRUE, the simplicial maps are inspected prior to passing to Simpers.
#' Namely, its checked that the faces of the each simplex being inserted exist prior to insertion, and if not they are
#' preprended to the insertion at the same time. Elementary collapses that don't collapse to a non-existent vertex throw an exception.
#' All operations are changed such that their vertices are listed in increasing order.
#' Simplex ids cannot be changed; it is up to the user to ensure the indices of the operation are consistent across all simplicial maps.
#' @examples
#' ## Elementary mode example
#' sp_dir <- system.file("extdata", package = "Simpers")
#' cmd <- sprintf(" -d %s -m %s", file.path(sp_dir, "Elementary_Mode", "iDC.txt"), file.path(sp_dir, "Elementary_Mode", "operations.txt"))
#' @import BH
#' @author R-wrapper written by Matt Piekenbrock. Source code for computing the persistence copyright 2015 by Fengtao Fan and Dayu Shi, developed by the Jyamiti group at The Ohio State University
#' @references
#' 1. Dey, Tamal K., Fengtao Fan, and Yusu Wang. "Computing topological persistence for simplicial maps." Proceedings of the thirtieth annual symposium on Computational geometry. ACM, 2014. \cr
#' 2. \href{http://web.cse.ohio-state.edu/~dey.8/SimPers/Simpers.html}{SimPers homepage}
#' @export
simpers <- function(x, sm, validate_maps = TRUE, elementary_mode=TRUE, dim=2L, n_maps=NULL, thresh=0, generators=elementary_mode, ...){
  if (!elementary_mode && is.null(n_maps)){ stop("'n_maps' must be specified in general mode.") }

  ## To generate temporary, random files
  random_file <- function() { paste0(tempfile(), as.hexmode(as.integer(Sys.time())), ".txt") }

  ## If x is a simplex tree, write it using BFS to ensure the requirement regarding the
  ## the boundary faces is met.
  if (is(x, "Rcpp_SimplexTree")){
    sc_file <- random_file()
    write(file = sc_file, x = sum(x$n_simplices))
    write_simplex <- function(simplex){ if (length(simplex)){  write(simplex, file = sc_file, append = TRUE) } }
    x$traverse(write_simplex, "bfs")
    x <- sc_file
  }
  stopifnot(is(x, "character"), file.exists(x))

  ## General mode only
  if (!missing(n_maps) && !is.null(n_maps)){ }

  ## Handle simplicial maps input
  stopifnot(is(sm, "character"))
  if (elementary_mode && validate_maps){
    initial_complex <- simplextree::simplex_tree()
    initial_complex$insert(lapply(strsplit(readLines(x)[-1], split = " "), as.integer))
    sm <- validate_sm(initial_complex, sm)
  }
  if (length(sm) > 1){
    stopifnot(is.character(sm), all(substr(sm, 0, 1) %in% c("i", "c", "#")))
    sm_fn <- random_file()
    writeLines(text = sm, con = sm_fn)
    sm <- sm_fn
  }
  stopifnot(is.character(sm), file.exists(sm) || dir.exists(sm))

  ## Where to save the output persistence diagram
  out_file <- random_file()

  ## Run simpers
  generator_file <- ""
  if (generators){ generator_file <- random_file()  }
  run_simpers(x, "", sm, "", out_file, generator_file, FALSE, FALSE, elementary_mode, generators, 0, 0, dim)

  ## Extract persistence diagram and return
  if (file.exists(out_file)){
    bar <- read.delim(out_file, sep = " ")
    if (generators && ncol(bar) == 4){
      dgm <- structure(bar, names=c("dim", "id", "birth", "death"), class="dgm")
      gen_params <- list(pattern = "^Generator for ID ([0-9]+):(.*)", x = readLines(generator_file))
      h1_ids <- as.integer(do.call(gsub, append(gen_params, list(replacement="\\1"))))
      cycles <- lapply(strsplit(trimws(do.call(gsub, append(gen_params, list(replacement="\\2")))), split = " "), as.integer)
      names(cycles) <- h1_ids
      dgm[["generators"]] <- cycles
      return(dgm)
    } else if (ncol(bar) == 3){
      return(structure(bar, names=c("dim", "birth", "death"), class="dgm"))
    }
    warning("Incompatible format detected. Returning delimited file output.")
    return(bar)
  }
  warning("Persistence file does not exist. Something wrong may have occurred.")
  return(NULL)
}


## Given an initial simplicial complex 'x' and a set of elementary simplicial maps 'sm', this
## function returns a possibly modified set of elementary simplicial maps that ensures the faces
## to every inserted simplex are inserted prior to insert
validate_sm <- function(x, sm){
  st <- simplextree::simplex_tree()
  st$deserialize(x$serialize()) ## copy complex
  si_maps <- if (is.vector(sm)){ sm } else { readLines(sm) }

  ## Read in elementary simplicial operations, and validate them.
  validated_sm <- unlist(lapply(si_maps, function(si_map){
    op <- strsplit(si_map, split = " ")[[1]]
    if (op[1] == "i"){
      simplex <- sort(as.integer(op[-1]))
      if (length(simplex) > 1){
        faces <- lapply(2:length(simplex), function(i){
          combs <- combn(length(simplex), i-1)
          apply(combs, 2, function(face_idx){
            face <- simplex[face_idx]
            face_exists <- st$find(face)
            if (!face_exists){ sprintf("i %s", paste0(face, collapse = " ")) } else { NULL }
          })
        })
        st$insert(simplex)
        si_map <- sprintf("i %s", paste0(simplex, collapse = " "))
        return(unlist(c(faces, si_map)))
      } else {
        st$insert(simplex)
        return(sprintf("i %s", paste0(simplex, collapse = " ")))
      }
    } else if (op[1] == "c"){
      from <- sort(as.integer(op[2:(which(op == "t")-1)]))
      to <- as.integer(tail(op, 1))
      if (!all(st$find(as.list(c(from, to))))){ stop(sprintf("Invalid collapse detected to non-existent vertex: %s", op)) }
      invisible(sapply(from, function(v_i){ st$collapse(v_i, to, to) }))
      return(sprintf("c %s t %d", paste0(from, collapse = " "), to))
    } else {
      return(si_map)
    }
  }))
  return(validated_sm)
}


#' @method format dgm
#' @export
format.dgm <- function(x){
  sprintf("Barcodes for dims: %s\n", paste0(as.character(unique(x[["dim"]])), collapse = ", "))
  # sprintf("Barcodes for dims: %s\n", paste0(as.character(unique(x[["dim"]])), collapse = ", "))
}

#' @method print dgm
#' @export
print.dgm <- function(x){
  writeLines(format(x))
}

#' @method plot dgm
#' @export
plot.dgm <- function(x, show_legend=TRUE){
  stopifnot(all(c("dim", "birth", "death") %in% names(x)))
  ggtda_installed <- requireNamespace("ggtda", quietly = TRUE)
  if (!ggtda_installed){
    warning("Diagram plotting requires the 'ggtda' package to be installed.\n
          Would you like to install it?")
    response <- readline("Install 'ggtda'? (y/n)")
    if (head(toupper(response), 1)  == "Y"){
      devtools::install_github('rrrlw/ggtda', vignettes = TRUE)
      loadNamespace(package = "ggtda")
    } else { return(NULL) }
  }
  barcodes <- data.frame(appear=x$birth, disappear=x$death, dim=as.character(x$dim))
  ggplot2::ggplot(barcodes, ggplot2::aes(start = appear, end = disappear, colour = dim, shape = dim)) +
    ggtda::geom_barcode() +
    ggtda::theme_tda()
  # plot.new()
  # all_times <- c(x[["birth"]], x[["death"]])
  # x_rng <- c(min(all_times), max(all_times[all_times != Inf]))
  # y_rng <- c(1, length(x[["birth"]]))
  # plot.window(xlim = x_rng, ylim = y_rng)
  # axis(1, at=x_rng)
  # n <- length(x$dim)
  # bd_rle <- rle(x[["dim"]])
  # len_cs <- cumsum(bd_rle$lengths)
  # # cc <- 1L
  # for (i in seq(n)){
  #   bd <- c(x$birth[i], x$death[i])
  #   y_h <- (n - i)+1
  #   rect(xleft = bd[1], xright=ifelse(bd[2]==Inf, max(y_rng), bd[2]), ybottom=y_h+0.25, ytop=y_h+0.50,
  #        col = palette()[min(which(i <= len_cs))], border = NA)
  #   # if (i %in% len_cs){
  #   #   hc_y <- max(y_rng)-len_cs[cc]+bd_rle$lengths[cc]/2
  #   #   mtext(side = 2, at=c(hc_y), text = parse(text=sprintf("H[%d]", cc-1L)), adj = 0) #min(x_rng)-0.05*max(x_rng)
  #   #   lines(x_rng, y=rep(y_h-0.1, 2), col = "darkblue")
  #   #   cc <- cc + 1L
  #   # }
  # }
  # if (!is.null(x[["id"]])){ mtext(as.character(x$id), side = 2, at = ((n-seq(n))+1)+0.25, adj = 0, padj=0.5) }
  # if (show_legend){
  #   legend_classes <- sapply(unique(x$dim), function(d) parse(text=sprintf("H[%d]", d)))
  #   legend("topright", legend = legend_classes, col = palette()[1:3], lty = 1, lwd = 4, cex=0.50, y.intersp=0.80)
  # }
  # title(main="Persistence barcodes")
}


# TODO: Map the id generators for H0 to the initial vertices!
# st <- simplextree::simplex_tree()
# st$insert(list(4, 7, 9, 12, 15))
# wut <- Simpers::simpers(st, sm = c("i 4 7", "# 1", "i 9 12", "# 2", "i 9 15", "# 3"))
# do.call(cbind, wut[c("id", "birth", "death")])
#
# st <- simplextree::simplex_tree()
# st$insert(list(c(4, 7), 9, 12, 15))
# wut <- Simpers::simpers(st, sm = c("i 9 12", "# 2", "i 9 15", "# 3"))
# do.call(cbind, wut[c("id", "birth", "death")])
#
#
# st <- simplextree::simplex_tree()
# wut <- Simpers::simpers(st, sm = c("i 4", "i 7", "i 9", "i 12", "i 15", "# 0", "i 4 7", "# 1", "i 9 12", "# 2", "i 9 15", "# 3"))
# do.call(cbind, wut[c("id", "birth", "death")])
#
#
# st <- simplextree::simplex_tree()
# st$insert(1:3)
# wut <- Simpers::simpers(st, sm = c("# 0", "i 4", "# 1", "i 5", "# 2", "i 4 5", "# 3"))

# Simpers:::run_simpers(strsplit(cmd, " ")[[1]])
#
# ## Make sure to add initial space!
# cmd <- sprintf(" -m %s", file.path(sp_dir, "operations_empty_domain_complex.txt"))
# Simpers:::run_simpers(strsplit(cmd, " ")[[1]])
