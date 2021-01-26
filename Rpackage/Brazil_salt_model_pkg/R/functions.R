## Brazil Salt Policy model: a decision support tool for primary prevention of NCDs
## Copyright (C) 2019 Chris Kypridemos

## Brazil Sodium Policy model is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have a copy of the GNU General Public License
## along with this program; if not, see <http://www.gnu.org/licenses/>
## or write to the Free Software Foundation, Inc., 51 Franklin Street,
## Fifth Floor, Boston, MA 02110-1301  USA.

# R function to fill columns in a data.table with random numbers
#' @export
gen_rn_cols <-
  function(dt,
           col_names,
           lower_bound = 0,
           upper_bound = 1) {
    stopifnot(is.data.table(dt))
    col_names <- paste0(col_names)
    N <- nrow(dt)
    for (j in col_names) {
      set(dt, NULL, j, my_runif(N, lower_bound, upper_bound))
    }
    invisible(dt)
  }

# delete output files
#' @export
delete_output_files <- function(x = output_dir()) {
  file.remove(list.files(
    path = x,
    full.names = TRUE,
    recursive = TRUE,
    all.files = TRUE
  ))
}


# dependencies <- function(x) {
#   for (j in x) {
#     # require returns T invisibly if it was able to load package
#     if (!base::require(j, character.only = TRUE)) {
#       # If package was not able to be loaded then re-install
#       utils::install.packages(j, dependencies = TRUE)
#       # Load package after installing
#       base::require(j, character.only = T)
#     }
#   }
# }


# truncated normal sampler that doesn't require you to throw out observations.
# 3x faster than truncnorm package.
# but truncnorm::rtruncnorm(1, 0, Inf, mean=0, sd = 0) gives NaN while this Inf
# because qnorn(1, x, x) = Inf
# rtruncnorm <-
#   function(n, a = -Inf, b = Inf, mean, sd) {
#     stats::qnorm(stats::runif(n, stats::pnorm(a, mean, sd), stats::pnorm(b, mean, sd)), mean, sd)
#   }

# define function for design$stochastic RR
# stochRRabov1 <-
#   function(n = .N, m, ci) {
#     # lognormal
#     rr <-
#       exp(rtruncnorm(n, 0, Inf, log(m), abs(log(m) - log(ci)) / 1.96))
#     rr[!is.finite(rr)] <- 1 # fix for rtruncnorm above
#     return(rr)
#   }

#' @export
stochRRtabl <- # need to run by id
  function(m, ci, stochastic = TRUE) {
    # lognormal
  	if (stochastic)
  		kk <- stats::runif(1)
  	else
  		kk <- 0.5
    rr <- exp(stats::qnorm(kk, log(m), abs(log(m) - log(ci)) / 1.96))
    rr[!is.finite(rr)] <- 1 # fix for rtruncnorm above
    return(rr)
  }

# stochRRbelow1 <-
#   function(n = .N, m, ci) {
#     # lognormal
#     rr <-
#       exp(rtruncnorm(n, -Inf, 0, log(m), abs(log(m) - log(ci)) / 1.96))
#     rr[!is.finite(rr)] <- 1 # fix for rtruncnorm above
#     return(rr)
#   }

# stochRRnorm <-
#   function(n = .N, m, ci) {
#     # normal distr
#     if (m < 1) {
#       a = 0
#       b = 1
#     } else {
#       a = 1
#       b = Inf
#     }
#     ifelse(m == ci,
#            rr <- rep(m, n),
#            rr <-
#              rtruncnorm(
#                n = n,
#                a = a,
#                b = b,
#                mean = m,
#                sd = abs(m - ci) / 1.96
#              ))
#     rr[!is.finite(rr)] <- 1 # fix for rtruncnorm above
#     return(rr)
#   }

# stochCOSTtabl <- # need to run by id
# 	function(shape1, shape2, scale, stochastic = TRUE)
# 	{
# 		# betapr
# 		if (stochastic)
# 			kk <- stats::runif(1)
# 		else
# 			kk <- 0.5
# 		cost <- extraDistr::qbetapr(kk, shape1, shape2, scale)
# 		return(cost)
# 	}

# function to clone dt for 2dmc
# clone_dt <-
#   function(DT, times = design$iteration_n) {
#     xx <- key(DT)
#     l <- rep(list(DT), times)
#     out <- setkeyv(rbindlist(l, idcol = T), xx)
#     return(out)
#   }

# Define function for sampling. Taken from sample man pages
#' @export
resample <-
  function(x, ...) {
    x <- na.omit(x)
    x[sample.int(length(x), ...)]
  }

# Define function for linear diffusion of policy effects over time.
#' @export
linear_diffusion <-
  function(effect,
           actual_time,
           diffusion_duration,
           previous_effect = 0)
  {
    previous_effect - bound(actual_time, 0, diffusion_duration) *
      (previous_effect - effect) /
      diffusion_duration
  }

# Define operator %!in%, meaning !%in%
#' @export
'%!in%' <- function(x, y)
  !('%in%'(x, y))


#' @export
group_pop <- function(dt) {
	tt <- setDT(expand.grid("agegroup" = unique(agegroup_fn(-100:100)),
													"sex" = c("men", "women"),
													"race" = c("hispanic", "white", "black", "other"),
													KEEP.OUT.ATTRS = FALSE,
													stringsAsFactors = TRUE), key = c("agegroup", "sex", "race"))
	tt[, group := paste0(agegroup, sex, race)]
	dt[tt, on = c("agegroup", "sex", "race"), group := i.group]
}

# Define outersect. Like setdiff but symmetrical. I.e. setdiff(a,b) is not the
# same as setdiff(b,a). outersect solve this by calculating both
#' @export
outersect <-
    function(x, y, ...) {
      big.vec <- c(x, y, ...)
      duplicates <- big.vec[duplicated(big.vec)]
      setdiff(big.vec, unique(duplicates))
    }

#' @export
agegroup_fn <-
  function(x) {
    breaks <- c(seq(20, 85, 5), Inf)
    labels <- c(
      "20-24",
      "25-29",
      "30-34",
      "35-39",
      "40-44",
      "45-49",
      "50-54",
      "55-59",
      "60-64",
      "65-69",
      "70-74",
      "75-79",
      "80-84",
      "85+"
    )
    ff <- data.table(age = 20:100)
    ff[, agegroup := cut(
      age,
      breaks = breaks,
      labels = labels,
      include.lowest = T,
      right = F,
      ordered_result = T)]

    if (is.data.table(x)) {
      x[ff, agegroup := i.agegroup, on = "age"]
      group_pop(x)
      #return(invisible(x)) # If active it needs to assign dt <- foo(dt)
    } else if (is.numeric(x)) {
      agegroup = cut(
        x,
        breaks = breaks,
        labels = labels,
        include.lowest = T,
        right = F,
        ordered_result = T
      )
      return(invisible(agegroup))
    } else return(print("only datatables and vectors are eligible inputs"))
  }


# Define function to bound a vector (numeric). fbound implemented in C++, fbound_inplace
#' @export
bound <-
  function(x,
           a = 0,
           b = 1,
           inplace = FALSE) {
    if (any(length(x) != length(a), length(x) != length(b))) {
      if (all(length(a) == 1L, length(b) == 1L)) {
        a <- rep(a, length(x))
        b <- rep(b, length(x))
      } else if(all(length(x) == length (a), length(b) == 1L)) {
        b <- rep(b, length(x))
      } else if(all(length(x) == length (b), length(a) == 1L)) {
        a <- rep(a, length(x))
      } else {
        stop("x, a, and b have different lengths")
      }
    }

    if (inplace) {
      if (typeof(x) == "integer") {
        a <- as.integer(a)
        b <- as.integer(b)
        return(fbound_inplace_int(x, a, b))
      } else if (typeof(x) == "double") {
        a <- as.numeric(a)
        b <- as.numeric(b)
        return(fbound_inplace(x, a, b))
      }
    } else {
      # if not inplace
      if (typeof(x) == "integer") {
        a <- as.integer(a)
        b <- as.integer(b)
        return(fbound_int(x, a, b))
      } else if (typeof(x) == "double") {
        a <- as.numeric(a)
        b <- as.numeric(b)
        return(fbound(x, a, b))
      }
    }
  }

#' @export
identical_elements <-
  function(x, tol = .Machine$double.eps ^ 0.5) {
    stopifnot(is.numeric(x))
    fequal(x, tol)
  }


# normalise a vector to 0,1 range
#' @export
normalise <-
  function(x, ...) {
    stopifnot(is.numeric(x))
    if (identical_elements(x))
      return(1)
    else
      return(fnormalise(x))
  }

# Add jitter and redistribute randomly those outside the range 0, 1
# jitter.constr <-
#   function(x, fortune = 0.1) {
#     x <- x + stats::runif(length(x), -fortune, fortune)
#     return(perc.rank(x))
#   }

# convert a DF (or DT) to matrix
#' @export
df2mat <-
  function(x) {
    stopifnot(is.data.frame(x))
    col_type <- lapply(x, typeof)
    if (!all(col_type %in% c("integer", "double")))
    {
      stop("Only works for integer and numeric matrices")
    }
    if ("double" %in% col_type)
    {
      x <- df2mat_numeric(x)
    } else {
      x <- df2mat_integer(x)
    }
    return(invisible(x))
  }

# shift by id
#' @export
shift_byid <-
  function(x, lag, replace, id, levels = NULL) {
    if (class(x) == "integer") {
      return(shift_byidInt(x, lag, replace, id))
    } else if (class(x) == "factor") {
      namx <- levels(x)
      out <- factor_(shift_byidInt(x, lag, replace, id), levels)
      levels(out) <- namx
      return(out)
    } else if (class(x) == "logical") {
      out <- factor_(shift_byidInt(x, lag, replace, id), levels)
      setattr(out, "levels", c("FALSE", "TRUE"))
      return(out)
    } else if (class(x) == "numeric") {
      return(shift_byidNum(x, lag, replace, id))
    } else
      stop("class of x not supported")
  }

#' @export
fit.betapr <- # NOT VECTORISED
  function(q, p = c(0.1, 0.5, 0.9), fnscale = 1) {
    if (length(q) == 1L)
      q <- c(q * 0.8, q, q * 1.2)
    ofn <- function(x)
      sum((q - extraDistr::qbetapr(p, x[1], x[2], x[3])) ^ 2)
    osol <-
      stats::optim(
        c(1, 1, 1),
        ofn,
        method = "L-BFGS-B",
        control = list(
          "fnscale" = fnscale,
          "maxit"   = 1e6,
          "ndeps"   = rep(1e-3, 3)
        ),
        lower = c(1, 1, 1),
        upper = c(Inf, Inf, Inf)
      )
    out <- as.list(osol$par)
    out <- stats::setNames(out, c("shape1", "shape2", "scale"))
    return(out)
  }

# From rriskDistributions to MLE fit beta to quantiles
is.error <- function(x) inherits(x, "try-error")

#' @export
get.beta.par <-
  function(p = c(0.025, 0.5, 0.975),
           q,
           show.output = TRUE,
           plot = TRUE,
           tol = 0.001,
           fit.weights = rep(1, length(p)),
           scaleX = c(0.1,  0.9),
           ...)
  {
    if (!is.numeric(p) | !is.numeric(q) | !is.numeric(fit.weights)) {
      stop(
        "INVALID INPUT, not numerical items in the input vectors 'p', 'q' and/or 'fit.weights'!",
        call. = FALSE
      )
    }
    if (prod(order(p) == seq(1:length(p))) == 0 | prod(order(q) ==
                                                       seq(1:length(q))) == 0) {
      stop("INVALID INPUT, the vector of probabilities/percentiles is not ordered!",
           call. = FALSE)
    }
    if (min(p) < 0 | max(p) > 1) {
      stop("INVALID INPUT, items of the probability vector should lie between 0 and 1!",
           call. = FALSE)
    }
    if (min(q) < 0 | max(q) > 1) {
      stop(
        "INVALID INPUT, percentiles are out of the domain (0, 1) => beta distribution couldn't be fitted!",
        call. = FALSE
      )
    }
    if (length(p) != length(q) | length(p) != length(fit.weights) |
        length(q) != length(fit.weights)) {
      stop(
        "INVALID INPUT, 'p', 'q' and 'fit.weights' are not of the same length! The vectors of quantiles, probabilities and weightings should be of the same length.",
        call. = FALSE
      )
    }
    if (length(q) < 2) {
      stop("INVALID INPUT, at least two quantiles must be known!",
           call. = FALSE)
    }
    if (!is.logical(show.output)) {
      stop("INVALID INPUT, the argument 'show.output' should be logical!",
           call. = FALSE)
    }
    if (!is.logical(plot)) {
      stop("INVALID INPUT, the argument 'plot' should be logical!",
           call. = FALSE)
    }
    if (!is.numeric(tol) | length(tol) != 1 | tol < 0) {
      stop(
        "INVALID INPUT, the argument 'tol' should be a single positive numerical value!",
        call. = FALSE
      )
    }
    fit.weights.original <- fit.weights
    fit.weights <- fit.weights / sum(fit.weights)
    minimize <- function(shape) {
      summand <- suppressWarnings(stats::pbeta(
        q = q,
        shape1 = shape[1],
        shape2 = shape[2]
      ) - p)
      summand <- summand * fit.weights
      sum(summand ^ 2)
    }
    fit <- c()
    fit$value <- tol + 1
    try1 <- try(fit <- stats::optim(
      par = c(0.1, 0.1),
      minimize,
      method = "L-BFGS-B",
      lower = 0.001,
      upper = 10000
    ),
    silent = TRUE)
    if (is.error(try1) || fit$value >= tol) {
      warning(
        "The fitting procedure 'L-BFGS-B' has failed (convergence error occurred or specified tolerance not achieved)!",
        call. = FALSE
      )
      fit <- c()
      fit$value <- tol + 1
      try2 <- try(fit <- stats::optim(minimize, method = "CG"),
                  silent = TRUE)
      if (is.error(try2) || fit$value >= tol) {
        warning(
          "The fitting procedure 'CG' has failed (convergence error occurred or specified tolerance not achieved)!",
          call. = FALSE
        )
        Par <- NA
      }
      else if (fit$value < tol) {
        message(
          "The fitting procedure 'CG' was successful!\n(Used this fallback optimization method because 'L-BFGS-B' has failed...)"
        )
        Par <- fit$par
        names(Par) <- c("shape1", "shape2")
        if (show.output)
          print(fit)
      }
    }
    else if (fit$value < tol) {
      message("The fitting procedure 'L-BFGS-B' was successful!")
      Par <- fit$par
      names(Par) <- c("shape1", "shape2")
      if (show.output)
        print(fit)
    }
    if (prod(!is.na(Par)) & plot) {
      main1 <- paste("shape1 = ", round(Par["shape1"], digits = 2))
      main2 <- paste("shape2 = ", round(Par["shape2"], digits = 2))
      main <- paste("Beta (", main1, ", ", main2, ")", sep = "")
      sub = paste("fit.weights = c(",
                  paste(fit.weights.original,
                        collapse = ", "),
                  ")",
                  sep = "")
      Support.lim <- c(
        stats::qbeta(
          p = min(p) * scaleX[1],
          shape1 = Par["shape1"],
          shape2 = Par["shape2"]
        ),
        stats::qbeta(
          p = (max(p) + (1 - max(p)) * scaleX[2]),
          shape1 = Par["shape1"],
          shape2 = Par["shape2"]
        )
      )
      Support <- seq(min(min(q), Support.lim[1]), max(max(q),
                                                      Support.lim[2]), length = 200)
      Probability <-
        stats::pbeta(Support, Par["shape1"], Par["shape2"])
      graphics::plot(
        Support,
        Probability,
        type = "l",
        xlim = range(Support.lim,
                     q),
        main = main,
        xlab = "Quantiles",
        sub = sub,
        ...
      )
      graphics::points(x = q, y = p, pch = 19, ...)
    }
    return(Par)
  }

# estimate beta params from mean and variance
#' @export
estim_beta_params <- function(mu, var) {
  # from https://stats.stackexchange.com/questions/12232/calculating-the-parameters-of-a-beta-distribution-using-the-mean-and-variance and wikipedia
  stopifnot(between(mu, 0, 1), var > 0)
  if (var >= (mu * (1 - mu))) var  <- mu * (1 - mu) * 0.9999
  alpha <- mu * ((mu * (1 - mu) / var) - 1) # if var < (mu * (1 - mu))
  beta <- (1 - mu) * ((mu * (1 - mu) / var) - 1) # var < (mu * (1 - mu))
  return(params = list(shape1 = alpha, shape2 = beta))
}
#plot(density(do.call(rbeta, c(list(n=1e5), estim_beta_params(design$chd_incidence_ratio, 0.005), list(ncp = 0)))))


# calculate perc ranks
# perc_rank <- function(x, n = .N) {
#   (frank(x,
#          na.last = F,
#          ties.method = "random") - 1) / (n - 1)
# }

# as.data.table.array - #1418 not merged yet to DT
# copied from https://github.com/Rdatatable/data.table/issues/1418
# original code by jangorecki

.onUnload <- function(libpath) {
  library.dynam.unload("BrazilSaltModelmisc", libpath)
}

`:=` = function(...)
  NULL # due to NSE notes in R CMD check

#' @import utils
#' @import data.table
NULL

#' @useDynLib BrazilSaltModelmisc
#' @importFrom Rcpp sourceCpp evalCpp
#' @importFrom methods as
#' @importFrom graphics abline legend lines par plot title
#' @importFrom stats density .checkMFClasses delete.response model.frame model.matrix plogis predict

NULL

#' Get Dropbox path
#'
#' `get_dropbox_path` returns the path of Dropbox. Works for both personal and business accounts
#'
#' This is an auxilliary function: It finds the Dropbox path in Windows, Linux, and OSX operating systems.
#'
#' @param pathtail A String vector (if not a string then it is converted to String).
#'    If present, it gets concatenated with the Dropbox path.
#'    See examples.
#' @param type A String scalar ("personal" or "business"). Which Dropbox path to return? The personal or the business one? It may be abbreviated.
#' @return Dropbax path as a String. If pathtail is present, it concatenates the Dropbox path with pathtail.
#' @export
#' @examples
#' \dontrun{
#' # Only work if Dropbox is installed.
#' get_dropbox_path() # Returns personal Dropbox path
#' get_dropbox_path(type = "business") # Returns business Dropbox path
#' get_dropbox_path("pathdownthetree") # Returns "Dropbox_path/pathdownthetree",
#'  # where Dropbox_path is the path to personal Dropbox
#' }
get_dropbox_path <-
  function(pathtail = character(0),
           type = c("personal", "business")) {
    if (!requireNamespace("jsonlite", quietly = TRUE))
      stop("Please install package jsonlite first.")
    type <- match.arg(type)
    if (type[[1]] == "personal") {
      if (.Platform$OS.type == "windows") {
        if (file.exists(paste0(Sys.getenv("APPDATA"), "/Dropbox/info.json"))) {
          # for older versions of Dropbox
          dropbox_path <-
            jsonlite::read_json(paste0(Sys.getenv("APPDATA"), "/Dropbox/info.json"))$personal$path
        }
        if (file.exists(paste0(Sys.getenv("LOCALAPPDATA"), "/Dropbox/info.json"))) {
          dropbox_path <-
            jsonlite::read_json(paste0(Sys.getenv("LOCALAPPDATA"), "/Dropbox/info.json"))$personal$path
        }
      } else {
        if (file.exists("~/.dropbox/info.json"))
          dropbox_path <-
            jsonlite::read_json("~/.dropbox/info.json")$personal$path
      }
    }
    if (type[[1]] == "business") {
      if (.Platform$OS.type == "windows") {
        if (file.exists(paste0(Sys.getenv("APPDATA"), "/Dropbox/info.json"))) {
          # for older versions of Dropbox
          dropbox_path <-
            jsonlite::read_json(paste0(Sys.getenv("APPDATA"), "/Dropbox/info.json"))$business$path
        }
        if (file.exists(paste0(Sys.getenv("LOCALAPPDATA"), "/Dropbox/info.json"))) {
          dropbox_path <-
            jsonlite::read_json(paste0(Sys.getenv("LOCALAPPDATA"), "/Dropbox/info.json"))$business$path
        }
      } else {
        if (file.exists("~/.dropbox/info.json"))
          dropbox_path <-
            jsonlite::read_json("~/.dropbox/info.json")$business$path
      }
    }
    if (is.null(dropbox_path))
      stop("Dropbox path cannot be located.")
    dropbox_path <-
      normalizePath(paste0(dropbox_path, "/", pathtail), mustWork = FALSE)
    return(dropbox_path)
  }




#' Get pCloud path
#'
#' `get_pcloud_path` returns the path of pCloud
#'
#' This is an auxilliary function: It finds the pCloud path in Windows, Linux, and OSX operating systems.
#'
#' @param pathtail A String vector (if not a string then it is converted to String).
#'    If present, it gets concatenated with the pCloud path.
#'    See examples.
#' @return pCloud path as a String. If pathtail is present, it concatenates the Dropbox path with pathtail.
#' @export
#' @examples
#' \dontrun{
#' # Only work if Dropbox is installed.
#' get_pcloud_path() # Returns pCloud path
#' get_pcloud_path("pathdownthetree") # Returns "pcloud_path_path/pathdownthetree",
#' # where pcloud_path is the path to pCloud
#' }
get_pcloud_path <- function(pathtail = character(0)) {
  if (.Platform$OS.type == "windows") {
    pcloud_path <- "p:\\"
  } else
    pcloud_path <- "~/pCloudDrive/"

  pcloud_path <-
    normalizePath(paste0(pcloud_path, "/", pathtail), mustWork = FALSE)
  return(pcloud_path)
}



#' Generate names for age-group bands
#'
#' `agegrp_name` generates names for age-group bands given lower and upper
#'   age limits, and band width
#'
#'
#' @param min_age A non-negative integer. The lower age limit for which
#'   names will be generated.
#' @param max_age A non-negative integer. The upper age limit for which
#'   names will be generated.
#' @param grp_width A positive integer. The band width of the age-groups.
#' @param grp_lessthan_1 A logical scalar. if \code{TRUE} and
#'  \code{min_age == 0}, then the first age-group name is "<1".
#' @param match_input A logical scalar. If \code{TRUE}, then the names
#'  are repeated to match every single year of age between \code{min_age}
#'  and \code{match_input_max_age}.
#' @param match_input_max_age a non-negative integer. See above.
#' @return A character vector of with the names for the age-groups.
#' @export
#' @examples
#' agegrp_name(20, 79, 5)
#' agegrp_name(20, 80, 5)
#' agegrp_name(0, 80, 10, TRUE)
#' agegrp_name(20, 30, 5, FALSE, TRUE)
#' agegrp_name(20, 30, 5, FALSE, TRUE)
#' agegrp_name(20, 30, 5, FALSE, TRUE, 32)
agegrp_name <-
  function(min_age = 0L,
           max_age = 85L,
           grp_width = 5L,
           grp_lessthan_1 = TRUE,
           match_input = FALSE,
           match_input_max_age = max_age) {
    stopifnot(
      min_age >= 0,
      max_age > 0,
      grp_width >= 1,
      max_age > min_age,
      match_input_max_age >= max_age
    )
    if (grp_width > 1) {
      x <- seq(min_age, max_age + 1L, grp_width)
      y <- shift(x, type = "lead") - 1L
      out <- paste0(sprintf("%02.0f", x), "-", sprintf("%02.0f", y))
      if ((tail(x, 1) - 1L) != max_age) {
        out[length(out)] <- paste0(x[length(x)], "+")
      } else {
        out <- head(out, length(out) - 1L)
      }
      if (grp_lessthan_1 && min_age == 0) {
        out <- c("<1", out)
        out[2] <- paste0("01", "-", sprintf("%02.0f", y[1]))
      }

      if (grp_lessthan_1 && min_age == 0) {
        if (match_input && ((tail(x, 1) - 1L) != max_age)) {
          out <- c(rep(out[2:((length(out) - 1L))], each = grp_width),
                   rep(out[length(out)], match_input_max_age - tail(x, 1) + 1L))
          out[1] <- "<1"
        }
        if (match_input && ((tail(x, 1) - 1L) == max_age)) {
          out <- rep(out[2:length(out)], each = grp_width)
          out[1] <- "<1"
        }
      } else {
        if (match_input && ((tail(x, 1) - 1L) != max_age)) {
          out <- c(rep(out[1:((length(out) - 1L))], each = grp_width),
                   rep(out[length(out)], match_input_max_age - tail(x, 1) + 1L))
        }
        if (match_input && ((tail(x, 1) - 1L) == max_age)) {
          out <- rep(out, each = grp_width)
        }
      }
    } else {
      out <- paste0(sprintf("%02.0f", seq(min_age, max_age, grp_width)))
    }
    return(out)
  }


#' Replace multiple values in a data.table column
#'
#' `replace_from_table` replace multiple values in a data.table column.
#'  . The values in \code{from} arguement are matched
#'    and replaced by those in \code{to} arguement.
#'    If \code{newcolname = NULL} the replace is by reference.
#'
#' @param dt A data.table to be changed by reference.
#' @param colname A string denoting the name of the column to be changed.
#' @param from A vector with values in \code{colname} to be replaced.
#' @param to A vector with values in \code{colname} to be replaced.
#' @param newcolname A string denoting the name of a new column
#'    to be created. If present, \code{colname} is not altered.
#'    If \code{newcolname = NULL}, \code{colname} is altered by reference
#' @return a data.table, invisibly.
#' @export
#' @examples
#' library(data.table)
#' library(BrazilSaltModelmisc)
#' dt <- data.table::data.table("a" = 1:5, "b" = seq(1, 2.2, 0.3),
#'  "d" = letters[1:5])
#' dt[, e := factor(a, labels = LETTERS[1:5])]
#' replace_from_table(data.table::copy(dt), "a", 1:3, 3L)[]
#' replace_from_table(data.table::copy(dt), "a", 3L, -11L)[]
#' replace_from_table(data.table::copy(dt), "a", 3L, -11L, "newcol")[]
#' replace_from_table(data.table::copy(dt), "b", 1.3, "a")[]
#' replace_from_table(data.table::copy(dt), "b", 1.3, "a", "newcol")[]
#' replace_from_table(data.table::copy(dt), "d", "a", "7")[]
#' replace_from_table(data.table::copy(dt), "d", "a", 7)[]
#' replace_from_table(data.table::copy(dt), "e", "B", "J")[]
replace_from_table <-
  function(dt,
           colname,
           from,
           to,
           newcolname = NULL) {
    old_ <- i.new_ <- NULL
    stopifnot(is.data.table(dt))
    stopifnot(length(colname) == 1L)
    stopifnot(is.null(newcolname) | length(newcolname) == 1L)
    stopifnot(colname %in% names(dt))
    stopifnot(length(from) >= length(to))
    # stopifnot(class(from) == dt[, class(get(colname))]) # not working for factors
    if (!is.null(newcolname) && newcolname %in% names(dt)) stop(
      "The new column name already exists in the data.table.")
    if (length(from) > length(to)) message("Note: many to few match.")

    colorder <- copy(names(dt))
    if (class(from) == class(to)) {
      reg <- data.table("old_" = from, "new_" = to)
      dt[, "old_" := get(colname)]
      dt[reg, on = "old_", old_ := i.new_]
      if (is.null(newcolname)) {
        dt[, (colname) := NULL]
        setnames(dt, "old_", colname)
        setcolorder(dt, colorder)
      } else {
        setnames(dt, "old_", newcolname)
      }
    } else {
      reg <- data.table("old_" = as(from, class(to)),
                        "new_" = to)
      dt[, "old_" := as(get(colname), class(to))]
      dt[reg, on = "old_", old_ := i.new_]
      if (is.null(newcolname)) {
        message(paste0(
          colname,
          " coerced to ",
          class(to),
          " to match target class."
        ))
        dt[, (colname) := NULL]
        setnames(dt, "old_", colname)
        setcolorder(dt, colorder)
      } else {
        setnames(dt, "old_", newcolname)
      }
    }
    return(invisible(dt))
  }



#' Generate age-group from age
#'
#' `to_agegrp` creates a new column
#'
#' @param dt A data.table with a column named \code{age}.
#' @param age_colname A string denoting the age column in \code{dt}.
#' @param colname A string denoting the name of the column that will be
#'   created for age-groups.
#' @param to_factor A logical. If \code{TRUE}, then the age-groups
#'   column is converted to factor.
#' @param ... Pass arguements to \code{\link{agegrp_name}}.
#' @return a data.table, invisibly.
#' @export
#' @examples
#' library(data.table)
#' library(BrazilSaltModelmisc)
#' to_agegrp(data.table(age = 0:99))[]
#' to_agegrp(data.table(age = 0:99), max_age = 80L)[]
#' to_agegrp(data.table(age = 0:99), grp_width = 10, max_age = 85)[]
to_agegrp <-
  function(dt,
           age_colname = "age",
           agegroup_colname = "agegroup",
           to_factor = TRUE,
           ...) {
    stopifnot(is.data.table(dt), age_colname %in% names(dt),
              length(age_colname) == 1L, length(agegroup_colname) == 1L,
              is.logical(to_factor))

    age_vec <- dt[, min(get(age_colname))]:dt[, max(get(age_colname))]
    replace_from_table(
      dt,
      colname = age_colname,
      from = age_vec,
      to = agegrp_name(
        min_age = min(age_vec),
        match_input = TRUE,
        match_input_max_age = max(age_vec),
        ...
      ),
      newcolname = agegroup_colname
    )
    if (to_factor) {
      dt[, (agegroup_colname) := factor(get(agegroup_colname))]
    }
    return(invisible(dt))
  }

# TODO add documentation
#' Clone a data.table
#'
#' `clone_dt` clones a data.table and binds the copies at the bottom of
#' the original data.table. It also creates an column named \code{`.id`}
#' to identify each iteration. The keys of the input data.table is retained.
#'
#' @export
clone_dt <-
  function(dt, times, idcol = TRUE) {
    xx <- key(dt)
    l <- rep(list(dt), times)
    out <- setkeyv(rbindlist(l, idcol = idcol), xx)
    return(invisible(out))
  }

# TODO add documentation
#' Calculate percentile rank
#'
#' `pctl_rank` calculates the percentile rank of a numeric vector
#'
#' @export
pctl_rank <- function(x, ties.method = c("average", "first", "random",
                                         "max", "min", "dense")) {
  stopifnot(is.numeric(x))
  ties.method <- match.arg(ties.method)
  n   <- length(x)
  out <- (frank(x,
                na.last = F,
                ties.method = ties.method) - 1) / (n - 1)
  return(out)
}


# TODO add documentation
#' Stochastic prediction from a gamlss object
#'
#' `validate_gamlss` returns a data.table with the observed and predicted
#'  variable. If \code{mc > 1} multiple predictions are drawn from the predicted
#'  distributions. Useful for plotting with ggplot
#'
#' @export
validate_gamlss <- function(dt, gamlss_obj, mc = 10L, orig_data = dt) {
  if (!requireNamespace("gamlss", quietly = TRUE))
    stop("Please install package gamlss first.")
  stopifnot("gamlss" %in% class(gamlss_obj), is.data.table(dt), mc >= 1,
            is.data.table(orig_data))
  nam_y <- as.character(gamlss_obj$call$formula[[2]])
  nam_var <- all.vars(gamlss_obj$call$formula[[3]])
  nam_dist <- paste0("r", gamlss_obj$family[[1]])
  nam_param <- gamlss_obj$parameters
  x <- copy(dt)
  x[, type := "Observed"]
  z <- copy(dt)
  z[, (nam_param) := gamlss::predictAll(gamlss_obj, type = "response",
                                        newdata = dt[, .SD, .SDcols = nam_var],
                                        data = orig_data[, .SD,
                                                         .SDcols = c(nam_var)])]
  z[, type := "Modelled"]
  z <- rbindlist(rep(list(z), mc))
  z[, (nam_y) := do.call(nam_dist, c(.N, .SD)), , .SDcols = nam_param]
  out <- rbind(x, z, use.names = TRUE, fill = TRUE)
  out[, (nam_param) := NULL]
}

# TODO add documentation
# If I name the function predict_gamlss BrazilSaltModelmisc:: is necessary when I
# call the function.
#' Prediction from a gamlss object in parallel
#'
#' `guess_gamlss` returns a data.table with the predicted
#'  variable. `dt` needs to have a column with percentiles named `rank_y`,
#'  where `y` the name of the predicted variable (i.e. bmi).
#'
#' @export
guess_gamlss <- function(dt, gamlss_obj, orig_data = gamlss_obj$data, nc = 1L) {
  if (!requireNamespace("gamlss", quietly = TRUE))
    stop("Please install package gamlss first.")
  stopifnot("gamlss" %in% class(gamlss_obj),
            is.data.table(dt), nc >= 1L,
            is.data.table(orig_data))
  nam_y <- as.character(gamlss_obj$call$formula[[2]])
  nam_var <- all.vars(gamlss_obj$call$formula[[3]])
  nam_dist <- paste0("q", gamlss_obj$family[[1]])
  nam_param <- gamlss_obj$parameters

  orig_data <- orig_data[, ..nam_var]
  dtu <- unique(dt[, ..nam_var]) # otherwise too slow
  dtu <- split(dtu, dtu$year)
  if("RevoUtilsMath" %in% (.packages())) tt <- getMKLthreads()
  if("RevoUtilsMath" %in% (.packages())) setMKLthreads(1L)
  dtu <- parallel::mclapply(dtu, function(x) {
    x[, (nam_param) := gamlss::predictAll(gamlss_obj,
                                          type = "response",
                                          newdata = .SD,
                                          data = orig_data)]
  },
  mc.preschedule = FALSE,
  mc.cores = nc)
  if("RevoUtilsMath" %in% (.packages())) setMKLthreads(tt)
  dtu <- rbindlist(dtu)
  # dtu[, (nam_param) := gamlss::predictAll(gamlss_obj,
  #                                        type = "response",
  #                                        newdata = .SD,
  #                                        data = orig_data)]
  dt[dtu, on = nam_var, (nam_param) := mget(paste0("i.", nam_param))]
  dt[, p := get(paste0("rank_", nam_y))]
  stopifnot(dt[, all(between(p, 0, 1, incbounds = FALSE))])
  dt[, (nam_y) := do.call(nam_dist, .SD), .SDcols = c("p", nam_param)]
  dt[, c("p", nam_param) := NULL]
}


# TODO add documentation
#' Prediction from a MASS:polr object in parallel
#'
#' `guess_polr` returns a data.table with the predicted
#'  variable. `dt` needs to have a column with percentiles named `rank_y`,
#'  where `y` the name of the predicted variable (i.e. active_days).
#'
#' @export
guess_polr <- function(dt, polr_obj) {
  if (!requireNamespace("MASS", quietly = TRUE))
    stop("Please install package MASS first.")
  if (!requireNamespace("matrixStats", quietly = TRUE))
    stop("Please install package matrixStats first.")
  stopifnot("polr" %in% class(polr_obj), is.data.table(dt))
  nam_y <- as.character(polr_obj$call$formula[[2]])
  nam_var <- all.vars(polr_obj$call$formula[[3]])
  #code adapted from method getAnywhere(predict.polr)
  Terms <- delete.response(polr_obj$terms)
  m <- model.frame(Terms, dt[, ..nam_var], na.action = function(x) x,
                   xlev = polr_obj$xlevels)
  if (!is.null(cl <- attr(Terms, "dataClasses")))
    .checkMFClasses(cl, m)
  X <- model.matrix(Terms, m, contrasts = polr_obj$contrasts)
  xint <- match("(Intercept)", colnames(X), nomatch = 0L)
  if (xint > 0L)
    X <- X[, -xint, drop = FALSE]
  n <- nrow(X)
  q <- length(polr_obj$zeta)
  eta <- drop(X %*% polr_obj$coefficients)
  cc <- plogis(matrix(polr_obj$zeta, n, q, byrow = TRUE) -
                 eta)
  dt[, p := get(paste0("rank_", nam_y))]
  dt[, (nam_y) := matrixStats::rowSums2(cc < p)]
  dt[, "p" := NULL]
}


# TODO add documentation
#' Deterministic prediction from a gamlss object
#'
#' `crossval_gamlss` returns the observed and predicted values of the dependent
#'  variable. Useful for cross-validation metrics.
#'
#' @export
crossval_gamlss <- function(dt, gamlss_obj, orig_data = dt, colnam = "rank") {
  stopifnot("gamlss" %in% class(gamlss_obj), is.data.table(dt),
            is.data.table(orig_data))
  out <- list()
  nam_y <- as.character(gamlss_obj$call$formula[[2]])
  nam_var <- all.vars(gamlss_obj$call$formula[[3]])
  nam_dist <- paste0("q", gamlss_obj$family[[1]])
  nam_param <- gamlss_obj$parameters
  out$observed <- dt[, get(nam_y)]
  z <- copy(dt)
  z[, (nam_param) := predictAll(gamlss_obj, type = "response",
                                newdata = dt[, .SD, .SDcols = nam_var],
                                data = orig_data[, .SD,
                                                 .SDcols = c(nam_y, nam_var)])]
  setnames(z, colnam, "p")
  z[p == 0, p := 0.0001]
  z[p == 1, p := 0.9999]
  z[, (nam_y) := do.call(nam_dist, .SD), .SDcols = c("p", nam_param)]
  out$predicted <- z[, get(nam_y)]
  return(out)
}

##' Generate Counts of Values in a Vector
##'
##' This function uses Rcpp sugar to implement a fast \code{table}, for
##' unique counts of a single vector. This implementation seeks to
##' produce identical output to \code{table(x, useNA="ifany")}. It is borrowed
##' from \code{Kmisc} package for convenience, since \code{Kmisc} is not in CRAN
##'  anymore. \code{Kmisc} is available at https://github.com/kevinushey/Kmisc

##'
##' The order of \code{NA}, \code{NaN} in the output may differ -- even
##' \R is inconsistent with the order that \code{NA} and \code{NaN} elements
##' are inserted.
##'
##' @param x A numeric, integer, character or logical vector, or a (potentially
##'   nested) list of such vectors. If \code{x} is a list, we recursively apply
##'   \code{counts} throughout elements in the list.
##' @export
##' @examples
##' x <- round( rnorm(1E2), 1 )
##' counts(x)
counts <- function(x) {
  if (is.list(x)) {
    output <- rapply(x, counts, how="list")
    return(output)
  } else {
    return(.Call('_BrazilSaltModelmisc_counts', x))
  }
}

# TODO add documentation
#' Obtain matching names corresponding to patterns
#'
#' `match_colnames_pattern` returns the matching names of the argument `dt`
#' (i.e. \code{names(dt)}) corresponding to the regular expression patterns
#' provided. The patterns must be supported by \code{\link{grep}}.
#' This is based on `data.table:::patterns`
#'
#' @export
match_colnames_pattern <- function(dt, ...) {
  p = unlist(list(...), use.names = FALSE)
  if (!is.character(p)) stop("Input patterns must be of type character.")
  cols = names(dt)
  cols[unlist(sapply(p, grep, cols))]
}

# TODO add documentation
#' Compare two distributions
#'
#' Summary statistics for the location/shape decomposition of the relative
#' distribution of the exposure: Modelled to Observed."
#'
#' @export
reldist_diagnostics <- function(comparison, reference, comparison_wt, reference_wt,
                                main, smooth = 0.35, discrete = FALSE) {
  if (!requireNamespace("reldist", quietly = TRUE)) {
    stop("Package \"reldist\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  opar <- par(no.readonly = TRUE)
  on.exit(par(opar))
  reference_dens  <- density(reference, weights = reference_wt)
  comparison_dens <- density(comparison, weights = comparison_wt)

  par(mfrow=c(2,2))
  plot(
    reference_dens,
    main = main,
    lty = 3,
    ylim = c(0, 1.1 * max(c(
      reference_dens$y, comparison_dens$y
    )))
  )
  lines(comparison_dens, col = "red", lty = 2)
  legend("topright", bg="transparent" ,
         legend=c("Comparison", "Reference"), box.lty = 0,
         col = c("red", "black"), lty = 2:3, cex = 0.8)
  g10 <- reldist::reldist(y=comparison, yo=reference,
                          smooth=smooth, ci=TRUE,
                          ywgt=comparison_wt, yowgt=reference_wt,
                          #ylim=c(0,1.5),
                          #yolabs=seq(-1,3,by=0.5),
                          bar="yes", quiet=FALSE, discrete = discrete,
                          xlab="proportion of the reference cohort")
  title(main="Overall relative density",cex=0.6)
  abline(h=1,lty=2)
  g1A <- reldist::reldist(y=comparison, yo=reference,
                          ywgt=comparison_wt, yowgt=reference_wt,
                          show="effect",
                          bar="yes", quiet=FALSE,
                          smooth=smooth, ci=TRUE, discrete = discrete,
                          #ylim=c(0,1.5),
                          #yolabs=seq(-1,3,by=0.5),
                          xlab="proportion of the reference cohort")
  title(main= "Effect of the median shift",cex=0.6)
  abline(h=1,lty=2)
  gA0 <- reldist::reldist(y=comparison, yo=reference,
                          smooth=smooth, ci=TRUE,
                          ywgt=comparison_wt, yowgt=reference_wt,
                          show="residual",
                          bar="yes", quiet=FALSE, discrete = discrete,
                          #ylim=c(0,1.5),
                          #yolabs=seq(-1,3,by=0.5),
                          xlab="proportion of the reference cohort")
  title(main="Median-adjusted relative density" ,cex=0.6)
  abline(h=1,lty=2)

  a1 <- reldist::rpy(y=comparison,yo=reference,
                     ywgt=comparison_wt,yowgt=reference_wt,pvalue=TRUE)
  a2 <- reldist::rpluy(y=comparison,yo=reference,
                       ywgt=comparison_wt,yowgt=reference_wt,pvalue=TRUE)
  a3 <- reldist::rpluy(y=comparison,yo=reference,
                       ywgt=comparison_wt,yowgt=reference_wt,pvalue=TRUE,
                       upper=TRUE)
  # p1 <- ifelse(a1[[4]]<0.001, "<0.001", format(a1[[4]], digits = 3))
  # p2 <- ifelse(a2[[4]]<0.001, "<0.001", format(a2[[4]], digits = 3))
  # p3 <- ifelse(a3[[4]]<0.001, "<0.001", format(a3[[4]], digits = 3))

  out <- data.table("Summary statistics" = c("Overall change entropy",
                                             "Median effect entropy",
                                             "Shape effect entropy",
                                             "Median polarization index",
                                             "Lower polarization index",
                                             "Upper polarization index"
  ),
  "Measure" = c(g10$entropy,
                reldist::entropy(g1A,g10),
                gA0$entropy,
                a1[[2]],
                a2[[2]],
                a3[[2]]
  ),
  "Lower 95% CI" = c(NA, NA, NA,
                     a1[[1]],
                     a2[[1]],
                     a3[[1]]
  ),
  "Upper 95% CI" = c(NA, NA, NA,
                     a1[[3]],
                     a2[[3]],
                     a3[[3]]
  ),
  "p-value" = c(NA, NA, NA,
                a1[[4]],
                a2[[4]],
                a3[[4]]
  )
  )
  return(out[])
}


#' Simplified loading and installing of packages
#'
#' This is a wrapper to \code{\link{require}} and \code{\link{install.packages}}.
#' Specifically, this will first try to load the package(s) and if not found
#' it will install then load and attach the packages. Additionally, if the
#' \code{update=TRUE} parameter is specified it will check the currently
#' installed package version with what is available on CRAN (or mirror) and
#' install the newer version.
#'
#' The function was originally created by Jason Bryer
#' \href{https://www.r-bloggers.com/function-to-simplify-loading-and-installing-packages/}{here}
#' and the source is available \href{https://gist.github.com/jbryer/9112634}{here}.
#' Note: I renamed the function to \code{dependencies} and adapted it to attach
#' instead of only load the packages.
#'
#' @param pkges a character vector with the names of the packages to load.
#' @param install if TRUE (default), any packages not already installed will be.
#' @param update if TRUE, this function will install a newer version of the
#'        package if available.
#' @param quiet if TRUE (default), package startup messages will be suppressed.
#' @param verbose if TRUE (default), diagnostic messages will be printed.
#' @param ... other parameters passed to \code{\link{require}},
#'            \code{\link{install.packages}}, and
#'            \code{\link{available.packages}}.
#' @return a data frame with four columns and rownames corresponding to the
#'         packages to be loaded. The four columns are: loaded (logical
#'         indicating whether the package was successfully loaded), installed
#'         (logical indicating that the package was installed or updated),
#'         loaded.version (the version string of the installed package), and
#'         available.version (the version string of the package currently
#'         available on CRAN). Note that this only reflects packages listed in
#'         the \code{pkges} parameter. Other packages may be loaded and/or
#'         installed as necessary by \code{install.packages} and \code{require}.
#'         If \code{verbose=FALSE} the data frame will be returned using
#'         \code{\link{invisible}}.
#' @export
#' @examples
#' \dontrun{
#' dependencies(c('devtools','lattice','ggplot2','psych'))
#' }
dependencies <-
  function(pkges,
           install = TRUE,
           update  = FALSE,
           quiet   = TRUE,
           verbose = TRUE,
           ...) {
    myrequire <- function(package, ...) {
      result <- FALSE
      if (quiet) {
        suppressMessages(suppressWarnings(result <- requireNamespace(package, ...)))
      } else {
        result <- suppressWarnings(requireNamespace(package, ...))
      }
      return(result)
    }
    mymessage <- function(msg) {
      if (verbose) {
        message(msg)
      }
    }

    installedpkgs <- installed.packages()
    availpkgs <- available.packages()[, c('Package', 'Version')]
    if (nrow(availpkgs) == 0) {
      warning(
        paste0(
          'There appear to be no packages available from the ',
          'repositories. Perhaps you are not connected to the ',
          'Internet?'
        )
      )
    }
    # It appears that hyphens (-) will be replaced with dots (.) in version
    # numbers by the packageVersion function
    availpkgs[, 'Version'] <- gsub('-', '.', availpkgs[, 'Version'])
    results <- data.frame(
      loaded = rep(FALSE, length(pkges)),
      installed = rep(FALSE, length(pkges)),
      loaded.version = rep(as.character(NA), length(pkges)),
      available.version = rep(as.character(NA), length(pkges)),
      stringsAsFactors = FALSE
    )
    row.names(results) <- pkges
    for (i in pkges) {
      loadedPkgs <- search()
      needInstall <- FALSE
      if (i %in% row.names(installedpkgs)) {
        v <- as.character(packageVersion(i))
        if (i %in% row.names(availpkgs)) {
          if (v != availpkgs[i, 'Version']) {
            if (!update) {
              mymessage(
                paste0(
                  'A different version of ',
                  i,
                  ' is available ',
                  '(current=',
                  v,
                  '; available=',
                  availpkgs[i, 'Version'],
                  ')'
                )
              )
            }
            needInstall <- update
          }
          results[i, ]$available.version <- availpkgs[i, 'Version']
        } else {
          mymessage(paste0(i, ' is not available on the repositories.'))
        }
      } else {
        if (i %in% row.names(availpkgs)) {
          needInstall <- TRUE & install
          results[i, ]$available.version <- availpkgs[i, 'Version']
        } else {
          warning(paste0(
            i,
            ' is not available on the repositories and ',
            'is not installed locally'
          ))
        }
      }
      if (needInstall | !myrequire(i)) {
        install.packages(pkgs = i, quiet = quiet)
        if (!myrequire(i, ...)) {
          warning(paste0('Error loading package: ', i))
        } else {
          results[i, ]$installed <- TRUE
          results[i, ]$loaded <- TRUE
          results[i, ]$loaded.version <- as.character(packageVersion(i))
        }
      } else {
        results[i, ]$loaded <- TRUE
        results[i, ]$loaded.version <- as.character(packageVersion(i))
      }
      loadedPkgs2 <- search()
      for (j in loadedPkgs2[!loadedPkgs2 %in% loadedPkgs]) {
        try(detach(j, character.only = TRUE), silent = TRUE)
      }
      library(i, character.only = TRUE)
    }
    # library(pkges, character.only	= TRUE)
    if (verbose) {
      return(results)
    } else {
      invisible(results)
    }
  }

# TODO add documentation
#' Scrambles rank trajectories of simulants
#'
#' Scrambles the rank trajectories using a continuous space random walk. The \code{jump} parameter defines the maximum distance of jump every year
#'
#' @export
scramble_trajectories <- function(x, pid, jump = 0.05) {
  if (all(x < 0 | x > 1))
    stop("Input needs to be between 0 and 1")
  if (is.unsorted(pid))
    stop("IDs must be sorted")
  if (jump >= 1)
    stop("Overlap needs to be <= 1")
  if (jump == 0) return(x) else return(fscramble_trajectories(x, pid, jump))
}


# tt <- data.table(x = runif(5e5), pid = 1:5e5)
# tt <- clone_dt(tt, 50)
# setkey(tt, pid, .id)
# tt[, y := scramble_trajectories(x, pid, 0.05)]
# tt[998:1003]
# print(tt[sample(.N, 1e4), ggplot2::qplot(x, y, alpha = I(1/20))])
# print(tt[.id == 50, hist(y)])
# print(tt[.id == 50 & y > 0.9, hist(y)])
# print(tt[.id == 50 & y < 0.1, hist(y)])
# print(tt[pid == 4,  plot(.id, y, ylim = c(0,1))])



# Ensures that when fwrites appent file colnames of file to be written, match those already in the file
#' @export
fwrite_safe <- function(x,
                        file = "",
                        append = TRUE,
                        ...) {
  if (append) {
    if (file.exists(file)) {
      col_names_disk <- names(fread(file, nrows = 0))
      col_names_file <- names(x)
      col_names <- outersect(col_names_disk, col_names_file)
      if (length(col_names) > 0)
        x[, (col_names) := NA]
      setcolorder(x, col_names_disk)
    }
  }
  fwrite(x, file, append, ...)
}

# from Mozaffarian NEJM
#' @export
salt_sbp_effect <-
  # if reduction in sodium use (-)
  function(salt_difference,
           age,
           sbp,
           black_race,
           n,
           stochastic) {
    if (stochastic) {
      # Y = the BP change in each trial in mm Hg, standardized to a sodium reduction of 100 mmol/d;
      y <-
        stats::rnorm(n, -3.735113, 0.7303861) +
        stats::rnorm(n, -0.1052782, 0.0294371) * (age - 50) +
        stats::rnorm(n, -1.873587, 0.8841412) * (sbp > 140) +
        stats::rnorm(n, -2.489173, 1.188258) * black_race
      y <-
        bound(-salt_difference * y / 5.85, Inf * sign(salt_difference), 0)
      # -salt_difference * y / 2300 instead of 5.85 to convert mmol of sodium to mg
    } else {
      y <-
        -3.735113 - 0.1052782 * (age - 50) - 1.873587 * (sbp > 140) - 2.489173 * black_race
      y <- -salt_difference * y / 5.85
      # -salt_difference * y / 2300 instead of 5.85 to convert mmol of sodium to mg
    }
    return(y)
  }

#' @export
prune_pop <- function(dt, scenario_name, design) {
  # generate col names
  col_nam <- c("pid", "new_simulant",
               design$strata_for_outputs,
               grep(
                 paste0("_", scenario_name, "$"),
                 grep(
                   "^prb_",
                   grep("^rn_", names(dt), value = TRUE, invert = TRUE),
                   value = TRUE,
                   invert = TRUE
                 ),
                 value = TRUE
               ))

  # keep cols/ rows of interest
  out <-
    dt[year >= design$init_year &
         between(age, design$ageL, design$ageH), .SD, .SDcols = (col_nam)]

  # delete deads of more than a year
  setnames(out, gsub(paste0("_", scenario_name), "", names(out)))
  setkey(out, pid, year) # order by pid, year
  out[, new_simulant := mk_new_simulant_markers(pid)]
  out[, del := identify_longdeads(all_cause_mrtl, new_simulant)]
  fn_env <- environment() # get environment in fn
  del_dt_rows(out, out$del, fn_env)
  out[, del := NULL]

  # bound prevalence based on design
  # TODO: Mutates .SD which is locked by default
  out[, lapply(.SD, bound, 0L, design$max_prvl_for_outputs[[1]], TRUE),
      .SDcols = patterns("_prvl$")]
  invisible(out)
}


#' @export
summarise_hlp <- function(dt, scenario_name, design) {
  strata_nam <- design$strata_for_outputs
  out <- vector("list", 2L)

  #create cross-join of strata cols
  out[[1]] <- unique(dt[, .SD, .SDcols = strata_nam])

  #expand for each duration year
  out[[1]] <-
    rbindlist(rep(list(out[[1]]), 1L + design$max_prvl_for_outputs), idcol = "duration")
  out[[1]][, duration := duration - 1L]

  # summarise population size (only works after prune_pop)
    out[[2]] <- dt[, .N, by = strata_nam]
    setnames(out[[2]], "N", "pop_size")
    out[[2]][, duration := 0L]
    out[[1]] <-
      merge(out[[1]],
            out[[2]],
            by = c(strata_nam, "duration"),
            all = TRUE)

  # summarise morbidity
  col_nam <- grep("_prvl$", names(dt), value = TRUE)

  for (j in col_nam) {
    tt <- dt[get(j) > 0, .SD, .SDcols = c(j, strata_nam)]
    setnames(tt, j, "duration")

    out[[2]] <- tt[, .N, by = c(strata_nam, "duration")]
    setnames(out[[2]], "N", j)
    out[[1]] <-
      merge(out[[1]],
            out[[2]],
            by = c(strata_nam, "duration"),
            all = TRUE)
  }

  # summarise mortality (always recorded as duration 0)
  col_nam <- setdiff(grep("_mrtl$", names(dt), value = TRUE), "all_cause_mrtl")

  for (j in col_nam) {
    out[[2]] <- dt[all_cause_mrtl == setdiff(unique(get(j)), 0L), .N, by = strata_nam]
    setnames(out[[2]], "N", j)
    out[[2]][, duration := 0L]
    out[[1]] <-
      merge(out[[1]],
            out[[2]],
            by = c(strata_nam, "duration"),
            all = TRUE)
  }

  # replace NA in all cols
  for (j in seq_len(ncol(out[[1]]))) set(out[[1]],which(is.na(out[[1]][[j]])),j,0)

  return(invisible(out[[1]]))
}

#' @export
summarise_pop <- function(dt, design, mc) { # mc: the Monte Carlo iteration identifier
  tt <- length(design$scenarios)
  out <- vector("list", tt)
  # names(out) <- design$scenarios

  for (i in seq_len(tt)) {
    out[[i]] <- summarise_hlp(prune_pop(dt, design$scenarios[[i]], design), design$scenarios[[i]], design)
    out[[i]][, scenario := design$scenarios[[i]]]
    out[[i]][, mc := mc]
  }
  out <- rbindlist(out, TRUE)
  setkeyv(out, c(design$strata_for_outputs, "duration"))

  invisible(out)
}

#' @export
summarise_lifecourse <- function(dt, design, mc) { # mc: the Monte Carlo iteration identifier
  tt <- length(design$scenarios)
  out <- vector("list", tt)
  # names(out) <- design$scenarios

  for (i in seq_len(tt)) {
    out[[i]] <- prune_pop(dt, design$scenarios[[i]], design)

    out[[i]][, lifecourse :=
               morbidity_resolve(
                 year,
                 # htn_prvl,
                 # t2dm_prvl,
                 chd_prvl,
                 stroke_prvl,
                 # lung_ca_prvl,
                 all_cause_mrtl,
                 design$friendly_disease_names,
                 design$friendly_fatal_diseases_names,
                 design$init_year
               )]

    del_col_nam <-
      setdiff(grep("_mrtl$", names(out[[i]]), value = TRUE), "all_cause_mrtl")
    out[[i]][, (del_col_nam) := NULL]

    out[[i]][, scenario := design$scenarios[[i]]]
    out[[i]][, mc := mc]
  }
  out <- rbindlist(out, TRUE)
  setkeyv(out, c("pid", "year"))

  invisible(out)
}

#' @export
del_dt_rows <- function(dt, indx_to_del, dt_env = .GlobalEnv) {
  stopifnot(is.data.table(dt), (is.integer(indx_to_del) | is.logical(indx_to_del)))

  dt_keys <- key(dt)
  if (is.integer(indx_to_del)) keep <- -indx_to_del
  if (is.logical(indx_to_del)) keep <- !indx_to_del

  name_of_dt <- deparse(substitute(dt))
  # dt_env <- pryr::where(name_of_dt) # to get dt envirnment
  dt_names <- copy(names(dt))
  dt_new <- dt[keep, dt_names[1L], with = F]
  set(dt, i = NULL, j = 1L, value = NULL)

  for (j in seq_len(ncol(dt))) {
    set(dt_new,
        i = NULL,
        j = dt_names[1L + j],
        value = dt[[1L]][keep])
    set(dt,
        i = NULL,
        j = 1L,
        value = NULL)
  }

  setkeyv(dt_new, dt_keys)
  assign(name_of_dt, value = dt_new, envir = dt_env)
}

