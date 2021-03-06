if (!exists("useful.cpp.load")) {
  sourceCpp(file = "usefulfunctions.cpp")
  useful.cpp.load <- TRUE
}

# do not foresee needing to use this on windows.
if (Sys.info()["sysname"] == "Linux") {
  shlib.ext <- ".so"
} else if (Sys.info()["sysname"] == "Darwin") {
  shlib.ext <- ".dylib"
} else if (Sys.info()["sysname"] == "Windows") {
  shlib.ext <- ".dll"
}

filename <- paste("/opt/OpenBLAS/lib/libopenblas", shlib.ext, sep = "")

if (file.exists(filename)) {
  require(inline)
  openblas.set.num.threads <- cfunction( signature(ipt="integer"),
                                         body = 'openblas_set_num_threads(*ipt);',
                                         otherdefs = c ('extern void openblas_set_num_threads(int);'),
                                         libargs = c ('-L/opt/OpenBLAS/lib -lopenblas'),
                                         language = "C",
                                         convention = ".C"
  )
} else {
  openblas.set.num.threads <- function(x) {
    return(NULL)
  }
}

# filename <- paste("~/packages/lib/libopenblas", shlib.ext, sep = "")
#
# if (file.exists(filename)) {
#   require(inline)
#   openblas.set.num.threads <- cfunction( signature(ipt="integer"),
#                                          body = 'openblas_set_num_threads(*ipt);',
#                                          otherdefs = c ('extern void openblas_set_num_threads(int);'),
#                                          libargs = c ('-L~/packages/lib -lopenblas'),
#                                          language = "C",
#                                          convention = ".C"
#   )
# } else {
#   openblas.set.num.threads <- function(x) {
#     return(NULL)
#   }
# }

###############################################
########## Common densities not in R ##########
###############################################
rtnorm <- function(mn, sd = 0.25, fudge = 0){
  upper <- pnorm(1 - fudge, mn, sd)
  lower <- pnorm(fudge, mn, sd)
  if (is.matrix(mn)) {
    U <- matrix(runif(prod(dim(mn)), lower, upper), dim(mn)[1], dim(mn)[2])
  }
  if (!is.matrix(mn)) {
    U <- runif(length(mn), lower, upper)
  }
  return(qnorm(U, mn, sd))
}

dtnorm <- function(y, mn, sd = 0.25, fudge = 0){
  upper <- pnorm(1 - fudge, mn, sd)
  lower <- pnorm(fudge, mn, sd)
  l <- dnorm(y, mn, sd, log = TRUE) - log(upper - lower)
  return(l)
}

dlognormal <- function(x, mu, sig){
  dnorm(log(x), log(mu), sig, log = T) - log(x)
}


###############################################
############# MH Update for MCMC ##############
###############################################
# update candidate standard deviation
mhUpdate <- function(acc, att, MH, nattempts = 50,
                     target.min = 0.3, target.max = 0.6,
                     lower = 0.8, higher = 1.2) {
  acc.rate     <- acc / att
  these.update <- att > nattempts
  these.low    <- (acc.rate < target.min) & these.update
  these.high   <- (acc.rate > target.max) & these.update

  MH[these.low]  <- MH[these.low] * lower
  MH[these.high] <- MH[these.high] * higher

  acc[these.update] <- 0
  att[these.update] <- 0

  results <- list(acc=acc, att=att, MH=MH)
  return(results)
}

###############################################
######### Common data transformations #########
###############################################
transform <- list(
  logit = function(x, lower = 0, upper = 1) {
    x <- (x - lower) / (upper - lower)
    return(log(x / (1 - x)))
  },
  inv.logit = function(x, lower = 0, upper = 1) {
    p <- 1 / (1 + exp(-x))
    p <- p * (upper - lower) + lower
    return(p)
  },
  probit = function(x, lower = 0, upper = 1) {
    x <- (x - lower) / (upper - lower)
    return(qnorm(x))
  },
  inv.probit = function(x, lower = 0, upper = 1) {
    p <- pnorm(x)
    p <- p * (upper - lower) + lower
    return(p)
  },
  log = function(x) log(x),
  exp = function(x) exp(x),
  copula = function(dens) {
    this.dens <- paste("p", dens, sep = "")
    function(x, ...) qnorm(do.call(this.dens, args = list(x, ...)))
  },
  inv.copula = function(dens) {
    this.dens <- paste("q", dens, sep = "")
    function(x, ...) do.call(this.dens, args = list(pnorm(x), ...))
  }
)

# Checks a function for use of global variables
# Returns TRUE if ok, FALSE if globals were found.
checkStrict <- function(f, silent=FALSE) {
  vars <- codetools::findGlobals(f)
  found <- !vapply(vars, exists, logical(1), envir=as.environment(2))
  names <- names(found)[found]

  if ((length(names) > 0)) {
    sum.nfncs <- 0
    for (i in 1:length(names)) {
      if(!is.function(eval(parse(text=names[i])))) {sum.nfncs <- sum.nfncs + 1}
    }
    if (sum.nfncs > 0) {
      warning("global variables used: ", paste(names(found)[found], collapse=', '))
      return(invisible(FALSE))
    }
  }

  !any(found)
}

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

###############################################
########### cross-validation setup ############
###############################################
get.cv.test.srs <- function(n, nfolds) {
  ## Returns a randomly selected set of cross-validation indices based on
  ## how many folds are selected.
  cat("Note: If you want to replicate cross-validation results, be sure to \n")
  cat("      set the seed before running get.cv.test.\n")
  random.cents <- sample(x = 1:n, size = n, replace = F)
  ntest <- ceiling(n / nfolds)
  cv.idx <- vector(mode = "list", length = nfolds)
  for (i in 1:nfolds) {
    start <- (i - 1) * ntest + 1
    if (i < nfolds) {
      end <- i * ntest
    } else {
      end <- n  # in case the last cv set has fewer sites in it
    }
    cv.idx[[i]] <- sort(random.cents[start:end])
  }

  return(cv.idx)
}

################################################################################
# cross-validation, stratified by site
################################################################################
get.cv.test.strat <- function(data, nfolds, idx = NULL) {
  ## Returns a randomly selected set of cross-validation indices based on
  ## how many folds are selected.
  ## data: matrix of the observations
  ## nfolds: how many folds you want
  ## idx: if you want to stratify over rows or colums
  ## returns cv.idx: a list (length = nfolds) of matrices the same size as data
  cat("Note: If you want to replicate cross-validation results, be sure to \n")
  cat("      set the seed before running get.cv.test. \n")

  samples <- data
  if (idx == 1) {
    n <- ncol(data)
    cat("Stratifying the cross-validation by rows. \n")
    for (i in 1:nrow(data)) {
      samples[i, ] <- sample(x = 1:n, size = n, replace = FALSE)
    }
  } else if (idx == 2) {
    n <- nrow(data)
    cat("Stratifying the cross-validation by column. \n")
    for (j in 1:ncol(data)) {
      samples[, j] <- sample(x = 1:n, size = n, replace = FALSE)
    }
  }

  ntest <- ceiling(n / nfolds)
  cv.idx <- vector(mode = "list", length = nfolds)
  for (fold in 1:nfolds) {
    start <- (fold - 1) * ntest + 1
    if (fold < nfolds) {
      end <- fold * ntest
    } else {
      end <- n  # in case the last cv set has fewer sites in it
    }

    # This should return a matrix of TF
    cv.temp <- matrix(FALSE, nrow(data), ncol(data))
    if (idx == 1) {
      for (i in 1:nrow(data)) {
        cv.temp[i, samples[i, start:end]] <- TRUE
      }
    } else if (idx == 2) {
      for (j in 1:ncol(data)) {
        cv.temp[samples[start:end, j], j] <- TRUE
      }
    }
    cv.idx[[fold]] <- cv.temp
  }

  return(cv.idx)
}

get.arr.idx <- function(idx, nrows) {
  # find the row, col indices from a vector index
  col <- ceiling(idx / nrows)
  row <- idx - (col - 1) * nrows
  return(c(row, col))
}

get.idx <- function(row, col, nrows) {
  # find the vector index from the row, col indices
  return((col - 1) * nrows + row)
}