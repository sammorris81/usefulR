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

################################################################################
# Common data transformations
################################################################################
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