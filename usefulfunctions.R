# do not foresee needing to use this on windows.
if (Sys.info()["sysname"] == "Linux") {
  shlib.ext <- ".so"
} else if (Sys.info()["sysname"] == "Darwin") {
  shlib.ext <- ".dylib"
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
    p <- exp(x) / (1 + exp(x))
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