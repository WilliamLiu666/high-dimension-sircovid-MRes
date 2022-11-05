load_transform <- function(path, region, assumptions) {
  e <- new.env()
  sys.source(file.path(path, "transform.R"), e)
  stopifnot(is.function(e$make_transform),
            is.function(e$apply_assumptions))
  make_transform <- e$make_transform
  base <- e$apply_assumptions(
    readRDS(file.path(path, "base.rds"))[[region]], assumptions)
  transform <- make_transform(base)
  list(base = base, transform = transform)
}


load_parameters <- function(path, region, assumptions, kernel_scaling) {
  parameters <- spimalot::spim_pars_pmcmc_load(path)

  info <- spimalot:::spim_pars_info(region, parameters$info)
  prior <- spimalot:::spim_pars_prior(region, info, parameters$prior)
  proposal <- spimalot:::spim_pars_proposal(region, info, parameters$proposal,
                                            kernel_scaling)
  
  transform <- load_transform(path, region, assumptions)

  ret <- list(region = region,
              assumptions = assumptions,
              info = info,
              prior = prior,
              proposal = proposal,
              transform = transform$transform,
              raw = parameters,
              base = transform$base)
  ret$mcmc <- build_parameters(ret)
  ret
}


build_parameters <- function(pars) {
  info <- pars$info
  prior <- pars$prior

  pars_mcmc_list <- Map(
    mcstate::pmcmc_parameter,
    name = info$name,
    initial = info$initial,
    min = info$min,
    max = info$max,
    discrete = info$discrete,
    prior = lapply(split(prior, prior$name), spimalot:::make_prior))

  ret <- mcstate::pmcmc_parameters$new(pars_mcmc_list,
                                       pars$proposal, pars$transform)

  ## Try and transform a single case and see if it works:
  ret$model(ret$initial())

  ret
}


pars_variable <- function(..., values. = list(...)) {
  stopifnot(is.list(values.) && !is.null(names(values.)))
  class(values.) <- "pars_variable"
  values.
}


pars_apply_assumptions <- function(x, assumptions) {
  f <- function(x) {
    if (inherits(x, "pars_variable")) {
    }
  }
  for (i in seq_along(x)) {
    el <- x[[i]]
    if (inherits(el, "pars_variable")) {
      if (!(assumptions %in% names(el))) {
        stop(sprintf("assumptions '%s' not found for '%s'",
                     assumptions, names(x)[[i]]))
      }
      x[[i]] <- el[[assumptions]]
    }
  }
  x
}


load_baseline <- function(path, assumptions) {
  baseline <- readRDS(path)
  if (!is.null(assumptions)) {
    baseline <- lapply(baseline, pars_apply_assumptions, assumptions)
  }
  baseline
}
