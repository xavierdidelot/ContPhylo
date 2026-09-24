#' Compute prior using stan
#' @param r Object of class ContPhylo
#' @param shape Shape of the gamma distribution for branch lengths
#' @param scale Scale of the gamma distribution for branch lengths
#' @returns Log of the prior
#' @export
prior_stan <- function(r,shape=1,scale=1) {
  a=r$coordinates[-1,1]
  b=r$coordinates[-1,2]
  n <- length(a) + 1L
  invisible(capture.output(suppressMessages({
  mod <- cmdstan_model(system.file('stan','model.stan',package='ContPhylo',mustWork = T))
  fit <- mod$sample(
    data=list(n=n,shape=shape,scale=scale),
    chains=1,parallel_chains=1,iter_warmup=0,iter_sampling=1,
    fixed_param=TRUE,init=list(list(a=a,b=b)),seed=1L,refresh=0)
  })))
  log_prior=as.numeric(fit$draws(variables='log_prior',format='matrix'))
  log_prior=log_prior+(n-1)*log(2)+lfactorial(n-1)-lfactorial(2*n-2)
  return(log_prior)
}
