#' #' Compute prior using stan
#' #' @param r Object of class ContPhylo
#' #' @param scale Scale of the exponential distribution for branch lengths
#' #' @returns Log of the prior
#' #' @export
#' prior_stan <- function(r,scale=1) {
#'   a=r$coordinates[-1,1]
#'   b=r$coordinates[-1,2]
#'   n <- length(a) + 1L
#'   invisible(capture.output(suppressMessages({
#'   mod <- cmdstan_model(system.file('stan','prior.stan',package='ContPhylo',mustWork = T))
#'   })))
#'   data_file <- tempfile(fileext = ".json")
#'   pars_file <- tempfile(fileext = ".json")
#'   output_file <- tempfile(fileext = ".csv")
#'   cmdstanr::write_stan_json(list(n=n,scale=scale), data_file)
#'   parameters <- list(
#'     a = array(a, dim = length(a)),
#'     b = array(b, dim = length(b))
#'   )
#'   cmdstanr::write_stan_json(parameters, pars_file)
#'   args <- c(
#'     "log_prob",
#'     paste0("constrained_params=", pars_file),
#'     paste0("jacobian=", 'true'),
#'     "data", paste0("file=", data_file),
#'     "output", paste0("file=", output_file),
#'     "refresh=0"
#'   )
#'     system2(mod$exe_file(), args = args,stdout=F,stderr=F)
#'   result <- utils::read.csv(output_file, comment.char = "#", check.names = FALSE)
#'   log_prior=unname(result$lp__[[1]])
#'   return(log_prior)
#' }
#'
#' #' Compute posterior using stan
#' #' @param mypml Object of class phangorn::pml containing the data
#' #' @param r Object of class ContPhylo
#' #' @param scale Scale of the exponential distribution for branch lengths
#' #' @returns Log of the prior
#' #' @export
#' posterior_stan <- function(mypml,r,scale=1) {
#'   a=r$coordinates[-1,1]
#'   b=r$coordinates[-1,2]
#'   n <- length(a) + 1L
#'   invisible(capture.output(suppressMessages({
#'     mod <- cmdstan_model(system.file('stan','posterior.stan',package='ContPhylo',mustWork = T))
#'   })))
#'   data_file <- tempfile(fileext = ".json")
#'   pars_file <- tempfile(fileext = ".json")
#'   output_file <- tempfile(fileext = ".csv")
#'   gendata=phyDat_to_stan(mypml$data)
#'   cmdstanr::write_stan_json(list(n=n,scale=scale,#root_distance=scale/2,root_sd=0.1,
#'     n_patterns=gendata$n_patterns,tip_state=gendata$tip_state,
#'     pattern_count=array(gendata$pattern_count,dim=length(gendata$pattern_count))), data_file)
#'   parameters <- list(
#'     a = array(a, dim = length(a)),
#'     b = array(b, dim = length(b))
#'   )
#'   cmdstanr::write_stan_json(parameters, pars_file)
#'   args <- c(
#'     "log_prob",
#'     paste0("constrained_params=", pars_file),
#'     paste0("jacobian=", 'true'),
#'     "data", paste0("file=", data_file),
#'     "output", paste0("file=", output_file),
#'     "refresh=0"
#'   )
#'   system2(mod$exe_file(), args = args,stdout=F,stderr=F)
#'   result <- utils::read.csv(output_file, comment.char = "#", check.names = FALSE)
#'   log_posterior=unname(result$lp__[[1]])
#'   return(log_posterior)
#' }
#'
#' #' HMC algorithm using stan
#' #' @param dat Object of class phyDat containing the genetic data
#' #' @param initTree Initial tree to use as start of MCMC
#' #' @param niter Number of iterations
#' #' @param thin Thining interval
#' #' @param scale Scale of the exponential distribution for branch lengths
#' #' @param showProgress Whether to show a progress bar
#' #' @returns HMC result
#' #' @export
#' hmc = function(dat=phyDatEmpty(10),initTree,niter=1e3,thin=max(1,round(niter/1e3)),scale=1,showProgress=T) {
#'   if (missing(initTree)) {
#'     n=length(dat)
#'     tree=rtree(n)
#'     tree$edge.length=tree$edge.length*scale
#'     tree$tip.label=names(dat)
#'   } else {
#'     n=Ntip(initTree)
#'     tree=initTree
#'   }
#'   r=phylo2contphylo(tree)
#'   a=r$coordinates[-1,1]
#'   b=r$coordinates[-1,2]
#'   root_distance=initTree$edge.length[which(initTree$edge[,2]==1)]/2
#'
#'   mod <- cmdstan_model(system.file('stan','posterior.stan',package='ContPhylo',mustWork = T))
#'   fit <- mod$sample(
#'     data = c(phyDat_to_stan(dat),list(scale = scale)),#,root_distance=root_distance,root_sd=1e-3)),
#'     adapt_delta = 0.99,
#'     chains = 1,
#'     parallel_chains = 1,
#'     iter_warmup = niter / 2,
#'     iter_sampling = niter / 2,
#'     init = list(list(a = a, b = b)),
#'     seed = 1L,
#'     refresh = ifelse(showProgress,thin,0))
#'   mat=fit$draws(format = 'matrix')
#'   lp=fit$lp()
#'   res=vector('list',niter/2/thin)
#'   l=list()
#'   l$tip.label=tree$tip.label
#'   class(l)<-'contphylo'
#'   for (i in seq(thin,niter/2,thin)) {
#'     l$coordinates=rbind(c(0,0),matrix(mat[i,2:(2*n-1)],ncol=2,byrow=F))
#'     res[[i/thin]]=contphylo2phylo(l)
#'     res[[i/thin]]$posterior=lp[i]
#'   }
#'   res[[1]]$fit=fit#interesting
#'   return(res)
#' }
