#' Function returning a phyDat with n sequences and zero site contribution
#' @param n Number of sequences
#' @returns Object of class phyDat
#' @export
phyDatEmpty=function(n=10) {
  x=matrix("n", nrow = n, dimnames = list(paste0("seq", seq_len(n)), NULL))
  dat=phangorn::phyDat(x, type = "DNA")
  attr(dat, "weight")=0#so that the site is ignored
  return(dat)
}
