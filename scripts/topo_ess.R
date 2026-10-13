#Computing the approximate topological ESS from rwty
topo_ess=function(trees) {
  filename=tempfile()
  f=file(filename,'w')
  for (i in 1:length(trees)) write.tree(trees[[i]],file=f,append=T)
  close(f)
  x=capture.output({
  r=rwty::load.trees(filename,type = "newick")
  v=rwty::topological.approx.ess(r,burnin=0)$approx.ess
  })
  return(v)
}
