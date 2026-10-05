#' Convert a contphylo into a phylo
#'
#' @param r Object of class `contphylo`.
#' @returns Object of class `phylo`.
#' @export
contphylo2phylo = function(r) {
  stopifnot(inherits(r, "contphylo"))
  type=r$type
  #Type of tree: 1 for ultrametric, 2 for rooted, 3 for unrooted.
  coords=r$coordinates
  n=nrow(coords)
  if (n<2) stop('Need at least two leaves.')
  if (type==3 && n<3) stop('Need at least three leaves for unrooted trees.')
  m=matrix(NA,n*2-1,4)#child1,child2,father,dist_to_father
  if (type==1) heights=rep(0,n*2-1)#for ultrametric case only
  for (i in 2:n) {
    if (type==3 && i==n) {
      graft=which(!is.na(m[,1])&is.na(m[,3]))
      graft_dist=0
    } else if (i==2) {
      graft=1
      graft_dist=abs(coords[i,1])
    } else {
      g=contour_traversal(m,coords[i,1])
      graft=g[1]
      graft_dist=g[2]
    }
    m[i,3]=n+i-1
    if (type==1) {
      m[i,4]=graft_dist+heights[graft]
      heights[n+i-1]=m[i,4]
    } else m[i,4]=exp(coords[i,2])
    m[n+i-1,1]=graft
    m[n+i-1,2]=i
    m[n+i-1,3]=m[graft,3]
    m[n+i-1,4]=m[graft,4]-graft_dist
    if (!is.na(m[m[graft,3],1])&&m[m[graft,3],1]==graft) m[m[graft,3],1]=n+i-1
    if (!is.na(m[m[graft,3],2])&&m[m[graft,3],2]==graft) m[m[graft,3],2]=n+i-1
    m[graft,3]=n+i-1
    m[graft,4]=graft_dist
  }
  root=which(is.na(m[,3]))
  if (root!=n+1) {
    #swap root and n+1
    row=m[n+1,]
    m[n+1,]=m[root,]
    m[root,]=row
    m[which(m[,1:3]==n+1, arr.ind=T)]=-1
    m[which(m[,1:3]==root,arr.ind=T)]=n+1
    m[which(m[,1:3]==-1  ,arr.ind=T)]=root
  }
  p=list()
  p$Nnode=n-1
  p$edge=cbind(m[,3],1:(n*2-1))
  p$edge=p$edge[-(n+1),]
  p$edge.length=m[-(n+1),4]
  p$tip.label=r$tip.label
  class(p)='phylo'
  if (type==3) p=unroot(p)
  return(p)
}

#' Perform contour traversal
contour_traversal = function(m,d) {
  root=which(!is.na(m[,1])&is.na(m[,3]))
  remainder=d
  from=root
  cur=m[root,1]
  if (d<0) return(c(root,-d))
  while (T) {
    if (cur==root && from==m[root,2]) break
    if (cur==root && from==m[root,1]) {
      from=root
      cur=m[root,2]
      next
    }
    if (remainder<m[cur,4] && from==m[cur,3]) {
      remainder=m[cur,4]-remainder
      return(c(cur,remainder))
    }
    if (is.na(m[cur,1]) && remainder<m[cur,4]*2) {
      remainder=remainder-m[cur,4]
      return(c(cur,remainder))
    }
    if (is.na(m[cur,1])) {
      remainder=remainder-2*m[cur,4]
      from=cur;cur=m[cur,3]
      next
    }
    if (from==m[cur,3]) {
      remainder=remainder-m[cur,4]
      from=cur
      cur=m[cur,1]
      next
    }
    if (from==m[cur,1]) {
      from=cur
      cur=m[cur,2]
      next
    }
    if (from==m[cur,2]) {
      if (remainder < m[cur,4]) return(c(cur,remainder))
      remainder=remainder-m[cur,4]
      from=cur
      cur=m[cur,3]
      next
    }
    stop('not happening')
    break
  }
  return(c(cur,remainder))
}

