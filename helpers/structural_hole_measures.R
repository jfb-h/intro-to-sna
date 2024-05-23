
# effective size

egonet_effsize <- function(egonet) {
  n <- network.size(egonet) - 1
  t <- network.edgecount(egonet) - n
  return(n - 2 * t / n)
}

# efficiency

egonet_efficiency <- function(egonet){
  size <- network.size(egonet) - 1
  effsize <- egonet_effsize(egonet)
  return(effsize/size)
}

# constraint

# this assumes that EGO is the first node in the network

dyadic_constraint <- function(egonet, alter) {
  pij <- egonet[1, alter] / sum(egonet[1,])
  pqj <- egonet[alter,-1] / sum(egonet[alter,])
  pqj[is.nan(pqj)] <- 0
  return((pij + sum(pij * pqj))^2)
}

egonet_constraint <- function(egonet) {
  alteri <- 2:network.size(egonet)
  dc <- sapply(alteri, function(a) dyadic_constraint(egonet, a))
  return(sum(dc))
}

# hierarchy

egonet_hierarchy <- function(egonet) {
  n <- network.size(egonet) - 1
  alteri <- 2:(n+1)
  dc <- sapply(alteri, function(a) dyadic_constraint(egonet, a))
  C <- sum(dc)
  dcnorm <- dc/(C/n)
  
  if (n == 0) return(NA) else if (n == 1) return(1) else {
    return(sum(dcnorm * log(dcnorm)) / (n*log(n)))
  }
} 

# apply to full network

apply_egonet <- function(net, FUN) {
  if (is.directed(net)) stop("Not implemented for directed graphs.")
  egonets <- lapply(ego.extract(net), FUN = network, directed = FALSE)
  return(sapply(egonets, FUN))
}

# run tests

# library(statnet)
# data(florentine)
# 
# net = flomarriage
# egonet = network(ego.extract(net, ego = 9)[[1]], directed = FALSE)
# gplot(egonet, gmode="graph", displaylabels = T)
# 

# egonet_effsize(egonet) == 6 - (4*0 + 2*1) / 6
# egonet_efficiency(egonet) == (6 - (4*0 + 2*1) / 6) / 6
# dyadic_constraint(egonet, 5) == (1/6 + 1/6*1/2)^2
# egonet_constraint(egonet) == 4 * 1/6^2 + 2 * (1/6 + 1/6 * 1/2)^2
# round(egonet_hierarchy(egonet), 3) == 0.045
# 
# apply_egonet(flomarriage, egonet_constraint)
# 

