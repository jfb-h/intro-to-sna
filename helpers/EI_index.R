

EI_index <- function(net, groups, out = "all"){
  
  if (is.directed(net)) {
    warning("Not implemented for directed networks. Network is symmetrized by weak rule.")
    net <- network(sna::symmetrize(net, rule = "weak"), directed = FALSE)
  }
  
  el <- as.edgelist(net, directed = FALSE)
  el <- cbind(el, NA)
  groups <- factor(groups)
  
  el[,3] <- ifelse(groups[el[,1]] == groups[el[,2]], 1, 0)
  
  # global EI-index
  
  I <- sum(el[,3])
  E <- nrow(el) - I
  global <- (E - I) / (E + I)

  if (out == "global") return(global)

  # individual EI-index 
  
  N <- network.size(net)
  df_v <- data.frame("E" = NA, "I" = NA, "EI_index" = NA)
  
  for (i in 1:network.size(net)) {
    deg <- degree(net, nodes = i, gmode = "graph")
    I <- sum(el[el[,1] == i | el[,2] == i, 3])
    E <- deg - I
    EI <- (E - I) / (E + I)
    df_v[i,] <- c(E, I, EI)
  }

  # group-level EI-index
  
  df_g <- data.frame("group" = levels(groups))
  df_g$N <- c(table(groups))
  df_g$E <- c(tapply(df_v$E, INDEX = groups, FUN = sum))
  df_g$I <- c(tapply(df_v$I, INDEX = groups, FUN = sum))
  df_g$EI <- c((df_g$E - df_g$I) / (df_g$E + df_g$I))

  # rescaled EI-index
  
  rescale <- function(x, maxold, minold, maxnew, minnew){
    res <- (maxnew-minnew)/(maxold-minold)*(x-maxold)+maxnew
    return(res)
  }
  
  S_i <- table(groups)
  L <- nrow(el)
  m <- S_i %*% t(S_i)
  I_star <- 1/2 * sum(diag(m) - S_i)
  E_star <- sum(m[upper.tri(m)])
  
  EI_expected <- (E_star - I_star)/(E_star + I_star)
  EI_max <- ifelse(E_star <= L, (E_star-(L-E_star))/(L), 1)
  EI_min <- ifelse(I_star <= L, ((L-I_star)-I_star)/(L), -1)
  
  rescaled <- rescale(global, EI_max, EI_min, 1, -1)

  results <- list()
  results$global <- global
  results$rescaled <- rescaled
  results$individual <- df_v
  results$group <- df_g
  return(results)
}


