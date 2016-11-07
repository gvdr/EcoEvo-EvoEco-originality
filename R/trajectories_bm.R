#' functions by Giulio Dalla Riva
#' email: gvdr@zoology.ubc.ca
#' trajectories_bm gets as input:
#' 
#' * `tree`, an object of class phylo
#' 
#' * `sigma`, the variance of the Brownian Motion steps
#' 
#' * `speciation_jump`, a scalar that amplifies (if > 1) or dumpen (if < 1)
#'  the divergence happening at the speciation event
#'  
#'  * `early_burst`, a power factor that dumpen the Brownian Motion as the
#'   age of the lineage decrease. If set to `FALSE` or 0 there's no dumping.
#'   
#'   * `trend_sd`, the standard variation of the random directional selection. If set to
#'   zero, the model is an adirectional brownian motion. The greater the variance,
#'   the larger the effect of directional selection. Each branch has its own,
#'   random, directional trend.
#'   
#'   * `set_parallel`, a logical parameter that determine whether plyr work
#'   on parallel or not. Set to `FALSE` by default.
#'   
#'   The function returns a data frame of all the lineage trait values across
#'   the history of the phylogeny.
trajectories_bm <- function(tree,sigma,
                            speciation_jump = 0,
                            early_burst = FALSE,
                            set_parallel = FALSE,
                            trend_sd = 0){
  
  if(trend_sd == 0){
    bm_path <- function(n_steps, sigma, speciation_jump,early_burst,...){
      if(is.numeric(early_burst)){
        sigmas <- sigma/(seq(n_steps-1))^(early_burst)
        path <- cumsum(rnorm(n = n_steps+1, mean = 0, sd = c(0,speciation_jump * sigma,sigmas)))
      } else {
        path <- cumsum(rnorm(n = n_steps+1, mean = 0, sd = c(0,speciation_jump * sigma,rep.int(sigma,n_steps-1))))
      }
      return(path)
    }
  } else {
    bm_path <- function(n_steps, sigma, speciation_jump,early_burst,...){
      if(is.numeric(early_burst)){
        sigmas <- sigma/(seq(n_steps-1))^(early_burst)
        tsd <- trend_sd * sigma
        path <- cumsum(rnorm(n = n_steps+1, mean = rnorm(1,sd = tsd), sd = c(0,speciation_jump * sigma,sigmas)))
      } else {
        tsd <- trend_sd * sigma
        path <- cumsum(rnorm(n = n_steps+1, mean = rnorm(1,sd = tsd), sd = c(0,speciation_jump * sigma,rep.int(sigma,n_steps-1))))
      }
      return(path)
    }
  }
  

  
  tree <- reorder.phylo(tree)
  
  Trajectories <- lapply(tree$edge.length, function(x) bm_path(x,sigma,speciation_jump,early_burst,trend_sd))
  
  final_state <- tree$edge.length + 1
  for (i in 1:nrow(tree$edge)) {
    end_node <- which(tree$edge[, 2] == tree$edge[i, 1])
    if (length(end_node) > 0) {
      Trajectories[[i]] <- Trajectories[[i]] + Trajectories[[end_node]][final_state[end_node]]
    } else {
      Trajectories[[i]] <- Trajectories[[i]] + Trajectories[[1]][1]
    }
  }
  
  H <- nodeHeights(tree)
  get_df <- function(i){
    data.frame(
      times = seq(H[i,1],H[i,2]),
      from = tree$edge[i,1],
      to = tree$edge[i,2],
      coord = Trajectories[[i]]
    )
  }
  
  df_traj <- plyr::ldply(seq_along(Trajectories),get_df, .parallel = set_parallel)
  return(df_traj)
}