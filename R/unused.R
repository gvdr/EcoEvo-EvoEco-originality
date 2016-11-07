
la_d_trj <- function(tree_list, # the list of trees
                     traj_parameters, # the array of parameters
                     sigma_simulations = 0.2
){
  
  trajectories <- ldply(tree_list,
                        traits,
                        sigma = sigma_simulations,
                        speciation_jump = traj_parameters$sj,
                        early_burst = traj_parameters$eb,
                        .id = "tree")
  
  return(trajectories)
}

la_dim_d_trj <- function(tree_list,
                         traj_parameters,
                         Dim_trait,
                         sigma_simulations = 0.2){
  
  trajectories <- ldply(tree_list, function(x)
    rdply(Dim_trait,
          traits(tree = x,
                 sigma = sigma_simulations,
                 speciation_jump = traj_parameters$sj,
                 early_burst = traj_parameters$eb),
          .id = "coord_n"
    ),
    .id = "tree")
  return(trajectories)
}

la_dim_trend_d_trj <- function(tree_list,
                               Dim_trait,
                               Trend_sd,
                               sigma_simulations = 0.2){
  
  trajectories <- llply(tree_list, function(x)
    rdply(Dim_trait,
          traits(tree = x,
                 sigma = sigma_simulations,
                 trend_sd = Trend_sd),
          .id = "coord_n"
    ))
  return(trajectories)
}


param_trj <- function(param_space,tree_list,...) {
  adply(param_space,
        1,
        function(x){
          la_d_trj(tree_list = tree_list,
                   list(sj = x[,1],
                        eb = x[,2]))
        }
  ) %>%
    as_tibble() %>%
    return()
}

param_dim_trj <- function(param_space,tree_list,...){
  adply(param_space,1,
        function(x){
          la_dim_d_trj(tree_list = tree_list,
                       list(sj = x[,1],eb = x[,2]),
                       Dim_trait = x[,3])
        }
  ) %>%
    as_tibble() %>%
    return()
}


tree_trait_gof <- function(all_data,mv = FALSE){ 
  
  if(mv){
    all_gof <- ldply(all_data,td_gof_mv) %>%
      as.data.frame()
  } else {
    all_gof <- ldply(all_data,td_gof) %>%
      as.tbl()
  }
  return(all_gof)
}