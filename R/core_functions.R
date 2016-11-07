one_tree_trajectory <- function(x, Dim_trait = 1, ...){
  trj <- rdply(Dim_trait,
               traits(tree = x, ...),
               .id = "coord_n"
  )
  return(trj)
}

d_trj <- function(tree_list,
                  Dim_trait = 1,
                  sigma = 0.2,
                  speciation_jump = 1,
                  early_burst = 0,
                  trend_sd = 0) {
  
  trajectories <- ldply(tree_list,
                        one_tree_trajectory,
                        Dim_trait = Dim_trait,
                        sigma = sigma,
                        speciation_jump = speciation_jump,
                        early_burst = early_burst,
                        trend_sd = trend_sd,
                        .id = "tree")
  
  return(trajectories)
}

get_eco_data <- function(trajectories,
                         variables = "tree",
                         across = "coord_n"){
  trajectories %>%
    group_by_(.dots = c(variables)) %>%
    widyr::pairwise_dist_(item = "Species",
                         feature = across,
                         value = "coord") %>%
    group_by_(.dots = c("item1",variables)) %>%
    summarise(eco_ave = mean(distance), # average Ecological Originality
              eco_min = min(distance), # minimum Ecological Originality
              eco_max = max(distance)) %>% # maximum Ecological Originality
    rename(Species = item1) %>%
    ungroup() %>%
    return()
}

consolidate <- function(eco_data, evo_data){
  all_data_trend <- eco_data %>%
    full_join(evo_data,.,by=c("tree","Species"))
  all_data_trend$tree <- as.integer(all_data_trend$tree)
  return(all_data_trend)
}

do_good_of_fit <- function(ecoevo_data, variables = "tree"){
  ecoevo_data %>%
    group_by_(.dots = variables) %>%
    do(good_of_fit(.)) %>%
    return()
}

traits <- function(tree,...){
  traj <- filter(trajectories_bm(tree,...),times == max(times))
  tips <- tree$tip
  traj <- mutate(traj, Species = tips)
  traj <- select(traj,-from,-to,-times)
  return(traj)
}

get_estimate <- function(model){
  model %>%
    tidy() %>%
    filter(grepl('eco', term)) %>%
    select_("estimate") %>%
    return()
}

good_of_fit <- function(uniqueness_dataframe){
  
  models <- uniqueness_dataframe %>%
    fit_with(lm,
             formulae(~w,
                      ave = ~eco_ave,
                      max = ~eco_max,
                      min = ~eco_min))
  
  return_data <- models %>%
    ldply(glance,
          .id = "ecological_uniqueness")
  
  models %>%
    ldply(get_estimate,
          .id = "ecological_uniqueness") %>%
    full_join(return_data,
              by = "ecological_uniqueness") %>%
    return()
}