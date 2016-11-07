#' A simple plotting function for plotting
#' one unidimensional trajectory.
plot_trajectory <- function(traject_data){
plot_x <- traject_data %>%
    ggplot(aes(x = times, 
               y = coord,
               group = to,
               colour = as.factor(to)
    )
    ) +
    geom_line() +
    scale_color_viridis(discrete=TRUE, guide = FALSE) +
    labs(y = "Trait value",
         x = "Time from root")  +
    theme_minimal()
  return(plot_x)
}