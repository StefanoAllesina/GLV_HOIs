library(tidyverse)
# you also need ggtern

pl_biomass <- function(dynamics, pars, eq){
  dt <- dynamics %>% group_by(time, run) %>% 
    summarise(abundance = sum(abundance), .groups = "drop")
  pl_b <- ggplot(dt, aes(x = time, y = abundance, colour = run)) + 
    geom_line()
}

pl_vstar <- function(dynamics, pars, eq){
  dt <- dynamics %>% 
    pivot_wider(names_from = variable, values_from = abundance)
  # compute vstar
  vstar <- numeric(0)
  for (i in 1:nrow(dt)){
    vstar <- c(vstar, Vstar_eq(pars, c(dt$`1`[i], dt$`2`[i], dt$`3`[i])))
  }
  dt <- dt %>% add_column(vstar = vstar)
  pl_v <- ggplot(dt, aes(x = time, y = vstar, colour = run)) + 
    geom_line()
}

pl_dynamics_tern <- function(dynamics, pars, eq){
  pl_tern <- ggtern::ggtern(eq, aes(x = p1, y  = p2, z = p3)) + 
    geom_point(size = 3, aes(shape = stability))
  dt <- dynamics %>% 
    pivot_wider(names_from = variable, values_from = abundance) %>% 
    mutate(p1 = `1` / (`1` + `2` + `3`), p2 = `2` / (`1` + `2` + `3`),
           p3 = `3` / (`1` + `2` + `3`))
  
  return(pl_tern + geom_line(data = dt, aes(colour = run)))
}