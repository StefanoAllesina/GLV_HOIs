library(tidyverse)

prune_equilibria <- function(eq){
  pruned <- list()
  for(eql in eq){
    if (all(is.infinite(eql) == FALSE)){
      if (all(Im(eql) == 0)){
        eql <- Re(eql)
        if (all(eql >= 0)){
          if (sum(eql) > 0){
            pruned[[length(pruned) + 1]] <- eql
          }
        }
      }
    }
  }
  return(pruned)
}

get_feasible_eq <- function(r, a1, a2, b, alpha){
  r <- as.complex(r)
  a1 <- as.complex(a1)
  a2 <- as.complex(a2)
  b <- as.complex(b)
  alpha <- as.complex(alpha)
  eq <- list(
    c(0, 0, -(r/(a1*alpha))),
    c(0, -(r/((a1 + a2)*alpha)), -(r/((a1 + a2)*alpha))),
    c(0, -(r/(a1*alpha)), 0),
    c(alpha*(a2 - a1)/(b*(alpha-1)), alpha*(a2 - a1)/(b*(alpha-1)), -(((a1^2 - a2^2)*alpha^2 + b * r *(1-alpha))/(a1 * b * alpha*(1-alpha))) ),
    c(alpha*(a2 - a1)/(b*(alpha-1)), -(((a1^2 - a2^2)*alpha^2 + b * r *(1-alpha))/(a1 * b * alpha*(1-alpha))), alpha*(a2 - a1)/(b*(alpha-1)) ),
    c(-(((a1^2 - a2^2)*alpha^2 + b * r *(1-alpha))/(a1 * b * alpha*(1-alpha))), alpha*(a2 - a1)/(b*(alpha-1)), alpha*(a2 - a1)/(b*(alpha-1))),
    c((((a1 + 2 * a2)*alpha - sqrt(((a1 + 2 * a2)*alpha)^2 + 4 * b * r *(alpha-1)))/(2 * b*(alpha-1))), (((a1 + 2 * a2)*alpha - sqrt(((a1 + 2 * a2)*alpha)^2 + 4 * b * r *(alpha-1)))/(2 * b*(alpha-1))), (((a1 + 2 * a2)*alpha - sqrt(((a1 + 2 * a2)*alpha)^2 + 4 * b * r *(alpha-1)))/(2 * b*(alpha-1)))),
    c((((a1 + 2 * a2)*alpha + sqrt(((a1 + 2 * a2)*alpha)^2 + 4 * b * r *(alpha-1)))/(2 * b*(alpha-1))), (((a1 + 2 * a2)*alpha + sqrt(((a1 + 2 * a2)*alpha)^2 + 4 * b * r *(alpha-1)))/(2 * b*(alpha-1))), (((a1 + 2 * a2)*alpha + sqrt(((a1 + 2 * a2)*alpha)^2 + 4 * b * r *(alpha-1)))/(2 * b*(alpha-1)))),
    c(-(r/((a1 + a2)*alpha)), 0, -(r/((a1 + a2)*alpha))),
    c(-(r/(a1*alpha)), 0, 0),
    c(-(r/((a1 + a2)*alpha)), -(r/((a1 + a2)*alpha)), 0),
    c(0, 0, 0)
  )
  
  if (alpha == 1.0){
    eq <- list(
      c(0, 0, -(r/(a1*alpha))),
      c(0, -(r/((a1 + a2)*alpha)), -(r/((a1 + a2)*alpha))),
      c(0, -(r/(a1*alpha)), 0),
      c(-r/(a1+2*a2), -r/(a1+2*a2), -r/(a1+2*a2)),
      c(-(r/((a1 + a2)*alpha)), 0, -(r/((a1 + a2)*alpha))),
      c(-(r/(a1*alpha)), 0, 0),
      c(-(r/((a1 + a2)*alpha)), -(r/((a1 + a2)*alpha)), 0),
      c(0, 0, 0)
    )
  }
  return(eq)
}


stability_eq <- function(pars, alpha, xstar){
  r <- pars[1]
  a1 <- pars[2]
  a2 <- pars[3]
  b <- pars[4]
  x1 <- xstar[1]
  x2 <- xstar[2]
  x3 <- xstar[3]
  M <- matrix(c(
    r + alpha*(2 * a1 * x1 + a2 * x2 + a2 * x3) + (1-alpha)*b * x2 * x3, 
    x1 * (alpha*a2 + (1-alpha)*b * x3), 
    x1 * (alpha*a2 + (1-alpha)*b * x2), 
    x2 * (alpha*a2 + (1-alpha)*b * x3), 
    r + alpha*(a2 * x1 + 2 * a1 * x2 + a2 * x3) + (1-alpha)*b * x1 * x3, 
    (alpha*a2 + (1-alpha)*b * x1) * x2, (alpha*a2 + (1-alpha)*b * x2) * x3, 
    (alpha*a2 + (1-alpha)*b * x1) * x3, 
    r + alpha*(a2 * x1 + a2 * x2 + 2 * a1 * x3) + (1-alpha)*b * x1 * x2
  ), 3, 3, byrow = TRUE)
  eM <- eigen(M, only.values = TRUE)$values
  if (max(Re(eM)) < 0) return(TRUE)
  return(FALSE)
}

Vstar_eq <- function(pars, alpha, xstar){
  r <- pars[1]
  a1 <- pars[2]
  a2 <- pars[3]
  b <- pars[4]
  x1 <- xstar[1]
  x2 <- xstar[2]
  x3 <- xstar[3]
  Vstar <- r * sum(xstar)
  Vstar <- Vstar + 0.5 * alpha * (a1 * sum(xstar^2) + 2 * a2 *(xstar[1] * xstar[2] + 
                                                         xstar[1] * xstar[3]+ 
                                                         xstar[2] * xstar[3]))
  Vstar <- Vstar + (1-alpha) * prod(xstar) * b
  return(Vstar)
}

get_equilibria <- function(pars, alpha){
  eq <- get_feasible_eq(pars[1], pars[2], pars[3], pars[4], alpha)
  eq <- prune_equilibria(eq)
  xstar <- matrix(0, 0, 3)
  vstar <- numeric(0)
  stability <- logical(0)
  for (eql in eq){
    xstar <- rbind(xstar, eql)
    stability <- c(stability, stability_eq(pars, alpha, eql))
    vstar <- c(vstar, Vstar_eq(pars, alpha, eql))
  }
  colnames(xstar) <- c("x1", "x2", "x3")
  xstar <- as_tibble(xstar)
  xstar <- xstar %>% add_column(stability = stability)
  xstar <- xstar %>% add_column(Vstar = vstar)
  xstar <- xstar %>% arrange(desc(Vstar), desc(x1), desc(x2), desc(x3))
  # add a label
  xstar <- xstar %>% mutate(label = row_number())
  # add total biomass and proportions
  xstar <- xstar %>% mutate(
    biom = x1 + x2 + x3,
    p1 = x1 / biom,
    p2 = x2 / biom,
    p3 = x3 / biom
  )
  return(xstar)
}
