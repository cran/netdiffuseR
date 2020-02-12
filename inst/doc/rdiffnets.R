## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----loading, message=FALSE, warning=FALSE------------------------------------
library(netdiffuseR)
library(sna)

## ----random-net---------------------------------------------------------------
set.seed(8826)

# Simulating a small world (10, 3) with pr = .3
net <- rgraph_ws(50, 3, .3)

# A bit more of rewiring
net <- rewire_graph(net, p=3, both.ends = TRUE)

## ----viz----------------------------------------------------------------------

# Visualizing with 
gplot(as.matrix(net))

## ----diffnet------------------------------------------------------------------
# Random diffusion with a fixed threshold of 1, simulating 5 time points
mydiffnet <- rdiffnet(
  seed.graph     = net,                    # The network we just created
  threshold.dist = 1,                      # Fixed threshold of 1
  t              = 5,                      # 5 time points
  rewire         = FALSE,                  # No rewire (defaults TRUE)
  exposure.args  = list(normalized=FALSE), # Exposure to be computed unnormalized 
                                           # so we use counts instead
  seed.nodes     = "random"                # Random set of initial adopters
  )

# Looking at the diffusion process
plot_diffnet(mydiffnet)

# Some summary stats
summary(mydiffnet)

## ----diffnet-multiple, warning=FALSE------------------------------------------
set.seed(871)
mydiffnet <- rdiffnet_multiple(
  statistic      = function(n) cumulative_adopt_count(n)["prop",],
  R              = 1000,
  stop.no.diff   = FALSE,                  # This option allows us to continue
                                           # The simulation process, even if there
                                           # is no adoption.
  seed.graph     = net,                    # The network we just created
  threshold.dist = 1,                      # Fixed threshold of 1
  t              = 5,                      # 5 time points
  rewire         = FALSE,                  # No rewire (defaults TRUE)
  exposure.args  = list(normalized=FALSE), # Exposure to be computed unnormalized 
                                           # so we use counts instead
  seed.nodes     = "random"                # Random set of initial adopters
  )

# Looking at the diffusion process
boxplot(
  t(mydiffnet),
  xlab = "Time",
  ylab = "Proportion of Adopters",
  main = "Simulation of 1,000 diffusion processes"
  )

