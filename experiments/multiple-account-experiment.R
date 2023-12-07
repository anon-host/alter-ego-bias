library(plyr)
library(dplyr)
library(expm)
library(reshape2)
library(bit)

source("alter-ego-bias/graphs/graph-utils.R")
source("alter-ego-bias/experiments/experiment-utils.R")

select.alter.egos <- function(alter.egos, treatment.assignments, setting) {
  egos <- which(alter.egos==1)
  treat <- which(treatment.assignments==1)
  ctrl <- which(treatment.assignments==0)
  if(setting == "random") {
    to.combine <- sample(egos, 2)
    return(to.combine)
  }
  else {
    treatment.egos <- intersect(egos, treat) 
    control.egos <- intersect(egos, ctrl) 

    if(length(treatment.egos) == 1){
      combine.treat <- as.integer(treatment.egos)
    }
    else{
      combine.treat <- as.integer(sample(treatment.egos, 1))
    }
    if(length(control.egos) == 1){
      combine.ctrl <- as.integer(control.egos)
    }
    else{
      combine.ctrl <- as.integer(sample(control.egos, 1))
    }
    return("selected"=c(combine.treat, combine.ctrl))
  }
}


multiple.account.experiment <- function(graph.params, clustering, ego.params, outcome.params, setting="greedy") { 
  # generate graph structure
  g <- generate.graph(graph.params)
  graph.properties <- get.graph.properties(g)
  graph.params$n <- graph.properties$n
  V(g)$name <- 1:graph.properties$n
  
  avg.degree <- mean(graph.properties$degrees)

  # generate graph clustering
  clusters <- generate.clusters(graph.properties$g, clustering)
  if(graph.params$graph.type=="stars"){
    clusters <- rep(0,graph.properties$n)
    clusters[1:graph.properties$n/2] <- 1
    clusters[(graph.properties$n/2+1):graph.properties$n] <- 2
  }
  else if(sum(clusters==1)==graph.properties$n) stop("Only one cluster found")

  # assign treatment 
  treatment <- treatment.assignment(graph.properties$g, clusters)
  if(graph.params$graph.type=="stars"){
    treatment.assignments <- clusters
    treatment.assignments[clusters==2] <- 0
  }
  else treatment.assignments <- treatment[clusters]

  # prepare outcome model parameters
  if(graph.params$graph.type=="facebook") { 
    noise <- FALSE
  } else { 
    noise <- TRUE
  }
  stochastic.vars <- get.stochastic.vars(graph.properties$n, 3, 0.1, noise)
  
  bias.behavior <- data.frame(index=numeric(), size.of.dom=logical(), method=character(), pt.uncovered=numeric(), alter.ego.influence=numeric(), ATE.true=numeric(), ATE.ego.gui=numeric(), gui.beta=numeric(), gui.gamma=numeric(), stringsAsFactors=FALSE)
  nonego.ATE <- as.numeric(calculate.ATE.various(0, graph.properties, matrix(0,1,graph.properties$n), outcome.params, ego.params, treatment.assignments, stochastic.vars, bias.behavior)$ATE.ego.gui[1])

  ego.params$setting <- setting
  ego.params$max <- TRUE
  ego.params$weighting <- "inf"
  ego.params$num.ego <- graph.properties$n/2
  if(graph.params$graph.type=="facebook") { 
    dominating.alter.egos.deg <- matrix(0, 1, graph.properties$n)
    if(ego.params$setting == "greedy"){
      alter.egos <- c(108, 3438, 1, 1685, 1913, 349, 415, 3981, 687, 699)
    }
    else{
      alter.egos <- sample(1:graph.properties$n,ego.params$num.ego,replace = FALSE)
    }
    dominating.alter.egos.deg[,sample(alter.egos, length(alter.egos))] <- 1
  }
  else if(graph.params$graph.type=="stars"){
    if(ego.params$setting == "greedy"){
      dominating.alter.egos.deg <- matrix(0, 1, graph.properties$n)
      dominating.alter.egos.deg[1] <- 1
      dominating.alter.egos.deg[graph.properties$n/2+1] <- 1
    }
    else{
      dominating.alter.egos.deg <- matrix(0, 1, graph.properties$n)
      alter.egos <- sample(1:graph.properties$n,ego.params$num.ego,replace = FALSE)
      dominating.alter.egos.deg[,sample(alter.egos, length(alter.egos))] <- 1
    }
  } 
  else { 
    alter.ego.list <- determine.alter.egos(graph.properties, ego.params)
    dominating.alter.egos.deg <- alter.ego.list
  }
  
  ego.params$max.dom.ego <- max(sum(dominating.alter.egos.deg), ego.params$max.dom.ego)

  ego.params$max <- FALSE
  egos.left <- clone(dominating.alter.egos.deg)
  alter.egos <- matrix(0,1,graph.properties$n)

  treat <- which(treatment.assignments==1)
  ctrl <- which(treatment.assignments==0)
  all.selected <- list()

  total.egos <- sum(egos.left)
  # cycle through increasing numbers of alter egos
  while( (sum(egos.left)>=2 | total.egos <= 10) & length(all.selected) <= graph.properties$n/4) { 
    egos <- which(egos.left==1)
    treat <- which(treatment.assignments==1)
    ctrl <- which(treatment.assignments==0)
    treatment.egos <- intersect(egos, treat) 
    control.egos <- intersect(egos, ctrl) 
     
    if((sum(treatment.egos)==0 | sum(control.egos)==0) & ego.params$setting=="greedy"){ # check if there's no nodes in treatment or control
      while(sum(treatment.egos)==0 | sum(control.egos)==0){
	      egos.left <- rep(0,graph.properties$n)
        egos.left[which(clone(dominating.alter.egos.deg)==0)] <- 1
        egos <- which(egos.left==1)
        treat <- which(treatment.assignments==1)
        ctrl <- which(treatment.assignments==0)
        treatment.egos <- intersect(egos, treat)
        control.egos <- intersect(egos, ctrl)
	      total.egos <- total.egos + sum(egos.left)
      }
    }
    selected <- select.alter.egos(egos.left, treatment.assignments, ego.params$setting)
    alter.egos[selected] <- 1
    egos.left[selected] <- 0
    all.selected <- append(all.selected, list(selected))

    ego.params$max.dom.ego <- ego.params$max.dom.ego-1
    ego.params$num.ego <- sum(alter.egos)

    if(graph.params$graph.type == "facebook"){
      if(length(all.selected) %% 20 == 0){
        bias.behavior <- calculate.ATE.various(length(all.selected), graph.properties, alter.egos, outcome.params, ego.params, treatment.assignments, stochastic.vars, bias.behavior, selected=all.selected)
      }
    }
    else{
      bias.behavior <- calculate.ATE.various(length(all.selected), graph.properties, alter.egos, outcome.params, ego.params, treatment.assignments, stochastic.vars, bias.behavior, selected=all.selected)
    }
  }

  ego.params$max.dom.ego <- max(sum(dominating.alter.egos.deg), ego.params$max.dom.ego)
  
  bias.behavior$index <- as.numeric(bias.behavior$index)
  bias.behavior$pt.uncovered <- as.numeric(bias.behavior$pt.uncovered)
  bias.behavior$ATE.true <- as.numeric(bias.behavior$ATE.true)
  bias.behavior$ATE.ego.gui <- as.numeric(bias.behavior$ATE.ego.gui)
  
  bias.behavior$pt.covered <- 1 - bias.behavior$pt.uncovered
  bias.behavior$nonego.ATE <- nonego.ATE
  bias.behavior$avg.degree <- avg.degree
  return(bias.behavior)
}
