library(plyr)
library(dplyr)
library(expm)
library(reshape2)
library(bit)

source("alter-ego-bias/graphs/graph-utils.R")

build.outcome.params <- function(lambda_0, lambda_1, lambda_2, sd.noise) { 
  outcome.params <- list()
  outcome.params$lambda_0 <- lambda_0
  outcome.params$lambda_1 <- lambda_1
  outcome.params$lambda_2 <- lambda_2
  outcome.params$sd.noise <- sd.noise
  
  return(outcome.params)
}

get.stochastic.vars <- function(num, steps, sdnoise, noise) { 
  stochastic.vars <- lapply(1:steps, function(x) { 
    if(noise) return(rnorm(num, 0, sdnoise))
    else return(matrix(0,num,1))
  })
  names(stochastic.vars) <- paste0("t", 1:steps)
  
  return(stochastic.vars)
}

treatment.assignment <- function(g, clusters, prob=0.5) { 
  return(rbinom(length(unique(clusters)), 1, prob))
}

determine.alter.egos <- function(graph.properties, ego.params) {
  alter.egos <- matrix(0, 1, graph.properties$n)
  if(ego.params$setting == "random") { 
    rand.order <- sample(1:graph.properties$n, graph.properties$n, replace = FALSE)
    if(ego.params$max) { 
      idx <- 1
      while(!check.dominating.set(graph.properties, alter.egos)) { 
        alter.egos[rand.order[idx]] <- 1
        idx <- idx + 1
      }  
    }
    else alter.egos[sample(1:graph.properties$n, ego.params$num.ego, replace=FALSE)] <- 1
  }
  if(ego.params$setting == "greedy") {
    dominating.set <- dominate.greedy.inf(graph.properties)
    if(ego.params$max) ego.params$num.ego <- length(dominating.set)
    alter.egos[,sample(dominating.set, ego.params$num.ego)] <- 1
  }
  return(alter.egos)
}


calculate.ATE.various <- function(idx, graph.properties, alter.egos, outcome.params, ego.params, treatment.assignments, stochastic.vars, bias.behavior, selected=NULL) { 
  tryCatch(
    expr = {
      uncovered.vertices <- 1 - alter.egos %*% graph.properties$adj - alter.egos
      pt.uncovered <- sum(uncovered.vertices == 1)/graph.properties$n
    },
    error = function(e){
      print(e)
      print("alter.egos")
      print(alter.egos)
    }
  )
  
  # calculate true outcome without alter egos
  ATE.true <- outcome.params$lambda_1 + outcome.params$lambda_2
  
  # calculate outcome with alter egos
  ego.params <- exposure.probs(ego.params, graph.properties, treatment.assignments, alter.egos)
  outcome.ego <- outcome.model(outcome.params, treatment.assignments, graph.properties, alter.egos, ego.params, stochastic.vars, selected)
  
  #print(treatment.assignments)
  # estimate ATE using the Gui framework
  lm.estimator.gui <- lam.I(graph.properties, treatment.assignments, outcome.ego)
  gui.beta <- lm.estimator.gui$coefficients[2]
  gui.gamma <- lm.estimator.gui$coefficients[3]
  ATE.ego.gui <- gui.beta + gui.gamma

  over.dom.max <- ifelse(ego.params$setting == "greedy", FALSE, ego.params$max.dom.ego < sum(alter.egos))
  if(idx == 0) over.dom.max <- FALSE
  ad.inf <- sum(ego.params$influence.as.ego[which(alter.egos==1)])/graph.properties$n
  
  method <- ifelse(ego.params$setting=="greedy", ego.params$weighting, ego.params$setting)
  if(is.null(ego.params$setting)) method <- "none"
  
  # determine bias in estimate due to alter egos
  bias.behavior[nrow(bias.behavior)+1,] <- c(idx, over.dom.max, method, pt.uncovered, ad.inf, ATE.true, ATE.ego.gui, gui.beta, gui.gamma)  
  return(bias.behavior)
}


lam.I <- function(graph.properties, treatment.assignments, outcome) {
  frac.treated <- as.numeric((graph.properties$adj %*% treatment.assignments) / graph.properties$degrees)
  outcome.model <- lm(outcome ~ treatment.assignments + frac.treated)
  
  return(outcome.model)
}


exposure.probs <- function(ego.params, graph.properties, treatment.assignments, alter.egos, lambda=0.1, p=2) { 
  nonego <- 1 - alter.egos
  treated.nonego <- treatment.assignments * nonego
  control.nonego <- (1 - treatment.assignments) * nonego
  treated.ego <- treatment.assignments * alter.egos
  control.ego <- alter.egos - treated.ego
  
  ego.params$empty <- as.vector(matrix(0, 1, graph.properties$n))
  ego.params$ego.exposure.neighbors <- as.vector(t(graph.properties$adj %*% t(alter.egos) / graph.properties$degrees))
  ego.params$ego.treat.exposure.neighbors <- as.vector(t(graph.properties$adj %*% t(treated.ego) / graph.properties$degrees))
  ego.params$ego.control.exposure.neighbors <- as.vector(t(graph.properties$adj %*% t(control.ego) / graph.properties$degrees))
  ego.params$nonego.treat.exposure.neighbors <- as.vector(t(graph.properties$adj %*% t(treated.nonego) / graph.properties$degrees))
  ego.params$nonego.control.exposure.neighbors <- as.vector(t(graph.properties$adj %*% t(control.nonego) / graph.properties$degrees))
  
  ego.params$treatment.exposure.neighbors <-as.vector( t(graph.properties$adj %*% treatment.assignments / graph.properties$degrees))
  ego.params$influence.as.ego <- colSums(graph.properties$transition)
  return(ego.params)
}

outcome.model <- function(outcome.params, treat, graph.properties, alter.egos, ego.params, stochastic.vars, selected=NULL) { 
  treated.ego <- treat * alter.egos
  control.ego <- alter.egos - treated.ego
  
  out.t0 <- matrix(0, 1, graph.properties$n)
  
  out.t1 <- outcome.params$lambda_0 + outcome.params$lambda_1 * treat + outcome.params$lambda_2 * rowSums(graph.properties$adj %*% diag(as.numeric(out.t0)) / graph.properties$degrees) + stochastic.vars$t1
  if(!is.null(selected)){
      for(sel in selected){
          out.t1[sel[2]] = out.t1[sel[1]]
      }
  }
  
  out.t2 <- outcome.params$lambda_0 + outcome.params$lambda_1 * treat + outcome.params$lambda_2 * rowSums(graph.properties$adj %*% diag(as.numeric(out.t1)) / graph.properties$degrees) + stochastic.vars$t2
  if(!is.null(selected)){
      for(sel in selected){
          out.t2[sel[2]] = out.t2[sel[1]]
      }
  }

  out.t3 <- outcome.params$lambda_0 + outcome.params$lambda_1 * treat + outcome.params$lambda_2 * rowSums(graph.properties$adj %*% diag(as.numeric(out.t2)) / graph.properties$degrees) + stochastic.vars$t3
  if(!is.null(selected)){
      for(sel in selected){
          out.t3[sel[2]] = out.t3[sel[1]]
      }
  }

  return(out.t3) 
}

add.graph.params <- function(bias.behavior, graph.params) { 
  params <- c("n", "graph.type", "degree", "p", "mu", "ncoms", "maxc", "minc")
  
  for(pa in params) { 
    bias.behavior[[pa]] <- ifelse(!is.null(graph.params[[pa]]), graph.params[[pa]], "")  
  }
  
  return(bias.behavior)
}

add.outcome.params <- function(bias.behavior, outcome.params) { 
  params <- c("lambda_0", "lambda_1", "lambda_2")
  
  for(pa in params) { 
    bias.behavior[[pa]] <- outcome.params[[pa]]
  }
  
  return(bias.behavior)
}

reduction.ego.model <- function(treated.ego, control.ego, outcome.params, clusters=NULL) { 
  n <- length(treated.ego)
  out <- matrix(0, 1, n)  
  out[1,which(treated.ego==1)] <- outcome.params$lambda_0 + outcome.params$lambda_1
  out[1,which(control.ego==1)] <- outcome.params$lambda_0 
  
  return(out) 
}



