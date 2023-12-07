source("alter-ego-bias/experiments/experiment-utils.R")
source("alter-ego-bias/experiments/multiple-account-experiment.R")

test.config <- function(idx, setting, configs, trials, all=FALSE) { 
  cat("Running", idx, "\n")
  print(configs[idx,])
  
  results <- data.frame(index=numeric(), size.of.dom=logical(), method=character(), 
                        pt.uncovered=numeric(), alter.ego.influence=numeric(), ATE.true=numeric(), 
                        variable=numeric(), value=numeric(), pt.covered=numeric(), n=numeric(), 
                        graph.type=character(), degree=numeric(), p=numeric(), 
                        mu=numeric(), ncoms=numeric(), maxc=numeric(), minc=numeric(), 
                        lambda_0=numeric(), lambda_1=numeric(), lambda_2=numeric(), stringsAsFactors=FALSE)
  
  graph.params <- build.graph.params(configs, idx)
  alter.ego.params <- list()
  alter.ego.params$model <- reduction.ego.model
  alter.ego.params$all <- all
  alter.ego.params$setting <- setting
  alter.ego.params$weighting <- "inf"
  outcome.params <- build.outcome.params(configs[idx,"lambda_0"], configs[idx,"lambda_1"], configs[idx,"lambda_2"], configs[idx,"sd.noise"])
  clustering <- "infomap"
  
  for(i in 1:trials) {
    graph.params$ind <- i
    
    cat("trial", i, "\n")
    bias.behavior.ATE <- multiple.account.experiment(graph.params, clustering, alter.ego.params, outcome.params, alter.ego.params$setting)
    bias.behavior.ATE$alter.ego.influence <- as.numeric(bias.behavior.ATE$alter.ego.influence)
    bias.behavior.ATE$gui.beta <- as.numeric(bias.behavior.ATE$gui.beta)
    bias.behavior.ATE$gui.gamma <- as.numeric(bias.behavior.ATE$gui.gamma)
    
    bias.behavior.ATE <- add.graph.params(bias.behavior.ATE, graph.params)
    bias.behavior.ATE <- add.outcome.params(bias.behavior.ATE, outcome.params)
    bias.behavior.ATE$graph.id <- configs[idx,"graph.no"]
    bias.behavior.ATE$ego.bias <- bias.behavior.ATE$nonego.ATE - bias.behavior.ATE$ATE.ego.gui
    
    results <- rbind(results, bias.behavior.ATE)
    write.csv(results, paste0("alter-ego-bias/results/new-",setting,"-results-", graph.params$graph.type, "-", outcome.params["lambda_1"], "-", outcome.params["lambda_2"], "-", i, ".csv"))
  }
}

args <- commandArgs(trailingOnly = TRUE)
idx <- as.integer(args[1])
configs <- read.csv("alter-ego-bias/experiments/configs/all_ego_configurations.csv")
setting <- args[2]
trials <- as.integer(args[3])
test.config(idx, setting, configs, trials, FALSE)
