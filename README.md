# The Effect of Alter Ego Accounts on A/B Tests for Social Networks
This repo contains the code for our work on A/B testing in networks with non-cooperative behavior. 

## Run the Experiments

[REPLICATION.md](REPLICATION.md) contains installation, configuration, and run instructions for experiments in the paper.

## Overview
Given a network graph and simulation parameters identified in a configuration file, 
* Cluster-randomized treatment assignment
* Non-cooperative participant network set construction
* Outcome simulation per n alter egos up to dominating set
* Model fit and outcome estimation with linear estimator
* Bias calculation in sample

#### Experiment Configs
Options: 
* Network structure for A/B test
* Alter ego network placement (random, greedy)
* Outcome simulation parameters

More information is provided in  [this README](experiments/configs/README.md).

#### Network Graph Resources
We provide some existing graph structures for the experiments:
* **Synthetic Corpus**:
We provide a corpus of graph types used for the simulation studies in the paper, generated under the same set of parameters. 
* **Real-World**:
We use networks released in the SNAP Library. 

Alternatively, run the scripts used to generate graphs in the corpus to generate new graphs. 
More details at [this README](REPLICATION.md).

