This markdown provides instructions for running the experiments in The Effect of Alter Ego Accounts on A/B Tests for Social Networks. 

## Setup 
1. Download R for your system: [https://cran.r-project.org/](https://cran.r-project.org/)

2. Clone the repo:
 `git clone git@github.com:KDL-umass/alter-ego-bias.git`

3. Set up R environment. Run install script and verify setup: `Rscript install/install.R`

## Configure Experiment Settings

1. Set up experiment configs. The experiments use the config located in experiments/configs/all_ego_configurations.csv, though you can make your own in the same folder. 

2. Create graphs. 
The forest-fire networks and small-world networks are automatically generated using igraph. SBM adjacency matrices and Facebook edge list are cached in the graphs/ folder.  

The Facebook network edge list can also be downloaded here: https://snap.stanford.edu/data/ego-Facebook.html. 

If you want to generate the SBMs from scratch, you can install the SBM benchmark algorithm scripts created by Fortunato et al. 
* Main page: https://www.santofortunato.net/resources#h.p_u6MEEWAKyhN0 (Download package 1)
* Download link: https://drive.google.com/file/d/0BwGBn8ta6pUrTG5tUGc1N2NjU2c/view?resourcekey=0-bIhuIu6lxYjgOvQMb2Pueg 

## Run the Experiments
Run the experiments using the following command:

`Rscript ./experiments/test-suit.R param1 param2 param3`

param1 is the index of the line in the configuration file to use. param2 the placement strategy, which is either "random" or "greedy" (worst-case). param3 is the number of trials to run. The following example command would run the third line of the config file (not including the header) with greedy placement for 10 trials:

`Rscript ./experiments/test-suit.R 3 greedy 10`