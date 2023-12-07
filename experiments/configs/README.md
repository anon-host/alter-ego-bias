# Configs
Brief description of all experiment config options.
* graph.type -- the type of graph (e.g. forest-fire, sbm, etc.)
* degree -- the degree of the nodes in lattice of a small-world network
* p -- the probability of randomly rewiring an edge in the lattice of a small-world network
* mu -- probability of adding an edge between a node inside and outside a community in an SBM
* fw -- forward-burning probability in a forest-fire model
* bw -- backward-burning probability in a forest-fire model
* size -- number of vertices in the graph
* graph.no -- index

See the "Additive Model of Peer Effects" section in the paper for a more full description of these parameters: 
* lambda_0 -- offset parameter 
* lambda_1 -- individual effect strength
* lambda_2 -- peer effect strength