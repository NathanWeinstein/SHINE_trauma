# The following are the scripts used to analyze our boolean model of 
# traumatic-shoch induced endotheliopathy, it is recommended to run them 
# in the order they appear below, if you run all of them, it will take 
# a few days for the script to finish running. shine_functions.R contains 
# most of the reusable tools used for this project. If you use an 
# operating system that is ot Windows you may need to modify the filepaths

# 1) find and classify the attractors "analyze_attractors.R":
# a) Loads the model in BoolNet format from "nets\\levels.bnet"
# b) The figures it creates are stored at the "Figures\\main" folder
# c) Stores the network, the attractors as the "attractors" data.frame
# and other data at "trap_spaces.RData"
# Caution! The order of the attractors can change from one run to the next.

source("analyze_attractors.R")

# 2) Explore the quiescent attractors robustness to noise

source("analyze_q_transitions.R")

# 3) Explore the SHINE attractors robustness to noise

source("analyze_shine_transitions.R")

# 4) Explore the transitions in other groups of attractors

source("analyze_other_transitions.R")

# 5) Analyze the effect of loss- and constitutive activity experiments 
# for each node in the network.
# Stores the resulting plots at "Figures\\Mutants"

source("all_mutants.R")

# 6) Analyze the effect of removing each edge in the network
# the modified networks can be found in the folder "nets"
# Stores the resulting plots at "Figures\\Edges"
# Stores the result at the "edges" data.frame"
# Saves the edges data.frame as "edges.csv"

source("analyze_edge_effect.R")

# 7) Analyze all the attractors using a SQUAD based ODE version of 
# the model

source("levels_ode.R")