## First specify the packages of interest
packages <- c("BoolNet", "ggplot2", "wesanderson", "forcats")

## Now load or install & load all
package_check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
    }
    library(x, character.only = TRUE)
  }
)

# A function that checks if an attractor has
# certain characteristics
# attr: the attractor
# conditions: a vector of integers (0s and 1s)
# names(conditions) contains the names of the variables 
# that have to remain active or inactive.

verify_attractor <- function(attr, conditions){
	for (i in rownames(attr)){
    state <- attr[i, ]
    for (n in names(conditions)){
      if(state[, n] != conditions[n]){
        return(FALSE)
      }
    }
	}
	return(TRUE)
}

show_nodes <- function(boolean_vector){
  cnames <- names(boolean_vector)
  s <- paste(cnames[boolean_vector], collapse=", ")
  print(noquote(s))
}

show_nodes_df <- function(boolean_vector){
  cnames <- colnames(boolean_vector)
  s <- paste(cnames[boolean_vector], collapse=", ")
  print(noquote(s))
}

sum_attractors <- function(group, network, attr_info){
  sum <- generateState(network, 
                specs = c(JPs = 0), 
                default = 0)
  n <- 0
  for(i in group){
    atri <- getAttractorSequence(attr_info, i)
    n  <- n + nrow(atri) 
    sum <- sum + colSums(atri)
  }
  sum <- sum / n
  return(sum)
}

analyze_attractors <- function(net, attr_info){
  behaviors <- data.frame(matrix(ncol = 8, nrow = 0))
  fields <- c(
    "i",
    "order",
    "jps",
    "bf", 
    "os", 
    "gd", 
    "ef",
    "pef"
  )
  for(i in 1:length(attr_info$attractors)){
    atri <- getAttractorSequence(attr_info, i)
    tsi <- sorting_hat(atri)
    order <- nrow(atri)
    obs <- c(
      "i" = i,
      "order" = order,
      "jps" = tsi$jp,
      "bf" = tsi$bf, 
      "os" = tsi$os, 
      "gd" = tsi$gd, 
      "ef" = tsi$ef,
      "pef" = tsi$pef
    )
    behaviors <- rbind(behaviors, obs)
  }
  colnames(behaviors) <- fields
  return(behaviors) 
}

# 
sorting_hat <- function(atri){
  ts <- list(jps = 0, bf = 0, os = 0, gd = 0, ef = 0, pef = 0)
  # MLC activity can even increase endothelial barrier function
  jp <- verify_attractor(atri, c(JPs = 1))
  ts$jps <- jp == 1
  spread <- verify_attractor(atri, c(MLC=0))
  ts$bf <- jp && spread 
  superoxide <- verify_attractor(atri, c(Superoxide=1))
  peroxynitrite <- verify_attractor(atri, c(ONOO=1))
  ts$os <- superoxide
  e <- verify_attractor(atri, c(ETS=1))
  ts$gd <- ts$os && e
  no_b1 <- verify_attractor(atri, c(NO_b1 = 1))
  no_b2 <- verify_attractor(atri, c(NO_b2 = 1))
  pgi2 <- verify_attractor(atri, c(PGI2=1))
  txa2 <- verify_attractor(atri, c(TXA2=1))
  ts$ef <- (pgi2 || !txa2) && (no_b1 || no_b2)
  ts$pef <- no_b1 || no_b2 
  return(ts)
}

estimate_trapspace_size <- function(net, initial_states){
  ts <- list(bf = 0, os = 0, gd = 0, ef = 0)
  for(i in 1:initial_states){
    # This creates a vector where each gene is randomly assigned 0 or 1
    randomState <- sample(0:1, length(net$genes), replace = TRUE)
    # Ensure the names are correct
    names(randomState) <- net$genes
    # Respect fixed gene states
    for(g in net$genes){
      if(net$fixed[g] >= 0){
        randomState[g] <- net$fixed[g]
      }
    }
    # Find the attractor
    attr_info <- getAttractors(
      net,
      type = "synchronous",
      method = "chosen",
      startStates = list(randomState)
    )
    atri <- getAttractorSequence(attr_info, 1)
    tsi <- sorting_hat(atri)
    ts$bf <- ts$bf + tsi$bf
    ts$os <- ts$os + tsi$os
    ts$gd <- ts$gd + tsi$gd
    ts$ef <- ts$ef + tsi$ef
  }
  ts$bf <- ts$bf/initial_states
  ts$os <- ts$os/initial_states
  ts$gd <- ts$gd/initial_states
  ts$ef <- ts$ef/initial_states
  return(ts)
}

analyze_net <- function(net, initial_states, reps){
  tssvs <- data.frame(matrix(ncol = 4, nrow = 0))
  for(i in 1:reps){
    p_tss <- estimate_trapspace_size(net, initial_states)
    p_tssv <- c(bf = p_tss$bf, os = p_tss$os, gd = p_tss$gd, ef = p_tss$ef)
    tssvs <- rbind(tssvs, p_tssv)
  }
  colnames(tssvs) <- c('bf', 'os', 'gd', 'ef')
  tssv <- apply(tssvs, 2, mean)
  tssv_sd <- apply(tssvs, 2, sd)
  return(list(tssv = tssv, tssv_sd = tssv_sd))
}

plot_mutant_behavior <- function(net, nodes, values, md, tssv, tssv_sd, save = FALSE){
  m_net <- fixGenes(net, nodes, values)
  na <- analyze_net(m_net, initial_states, reps)
  m_tssv <- na$tssv
  m_tssv_sd <- na$tssv_sd
  # Make the plot
  observation <- c('bf', 'os', 'gd', 'ef', 'bf', 'os', 'gd', 'ef')
  level <- c('Control', 'Control', 'Control', 'Control', md, md, md, md)
  level <- factor(level, levels = c('Control', md))
  fraction <- c(tssv, m_tssv)
  sd <- c(tssv_sd, m_tssv_sd)
  data <- data.frame(observation, level, fraction, sd)
  m_barplot <- ggplot(data = data, aes(x=observation, y=fraction, fill = level)) +
    geom_bar(stat="identity" , position=position_dodge()) +
    scale_y_continuous(limit = c(0, 1)) +
    geom_errorbar(aes(ymin=fraction-sd, ymax=fraction+sd), width=.2,
                  position=position_dodge(.9)) +
    theme_linedraw(base_size = 15) +
    scale_fill_manual(values=wes_palette(n=2, name="Cavalcanti1"))
  file_path <- paste('Figures\\Mutants\\', md, '.jpg', sep='')
  if(save){
    ggsave(file_path, width = 5, height = 3, dpi = 300, units = "in")
  }
  return(list(tssv = m_tssv, tssv_sd = m_tssv_sd))
}

plot_edge_behavior <- function(m_net, md, tssv, tssv_sd, save = TRUE, initial_states=1000, reps=10){
  na <- analyze_net(m_net, initial_states, reps)
  m_tssv <- na$tssv
  m_tssv_sd <- na$tssv_sd
  # Make the plot
  observation <- c('bf', 'os', 'gd', 'ef', 'bf', 'os', 'gd', 'ef')
  level <- c('Control', 'Control', 'Control', 'Control', md, md, md, md)
  level <- factor(level, levels = c('Control', md))
  fraction <- c(tssv, m_tssv)
  sd <- c(tssv_sd, m_tssv_sd)
  data <- data.frame(observation, level, fraction, sd)
  m_barplot <- ggplot(data = data, aes(x=observation, y=fraction, fill = level)) +
    geom_bar(stat="identity" , position=position_dodge()) +
    scale_y_continuous(limit = c(0, 1)) +
    geom_errorbar(aes(ymin=fraction-sd, ymax=fraction+sd), width=.2,
                  position=position_dodge(.9)) +
    theme_linedraw(base_size = 15) +
    scale_fill_manual(values=wes_palette(n=2, name="Cavalcanti1"))
  file_path <- paste('Figures\\Edges\\', md, '.jpg', sep='')
  if(save){
    ggsave(file_path, width = 5, height = 3, dpi = 300, units = "in")
  }
  return(list(tssv = m_tssv, tssv_sd = m_tssv_sd))
}

# The followin function expects the results to be a list
# the first element of the list are the results for the wild level
# the second element contains the results for a mutant with node set to 0
# the third element contains the results for a mutant with node set to 1
# the fourth element, if present, contains the results for a mutant with node set to 2
# node is the name of the node in the network
plot_node_behavior <- function(node, results){
  n <- length(results)
  observation <- rep(c('bf', 'os', 'gd', 'ef'),times=n)
  level <- c('Control')
  for(i in 0:(n-2)){
    level <- c(level, paste(node, i, sep='_'))
  }
  desired_order <- level
  level <- rep(level, each=length(c('bf', 'os', 'gd', 'ef')))
  level <- factor(level, levels = desired_order)
  fraction <- c()
  sd <- c()
  for(i in 1:n){
    ri <- results[[i]]
    fraction <- c(fraction, ri$tssv)
    sd <- c(sd, ri$tssv_sd)
  }
  data <- data.frame(observation, level, fraction, sd)
 # data <- data %>%
 #   mutate(level = factor(level, levels = desired_order))
  m_barplot <- ggplot(data = data, aes(x=observation, y=fraction, fill = level)) +
    geom_bar(stat="identity" , position=position_dodge()) +
    scale_y_continuous(limit = c(0, 1)) +
    geom_errorbar(aes(ymin=fraction-sd, ymax=fraction+sd), width=.2,
                  position=position_dodge(.9)) +
    theme_linedraw(base_size = 15) +
    scale_fill_manual(values=wes_palette(n=n, name="Cavalcanti1"))
  file_path <- paste('Figures\\', node, '.jpg', sep='')
  ggsave(file_path, width = 5, height = 3, dpi = 300, units = "in")
  return(data)
}

perturb_attractor <- function(
    attr_info, # getAttractors(...)
    shine, # a network
    i_attr, # number of the attractor in attr_info
    max_flipped, # the maximum number of molecule activities affected
    reps)
{
  landscape <- list()
  init <- getAttractorSequence(attr_info, i_attr)
  landscape$initial <- init[1,]
  genes <- attr_info$stateInfo$genes
  ngenes <- length(genes)
  set.seed(123)
  results <- data.frame(matrix(ncol = 6, nrow = 0))
  fields <- c(
    "modified",
    "n",
    "jps",
    "bf", 
    "os", 
    "gd", 
    "ef",
    "pef"
  )
  for(i in 1:reps){
    # Make a copy of the initial state
    cp <- init[1,]
    # Random number between 1 and max_flipped
    nf <- floor(runif(1, min = 2, max = max_flipped + 1))
    # Which molecules to change
    chosen <- sample(1:ngenes, size = nf, replace = FALSE)
    chosen <- sort(chosen)
    sg <- genes[chosen]
    for(m in sg){
      cp[m] <- (cp[m] + 1) %% 2
    }
    startstate <- unlist(cp[1, ])
    attr_qatr1 <- getAttractors(
      shine,
      type = "synchronous",
      method = "chosen",
      startStates = list(startstate)
    )
    atri <- getAttractorSequence(attr_qatr1, 1)
    tsi <- sorting_hat(atri) # this classifies the reached attractor
    psg <- paste(sg, collapse = ', ') # A string with the molecules that will be affected
    obs <- c(
      "modified" = psg,
      "n" = nf,
      "jps" = tsi$jps,
      "bf" = tsi$bf, 
      "os" = tsi$os, 
      "gd" = tsi$gd, 
      "ef" = tsi$ef,
      "pef" = tsi$pef
    )
    results <- rbind(results, obs)
  }
  colnames(results) <- fields
  results$n <- as.numeric(results$n)
  results$jps <- as.logical(results$jps)
  results$bf <- as.logical(results$bf)
  results$os <- as.logical(results$os)
  results$gd <- as.logical(results$gd)
  results$ef <- as.logical(results$ef)
  results$pef <- as.logical(results$pef)
  landscape$results <- results
  return(landscape)
}

combinations_of_n_perturbations <- function(
    attr_info, # getAttractors(...)
    shine, # a network
    attrs, # list of attractors in attr_info (numeric)
    flipped, # the number of molecule activities affected smaller than 7
    tb # list that contains target bf os gd and ef activity
    )
{
  # Generate a  the possible combinations
  genes <- attr_info$stateInfo$genes
  ngenes <- length(genes)
  combinations <- combn(x = ngenes, m = flipped)
  scores <- c()
  # Iterate through the combinatons
  for(i in 1:dim(combinations)[2]){
    comb_i <- combinations[,i]
    score_i <- 0
    sg <- genes[comb_i] # a vector of gene names
    psg <- paste(sg, collapse = ', ') # A string with the molecules that will be affected
    for(atr in attrs){
      init <- getAttractorSequence(attr_info, atr)
      cp <- init[1,]
      for(m in sg){
        cp[m] <- (cp[m] + 1) %% 2
      }
      startstate <- unlist(cp[1, ])
      attr_qatr1 <- getAttractors(
        shine,
        type = "synchronous",
        method = "chosen",
        startStates = list(startstate)
      )
      atri <- getAttractorSequence(attr_qatr1, 1)
      tsi <- sorting_hat(atri) # this classifies the reached attractor
      if(tsi$bf == tb$bf & tsi$os == tb$os & tsi$gd == tb$gd & tsi$ef == tb$ef){
        score_i <- score_i + 1
      }
    }
    scores[psg] <- score_i
  }
  return(scores)
}
