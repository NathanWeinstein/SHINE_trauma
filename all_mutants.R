start_time <- Sys.time()
if (!dir.exists("Figures\\Mutants")) {
  dir.create("Figures\\Mutants", recursive = TRUE)
}
source("analyze_mutation_effct.R")
source("analyze_mutation_effct_q.R")
source("analyze_mutation_effct_context.R")
source("analyze_mutation_effct_T.R")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)