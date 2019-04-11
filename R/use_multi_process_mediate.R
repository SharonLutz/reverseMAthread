library(tictoc)
library(parallel)
source("multi_process_mediate.R")

create_env_list <- function(num_envs) {
  env_list <- vector("list", num_envs)
  for(i in 1:num_envs){env_list[[i]] = MediateVariables()}
  return (env_list)
}

test_n_cores_k_jobs <- function(n, k){
  num_cores <- n
  num_jobs <- k
  sim_args <- create_env_list(k)
  options(mediate.cores = num_cores)
  
  tic(paste("all_mediate_jobs", "n", n, "k", k))
  result_list <- do_mediate(sim_args)
  toc()
}
##num_cores <- detectCores() - 1
#num_jobs <- (detectCores() - 1)*4

test_single_core <- function(n_times){
  tic(paste("just one env, direct function call", n_times, "times"))
  for(i in 1:n_times){
    simulate_and_mediate(MediateVariables())
  }
  toc()
}

cores_to_test <- 5

test_single_core(cores_to_test)
test_n_cores_k_jobs(cores_to_test, cores_to_test)
#test_n_cores_k_jobs(cores_to_test, cores_to_test*5)
