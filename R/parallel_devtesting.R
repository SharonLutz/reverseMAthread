#' @include multi_process_mediate.R

test_single_core <- function(n_times){
  tictoc::tic(paste("just one env, direct function call", n_times, "times"))
  for(i in 1:n_times){
    simulate_and_mediate()
  }
  tictoc::toc()
}

test_n_cores_k_jobs <- function(n, k){
  g_env = globalenv()
  med_vars = MediateDataGenerationParameters(nSim=k)
  sim_args = generateDataMatrix(med_vars)
  g_env[["med_vars"]] = med_vars
  g_env[["sim_args"]] = sim_args
  options(mediate.cores = n)
  tictoc::tic(paste("all_mediate_jobs", "n", n, "k", k, "matrix length:", length(sim_args)))
  result_list = mediate_parallel(sim_args, nSimImai=10000)
  attr(result_list, "dim") <- dim(sim_args)
  tictoc::toc()
  return(result_list)
}

run_test = function(cores_to_test = 3){
  #test_single_core(cores_to_test)
  result = test_n_cores_k_jobs(cores_to_test, cores_to_test)
  #test_n_cores_k_jobs(cores_to_test, cores_to_test*5)
}
