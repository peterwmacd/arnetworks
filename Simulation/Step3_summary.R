res_list <- list()
p = 50
n = 100
for (ix1 in c(1)){
  load( paste0("inputs", "_ix", ix1, "_p_", p,"_n_", n,".RData"))
  
  out = readRDS(paste0("out_all_ix",ix1, ".rds"))
  
  analyze_single <- function(out, extract_fun) {
    results <- lapply(out, function(sim_res) {
      lapply(sim_res[1:9], extract_fun)  # assuming 9 initializations per replicate
    })
    results_flat <- do.call(c, lapply(results, function(x) x))
    results_matrix <- na.omit(do.call(rbind, results_flat))
    
    res_vec <- numeric(ncol(results_matrix) * 3)
    idx <- 1
    for (j in 1:ncol(results_matrix)) {
      sd_val <- sd(results_matrix[, j], na.rm = TRUE)
      res_vec[idx]     <- mean(abs(results_matrix[, j]) >= qnorm(0.975))
      res_vec[idx + 1] <- 1 / sd_val
      res_vec[idx + 2] <- mean((1 / sd_val) * abs(results_matrix[, j]) >= qnorm(0.975))
      idx <- idx + 3
    }
    return(res_vec)
  }
  
  # ---------- Initial Estimation Error ----------
  # --- xi ---
  xi_per_rep <- sapply(out, function(sim_res) {
    rmaes <- sapply(sim_res[1:9], function(res) if (!is.null(res)) mean(abs(res$xi_init - xi) / xi[1]) else NA)
    mean(rmaes, na.rm = TRUE)
  })
  xi_rmae_mean <- mean(xi_per_rep, na.rm = TRUE)
  xi_rmae_sd <- sd(xi_per_rep, na.rm = TRUE)
  
  # --- Eta ---
  eta_per_rep <- sapply(out, function(sim_res) {
    rmaes <- sapply(sim_res[1:9], function(res) if (!is.null(res)) mean(abs(res$eta_init - eta) / eta[1]) else NA)
    mean(rmaes, na.rm = TRUE)
  })
  eta_rmae_mean <- mean(eta_per_rep, na.rm = TRUE)
  eta_rmae_sd <- sd(eta_per_rep, na.rm = TRUE)
  
  # --- a ---
  a_per_rep <- sapply(out, function(sim_res) {
    rmaes <- sapply(sim_res[1:9], function(res) if (!is.null(res)) abs(res$gVal_init[1] - a) / a else NA)
    mean(rmaes, na.rm = TRUE)
  })
  a_rmae_mean <- mean(a_per_rep, na.rm = TRUE)
  a_rmae_sd <- sd(a_per_rep, na.rm = TRUE)
  
  # --- b ---
  b_per_rep <- sapply(out, function(sim_res) {
    rmaes <- sapply(sim_res[1:9], function(res) if (!is.null(res)) abs(res$gVal_init[2] - b) / b else NA)
    mean(rmaes, na.rm = TRUE)
  })
  b_rmae_mean <- mean(b_per_rep, na.rm = TRUE)
  b_rmae_sd <- sd(b_per_rep, na.rm = TRUE)
  
  
  
  

  
  
  # ---------- Final Estimation Error ----------
  xi_per_rep <- sapply(out, function(sim_res) {
    rmaes <- sapply(sim_res[1:9], function(res) if (!is.null(res)) mean(abs(res$xi - xi) / xi[1]) else NA)
    mean(rmaes, na.rm = TRUE)
  })
  xi_rmae_mean2 <- mean(xi_per_rep, na.rm = TRUE)
  xi_rmae_sd2 <- sd(xi_per_rep, na.rm = TRUE)
  
  # --- Eta ---
  eta_per_rep <- sapply(out, function(sim_res) {
    rmaes <- sapply(sim_res[1:9], function(res) if (!is.null(res)) mean(abs(res$eta - eta) / eta[1]) else NA)
    mean(rmaes, na.rm = TRUE)
  })
  eta_rmae_mean2 <- mean(eta_per_rep, na.rm = TRUE)
  eta_rmae_sd2 <- sd(eta_per_rep, na.rm = TRUE)
  
  # --- a ---
  a_per_rep <- sapply(out, function(sim_res) {
    rmaes <- sapply(sim_res[1:9], function(res) if (!is.null(res)) abs(res$gVal[1] - a) / a else NA)
    mean(rmaes, na.rm = TRUE)
  })
  a_rmae_mean2 <- mean(a_per_rep, na.rm = TRUE)
  a_rmae_sd2 <- sd(a_per_rep, na.rm = TRUE)
  
  # --- b ---
  b_per_rep <- sapply(out, function(sim_res) {
    rmaes <- sapply(sim_res[1:9], function(res) if (!is.null(res)) abs(res$gVal[2] - b) / b else NA)
    mean(rmaes, na.rm = TRUE)
  })
  b_rmae_mean2 <- mean(b_per_rep, na.rm = TRUE)
  b_rmae_sd2 <- sd(b_per_rep, na.rm = TRUE)
  
  
  # ---------- Inference via analyze_single ----------
  i = p
  xi_vec <- analyze_single(out, function(res) {
    1/ res$se_estimates[2 + i] * (res$xi[i] - xi[i])
  })
  eta_vec<- analyze_single(out, function(res) {
    1/ res$se_estimates[2 + p + i] * (res$eta[i] - eta[i])
  })
  
  ab_vec <- analyze_single(out, function(res) {
    1/ res$se_estimates[1:2] * (res$gVal - c(a, b))
  })
  
  # ---------- Store rows ----------
  res_list[[length(res_list)+1]] <- data.frame(
    param = "xi", 
    initial = xi_rmae_mean,
    sd = xi_rmae_sd,
    est = xi_rmae_mean2,
    sd2 = xi_rmae_sd2,
    cv = 1-xi_vec[1]
  )
  res_list[[length(res_list)+1]] <- data.frame(
    param = "eta", 
    initial = eta_rmae_mean,
    sd = eta_rmae_sd,
    est = eta_rmae_mean2,
    sd2 = eta_rmae_sd2,
    cv = 1-eta_vec[1]
  )
  res_list[[length(res_list)+1]] <- data.frame(
    param = "a", 
    initial = a_rmae_mean,
    sd = a_rmae_sd,
    est = a_rmae_mean2,
    sd2 = a_rmae_sd2,
    cv = 1- ab_vec[1]
  )
  res_list[[length(res_list)+1]] <- data.frame(
    param = "b", 
    initial = b_rmae_mean,
    sd = b_rmae_sd,
    est = b_rmae_mean2,
    sd2 = b_rmae_sd2,
    cv = 1- ab_vec[4]
  )
  
  # ---------- Final table ----------
  resTable_long <- do.call(rbind, res_list)
  
}

library(dplyr)

resTable <- resTable_long %>%
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
  mutate(
    sd  = ifelse(sd  < 0.001, 0.001, sd),
    sd2 = ifelse(sd2 < 0.001, 0.001, sd2)
  )
write.csv(resTable, file = paste0("resTable_p_",p,"_n_",n,".csv"))

