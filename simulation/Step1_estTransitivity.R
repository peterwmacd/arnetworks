args <- commandArgs(trailingOnly = TRUE)
start_try <- as.integer(args[1])
end_try <- as.integer(args[2])
tag <- args[3]

p = 50
n = 100
load(paste0("Setting/inputs_", tag, "_p_", p,"_n_", n,".RData"))

library(arnetworks)
library(doParallel)

cl <- makeCluster(8)
registerDoParallel(cl)

out <- foreach(try = start_try:end_try, .errorhandling = "remove") %dopar% {
  set.seed(try + 1000)
  data <- simulateTransitivity(p, n, xi, eta, a, b)
  X <- data$X
  U <- data$U[, , 1:(n - 1)]
  V <- data$V[, , 1:(n - 1)]
  rm(data)
  
  full_res <- vector("list", length = nrow(initset))
  for (ix2 in 1:nrow(initset)) {
    xiE <- rep(initset[ix2, 1], p)
    etaE   <- rep(initset[ix2, 2], p)
    
    full_res[[ix2]] <- estTransitivity(X, U, V, initXi = xiE, initEta = etaE,
                            tauSeq_a = tauSeq_a,
                            tauSeq_b = tauSeq_b,
                            tauSeq_xi = tauSeq_xi,
                            tauSeq_eta = tauSeq_eta,
                            verbose=FALSE, doInference = TRUE)
  }
  full_res
}

stopCluster(cl)
saveRDS(out, file = paste0("out_", tag, "_", start_try, "_", end_try, ".rds"))
