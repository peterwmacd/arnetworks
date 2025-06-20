tag <- "ix1"
files <- list.files(pattern = paste0("^out_", tag, "_.*\\.rds$"))
all_out <- lapply(files, readRDS)
out_combined <- do.call(c, all_out)
saveRDS(out_combined, file = paste0("out_all_", tag, ".rds"))