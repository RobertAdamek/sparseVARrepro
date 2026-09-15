library(ggplot2)
library(ggh4x)
library(reshape2)
setwd(this.path::here())

sim_result_files <- list.files("../Results")
dgps <- unique(substr(sim_result_files, 1, 4))

level <- 0.05

clean_name<-function(meth){
  meth<-gsub(pattern="-L1-unpen-own", replacement="", meth)
  meth<-gsub(pattern="VAR-GP-unpen-own", replacement="GPI", meth)
  meth<-gsub(pattern="-11", replacement="", meth)
  meth<-gsub(pattern="04", replacement="0.4", meth)
  meth<-gsub(pattern="08", replacement="0.8", meth)
  return(meth)
}

extract_pars <- function(x) {
  N <- substr(x, regexpr("N = ", x), regexpr("n = ", x) - 3)
  n <- substr(x, regexpr("n = ", x) + 4, regexpr("prop = ", x) - 3)
  prop <- substr(x, regexpr("prop = ", x) + 7, regexpr("mu = ", x) - 3)
  mu <- substr(x, regexpr("mu = ", x) + 5, regexpr(")", x) - 1)
  sp <- ifelse(mu == 0, "Size", paste0("Pw: µ = ", mu, ", p = ", prop))
  return(c(sp = sp, N = N, n = n))
}

for (dgp in dgps) {
  results <- sim_result_files[substr(sim_result_files, 1, 4) == dgp & grepl(".RData", sim_result_files)]
  nr_results <- length(results)
  results <- c(results[nr_results], results[-nr_results])
  load(paste0("../Results/", results[1]))
  methods <- clean_name(colnames(reject))
  nr_rows <- nrow(reject)
  table <- data.frame(matrix(nrow = nr_rows * nr_results, ncol = 3 + length(methods)))
  colnames(table) <- c("size_power", "N", "T", methods)
  for (i in 1:nr_results) {
    load(paste0("../Results/", results[i]))
    table[(i - 1) * nr_rows + 1:nr_rows, 3 + 1:length(methods)] <- reject[, , dimnames(reject)[[3]] == as.character(level)]
    table[(i - 1) * nr_rows + 1:nr_rows, 1:3] <- extract_pars(rownames(reject))
  }
  
  df <- melt(table, id.vars = 1:3)
  df$size_power <- factor(df$size_power, levels = unique(df$size_power))
  df$N <- factor(df$N, levels = unique(df$N))
  df$T <- factor(df$T, levels = unique(df$T))
  colnames(df)[4:5] <- c("Method", "Rejection Frequency")
  
  ggplot(df, aes(x = T, y = `Rejection Frequency`, fill = Method)) +
    geom_col(position = "dodge") +
    scale_fill_brewer(palette = "Paired") +
    facet_grid(vars(size_power), vars(N), scales = "free_y") +
    theme(legend.position = "bottom") +
    geom_hline(data = data.frame(size = level, size_power = factor("Size", levels = unique(df$size_power))),
      aes(yintercept = size)) +
    scale_y_continuous(limits = c(0, 1)) +
    facetted_pos_scales(
      y = list(
        size_power == "Size" ~ scale_y_continuous(limits = c(0, 0.2))
      )
    )
  ggsave(paste0("../Results/", dgp, "_size_power.pdf"), width = 16, height = 4 * length(unique(table$N)))
}
