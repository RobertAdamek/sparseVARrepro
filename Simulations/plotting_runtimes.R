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

method_colors <- c("VAR-BIC" = "steelblue1", "VAR-TF"="steelblue",     
                   "MBB" = "sienna1", "BWB" = "sienna3", "DWB"="sienna4",
                   "GPI-BIC" = "springgreen3", "GPI-TF" = "springgreen4",
                   "VAR-oracle" = "purple4")
which_methods_plot <- c("VAR-TF", "GPI-TF")
for (dgp in dgps) {
  results <- sim_result_files[substr(sim_result_files, 1, 4) == dgp & grepl(".RData", sim_result_files)]
  nr_results <- length(results)
  results <- c(results[nr_results], results[-nr_results])
  load(paste0("../Results/", results[1]))
  methods <- clean_name(colnames(timing))
  nr_rows <- nrow(timing)
  table <- data.frame(matrix(nrow = nr_rows * nr_results, ncol = 3 + length(which_methods_plot)))
  colnames(table) <- c("size_power", "N", "T", which_methods_plot)
  for (i in 1:nr_results) {
    load(paste0("../Results/", results[i]))
    median_runtimes <- log(-1*apply(timing, c(1,2), median))
    colnames(median_runtimes) <- c(methods)
    table[(i - 1) * nr_rows + 1:nr_rows, 3 + 1:length(which_methods_plot)] <- median_runtimes[,which_methods_plot]
    table[(i - 1) * nr_rows + 1:nr_rows, 1:3] <- extract_pars(rownames(median_runtimes))
  }
  
  df <- melt(table, id.vars = 1:3)
  df$size_power <- factor(df$size_power, levels = unique(df$size_power))
  df$N <- factor(df$N, levels = unique(df$N))
  df$T <- factor(df$T, levels = unique(df$T))
  colnames(df)[4:5] <- c("Method", "Runtime")
  
  ggplot(df, aes(x = T, y = Runtime, fill = Method)) +
    #geom_line() +
    geom_col(position = "dodge") +
    #scale_fill_brewer(palette = "Paired") +
    scale_fill_manual(values = method_colors)+
    facet_grid(vars(size_power), vars(N), scales = "free_y") +
    theme(legend.position = "bottom") +
    ylab("Median runtime in log(seconds)")
    #geom_hline(data = data.frame(size = level, size_power = factor("Size", levels = unique(df$size_power))),
    #           aes(yintercept = size)) +
    #scale_y_continuous(limits = c(0, 1)) #+
    #facetted_pos_scales(
    #  y = list(
    #    size_power == "Size" ~ scale_y_continuous(limits = c(0, 0.2))
    #  )
    #)
  ggsave(paste0("../Results/", dgp, "_runtime.pdf"), width = 16, height = 4 * length(unique(table$N)))
}
