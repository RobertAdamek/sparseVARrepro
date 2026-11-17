library(ggplot2)
library(ggh4x)
library(reshape2)
setwd(this.path::here())

sim_result_files <- list.files("../Results")
sim_result_files <- sim_result_files[sim_result_files!="power_curve_results"]
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

color_switch <- c("VAR-TF" = "steelblue", "GPI-TF" = "springgreen4")
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
    mean_runtimes <- -1*apply(timing, c(1,2), mean)
    colnames(mean_runtimes) <- c(methods)
    table[(i - 1) * nr_rows + 1:nr_rows, 3 + 1:length(which_methods_plot)] <- mean_runtimes[,which_methods_plot]
    table[(i - 1) * nr_rows + 1:nr_rows, 1:3] <- extract_pars(rownames(mean_runtimes))
  }
  table$Relative <- table$`GPI-TF` / table$`VAR-TF`
  table$`GPI-TF` <- table$`VAR-TF` <- NULL
  table$which_faster <- ifelse(table$Relative >= 1, "VAR-TF", "GPI-TF")

  df <- melt(table, id.vars = c(1:3, 5))
  df$size_power <- factor(df$size_power, levels = unique(df$size_power))
  df$N <- factor(df$N, levels = unique(df$N))
  df$T <- factor(df$T, levels = unique(df$T))
  colnames(df)[5:6] <- c("Method", "Runtime")
  y_facet_size <- 9*18/max(nchar(as.character(df$size_power)))
  ggplot(df, aes(x = T, y = Runtime, group = 1, colour = which_faster)) +
    scale_colour_manual(values = color_switch)+
    facet_grid(vars(size_power), vars(N), scales = "free_y") +
    #theme(strip.text.y = element_text(size = y_facet_size)) +
    #theme(legend.position = "bottom") +
    #guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(legend.position = "none") +
    geom_hline(aes(yintercept = 1)) +
    geom_point() +
    ylab(expression("Relative runtime  "*over("GPI-TF", "VAR-TF"))) +
    scale_y_continuous(limits = c(0.9, 1.45))
    #facetted_pos_scales(
    #  y = list(
    #    size_power == "Size" ~ scale_y_continuous(limits = c(0, 0.2))
    #  )
    #)
  ggsave(paste0("../Results/", dgp, "_runtime.pdf"), width = 8, height = 2 * length(unique(table$size_power)))
}
