setwd(this.path::here())

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

level <- 0.05

methods <- clean_name(colnames(timing))
nr_rows <- nrow(timing)
table <- data.frame(matrix(nrow = nr_rows * length(mus), ncol = 3 + length(methods)))
colnames(table) <- c("mu", "N", "T", methods)
mus <- seq(from = -0.25, to = 0.25, length.out = 51)
for(m in 1:length(mus)){
  load(paste0("dgp2_power_curve_mu", m, ".RData"))
  table[(m - 1) * nr_rows + 1:nr_rows, 3 + 1:length(methods)] <- -1*apply(X = timing, MARGIN = 1:2, FUN = mean)
  table[(m - 1) * nr_rows + 1:nr_rows, 1:3] <- extract_pars(rownames(timing))
}
table[,1] <- rep(mus, each = nr_rows)
df <- melt(table, id.vars = 1:3)
#df$mu <- factor(df$mu, levels = unique(df$mu))
df$N <- factor(df$N, levels = unique(df$N))
df$T <- paste0("T = ",df$T)
df$T <- factor(df$T, levels = unique(df$T))
colnames(df)[4:5] <- c("Method", "Runtime")

ggplot(df, aes(x = mu, y = Runtime, color = Method)) +
  #geom_col(position = "dodge") +
  #scale_fill_brewer(palette = "Paired") +
  geom_line() +
  scale_color_manual(values = method_colors)+
  facet_grid(vars(N), vars(T), scales = "free_y") +
  theme(legend.position = "bottom") +
  ylab("Runtime (seconds)")
  #geom_hline(data = data.frame(size = level, size_power = factor("Size", levels = unique(df$mu))),
  #           aes(yintercept = size)) +
  #scale_y_continuous(limits = c(0, 1))
#facetted_pos_scales(
#  y = list(
#    size_power == "Size" ~ scale_y_continuous(limits = c(0, 0.2))
#  )
#)
ggsave(paste0("power_curve_runtimes.pdf"), width = 4 * length(unique(table$T)), height = 4 * length(unique(table$N)))
