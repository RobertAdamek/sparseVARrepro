setwd(this.path::here())
setwd("../Results")
library(ggplot2)
library(ggpubr)

sim_file_names <- NULL
for(f in list.files()){
  if(grepl(".RData", f, fixed=TRUE)){
    sim_file_names <- c(sim_file_names, gsub(".RData", "", f, fixed = TRUE))
  }
}
sim_results <- vector("list", length(sim_file_names))
names(sim_results) <- sim_file_names
for(name in sim_file_names){
  load(paste0(name,".RData"))
  sim_results[[name]] <- list("reject" = reject, "timing" = timing, "tuning" = tuning)
  rm(reject)
  rm(timing)
  rm(tuning)
}


plot_DGPs<-vector("list", length(sim_file_names))
names(plot_DGPs)<-sim_file_names
for(name in sim_file_names){
  timing <- sim_results[[name]][["timing"]]
  pars <- dimnames(timing)$pars
  boot <- dimnames(timing)$boot
  sim <- dimnames(timing)$sim
  p <- vector("list", length(pars))
  names(p) <- pars
  for(setup in pars){
    box_df<-data.frame(matrix(0, nrow = length(boot)*length(sim), ncol = 2))
    names(box_df) = c("time", "method")
    for(i in 1:length(boot)){
      box_df$time[(1+(i-1)*length(sim)):(i*length(sim))] <- -timing[setup, i,]
      box_df$method[(1+(i-1)*length(sim)):(i*length(sim))] <- rep(gsub("unpen-own-", "", boot[i]), length(sim))
    }
    p[[setup]] <- ggplot(box_df, aes(method, time)) +
      geom_boxplot() +
      stat_summary(fun=mean, geom="point", shape=20, size=3, color="red", fill="red") +
      ggtitle(setup) +
      theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust = 1, size = 15),
            title = element_text(size = 15))
  }
  plot_DGPs[[name]] <- annotate_figure(ggarrange(plotlist = p, nrow = 4, ncol = 4), top = text_grob(name, face = "bold", size = 30))
}
for(name in sim_file_names){
  ggsave(paste0("timing_",name,".pdf"), plot = plot_DGPs[[name]], width = 30, height = 40)
}
