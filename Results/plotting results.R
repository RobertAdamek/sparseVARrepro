setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(ggplot2)
library(ggpubr)
library(data.table)
library(abind)
library(viridis)
 
clean_name<-function(meth){
  meth<-gsub(pattern="-L1-unpen-own", replacement="", meth)
  meth<-gsub(pattern="-11", replacement="", meth)
  meth<-gsub(pattern="04", replacement="0.4", meth)
  meth<-gsub(pattern="08", replacement="0.8", meth)
  return(meth)
}

Ns<-c(20, 40, 100, 200)
Ts<-c(50, 100, 200, 500)
DGPs<-c(0:10)
means<-c("0013","00175","0035","01","0175","025","05","075","1","10"); means_val<-c(0.013,0.0175,0.035,0.1,0.175,0.25,0.5,0.75,1,10)
proportions<-c("01","05","09"); proportions_val<-c(0.1,0.5,0.9)

sim_result_files<-list.files(); sim_result_files<-sim_result_files[which(grepl(".RData", sim_result_files, fixed=TRUE))]

# merge size results for the same DGP
for(dgp in DGPs){
  if(!any(grepl(paste0("dgp",dgp), sim_result_files, fixed=TRUE))){
    next
  }
  plot_ind<-which(as.logical(grepl(paste0("dgp",dgp,"_size"), sim_result_files, fixed=TRUE) * !grepl(paste0("merged"), sim_result_files, fixed=TRUE)))
  for(i in 1:length(plot_ind)){
    load(sim_result_files[plot_ind[i]])
    size<-reject
    rm(reject)
    if(i==1){
      merged<-size
    }else{
      merged<-abind(merged,size,along=2)
    }
  }
  merged<-merged[,order(dimnames(merged)[[2]]),]
  save(merged, file=paste0("dgp",dgp,"_size_merged.RData"))
  rm(merged)
}
sim_result_files<-list.files(); sim_result_files<-sim_result_files[which(grepl(".RData", sim_result_files, fixed=TRUE))]
plot_methods<-c("VAR-BIC", "VAR-TF", "VAR-oracle", "BWB", "MBB")

plot_scale<-70
# plot all ----------------------------------------------------------------

for(dgp in DGPs){
  if(!any(grepl(paste0("dgp",dgp), sim_result_files, fixed=TRUE))){
    next
  }
  if(dgp==0){
    tit<-"Factor model with sparse idiosyncratic components"
  }else if(dgp==1){
    tit<-"Block-diagonal VAR(1)"
  }else if(dgp==2){
    tit<-"Diagonal VAR(1)"
  }else if(dgp==5){
    tit<-"Weakly sparse VAR(1), strong dependence"
  }else if(dgp==9){
    tit<-"Weakly sparse VAR(1)"
  }else if(dgp==10){
    tit<-"Example 1 DGP"
  }
  # size plot
  plot_ind<-which(grepl(paste0("dgp",dgp,"_size_merged"), sim_result_files, fixed=TRUE))
  for(i in 1:length(plot_ind)){
    load(sim_result_files[plot_ind[i]])
    Methods<-clean_name(dimnames(merged)[[2]])
    size<-merged
    rm(merged)
    p<-list()
    count<-0
    # plots for the different Ns
    for(ns in 1:length(Ns)){
      count<-count+1
      D_mat<-matrix(0, nrow=length(Ts), ncol=1+length(Methods))
      colnames(D_mat)<-c("T", Methods)
      D_mat[,1]<-Ts
      D<-data.table(D_mat) 
      for(ms in (1:length(Methods))){
        D[,1+ms]<- size[(ns-1)*length(Ts)+1:length(Ts),ms,2]
      }
      # remove methods that I don't want to plot
      cols_select<-c("T",plot_methods); if(dgp==0){cols_select<-cols_select[-which(cols_select=="VAR-oracle")]}
      D<-D[,..cols_select]
      D_melt<-melt.data.table(data=D, id.vars="T", value.name="rejection", variable.name = "Method")
      D_melt$T<-as.factor(D_melt$T)
      p[[count]]<-ggplot(data=D_melt, mapping=aes(x=T, y=rejection,  fill=Method))+#color=Method,
        coord_cartesian(ylim = c(0, 0.5))+
        #geom_line()+
        #geom_segment()+
        #geom_point()+
        geom_col(position="dodge")+
        #scale_x_continuous(breaks=Ts)+
        geom_hline(yintercept = 0.05, color="black")+
        #scale_color_manual(values=c("VAR-BIC"="blue", "VAR-TF"="violet", "VAR-TF-0.4"="red","VAR-TF-0.8"="purple","MBB"="green", "BWB"="darkgreen", "VAR-oracle"="darkorange"))+
        scale_fill_manual(values=c("VAR-BIC"=rgb(red=68, green=1, blue=83, alpha=255, maxColorValue = 255), 
                                   "VAR-TF"=rgb(red=58, green=82, blue=138, alpha=255, maxColorValue = 255), 
                                   "VAR-oracle"=rgb(red=32, green=144, blue=139, alpha=255, maxColorValue = 255), 
                                   "VAR-TF-0.4"="red",
                                   "VAR-TF-0.8"="pink",
                                   "BWB"=rgb(red=93, green=199, blue=98, alpha=255, maxColorValue = 255), 
                                   "MBB"=rgb(red=252, green=231, blue=36, alpha=255, maxColorValue = 255)))+
        ggtitle(paste0("N=",Ns[ns]))
    }
    # combine horizontally
    p_combined<-ggarrange(plotlist=p, nrow=1, ncol=length(Ns), common.legend = TRUE, legend="none")
      p_combined_legend<-ggarrange(plotlist=p, nrow=1, ncol=length(Ns), common.legend = TRUE, legend="bottom")
    p_combined_title<-annotate_figure(p_combined, top=text_grob(paste0("Size")))
      p_combined_title_legend<-annotate_figure(p_combined_legend, top=text_grob(paste0("Size")))
    p_combined_title_DGP<-annotate_figure(p_combined_title, top=text_grob(paste0(tit), face="bold"))
      p_combined_title_DGP_legend<-annotate_figure(p_combined_title_legend, top=text_grob(paste0(tit), face="bold"))
  }
  size_for_later<-ggarrange(p_combined_title, nrow=1, ncol=1, common.legend = TRUE, legend="none")
  #power plot
  if(!any(grepl(paste0("dgp",dgp,"_power"), sim_result_files, fixed=TRUE))){
    next
  }
  plot_ind<-which(grepl(paste0("dgp",dgp,"_power"), sim_result_files, fixed=TRUE))
  q<-list()
  for(i in 1:length(plot_ind)){
    load(sim_result_files[plot_ind[i]])
    Methods<-clean_name(dimnames(reject)[[2]])
    power<-reject
    rm(reject)
    p<-list()
    for(j in 1:length(means)){
      m<-means[j]
      if(grepl(paste0("_mu",m,"_"), sim_result_files[plot_ind[i]], fixed=TRUE)){
        mu<-means_val[j]
      }
    }
    for(k in 1:length(proportions)){
      pr<-proportions[k]
      if(grepl(paste0("_prop",pr,"."), sim_result_files[plot_ind[i]], fixed=TRUE)){
        prop<-proportions_val[k]
      }
    }
    count<-0
    # plots for the different Ns
    for(ns in 1:length(Ns)){
      count<-count+1
      D_mat<-matrix(0, nrow=length(Ts), ncol=1+length(Methods))
      colnames(D_mat)<-c("T", Methods)
      D_mat[,1]<-Ts
      D<-data.table(D_mat)
      for(ms in (1:length(Methods))){
        D[,1+ms]<- power[(ns-1)*length(Ts)+1:length(Ts),ms,2]
      }
      # remove methods that I don't want to plot
      cols_select<-c("T",plot_methods); if(dgp==0){cols_select<-cols_select[-which(cols_select=="VAR-oracle")]}
      D<-D[,..cols_select]
      D_melt<-melt.data.table(data=D, id.vars="T", value.name="rejection", variable.name = "Method")
      D_melt$T<-as.factor(D_melt$T)
      
    p[[count]]<-ggplot(data=D_melt, mapping=aes(x=T, y=rejection,  fill=Method))+
      coord_cartesian(ylim = c(0, 1))+
      geom_col(position="dodge")+
      scale_fill_manual(values=c("VAR-BIC"=rgb(red=68, green=1, blue=83, alpha=255, maxColorValue = 255), 
                                 "VAR-TF"=rgb(red=58, green=82, blue=138, alpha=255, maxColorValue = 255), 
                                 "VAR-oracle"=rgb(red=32, green=144, blue=139, alpha=255, maxColorValue = 255), 
                                 "VAR-TF-0.4"="red",
                                 "VAR-TF-0.8"="pink",
                                 "BWB"=rgb(red=93, green=199, blue=98, alpha=255, maxColorValue = 255), 
                                 "MBB"=rgb(red=252, green=231, blue=36, alpha=255, maxColorValue = 255)))+
      ggtitle(paste0("N=",Ns[ns]))+
      guides(color = guide_legend(override.aes = list(shape = NA ) ) )+ 
      theme(legend.position="bottom")
    }
    # combine horizontally

    leg<-(i==length(plot_ind))
    p_combined<-ggarrange(plotlist=p, nrow=1, ncol=length(Ns), common.legend = TRUE, legend="none")
    p_combined_title<-annotate_figure(p_combined, top=text_grob(paste0("Power - means: ", mu, ", proportion: ", prop)))
    q[[i]]<-p_combined_title
  }
  powers_for_later<-q
  legend_for_later<-get_legend(p[[1]])
  # combine vertically
  all_rows<-ggarrange(plotlist=q, nrow=length(plot_ind), ncol=1, common.legend = TRUE, legend="bottom", legend.grob=get_legend(p[[1]]))#
  all_rows_title<-annotate_figure(all_rows, top=text_grob(paste0(tit), face="bold"))

  
  # power+size
  size_power_list<-c(list(size_for_later),powers_for_later)
  size_power<-ggarrange(plotlist=size_power_list, nrow=length(size_power_list), ncol=1, common.legend = TRUE, legend="bottom", legend.grob=legend_for_later)
  ggsave(paste0("dgp",dgp,"_size_power.pdf"), plot=size_power, width=length(Ns)*plot_scale, height=length(size_power_list)*plot_scale, units="mm")
}

