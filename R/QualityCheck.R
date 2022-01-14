#' Quality check
#'
#' This function is for automatically generating quality check plots, including intra-group CV plot and protein intensity rank plot
#'
#' @param quan Protein quantification result
#' @param group a matrix contains file names and their group information
#'
#' @export


# Quality Illstration
Quality_check <- function(quan, group){
  quan <- as.data.frame(quan)
  quan$PepNum <- NULL
  output <- list()
  quan <- 2^quan
  RSD <- numeric()
  group_name <- unique(group$GroupID)
  for (i in 1:length(group_name)) {
    quan_val <- quan[,group$GroupID == group_name[i]]
    SD <- apply(quan_val,1,function(x) sd(x, na.rm = TRUE))
    AVE <- apply(quan_val,1,function(x) mean(x, na.rm = TRUE))
    RSD_ind <- SD/AVE
    RSD <- cbind(RSD,RSD_ind)
  }
  RSD <- as.data.frame(RSD)
  names(RSD) <- group_name
  RSD_plot <- suppressMessages(melt(RSD))
  names(RSD_plot) <- c('Group','CV')
  RSD_plot$CV <- 100*RSD_plot$CV

  IGCV <- ggplot(RSD_plot,aes(x=Group,y=CV,fill=Group)) + geom_boxplot() +
    scale_fill_brewer(palette="Blues") + scale_x_discrete(limits=unique(RSD_plot$Group)) +
    ylab('CV %') +
    scale_y_continuous(limits = c(0, 150), breaks = seq(0, 150, by = 20)) +
    geom_hline(yintercept=15, linetype="dashed", color = "red", size=1) +
    theme(legend.position="none",panel.grid.minor=element_blank(),
          panel.background=element_rect(fill=NA, colour=NA),
          panel.border=element_rect(size=1,fill=NA,colour="black"))+
    theme(axis.text=element_text(size=14, family="Helvetica", colour="black"),
          axis.title=element_text(size=14, family="Helvetica",colour="black"),
          legend.text=element_text(size=14, family="Helvetica",colour="black"),
          legend.title=element_text(size=14, family="Helvetica", colour="black"))
  #Plotting the intensity curve of mean protein intensities
  MeanInt <- apply(as.matrix(quan),1,function(x) log2(mean(x, na.rm = TRUE)))
  MeanInt <- sort(MeanInt,decreasing=TRUE)
  ProNum <- 1:length(MeanInt)
  Int_plot <- as.data.frame(matrix(c(MeanInt,ProNum),ncol=2,byrow=FALSE))
  colnames(Int_plot) <- c("Log2Intensities","ProteinNumber")

  PRank <- ggplot(Int_plot,aes(x=ProteinNumber,y=Log2Intensities)) +
    geom_point(alpha=0.3,shape=21) +
    theme(legend.position="none",panel.grid.minor=element_blank(),
          panel.background=element_rect(fill=NA, colour=NA),
          panel.border=element_rect(size=1,fill=NA,colour="black")) +
    theme(axis.text=element_text(size=14, family="Helvetica", colour="black"),
          axis.title=element_text(size=14, family="Helvetica",colour="black"),
          legend.text=element_text(size=14, family="Helvetica",colour="black"),
          legend.title=element_text(size=14, family="Helvetica", colour="black"))
  output$CVplot <- IGCV
  output$PRplot <- PRank
  return(output)
}
