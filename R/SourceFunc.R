#' UHR-IonStar Source Functions
#'
#' These functions are the application functions in the R shiny web app 'UHR-IonStar'.
#'
#'
#' @export

### Data Processing, Ratio_p
Ratio_p <- function(quan, group, decoy, control_group=NA, test_type="t", equalVar = TRUE){
  withProgress(message = 'Protein Quantification', value = 0, {
    incProgress(0.5, detail = 'statistical testing...')
  outputs <- list()
  if ('PepNum' %in% names(quan)){
    quan$PepNum <- NULL
  }

  ##Data reordering and cleanup (keep the order in quan same with group)
  group <- group[order(group$GroupID, group$Rawfiles),]
  group_order <- match(as.character(group$Rawfiles), names(quan))
  quan <- quan[,group_order]
  #Remove decoy protein entries
  quan <- quan[!grepl(decoy, rownames(quan)),]
  quan_results <- quan
  quan <- 2^as.matrix(quan)
  ave_df <- as.data.frame(avearrays(quan, ID=group$GroupID, weights = NULL)) # calculate mean for each group

  #Assign the control group, substitute with input from Shiny
  #Calculate ratios between case and control conditions
  con <- match(control_group,unique(group$GroupID))

  if(is.na(con)){
    con <- 1
    control_group <- unique(as.character(group$GroupID))[1]
  }
  control_group <<- control_group # turn to global variable
  ave_ratio <- ave_df/ave_df[,con]
  group_char <- as.character(unique(group$GroupID))
  group_num <- length(group_char)
  r_chara <- rep(NA,group_num)
  for (i in 1:group_num) {
    r_chara[i] <- paste("Ratio_",group_char[i],".",control_group,sep="")
  }

  #Paired t.test between case and control conditions
  quan_test <- quan
  colnames(quan_test) <- group$GroupID
  if(test_type=="Original t-test"){
    p_chara <- rep(NA,group_num)
    for (i in 1:group_num){
      pvals <- apply(quan_test, 1, function(x)
        t.test(x[which(colnames(quan_test)== control_group)],
               x[which(colnames(quan_test)== group_char[i])],var.equal = as.logical(equalVar))$p.value)
      p_chara[i] <- paste("pVal_",group_char[i],".",control_group,sep="")
      ave_ratio <- cbind(ave_ratio,pvals)
    }
  }
  if(test_type=="Paired t-test"){
    for (i in 1:group_num) {
      p_chara <- rep(NA,group_num)
      pvals <- apply(quan_test, 1, function(x)
        t.test(x[which(colnames(quan_test)== control_group)],
               x[which(colnames(quan_test)== group_char[i])],var.equal = as.logical(equalVar), paired = TRUE)$p.value)
      p_chara[i] <- paste("pVal_",group_char[i],".",control_group,sep="")
      ave_ratio <- cbind(ave_ratio,pvals)
    }
  }
  if(test_type=="Wilcoxon signed-rank test"){
    p_chara <- rep(NA,group_num)
    for (i in 1:group_num) {
      suppressWarnings(pvals <- apply(quan_test, 1, function(x)
        wilcox.test(x[which(colnames(quan_test)== control_group)],
                    x[which(colnames(quan_test)== group_char[i])])$p.value))
      p_chara[i] <- paste("pVal_",group_char[i],".",control_group,sep="")
      ave_ratio <- cbind(ave_ratio,pvals)
    }
  }
  if(test_type=="Kruskal–Wallis ANOVA"){
    pvals <- apply(quan_test, 1, function(x)
      kruskal.test(x, colnames(quan_test))$p.value)
    p_chara <- "pVal_KruskalWallis"
    ave_ratio <- cbind(ave_ratio,pvals)
  }
  colnames(ave_ratio) <- c(r_chara,p_chara)
  tmp_col <- gsub('\\..*','',colnames(ave_ratio))
  tmp_col <- gsub('.*_','',tmp_col)
  ave_ratio <- ave_ratio[,!(tmp_col %in% control_group)]

  outputs$quan_results <- quan_results
  outputs$group_ordered <- group
  outputs$ave_ratio <- ave_ratio
  outputs$all <- cbind(quan_results,ave_ratio)
  outputs$control_group <- control_group
  outputs$group_char <- group_char
  outputs$notice <- "Processing completed."
  outputs$test_type <- test_type

  return(outputs)
})
}

### Data Visualization, plot_for_all
plot_for_all <- function(quan, group, plot_selection, Group1, Group2){
  quan <- 2^quan
  ave_df <- as.data.frame(avearrays(quan, ID=group$GroupID, weights = NULL))
  if(plot_selection == "Intra-group CV boxplot"){
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

    selected_plot <- ggplot(RSD_plot,aes(x=Group,y=CV,fill=Group)) + geom_boxplot() +
      scale_fill_brewer(palette="Blues") + scale_x_discrete(limits=unique(RSD_plot$Group)) +
      ylab('CV %') +
      scale_y_continuous(limits = c(0, 150), breaks = seq(0, 150, by = 20)) +
      theme(legend.position="none",panel.grid.minor=element_blank(),
            panel.background=element_rect(fill=NA, colour=NA),
            panel.border=element_rect(size=1,fill=NA,colour="black"))+
      theme(axis.text=element_text(size=14, family="Helvetica", colour="black"),
            axis.title=element_text(size=14, family="Helvetica",colour="black"),
            legend.text=element_text(size=14, family="Helvetica",colour="black"),
            legend.title=element_text(size=14, family="Helvetica", colour="black"))
  }
  if(plot_selection == "Protein intensity-rank curve"){
    #Plotting the intensity curve of mean protein intensities
    MeanInt <- apply(as.matrix(quan),1,function(x) log2(mean(x, na.rm = TRUE)))
    MeanInt <- sort(MeanInt,decreasing=TRUE)
    ProNum <- 1:length(MeanInt)
    Int_plot <- as.data.frame(matrix(c(MeanInt,ProNum),ncol=2,byrow=FALSE))
    colnames(Int_plot) <- c("Log2Intensities","ProteinNumber")

    selected_plot <- ggplot(Int_plot,aes(x=ProteinNumber,y=Log2Intensities)) +
      geom_point(alpha=0.3,shape=21) +
      theme(legend.position="none",panel.grid.minor=element_blank(),
            panel.background=element_rect(fill=NA, colour=NA),
            panel.border=element_rect(size=1,fill=NA,colour="black")) +
      theme(axis.text=element_text(size=14, family="Helvetica", colour="black"),
            axis.title=element_text(size=14, family="Helvetica",colour="black"),
            legend.text=element_text(size=14, family="Helvetica",colour="black"),
            legend.title=element_text(size=14, family="Helvetica", colour="black"))

  }
  if(plot_selection == "Inter-group correlation plot"){
    con2 <- match(c(Group1,Group2),group$GroupID)
    if(NA %in% con2){
      Group1 <- colnames(ave_df)[1]
      Group2 <- colnames(ave_df)[2]
    }

    Corr_plot <- ave_df[,match(c(Group1,Group2),colnames(ave_df))]
    Corr_plot <- as.data.frame(matrix(apply(Corr_plot,c(1,2),log2),ncol=2,byrow=FALSE))
    #Data processing and plotting
    colnames(Corr_plot) <- c("G1","G2")
    corr <- round(cor(Corr_plot$G1,Corr_plot$G2),digits=3)
    corr_label <- paste("R Square=",corr,sep="")
    selected_plot <- ggplot(Corr_plot,aes(x=G1,y=G2)) + geom_point(alpha=0.3,shape=21) +
      theme(legend.position="none",panel.grid.minor=element_blank(),
            panel.background=element_rect(fill=NA, colour=NA),
            panel.border=element_rect(size=1,fill=NA,colour="black")) +
      theme(axis.text=element_text(size=14, family="Helvetica", colour="black"),
            axis.title=element_text(size=14, family="Helvetica",colour="black"),
            legend.text=element_text(size=14, family="Helvetica",colour="black"),
            legend.title=element_text(size=14, family="Helvetica", colour="black")) +
      labs(x=Group1,y=Group2) + annotate("text", x = 35, y = 22, label = corr_label)
  }
  if(plot_selection == "Pearson correlation matrix plot"){
    #Pearson correlation matrix
    Correlation <- cor(quan)
    cor_min <- min(Correlation)
    selected_plot <- ggcorrplot(Correlation,colors = c("royalblue3", "gray90", "brown3"),
                                type = "lower", method = "circle") +
      theme(axis.text=element_text(size=14, family="Helvetica", colour="black"),
            axis.title=element_text(size=14, family="Helvetica",colour="black"),
            legend.text=element_text(size=14, family="Helvetica",colour="black"),
            legend.title=element_text(size=14, family="Helvetica", colour="black")) +
      scale_fill_gradientn(name="Correlation", colors=c('#f6ff50','#5eff50','#50fff6','#5055ff'), limits=c(cor_min,1),
                           na.value = 'gray90')
  }
  if(plot_selection == "Principal Component Analysis (PCA)"){
    quan <- log2(quan)
    PCA_data <- as.data.frame(t(quan))
    PCA_results <- prcomp(PCA_data,center=TRUE,scale.=TRUE)

    PCA_level <- group$GroupID
    selected_plot <- ggbiplot(PCA_results, obs.scale=1, var.scale=1, groups=PCA_level,
                              circle=FALSE, var.axes=FALSE) +
      geom_point(aes(colour=PCA_level),size=5) + scale_colour_brewer(name="Groups",palette="Set2") +
      #geom_hline(yintercept=0) + geom_vline(xintercept=0) +
      geom_label_repel(aes(label=group$Rawfiles),color="black",size=3,fill=NA, label.size=NA) +
      theme(panel.grid.minor=element_blank(),
            #panel.grid.major=element_line(colour="gray80",linetype="dashed"),
            panel.background=element_rect(fill=NA, colour=NA),
            panel.border=element_rect(size=1,fill=NA,colour="black")) +
      theme(axis.text=element_text(size=14, family="Helvetica", colour="black"),
            axis.title=element_text(size=14, family="Helvetica",colour="black"),
            legend.text=element_text(size=14, family="Helvetica",colour="black"),
            legend.title=element_text(size=14, family="Helvetica", colour="black")) +
      coord_fixed(ratio=1.05)

  }
  return(selected_plot)

}

### Data Visualization, plot_ratio_distribution
plot_ratio_distribution <- function(ave_ratio, group, selected, dist_corr){
  outputs <- list()
  group_char <- unique(group$GroupID)
  ratio_matrix <- as.matrix(ave_ratio[,1:(ncol(ave_ratio)/2)])
  tmp <- as.vector(match(selected, group_char[group_char != control_group]))
  ratio_matrix <- as.data.frame(ratio_matrix[,tmp])
  if(length(tmp)==1){
    colnames(ratio_matrix) <- paste("Ratio_", selected,"/",control_group, sep="")
  }

  #ratio_matrix <- log2(ratio_matrix)
  if(!is.null(dist_corr)){
    if(dist_corr == "Correction"){
      for(i in 1:ncol(ratio_matrix)){
        dens <- density(ratio_matrix[,i])
        correc <- dens$x[which(dens$y == max(dens$y))]
        ratio_matrix[,i] <- ratio_matrix[,i] - correc
      }
    }
  }
  outputs$correction_ratio <- ratio_matrix

  ratio_plot <- melt(ratio_matrix)
  axis_limit <- ceiling(max(abs(ratio_plot$value)))

  Ratio_denplot <- ggplot(ratio_plot,aes(x=value,colour=variable)) +
    geom_density(size=0.6) + scale_x_continuous(limits=c(-axis_limit,axis_limit)) +
    scale_colour_brewer(name="Groups",palette="Set2") +
    theme(panel.grid.minor=element_blank(), panel.background=element_rect(fill=NA, colour=NA),
          panel.border=element_rect(size=1,fill=NA,colour="black")) +
    theme(axis.text=element_text(size=12, family="Helvetica", colour="black"),
          axis.title=element_text(size=12, family="Helvetica",colour="black"),
          legend.text=element_text(size=12, family="Helvetica",colour="black"),
          legend.title=element_text(size=12, family="Helvetica", colour="black")) +
    labs(x="Log2Ratio",y="Density") + geom_vline(xintercept=0,linetype="dashed",colour="gray80")

  outputs$plotting <- Ratio_denplot
  return(outputs)

}

### Discovery of Changes, sign_changes
sign_changes <- function(quan_all, test_type, chosen_group, R_cutoff=1.4, p_cutoff=0.05){
  outputs <- list()
  if(test_type %in% c("Original t-test","Paired t-test","Wilcoxon signed-rank test")){
    R_col <- paste("Ratio_",as.character(chosen_group),".",control_group,sep="")
    p_col <- paste("pVal_",as.character(chosen_group),".",control_group,sep="")
    quan_select <- quan_all[,c(match(c(R_col,p_col),colnames(quan_all)))]
  }else{
    quan_select <- quan_all[,c((ncol(quan_all)-1), ncol(quan_all))]
  }
  names(quan_select) <- c("Ratio","pVal")
  R_upper <- log2(as.numeric(R_cutoff))
  R_lower <- log2(1/as.numeric(R_cutoff))
  quan_up <- subset(quan_select, log2(Ratio) >= R_upper & pVal<=as.numeric(p_cutoff))
  quan_down <- subset(quan_select, log2(Ratio) <= R_lower & pVal<=as.numeric(p_cutoff))
  quan_changed <- rbind(quan_up,quan_down)
  quan_select$match <- rownames(quan_select) %in% rownames(quan_changed)
  quan_unchanged <- quan_select[quan_select$match == FALSE,]
  quan_unchanged$match <- NULL

  outputs$quan_up <- quan_up
  outputs$quan_down <- quan_down
  outputs$quan_changed <- quan_changed
  outputs$quan_unchanged <- quan_unchanged
  outputs$R_cutoff <- R_cutoff
  outputs$p_cutoff <- p_cutoff
  return(outputs)

}

### Discovery of Changes, plot_for_all2
plot_for_all2 <- function(sig_list, quan, plot_selection,TOPx = 50){
  quan_up <- sig_list$quan_up
  quan_down <- sig_list$quan_down
  quan_unchanged <- sig_list$quan_unchanged
  quan_up$Ratio <- apply(as.matrix(quan_up$Ratio),1, log2)
  quan_up$pVal <- apply(as.matrix(quan_up$pVal),1,function(x) -log10(x))
  quan_down$Ratio <- apply(as.matrix(quan_down$Ratio),1, log2)
  quan_down$pVal <- apply(as.matrix(quan_down$pVal),1,function(x) -log10(x))
  quan_unchanged$Ratio <- apply(as.matrix(quan_unchanged$Ratio),1, log2)
  quan_unchanged$pVal <- apply(as.matrix(quan_unchanged$pVal),1,function(x) -log10(x))
  quan_up$Group <- "Up-changed"
  quan_down$Group <- "Down-changed"
  quan_unchanged$Group <- "Unchanged"
  quan_vol <- rbind(quan_up,quan_down,quan_unchanged)
  quan_vol$Group <- factor(quan_vol$Group, levels = c("Up-changed","Down-changed","Unchanged"))
  if(plot_selection == "Valcano Plot"){
    # Volcano plot of data distribution
    cols <- c("red","green","black")
    x_limit <- round(max(abs(quan_vol$Ratio)))
    y_limit <- round(max(abs(quan_vol$pVal)))
    plotting <- ggplot(quan_vol,aes(x=Ratio,y=pVal,fill=Group)) + geom_point(alpha=0.4,shape=21) +
      scale_fill_manual(values=cols) + scale_x_continuous(limits=c(-x_limit,x_limit)) + scale_y_continuous(limits=c(0,y_limit)) +
      geom_hline(yintercept=-log10(as.numeric(sig_list$p_cutoff)), linetype="dashed",
                 color = "red", size=0.5)+
      geom_vline(xintercept=-log2(as.numeric(sig_list$R_cutoff)), linetype="dashed",
                 color = "red", size=0.5)+
      geom_vline(xintercept=log2(as.numeric(sig_list$R_cutoff)), linetype="dashed",
                 color = "red", size=0.5)+
      theme(panel.grid.minor=element_blank(),
            panel.grid.major=element_line(colour="gray80",linetype="dashed"),
            panel.background=element_rect(fill=NA, colour=NA),
            panel.border=element_rect(size=1,fill=NA,colour="black")) +
      theme(axis.text=element_text(size=12, family="Helvetica", colour="black"),
            axis.title=element_text(size=12, family="Helvetica",colour="black"),
            legend.text=element_text(size=12, family="Helvetica",colour="black"),
            legend.title=element_text(size=12, family="Helvetica", colour="black"))+
      labs(x="Log2(Ratio)",y="-Log10(p-value)")
    # Vol_plot <- ggplot(quan_vol,aes(x=Ratio,y=pVal,fill=Group)) + geom_point(alpha=0.3,shape=21) + scale_fill_manual(values=cols) + scale_x_continuous(limits=c(-1.2*x_limit,1.2*x_limit))
  }
  if(plot_selection == "Intensity curve with changes highlighted"){
    # Plotting the intensity curve highlighting significantly changed proteins
    # NOTE: the highlighted changed proteins are based on the results above.
    # When changing the comparison groups the output will be different as well
    quan_new <- quan
    quan_new$Groups <- rownames(quan) %in% rownames(quan_vol[quan_vol$Group != "Unchanged",])
    quan_new$Log2Intensities <-apply(2^quan,1,function(x) log2(mean(x, na.rm = TRUE)))
    quan_new <- quan_new[order(quan_new$Log2Intensities, decreasing = TRUE),]
    quan_new$ProteinNumber <- 1:nrow(quan_new)
    quan_new <- quan_new[,c(ncol(quan_new)-2,ncol(quan_new)-1, ncol(quan_new))]
    quan_new$Groups[quan_new$Groups == FALSE] <- "Unchanged"
    quan_new$Groups[quan_new$Groups == TRUE] <- "Changed"
    plotting <- ggplot(quan_new,aes(x=ProteinNumber,y=Log2Intensities)) +
      geom_point(alpha=1,aes(colour=Groups, shape=Groups, size=Groups)) +
      scale_colour_manual(values=c("red","black")) +
      scale_shape_manual(values=c(17,1))+
      scale_size_manual(values=c(2,2))+
      theme(panel.grid.minor=element_blank(), panel.background=element_rect(fill=NA, colour=NA),
            panel.border=element_rect(size=1,fill=NA,colour="black")) +
      theme(axis.text=element_text(size=12, family="Helvetica", colour="black"),
            axis.title=element_text(size=12, family="Helvetica",colour="black"),
            legend.text=element_text(size=12, family="Helvetica",colour="black"),
            legend.title=element_text(size=12, family="Helvetica", colour="black"))
  }
  if(plot_selection == "Plot for Up-regulated"){
    #Plotting the TOPx Up/Down-regulated proteins
    # TOPx <- "20" #assign how many proteins to show
    #quan_up$Fold <- apply(as.matrix(quan_up$Ratio),1,function(x) 2^x)
    #quan_down$Fold <- apply(as.matrix(quan_down$Ratio),1,function(x) 2^(-x))
    quan_up$Fold <- as.matrix(quan_up$Ratio)
    quan_up$acName <- substr(rownames(quan_up), start = 1, stop = 20)
    quan_up <- quan_up[order(quan_up$Fold,decreasing=TRUE),]
    Up_plot <- quan_up[1:min(nrow(quan_up),as.numeric(TOPx)),]
    plotting <- ggplot(Up_plot,aes(x=reorder(acName,-Fold),y=Fold)) +
      geom_bar(stat="identity",colour="black",fill="brown2",alpha=0.7) +
      theme(legend.position="none",panel.grid.minor=element_blank(),
            panel.background=element_rect(fill=NA, colour=NA),
            panel.border=element_rect(size=1,fill=NA,colour="black")) +
      theme(axis.text=element_text(size=10, family="Helvetica", colour="black"),
            axis.text.x=element_text(angle=90,hjust=1),
            axis.title=element_text(size=12, family="Helvetica",colour="black"),
            legend.text=element_text(size=12, family="Helvetica",colour="black"),
            legend.title=element_text(size=12, family="Helvetica", colour="black")) +
      labs(x="",y="Log2 Fold change")
  }
  if(plot_selection == "Plot for Down-regulated"){
    quan_down$Fold <-as.matrix(quan_down$Ratio)
    quan_down$acName <- substr(rownames(quan_down), start = 1, stop = 15)
    quan_down <- quan_down[order(quan_down$Fold,decreasing=FALSE),]
    Down_plot <- quan_down[1:min(nrow(quan_down),as.numeric(TOPx)),]
    plotting <- ggplot(Down_plot,aes(x=reorder(acName,Fold),y=Fold)) +
      geom_bar(stat="identity",colour="black",fill="springgreen2",alpha=0.7) +
      theme(legend.position="none",panel.grid.minor=element_blank(),
            panel.background=element_rect(fill=NA, colour=NA),
            panel.border=element_rect(size=1,fill=NA,colour="black")) +
      theme(axis.text=element_text(size=10, family="Helvetica", colour="black"),
            axis.text.x=element_text(angle=90,hjust=1),
            axis.title=element_text(size=12, family="Helvetica",colour="black"),
            legend.text=element_text(size=12, family="Helvetica",colour="black"),
            legend.title=element_text(size=12, family="Helvetica", colour="black")) +
      labs(x="",y="Log2 Fold change")
  }
  if(plot_selection == "Gene Ontology (Accession)"){
    annotation="GOTERM_BP_DIRECT"
    pvalueCutoff=1
    qvalueCutoff=1
    minGSSize=5

    quan_sum <- rbind(quan_up,quan_down)
    rownames(quan_sum) <- gsub("\\:.*","",rownames(quan_sum))

    sum_list <- rownames(quan_sum)
    sum_David <- enrichDAVID(gene=sum_list, idType="UNIPROT_ACCESSION",annotation=annotation,david.user="a",
                             pvalueCutoff=pvalueCutoff,qvalueCutoff=qvalueCutoff,minGSSize=minGSSize)
    sum_David_results <- as.data.frame(summary(sum_David))
    ##Just for plotting, should be customized.
    David_plot <- subset(sum_David_results,pvalue<0.05 & qvalue<0.2 & Count>=10)
    Drow <- 20
    David_plot2 <- David_plot[1:Drow,]
    colnames(David_plot2)[2] <- "GOTerm"
    plotting <- ggplot(David_plot2,aes(x=reorder(GOTerm,Count),y=Count,fill=qvalue)) +
      geom_bar(stat="identity") + coord_flip() +
      scale_fill_distiller(name="p-value",palette="Spectral") +
      theme(panel.grid.minor=element_blank(),
            panel.background=element_rect(fill=NA, colour=NA),
            panel.border=element_rect(size=1,fill=NA,colour="black")) +
      theme(axis.text=element_text(size=10, family="Helvetica", colour="black"),
            axis.title=element_text(size=12, family="Helvetica",colour="black"),
            legend.text=element_text(size=12, family="Helvetica",colour="black"),
            legend.title=element_text(size=12, family="Helvetica", colour="black")) +
      labs(title="GO - Biological Processes",x="",y="Protein number")
  }

  if(plot_selection == "Gene Ontology (ProteinID)"){
    annotation="GOTERM_BP_DIRECT"
    pvalueCutoff=1
    qvalueCutoff=1
    minGSSize=5

    quan_sum <- rbind(quan_up,quan_down)
    rownames(quan_sum) <- gsub("\\:.*","",rownames(quan_sum))

    sum_list <- as.character(rownames(quan_sum))
    sum_David <- enrichDAVID(gene=sum_list, idType="UNIPROT_ID",annotation=annotation,david.user="a",
                             pvalueCutoff=pvalueCutoff,qvalueCutoff=qvalueCutoff,minGSSize=minGSSize)
    sum_David_results <- as.data.frame(summary(sum_David))
    ##Just for plotting, should be customized.
    David_plot <- subset(sum_David_results,pvalue<0.05 & qvalue<0.2 & Count>=10)
    Drow <- 20
    David_plot2 <- David_plot[1:Drow,]
    colnames(David_plot2)[2] <- "GOTerm"
    plotting <- ggplot(David_plot2,aes(x=reorder(GOTerm,Count),y=Count,fill=qvalue)) +
      geom_bar(stat="identity") + coord_flip() +
      scale_fill_distiller(name="p-value",palette="Spectral") +
      theme(panel.grid.minor=element_blank(),
            panel.background=element_rect(fill=NA, colour=NA),
            panel.border=element_rect(size=1,fill=NA,colour="black")) +
      theme(axis.text=element_text(size=10, family="Helvetica", colour="black"),
            axis.title=element_text(size=12, family="Helvetica",colour="black"),
            legend.text=element_text(size=12, family="Helvetica",colour="black"),
            legend.title=element_text(size=12, family="Helvetica", colour="black")) +
      labs(title="GO - Biological Processes",x="",y="Protein number")
  }
  return(plotting)
}
