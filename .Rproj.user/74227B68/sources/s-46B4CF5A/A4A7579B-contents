#' UHR-IonStar Source Functions
#'
#' These functions are the application functions in the R shiny web app 'UHR-IonStar'.
#'
#'
#' @export

### IonStarStat, Frame Generation
FrameGen <- function(db,spectrum,col_filename,col_scannum,accession, pepseq){

  withProgress(message = 'Frame Generation', value = 0, {
    incProgress(0.2, detail = 'frame extraction...')
    output <- list()
    col_filename <- as.numeric(col_filename)
    col_scannum <- as.numeric(col_scannum)
    col_framelist <- as.numeric(c(accession,pepseq))
    sql <- "select frames.frameid,(frames.timestart+frames.timestop)/2 as Time,(frames.mzstart+frames.mzstop)/2 as MZ, ms2scans.ms2scan, ms2scans.RawFileID, FrameAttrib_1.* from Frames
  INNER JOIN MS2Scans on MS2Scans.FrameID=Frames.frameid
  INNER JOIN FrameAttrib_1 ON FrameAttrib_1.FrameID=Frames.FrameID
  order by Frames.Frameid"
    frame <- dbGetQuery(db, sql)
    RawFiles <- dbGetQuery(db, "select * from RawFiles")
    dbDisconnect(db)
    raw <- RawFiles[match(frame[,"RawFileID"], RawFiles[,"RFID"]), "RawFile"]
    raw <- sub(".raw", "", raw)
    frame <- cbind(frame, RAWFILE=raw)
    frame[,"RAWFILE"]<-toupper(sub(".RAW","",frame[,"RAWFILE"]))
    nframe <- toupper(paste(frame[,"RAWFILE"], frame[,"MS2Scan"], sep="_"))

    incProgress(0.8, detail = 'merging with spectrum report...')
    spectrum$SpecFile <- apply(as.matrix(spectrum$SpecFile),1,function(x) gsub(".raw.mzXML","",x))
    nspec <- toupper(paste(spectrum[,col_filename], spectrum[,col_scannum], sep="_"))
    idx <- match(nspec, nframe)
    message1 <- paste("Total number of spectra with matching frames:",sum(!is.na(idx)))
    message2 <- paste("Total number of spectra without matching frames:",sum(is.na(idx)))
    print(message1)
    print(message2)
    FS <- cbind(spectrum[!is.na(idx),], frame[na.omit(idx),])
    #FS[,col_framelist[1]]<- sub("\\,.*", "", FS[,col_framelist[1]])
    FS[,col_framelist[1]]<- toupper(sub("sp\\|", "", FS[,col_framelist[1]]))
    FS[,col_framelist[1]]<- toupper(sub("\\|", ":", FS[,col_framelist[1]]))
    FS[,col_framelist[2]]<- toupper(FS[,col_framelist[2]])
    FS$ID<-paste(FS[,col_framelist[1]],FS[,col_framelist[2]],FS[,"FrameID"],sep="|")
    FS2<-FS[!duplicated(FS$ID),]
    row.names(FS2)<-FS2$ID

    sub_start <- ncol(spectrum)+12
    sub_end <- ncol(spectrum)+11+length(unique(frame[,"RAWFILE"]))
    sub_data<-FS2[,c(col_framelist,(ncol(spectrum)+1),sub_start:sub_end)]

    sampleid_list <- colnames(sub_data[,4:ncol(sub_data)])
    sampleid_mat <- as.data.frame(matrix(c(sampleid_list,rep("NA",times=length(sampleid_list))),ncol=2,nrow=length(sampleid_list)))
    colnames(sampleid_mat) <- c("Rawfiles","GroupID")
  })

  output$framelist <- sub_data
  output$sampleid <- sampleid_mat
  output$notice <- 'Annotated Frames and Sample List have been generated.'
  output$message1 <- message1
  output$message2 <- message2
  return(output)
}

### IonStarStat, Protein Quantification
PSummerize <- function(eset){
  if (class(eset) != "ExpressionSet") {
    stop("Input data is not an ExpressionSet!")
  }
  df <- exprs(eset)
  fn <- paste(featureData(eset)$Protein, featureData(eset)$Peptide,
              sep = "|")
  PData <- tapply(1:nrow(df), fn, function(x) colSums(df[x,
                                                         , drop = FALSE]))
  PData <- do.call(rbind, PData)
  featureData <- annotatedDataFrameFrom(PData, byrow = TRUE)
  pp <- do.call(rbind, strsplit(rownames(PData), split = "\\|"))
  featureData$Protein <- pp[, 1]
  featureData$Peptide <- pp[, 2]
  eset_2 <- ExpressionSet(PData)
  phenoData(eset_2) <- phenoData(eset)
  eset <- eset_2
  featureData(eset) <- featureData
  col.names <- colnames(exprs(eset))
  row.names <- rownames(exprs(eset))
  exprs(eset) <- log2(exprs(eset) + 1)
  colnames(exprs(eset)) <- col.names
  rownames(exprs(eset)) <- row.names
  return(eset)
}

IonStarFunc <- function(raw, cond, norm_method = 'TIC', OutlierRM = 'A', quan_method = 'sum'){
  outputs <- list()
  raw$shared <- NULL
  # remove several rows with missing data
  raw$Label <- apply(raw[,4:ncol(raw)], 1, function(x) sum(x<2^4))
  raw <- subset(raw, raw$Label == 0)
  raw$Label <- NULL
  condition <- cond[match(colnames(raw)[-c(1:3)], cond[,1]),2]

  pdata <- newProDataSet(proData=raw, condition=condition)
  if(norm_method != 'TIC' & norm_method != 'quantiles'){
    ndata <- pnormalize(pdata, summarize=TRUE, method=norm_method)
  }else{
    ndata <- PSummerize(pdata)
  }
  if(!is.null(OutlierRM)){
    if(length(unique(condition)) <= 2){
      ndata<- OutlierPeptideRM(ndata,condition,variance=0.7,critM1=1/3,critM2=1/4,ratio=FALSE)
    }else{
      ndata<- OutlierPeptideRM(ndata,condition,variance=0.7,critM1=1/3,critM2=1/4,ratio=TRUE)
    }
  }
  quan <- ProteinQuan(eset=ndata, method=quan_method)
  outputs$quan <- quan
  outputs$pep <- exprs(ndata)
  return(outputs)
}

ProQuant <- function(raw, cond, database, rm, norm_method, quan_method, outlier){

  withProgress(message = 'Protein Quantification', value = 0, {
    incProgress(0.2, detail = 'frame cleanup...')
    outputs <- list()
    raw$shared <- NULL
    # remove several rows with missing data
    raw$Label <- apply(raw[,4:ncol(raw)], 1, function(x) sum(x<2^4))
    raw <- subset(raw, raw$Label == 0)
    raw$Label <- NULL
    condition <- cond[match(colnames(raw)[-c(1:3)], cond[,1]),2]
    incProgress(0.4, detail = 'normalization and outlier rejection...')
    pdata <- newProDataSet(proData=raw, condition=condition)
    if(norm_method != 'none'){
      ndata <- pnormalize(pdata, summarize=TRUE, method=norm_method)
    }else{
      ndata <- PSummerize(pdata)
    }
    if(!is.null(outlier)){
      if(length(unique(condition)) <= 2){
        ndata<- OutlierPeptideRM(ndata,condition,variance=0.7,critM1=1/3,critM2=1/4,ratio=FALSE)
      }else{
        ndata<- OutlierPeptideRM(ndata,condition,variance=0.7,critM1=1/3,critM2=1/4,ratio=TRUE)
      }
    }
    ndata<-SharedPeptideRM(ndata)
    incProgress(0.6, detail = 'protein quantification...')
    quan <- ProteinQuan(eset=ndata, method=quan_method)
    if(is.null(rm)){
      outputs$quan <- quan
      outputs$pep <- exprs(ndata)
    }else{
      incProgress(0.8, detail = 'data deconvolution...')
      peptides_checked <- Deconvolution_check(as.data.frame(exprs(ndata)), database = database)
      peptides_double_check <- Deconvolution_double_check(peptides_checked, database)
      final_out <- peptide_output(peptides_double_check, cond)

      outputs$quan <- rbind(final_out$unique_protein, final_out$shared_protein)
      outputs$pep <- rbind(final_out$unique_peptide, final_out$shared_peptide)
    }
  })
  outputs$notice <- "Protein Quantification Completed."
  return(outputs)
}

### Data Processing, Ratio_p
Ratio_p <- function(quan, group, decoy, control_group=NA, is_t_test="t", pro_or_pep){
  outputs <- list()
  if ('PepNum' %in% names(quan)){
    quan$PepNum <- NULL
  }
  # The first character in colnames cannot be a "number".
  #names(quan) <- gsub("X","",names(quan))
  ##Data reordering and cleanup
  quan_order <- quan[,-1]
  group_names <<- as.character(group[,1]) # for PCA plot
  sample_order <- match(colnames(quan_order),group_names)
  quan_results <- cbind(quan[,1],quan_order[,sample_order])


  #Remove decoy protein entries
  decoy_id <- decoy
  quan_results <- quan_results[grep(
    decoy_id, quan_results[, 1], fixed = FALSE, invert = TRUE),]

  #Remove protein entries with missing data
  quan_results$mv.num <- apply(
    quan_results[, -1], 1, function(x) sum(x<4))

  quan_results <- subset(quan_results, mv.num==0)
  quan_results <- quan_results[,1:(ncol(quan_results)-1)]

  ##Ratio & p-value calculation
  sample_num <<- nrow(group)
  Group_id <<- as.character(group[,2]) # for PCA plot
  group_num <<- length(unique(Group_id))
  rep_num <<- sample_num/group_num
  quan_matrix <- quan_results[,-1]
  sample_name <- as.character(group[,2])
  group_char <<- unique(sample_name)
  rownames(quan_matrix) <- quan_results[,1]
  colnames(quan_matrix) <- sample_name
  quan_matrix <<- quan_matrix #update to global variable
  quan_matrix2 <- matrix(sapply(quan_matrix,function(x) 2^x),nrow=nrow(quan_matrix))
  rownames(quan_matrix2) <- rownames(quan_matrix)
  colnames(quan_matrix2) <- colnames(quan_matrix)
  quan_matrix2 <<- as.data.frame(quan_matrix2)
  #Assign the control group, substitute with input from Shiny
  #Calculate ratios between case and control conditions

  ave_df <<- as.data.frame(avearrays(quan_matrix2, ID=colnames(quan_matrix), weights = NULL))

  con <- match(control_group,group_char)
  if(is.na(con)){
    con <- 1
    control_group <- as.character(unique(as.character(group[,2])))[1]
  }
  control_group <<- control_group # turn to global variable
  ave_ratio <- ave_df/ave_df[,con] # plot ratio distribution
  r_chara <- rep(NA,group_num)
  for (i in 1:group_num) {
    r_chara[i] <- paste("Ratio_",group_char[i],"/",control_group,sep="")

  }
  withProgress(message = 'Data Processing', value = 0, {
    incProgress(0.5, detail = 'statistical tests...')
    #Paired t.test between case and control conditions
    if(is_t_test=='t'){
      p_chara <- rep(NA,group_num)
      tmp <- as.character(unique(Group_id))
      for (i in 1:group_num) {
        pvals <- apply(quan_matrix, 1, function(x)
          t.test(x[which(names(quan_matrix)== control_group)],
                 x[which(names(quan_matrix)== tmp[i])],var.equal = TRUE)$p.value)
        p_chara[i] <- paste("pVal_",group_char[i],"/",control_group,sep="")
        ave_ratio <- cbind(ave_ratio,pvals)
      }
    }

    if(is_t_test=='paired-t'){
      p_chara <- rep(NA,group_num)
      tmp <- as.character(unique(Group_id))
      for (i in 1:group_num) {
        pvals <- apply(quan_matrix, 1, function(x)
          t.test(x[which(names(quan_matrix)== control_group)],
                 x[which(names(quan_matrix)== tmp[i])],var.equal = TRUE, paired = TRUE)$p.value)
        p_chara[i] <- paste("pVal_",group_char[i],"/",control_group,sep="")
        ave_ratio <- cbind(ave_ratio,pvals)
      }
    }


    if(pro_or_pep == 'pro'){
      colnames(ave_ratio) <- c(r_chara,p_chara)
      ave_ratio <<- ave_ratio[,-c(con,(con+group_num))]
      colnames(quan_results)[1] <- "ProteinAC"
      quan_results <<- cbind(quan_results,ave_ratio)
      quan_all <- merge(quan_results,ave_ratio, by.x = 'ProteinAC', by.y = 'row.names')
      tmp <- as.data.frame(matrix(rep(NA,2*nrow(quan_all)), ncol = 2))
      tmp[,1] <- gsub("\\:.*","", quan_all[,1])
      tmp[,2] <- gsub("^.*?:","", quan_all[,1])
      names(tmp) <- c('Accession','ProteinID')
      quan_all2 <- cbind(tmp,quan_all[,-1])
    }
    if(pro_or_pep == 'pep'){
      colnames(ave_ratio) <- c(r_chara,p_chara)
      ave_ratio <<- ave_ratio[,-c(con,(con+group_num))]
      colnames(quan_results)[1] <- "PeptideInfo"
      quan_results <<- cbind(quan_results,ave_ratio)
      quan_all <- merge(quan_results,ave_ratio, by.x = 'PeptideInfo', by.y = 'row.names')
      tmp_ProteinAC <- sub('\\|.*','',quan_all[,1])
      tmp <- as.data.frame(matrix(rep(NA,3*nrow(quan_all)), ncol = 3))
      tmp[,1] <- gsub("\\:.*","", tmp_ProteinAC)
      tmp[,2] <- gsub("^.*?:","", tmp_ProteinAC)
      tmp[,3] <- sub('.*\\|','',quan_all[,1])
      names(tmp) <- c('Accession','ProteinID','Peptide_Sequence')
      quan_all2 <- cbind(tmp,quan_all[,-1])
    }
  })
  outputs$quan_results <- quan_results
  outputs$ave_ratio <- ave_ratio
  outputs$all <- quan_all2
  outputs$control_group <- control_group
  outputs$group_char <- group_char
  outputs$notice <- "Processing completed."
  return(outputs)

}

### Data Visualization, plot_for_all
plot_for_all <- function(plot_selection, Group1, Group2,corr){
  if(plot_selection == "Intra-group CV boxplot"){

    RSD <- numeric()
    for (m in 1:group_num) {
      quan_val <- quan_matrix2[,(rep_num*(m-1)+1):(rep_num*m)]
      SD <- apply(quan_val,1,function(x) sd(x))
      RSD_ind <- SD/ave_df[,m]*100
      RSD <- c(RSD,RSD_ind)
    }
    RSD_save <- as.data.frame(matrix(RSD,ncol=group_num,byrow=FALSE))
    c_chara <- character()
    for (n in 1:group_num) {
      c_chara2 <- paste("CV.",group_char[n],sep="")
      c_chara <- c(c_chara,c_chara2)
    }
    colnames(RSD_save) <- c_chara
    rownames(RSD_save) <- row.names(quan_matrix)
    Group_plot <- rep(group_char,each=nrow(quan_matrix))
    RSD_plot <- as.data.frame(matrix(c(RSD,Group_plot),ncol=2,byrow=FALSE))
    colnames(RSD_plot) <- c("CV","Group")
    RSD_plot$CV <- apply(as.matrix(RSD_plot$CV),1,as.numeric)


    selected_plot <- ggplot(RSD_plot,aes(x=Group,y=CV,fill=Group)) + geom_boxplot() +
      scale_fill_brewer(palette="Blues") + scale_x_discrete(limits=unique(RSD_plot$Group)) +
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
    MeanInt <- apply(as.matrix(quan_matrix2),1,function(x) log2(mean(x)))
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
    #Plotting correlation between two replicates/groups

    #Inter-replicate
    #Group1 <- "A1"
    #Group2 <- "C2"
    #Corr_plot <-quan_results[,match(c(Group1,Group2),colnames(quan_results))]

    #Inter-group
    con2 <- match(c(Group1,Group2),colnames(ave_df))
    if(NA %in% con2){
      Group1 <- colnames(ave_df)[1]
      Group2 <- colnames(ave_df)[2]
    }

    Corr_plot <- ave_df[,match(c(Group1,Group2),colnames(ave_df))]
    Corr_plot <- as.data.frame(matrix(apply(Corr_plot,c(1,2),log2),ncol=2,byrow=FALSE))
    #Data processing and plotting
    colnames(Corr_plot) <- c("G1","G2")
    corr <- round(cor(Corr_plot$G1,Corr_plot$G2),digits=3)
    corr_label <- paste("R^2=",corr,sep="")
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
    Pearson <- quan_results[,2:(sample_num+1)]
    rownames(Pearson) <- quan_results$ProteinAC
    Correlation <- cor(Pearson)
    selected_plot <- ggcorrplot(Correlation,colors = c("royalblue3", "gray90", "brown3"),
                                type = "lower", method = "circle") +
      theme(axis.text=element_text(size=14, family="Helvetica", colour="black"),
            axis.title=element_text(size=14, family="Helvetica",colour="black"),
            legend.text=element_text(size=14, family="Helvetica",colour="black"),
            legend.title=element_text(size=14, family="Helvetica", colour="black")) +
      scale_fill_gradientn(name="Correlation", colors=c('#f6ff50','#5eff50','#50fff6','#5055ff'), limits=c(as.numeric(corr),1),
                           na.value = 'gray90')
    #cor_plot <- cor_plot + scale_fill_gradientn(colors=c('#fff050','#FF5733','#C70039','#900C3F'),limits=c(0.93,1),na.value = 'gray90')

  }
  if(plot_selection == "Principal Component Analysis (PCA)"){
    PCA_data <- quan_results[,match(group_names,colnames(quan_results))]
    PCA_data <- as.data.frame(t(PCA_data))
    colnames(PCA_data) <- quan_results[,1]
    PCA_results <- prcomp(PCA_data,center=TRUE,scale.=TRUE)

    #tmp <- PCA_results$x[,1]
    #PCA_results$x[,1] <- PCA_results$x[,2]
    #PCA_results$x[,2] <- tmp

    PCA_level <- Group_id
    selected_plot <- ggbiplot(PCA_results, obs.scale=1, var.scale=1, groups=PCA_level,
                              circle=FALSE, var.axes=FALSE) +
      geom_point(aes(colour=PCA_level),size=5) + scale_colour_brewer(name="Groups",palette="Set2") +
      #geom_hline(yintercept=0) + geom_vline(xintercept=0) +
      geom_label_repel(aes(label=group_names),color="black",size=3,fill=NA, label.size=NA) +
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
plot_ratio_distribution <- function(selected, dist_corr){
  outputs <- list()
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
      #      dens <- density(ratio_matrix$value)
      #      correc <- dens$x[which(dens$y == max(dens$y))]
      #      ratio_plot$value <- ratio_plot$value - correc
    }
  }
  corr_ratio <<- ratio_matrix
  ratio_plot <- melt(ratio_matrix)
  #  ratio_plot$value <- apply(as.matrix(ratio_plot$value),1,log2)

  #  if(!is.null(dist_corr)){
  #    if(dist_corr == "Correction"){
  #    dens <- density(ratio_plot$value)
  #    correc <- dens$x[which(dens$y == max(dens$y))]
  #    ratio_plot$value <- ratio_plot$value - correc
  #    }
  #  }
  axis_limit <- ceiling(max(max(ratio_plot$value),abs(min(ratio_plot$value))))

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
  outputs$correction_ratio <- corr_ratio
  return(outputs)

}

### Discovery of Changes, sign_changes
sign_changes <- function(chosen_group, R_cutoff=1.4, p_cutoff=0.05){
  outputs <- list()
  R_col <- paste("Ratio_",chosen_group,"/",control_group,sep="")
  p_col <- paste("pVal_",chosen_group,"/",control_group,sep="")
  quan_select <- quan_results[,c(1,match(c(R_col,p_col),colnames(quan_results)))]
  colnames(quan_select) <- c("ProteinAC","Ratio","pVal")
  # if(pro_or_pep == 'pep'){
  #  quan_select$ProteinAC <- sub('.*\\|','',quan_select$ProteinAC)
  #}
  R_upper <- as.numeric(R_cutoff)
  R_lower <- 1/as.numeric(R_cutoff)
  quan_up <<- subset(quan_select, Ratio>=R_upper & pVal<=as.numeric(p_cutoff))
  quan_down <<- subset(quan_select, Ratio<=R_lower & pVal<=as.numeric(p_cutoff))
  quan_changed <<- rbind(quan_up,quan_down)
  quan_select$match <- quan_select$ProteinAC %in% quan_changed$ProteinAC
  quan_unchanged <<- subset(quan_select, match=="FALSE")
  outputs$quan_up <- quan_up
  outputs$quan_down <- quan_down
  outputs$quan_changed <- quan_changed
  outputs$quan_unchanged <- quan_unchanged
  return(outputs)

}

### Discovery of Changes, plot_for_all2
plot_for_all2 <- function(plot_selection,TOPx = 50){
  if(plot_selection == "Valcano Plot"){
    # Volcano plot of data distribution
    quan_up$Ratio <- apply(as.matrix(quan_up$Ratio),1, log2)
    quan_up$pVal <- apply(as.matrix(quan_up$pVal),1,function(x) -log10(x))
    quan_down$Ratio <- apply(as.matrix(quan_down$Ratio),1, log2)
    quan_down$pVal <- apply(as.matrix(quan_down$pVal),1,function(x) -log10(x))
    quan_unchanged$Ratio <- apply(as.matrix(quan_unchanged$Ratio),1, log2)
    quan_unchanged$pVal <- apply(as.matrix(quan_unchanged$pVal),1,function(x) -log10(x))
    quan_unchanged <- quan_unchanged[,1:3]
    quan_up$Group <- "Up-changed"
    quan_down$Group <- "Down-changed"
    quan_unchanged$Group <- "Unchanged"
    quan_vol <- rbind(quan_up,quan_down,quan_unchanged)
    up_col <- "red"
    down_col <- "green"
    cols <- c(down_col,"black",up_col)
    x_limit <- max(max(quan_vol$Ratio),min(quan_vol$Ratio))
    plotting <- ggplot(quan_vol,aes(x=Ratio,y=pVal,fill=Group)) + geom_point(alpha=0.4,shape=21) +
      scale_fill_manual(values=cols) + scale_x_continuous(limits=c(-6,6)) + scale_y_continuous(limits=c(0,3.5)) +
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
    quan_matrix2$Groups <- rownames(quan_matrix2) %in% quan_changed$ProteinAC
    quan_matrix2$Log2Intensities <-apply(as.matrix(quan_matrix2),1,function(x) log2(mean(x)))
    quan_matrix3 <- quan_matrix2[order(quan_matrix2$Log2Intensities,decreasing=TRUE),]
    quan_matrix3$ProteinNumber <- 1:nrow(quan_matrix3)
    quan_matrix3 <- quan_matrix3[,c(ncol(quan_matrix3)-2,ncol(quan_matrix3)-1, ncol(quan_matrix3))]
    quan_matrix3[which(quan_matrix3$Groups == FALSE),1] <- "Unchanged"
    quan_matrix3[which(quan_matrix3$Groups == TRUE),1] <- "Changed"
    plotting <- ggplot(quan_matrix3,aes(x=ProteinNumber,y=Log2Intensities)) +
      geom_point(alpha=1,aes(colour=Groups, shape=Groups, size=Groups)) +
      scale_colour_manual(values=c("red","black")) +
      scale_shape_manual(values=c(17,1))+
      scale_size_manual(values=c(4,2))+
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
    quan_up$acName <- substr(quan_up$ProteinAC, start = 1, stop = 15)
    quan_up <- quan_up[order(quan_up$Fold,decreasing=TRUE),]
    Up_plot <- quan_up[1:min(nrow(quan_up),as.numeric(TOPx)),]
    plotting <- ggplot(Up_plot,aes(x=reorder(acName,-Fold),y=Fold)) +
      geom_bar(stat="identity",colour="black",fill="brown2",alpha=0.7) +
      theme(legend.position="none",panel.grid.minor=element_blank(),
            panel.background=element_rect(fill=NA, colour=NA),
            panel.border=element_rect(size=1,fill=NA,colour="black")) +
      theme(axis.text=element_text(size=12, family="Helvetica", colour="black"),
            axis.text.x=element_text(angle=90,hjust=1),
            axis.title=element_text(size=12, family="Helvetica",colour="black"),
            legend.text=element_text(size=12, family="Helvetica",colour="black"),
            legend.title=element_text(size=12, family="Helvetica", colour="black")) +
      labs(x="",y="Fold change")
  }
  if(plot_selection == "Plot for Down-regulated"){
    quan_down$Fold <-as.matrix(quan_down$Ratio)
    quan_down$acName <- substr(quan_down$ProteinAC, start = 1, stop = 15)
    quan_down <- quan_down[order(quan_down$Fold,decreasing=FALSE),]
    Down_plot <- quan_down[1:min(nrow(quan_down),as.numeric(TOPx)),]
    plotting <- ggplot(Down_plot,aes(x=reorder(acName,Fold),y=Fold)) +
      geom_bar(stat="identity",colour="black",fill="springgreen2",alpha=0.7) +
      theme(legend.position="none",panel.grid.minor=element_blank(),
            panel.background=element_rect(fill=NA, colour=NA),
            panel.border=element_rect(size=1,fill=NA,colour="black")) +
      theme(axis.text=element_text(size=12, family="Helvetica", colour="black"),
            axis.text.x=element_text(angle=90,hjust=1),
            axis.title=element_text(size=12, family="Helvetica",colour="black"),
            legend.text=element_text(size=12, family="Helvetica",colour="black"),
            legend.title=element_text(size=12, family="Helvetica", colour="black")) +
      labs(x="",y="Fold change")
  }
  if(plot_selection == "Gene Ontology (Accession)"){
    annotation="GOTERM_BP_DIRECT"
    pvalueCutoff=1
    qvalueCutoff=1
    minGSSize=5

    quan_sum <<- rbind(quan_up,quan_down)
    quan_sum$ProteinAC <- gsub("\\:.*","",quan_sum$ProteinAC)

    sum_list <- as.character(quan_sum$ProteinAC)
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
    quan_sum$ProteinAC <- gsub(".*\\:","",quan_sum$ProteinAC)

    sum_list <- as.character(quan_sum$ProteinAC)
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

### Re-verification, ProPepCombine
ProPepCombine <- function(subpro,subpep){
  outputs <- list()
  subpro$Accession <- as.character(subpro$Accession)
  subpep$Accession <- as.character(subpep$Accession)
  subpro$ProteinID <- as.character(subpro$ProteinID)
  subpep$ProteinID <- as.character(subpep$ProteinID)
  subpep$Peptide_Sequence <- as.character(subpep$Peptide_Sequence)

  TRpro <- which(subpro$Accession == 'TR')
  for (i in TRpro) {
    subpro$Accession[i] <- paste(subpro$Accession[i],subpro$ProteinID[i], sep = ':')
  }
  TRpep <- which(subpep$Accession == 'TR')
  for (i in TRpep) {
    subpep$Accession[i] <- paste(subpep$Accession[i],subpep$ProteinID[i], sep = ':')
  }

  subpep$ProteinID <- NULL
  names(subpep) <- names(subpro)


  Output_df <- rbind(subpro[1,],subpep[which(subpep$Accession == subpro$Accession[1]),])
  for (j in 2:nrow(subpro)) {
    tmp <- which(subpep$Accession == subpro$Accession[j])
    tmp_output <- rbind(subpro[j,],subpep[tmp,])
    Output_df <- rbind(Output_df,tmp_output)
  }

  name_tmp <- names(Output_df)
  name_tmp[which(name_tmp == 'ProteinID')] <- 'ProteinID/PepSeq'
  names(Output_df) <- name_tmp

  outputs$df <- Output_df
  outputs$notice <- 'Integration complicated.'
  return(outputs)
}

### Re-verification, CustomSearchPlot
CustomSearchPlot <- function(df,Accession){
  output <- list()
  subtest_df <- df[,c(1,2,grep('Ratio',names(df)))]
  tmp_df <- subtest_df[which(subtest_df$Accession == Accession),]
  tmp_df <- tmp_df[,-1]
  tmp_df2 <- melt(tmp_df)
  names(tmp_df2)<- c('Name','group','RatioValue')
  tmp_df2$group <- as.character(tmp_df2$group)
  Search_plot <- ggplot(tmp_df2,aes(x = group, y = RatioValue, fill = Name)) +
    geom_bar(stat="identity",position=position_dodge()) +
    theme(axis.text=element_text(size=14, family="Helvetica", colour="black"),
          axis.title=element_text(size=14, family="Helvetica",colour="black"),
          legend.text=element_text(size=14, family="Helvetica",colour="black"),
          legend.title=element_text(size=14, family="Helvetica", colour="black"))

  Search_plot_nolegend <- ggplot(tmp_df2,aes(x = group, y = RatioValue, fill = Name)) +
    geom_bar(stat="identity",position=position_dodge()) +
    theme(legend.position="none",
          axis.text=element_text(size=14, family="Helvetica", colour="black"),
          axis.title=element_text(size=14, family="Helvetica",colour="black"),
          legend.text=element_text(size=14, family="Helvetica",colour="black"),
          legend.title=element_text(size=14, family="Helvetica", colour="black"))
  if(nrow(tmp_df) > 6){
    output$plot1 <- Search_plot
    output$plot2 <- Search_plot_nolegend
  }
  else{
    output$plot1 <- Search_plot
    output$plot2 <- Search_plot
  }
  return(output)
}


