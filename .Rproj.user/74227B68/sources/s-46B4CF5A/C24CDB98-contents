#' Protein Deconvolution Functions
#'
#' These functions can separate the protein groups to two categories: species-unique and species-shared protein groups
#'
#' @export


# one function in 'Deconvolution_check'
shared_checking <- function(cha, database){
  tmp <- grep(cha, database, value = TRUE)
  tmp_names <- gsub('sp\\|','',names(tmp))
  if(length(tmp) > 1){
    tmp2 <- gsub('.*\\_','',names(tmp))
    if(length(unique(tmp2)) > 1){
      return(c('Species-Shared',paste(tmp_names, collapse = ','), length(tmp)))
    }
    else{
      return(c('Gene-Shared',paste(tmp_names, collapse = ','), length(tmp)))
    }
  }
  else{
    return(c('',tmp_names, length(tmp)))
  }
}

# one function in 'Deconvolution_double_check'
double_check <- function(cha, database){
  tmpR <- grep(paste('R',cha, sep = ''), database, value = TRUE)
  tmpK <- grep(paste('K',cha, sep = ''), database, value = TRUE)
  tmpM <- grep(paste('M',cha, sep = ''), database, value = TRUE)
  tmpM <- tmpM[gregexpr(paste('M',cha, sep = ''), tmpM) == 1]
  tmp_namesM <- gsub('sp\\|','',names(tmpM))
  tmp_namesR <- gsub('sp\\|','',names(tmpR))
  tmp_namesK <- gsub('sp\\|','',names(tmpK))
  output <- base::union(base::union(tmp_namesK, tmp_namesR), tmp_namesM)
  len_output <- length(output)
  # if(output == 0){
  #   tmpM <- grep(paste('M',cha, sep = ''), database, value = TRUE)
  #   tmpM <- tmpM[gregexpr(paste('M',cha, sep = ''), tmpM) == 1]
  #   tmp_namesM <- gsub('sp\\|','',names(tmpM))
  #   return(c(paste(tmp_namesM, collapse = ','), length(tmpM)))
  # }else{
  return(c(paste(output, collapse = ','), len_output))
}

Deconvolution_check <- function(peptides, database){
  ###### shared peptide examination ######
  peptides$Accession <- gsub('\\|.*','',rownames(peptides))
  peptides$Pepseq <- gsub('.*\\|','',rownames(peptides))

  pep_num <- as.data.frame(table(peptides$Pepseq))
  names(pep_num) <- c('Pepseq','Pep_Num')
  peptides$Pep_Num <- pep_num$Pep_Num[match(peptides$Pepseq, pep_num$Pepseq)]
  peptides$Protein_Group_Num <- str_count(peptides$Accession, pattern = ',') + 1
  pep_unique <- peptides$Pepseq[peptides$Pep_Num == 1]

  Sys.time()
  shared_status <- lapply(pep_unique, function(x) shared_checking(x, database = database))
  Sys.time()
  output <- as.data.frame(do.call(rbind, shared_status))
  names(output) <- c('Shared_status','Rematched_Accession','Rematched_Protein_Num')
  output_shared <- data.frame(Pepseq = pep_unique, output)
  output_all <-merge(peptides, output_shared, by='Pepseq', sort = FALSE, all.x = TRUE)
  # keep the column order
  output_all <- output_all[,c(2:(ncol(peptides)-2), 1, (ncol(output_all)-4):ncol(output_all))]
  output_all$Rematched_Accession <- gsub('\\|',':',output_all$Rematched_Accession)
  return(output_all)
}

###### protein examination ######
Deconvolution_double_check <- function(peptides, database){
  # remove peptides assigned to multiple protein groups
  peptides <- peptides[!(peptides$Pep_Num > 1),]
  # remove decoys
  peptides <- peptides[!(peptides$Rematched_Accession == '0'),]
  # subset peptides that rematched_protein_num > protein_group_num
  peptides_subset <- peptides[peptides$Protein_Group_Num < peptides$Rematched_Protein_Num,]
  recheck_pep <- peptides_subset$Pepseq
  Sys.time()
  test <- lapply(recheck_pep, function(x) double_check(x, database))
  Sys.time()
  test_out <-  as.data.frame(do.call(rbind, test))
  names(test_out) <- c('Digestion_recheck_accession', 'Digestion_recheck_num')
  peptides_subset <- cbind(peptides_subset, test_out)

  peptides_subset$Digestion_recheck_accession <- gsub('\\|', ':', peptides_subset$Digestion_recheck_accession)
  peptides$Rematched_Accession[match(peptides_subset$Pepseq, peptides$Pepseq)] <- peptides_subset$Digestion_recheck_accession
  peptides$Rematched_Protein_Num[match(peptides_subset$Pepseq, peptides$Pepseq)] <- peptides_subset$Digestion_recheck_num

  # share status update
  for (i in 1:nrow(peptides)) {
    tmp <- peptides$Rematched_Accession[i]
    tmp2 <- gsub('.*\\_','',unlist(strsplit(tmp, ',')))
    if(length(unique(tmp2)) > 1){
      peptides$Shared_status[i] <- 'Species-Shared'
    }else{
      if(length(tmp2)>1){
        peptides$Shared_status[i] <- 'Gene-Shared'
      }else{
        peptides$Shared_status[i] <- ''
      }
    }
  }
  return(peptides)
}

# IonStar streamline after shared peptide removal
peptide_output <- function(peptides, cond){
  outputs <- list()
  # peptide outputs
  # 1. Unique protein groups
  peptide_out <- peptides[peptides$Shared_status != 'Species-Shared',]
  peptide_out_id <- peptide_out[,c('Rematched_Accession','Pepseq')]
  peptide_out_frame <- 1:nrow(peptide_out)

  peptide_out <- cbind(peptide_out_id, peptide_out_frame, peptide_out[,1:(ncol(peptide_out)-7)])
  names(peptide_out)[1:3] <- c('Accession','Peptide.Sequence','FrameID')

  raw <- peptide_out
  condition <- cond[match(colnames(raw)[-c(1:3)], cond[,1]),2]

  pdata <- newProDataSet(proData=raw, condition=condition)
  quan <- ProteinQuan(pdata, method="sum")

  outputs$unique_protein <- cbind(as.data.frame(quan), Species_Shared = FALSE)
  pep_tmp <- exprs(pdata)
  rownames(pep_tmp) <- gsub('\\|[0-9].*', '', rownames(pep_tmp))
  outputs$unique_peptide <- cbind(as.data.frame(pep_tmp), Species_Shared = FALSE)

  # 2. Species-mixed protein groups
  peptide_out2 <- peptides[peptides$Shared_status == 'Species-Shared',]
  peptide_out_id <- peptide_out2[,c('Rematched_Accession','Pepseq')]
  peptide_out_frame <- 1:nrow(peptide_out2)

  peptide_out2 <- cbind(peptide_out_id, peptide_out_frame, peptide_out2[,1:(ncol(peptide_out2)-7)])
  names(peptide_out2)[1:3] <- c('Accession','Peptide.Sequence','FrameID')

  raw <- peptide_out2
  condition <- cond[match(colnames(raw)[-c(1:3)], cond[,1]),2]

  pdata <- newProDataSet(proData=raw, condition=condition)
  quan <- ProteinQuan(pdata, method="sum")

  outputs$shared_protein <- cbind(as.data.frame(quan), Species_Shared = TRUE)
  pep_tmp <- exprs(pdata)
  rownames(pep_tmp) <- gsub('\\|[0-9].*', '', rownames(pep_tmp))
  outputs$shared_peptide <- cbind(as.data.frame(pep_tmp), Species_Shared = TRUE)

  return(outputs)
}
