#' IonStar frame-to-protein pipeline
#'
#' These functions are IonStar frame-to-protein pipeline and its modules
#'
#' @param raw annotated frame file
#' @param cond a matrix contains file names and their group information
#' @param log2form choose if performing frame normalization on log2 level or not
#' @param norm_method choose normalization method ("quantiles" or "TIC")
#' @param quan_method choose peptide aggregation method ("sum" or "fit")
#' @return Protein and peptide quantification results
#' @export

IonStar_DPpipeline <- function(raw, cond, log2form = TRUE, norm_method, TIC_type = 'old',
                               frameOR = NULL, pepOR = NULL, OR_method = "ratio", quan_method = 'sum'){
  #raw$Label <- apply(raw[,4:ncol(raw)], 1, function(x) sum(x<2^4))
  #raw <- subset(raw, raw$Label == 0)
  #raw$Label <- NULL
  condition <- cond[match(colnames(raw)[-c(1:3)], cond[,1]),2]
  withProgress(message = 'Protein Quantification', value = 0, {
    incProgress(0.2, detail = 'Redundant Frame Removal...')
    fdata <- Frame_redundant_RM(raw)
    ndata <- Frame_norm(fdata, log2form = log2form, method = norm_method, TICtype = TIC_type)
    incProgress(0.4, detail = 'Low Quality Data Removal...')
    if(!is.null(frameOR)) ndata <- RevisedOutR(ndata, frame = TRUE, condition, variance = 0.7,critM1 = 1/3,critM2 = 1/4, method = OR_method)
    pdata <- Frame_aggre(ndata)
    if(!is.null(pepOR)) pdata <- RevisedOutR(pdata, frame = FALSE, condition, variance = 0.7,critM1 = 1/3,critM2 = 1/4, method = OR_method)
    incProgress(0.6, detail = 'Peptide Aggregation...')
    quan <- Pep_aggre(pdata, method = quan_method, log2 = TRUE)
  })
  output <- list()
  output$protein <- quan
  output$peptide <- pdata
  output$group_info <- data.frame(Rawfiles = colnames(pdata), GroupID = condition)
  output$notice <- 'Protein Quantification is finished!'
  return(output)
}
#' @export
Frame_redundant_RM <- function(df){
  df <- unique(df)
  message(paste("Input", length(unique(df[, 1])), "proteins."))
  pf <- unique(df[, c(2, 3)])
  dupFrame <- unique(pf[duplicated(pf[, 2]), 2])
  message(paste(length(dupFrame), "duplicated frames founded."))
  df <- df[!(df[, 3] %in% dupFrame), ]
  message(paste(length(unique(df[, 1])), "proteins left after filtering."))
  fData <- df[,1:3]
  df <- as.matrix(df[, -c(1:3)])
  rownames(df) <- paste(fData[, 1], fData[, 2], fData[,3], sep = ";")
  return(df)
}

#' @export
Frame_norm <- function(mat, log2form = FALSE, method, TICtype = 'old'){
  col.names <- colnames(mat)
  row.names <- rownames(mat)
  if(log2form) mat <- log2(mat + 1)
  if ((!is.null(method)) && (method == "quantiles")) {
    eval(parse(text = paste("nmethod <- ", paste("normalize",
                                                 method, sep = "."))))
    mat <- nmethod(mat)
  }
  if ((!is.null(method)) && (method == "TIC")) {
    if(log2form){
      if(TICtype == 'old'){
        average <- apply(mat, 2, mean, na.rm = TRUE)
        mean_all <- mean(average)
        correction <- average - mean_all
        mat <- t(apply(mat, 1, function(x) x - correction))
      }
      if(TICtype == 'new'){
        average <- apply(mat, 2, mean, na.rm = TRUE)
        correction <- max(average) - average
        mat2 <- t(apply(mat, 1, function(x) x + correction))
      }
    }else{
      total <- apply(mat, 2, sum, na.rm = TRUE)
      ratio <- max(total)/total
      mat <- t(apply(mat, 1, function(x) x * ratio))
    }
  }
  colnames(mat) <- col.names
  rownames(mat) <- row.names
  if(!log2form) mat <- log2(mat + 1)
  return(mat)
}

#' @export
RevisedOutR <- function(mat, frame = TRUE, condition, variance = 0.7, critM1 = 1/3, critM2 = 1/4, method = "center") {
  fData <- strsplit(rownames(mat),'\\;')
  fData <- data.frame(do.call(rbind,fData))
  if(frame){
    names(fData) <- c('Protein','Peptide','FrameID')
  }else{
    names(fData) <- c('Protein','Peptide')
  }
  mat0 <- data.frame(Protein = factor(fData$Protein), Peptide = factor(fData$Peptide), mat)
  mean <- t(apply(mat, 1, function(x) {
    tapply(x, condition, FUN = mean, simplify = TRUE)
  }))
  if (is.null(method) | method == 'none'){
    input <- data.frame(Protein = factor(fData$Protein),
                        Peptide = factor(fData$Peptide), mean)
  }
  if (method == "ratio") {
    ratio <- mean[, -1] - mean[, 1]
    input <- data.frame(Protein = factor(fData$Protein),
                        Peptide = factor(fData$Peptide), ratio)
  }
  if(method == "center"){
    frameMeans <- apply(mean, 1, mean, na.rm = TRUE)
    mean = mean - frameMeans
    input <- data.frame(Protein = factor(fData$Protein),
                        Peptide = factor(fData$Peptide), mean)
  }
  LIST <- split(input, input$Protein)
  outlier_rm <- function(list) {
    test <- as.data.frame(list)[, -c(1:2)]
    if (nrow(test) >= 3) {
      A <- pcout(test, makeplot = FALSE, explvar = variance,
                 crit.M1 = critM1, crit.c1 = 2.5, crit.M2 = critM2,
                 crit.c2 = 0.99, cs = 0.25, outbound = 0.25)$wfinal01
      return(A)
    }
    else {
      A <- rep(1, nrow(test))
      names(A) <- row.names(test)
      return(A)
    }
  }

  analysis <- lapply(LIST, outlier_rm)
  outlier_result <- data.frame((unlist(analysis)))
  row.names(outlier_result) <- sub(".*\\.", "", row.names(outlier_result))
  colnames(outlier_result) <- "outlier"
  message(paste(table(outlier_result)[1], "outliers were removed;",
                table(outlier_result)[2], "items left after outlier removal."))
  #new_df <- cbind(df0, outlier_result)
  new_df <- merge(mat0, outlier_result, by = 0)
  rownames(new_df) <- new_df$Row.names
  new_df$Row.names <- NULL
  new_df <- new_df[new_df$outlier == 1, !names(new_df) %in%
                     c("outlier")]
  return(as.matrix(new_df[,-c(1:2)]))
}

#' @export
Frame_aggre <- function(mat){
  mat <- 2^mat
  fData <- strsplit(rownames(mat),'\\;')
  fData <- data.frame(do.call(rbind,fData))
  names(fData) <- c('Protein','Peptide','FrameID')
  fn <- paste(fData$Protein, fData$Peptide,sep = ";")
  PData <- tapply(1:nrow(mat), fn, function(x) colSums(mat[x,
                                                           , drop = FALSE]))
  PData <- data.frame(do.call(rbind, PData))
  return(log2(as.matrix(PData)))
}

#' @export
Pep_aggre <- function(mat, method = "sum", log2 = TRUE){
  SID <- colnames(mat)
  fData <- rownames(mat)
  df1 <- data.frame(Protein = factor(gsub(';.*','',fData)), Peptide = factor(gsub('.*;','',fData)),
                    mat)
  dfm <- melt(df1, c("Protein", "Peptide"), SID)
  quan <- lapply(split(dfm, dfm$Protein), function(x) qmodel(dfm1 = x,
                                                             method = method, log2 = log2))
  quan <- do.call(rbind, quan)
  quan <- quan[, c("PepNum", SID)]
  return(quan)
}

#' @export
qmodel <- function (dfm1, method = "sum", log2 = TRUE){
  if (method == "fit" & log2 == TRUE) {
    lm1 <- NA
    qres <- rep(NA, length(unique(dfm1$variable)))
    if (any(duplicated(dfm1$variable))) {
      if ("Frame" %in% colnames(dfm1)) {
        try(lm1 <- MCMCglmm(value ~ 0 + variable, random = ~Frame +
                              Frame:Peptide, data = dfm1, verbose = FALSE),
            silent = TRUE)
        qres <- summary(lm1)$solutions[, 1]
      }
      else {
        try(lm1 <- MCMCglmm(value ~ 0 + variable, random = ~Peptide,
                            data = dfm1, verbose = FALSE), silent = TRUE)
        qres <- summary(lm1)$solutions[, 1]
      }
    }
    else {
      try(lm1 <- glm(value ~ 0 + variable, data = dfm1),
          silent = TRUE)
      qres <- coef(lm1)
    }
    names(qres) <- sub("variable", "", names(qres))
  }
  else {
    qres <- tapply(2^dfm1$value, as.character(dfm1$variable),
                   FUN = method)
    if (log2 == TRUE) {
      qres <- log2(qres)
    }
  }
  PepNum <- length(unique(dfm1$Peptide))
  qres <- c(PepNum = PepNum, qres)
  return(qres)
}
