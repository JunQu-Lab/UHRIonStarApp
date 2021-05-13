#' Peptide Abundance Aggregating to Protein
#'
#' These are former functions from package 'IonStarStat' involved in peptide aggregation to protein level.
#'
#'
#' @param eset A ExpressionSet containing peptide abundance information
#' @param method Peptide aggregating method. 'sum': sum abundance; 'fit': general linear mixed model.
#' @param log2 Whether the data is with log2 form.
#' @return A 'ExpressionSet' containing protein abundance information
#' @export

ProteinQuan <- function (eset, method = "sum", log2 = TRUE)
{
  df <- exprs(eset)
  SID <- sampleNames(eset)
  fData <- featureData(eset)
  if ("Frame" %in% varLabels(fData)) {
    df1 <- data.frame(Protein = factor(fData$Protein), Peptide = factor(fData$Peptide),
                      Frame = factor(fData$Frame), exprs(eset))
    dfm <- melt(df1, c("Protein", "Peptide", "Frame"), SID)
  }
  else {
    df1 <- data.frame(Protein = factor(fData$Protein), Peptide = factor(fData$Peptide),
                      exprs(eset))
    dfm <- melt(df1, c("Protein", "Peptide"), SID)
  }
  qres <- lapply(split(dfm, dfm$Protein), function(x) qmodel(dfm1 = x,
                                                             method = method, log2 = log2))
  qres <- do.call(rbind, qres)
  qres <- qres[, c("PepNum", sampleNames(eset))]
  return(qres)
}

qmodel <- function (dfm1, method = "sum", log2 = TRUE)
{
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
