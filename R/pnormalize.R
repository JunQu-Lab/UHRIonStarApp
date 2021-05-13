#' Normalization for Peptide Data
#'
#' This is a former function from package 'IonStarStat'. This function performs data normalization at the frame level.
#'
#'
#' @param eset A ExpressionSet containing peptide abundance information
#' @param summarize Summarize
#' @param method TIC: total ion current normalization; quantiles: quantile normalization.
#' @return A 'ExpressionSet' containing pepitde information
#' @export

pnormalize <- function (eset, summarize = FALSE, method)
{
  if (class(eset) != "ExpressionSet") {
    stop("Input data is not an ExpressionSet!")
  }
  if (summarize == TRUE) {
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
  }
  col.names <- colnames(exprs(eset))
  row.names <- rownames(exprs(eset))
  exprs(eset) <- log2(exprs(eset) + 1)
  if ((!is.null(method)) && (method == "quantiles")) {
    eval(parse(text = paste("nmethod <- ", paste("normalize",
                                                 method, sep = "."))))
    exprs(eset) <- nmethod(exprs(eset))
  }
  if ((!is.null(method)) && (method == "TIC")) {
    average <- apply(exprs(eset), 2, mean)
    mean_all <- mean(average)
    correction <- average - mean_all
    exprs(eset) <- t(apply(exprs(eset), 1, function(x) {
      x - correction
    }))
  }
  colnames(exprs(eset)) <- col.names
  rownames(exprs(eset)) <- row.names
  return(eset)
}
