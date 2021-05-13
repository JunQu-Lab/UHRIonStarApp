#' Outlier Rejection for Peptide Data
#'
#' This is a former function from package 'IonStarStat'.
#'
#'
#' @param eset A ExpressionSet containing peptide abundance information
#' @param condition The group information for each file (e.g., c('A','A','A','B','B','B'))
#' @param variance Variance
#' @param ratio Ratio
#' @return A 'ExpressionSet' containing pepitde information
#' @export


OutlierPeptideRM <- function (eset, condition, variance, critM1, critM2, ratio = FALSE)
{
  if (class(eset) != "ExpressionSet") {
    stop("Input data is not an ExpressionSet!")
  }
  df <- exprs(eset)
  SID <- sampleNames(eset)
  fData <- featureData(eset)
  df0 <- data.frame(Protein = factor(fData$Protein), Peptide = factor(fData$Peptide),
                    df)
  mean <- t(apply(df, 1, function(x) {
    tapply(x, condition, FUN = mean, simplify = TRUE)
  }))
  input <- data.frame(Protein = factor(fData$Protein), Peptide = factor(fData$Peptide),
                      mean)
  if (ratio == TRUE) {
    ratio <- mean[, -1] - mean[, 1]
    input <- data.frame(Protein = factor(fData$Protein),
                        Peptide = factor(fData$Peptide), ratio)
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
                table(outlier_result)[2], "peptides left after outlier removal."))
  table(outlier_result)
  new_df <- cbind(df0, outlier_result)
  new_df <- new_df[new_df$outlier == 1, !names(new_df) %in%
                     c("outlier")]
  featureData <- annotatedDataFrameFrom(as.matrix(new_df[-c(1:2)]),
                                        byrow = TRUE)
  featureData$Protein <- new_df[, 1]
  featureData$Peptide <- new_df[, 2]
  eset <- ExpressionSet(assayData = as.matrix(new_df[, -c(1:2)]),
                        featureData = featureData)
  return(eset)
}
