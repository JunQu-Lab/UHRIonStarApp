#' Shared Frames Removal
#'
#' This are former functions from package 'IonStarStat'. This function removes the frames that are assigned to
#' multiple protein groups.
#'
#' @param eset A ExpressionSet containing peptide abundance information
#' @return A 'ExpressionSet' containing pepitde information
#' @export


SharedPeptideRM <- function (eset)
{
  pepdata <- exprs(eset)
  text <- row.names(pepdata)
  str <- strsplit(text, "|", fixed = TRUE)
  str <- do.call(rbind.data.frame, str)
  stringa <- str[, c(1:2)]
  colnames(stringa) <- c("ProteinAC", "Peptideseq")
  newdata <- cbind(stringa, pepdata)
  row.names(newdata) <- NULL
  n_occur <- data.frame(table(newdata$Peptideseq))
  newdata2 <- newdata[newdata$Peptideseq %in% n_occur$Var1[n_occur$Freq ==
                                                             1], ]
  row.names(newdata2) <- paste(newdata2[, 1], newdata2[, 2],
                               sep = "|")
  featureData <- annotatedDataFrameFrom(as.matrix(newdata2),
                                        byrow = TRUE)
  featureData$Protein <- newdata2[, 1]
  featureData$Peptide <- newdata2[, 2]
  eset <- ExpressionSet(assayData = as.matrix(newdata2[, -c(1:2)]),
                        featureData = featureData)
}
