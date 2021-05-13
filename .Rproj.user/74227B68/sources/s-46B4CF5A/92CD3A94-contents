#' Format the Frame Dataset
#'
#' This is a former function from package 'IonStarStat'. This function firstly removes frames
#' that are assigned to different protein groups, then format the dataset to 'ExpressionSet'.
#'
#' @param proData A frame data set
#' @param condition The group information for each file (e.g., c('A','A','A','B','B','B'))
#' @param phenoData Former-designed phenotype (default: NULL)
#' @param featureData Former-designed feature (default: NULL)
#' @return A 'ExpressionSet' containing frame (peptide) information
#' @export


newProDataSet <- function (proData, condition, phenoData = NULL, featureData = NULL)
{
  proData <- unique(proData)
  message(paste("Input", length(unique(proData[, 1])), "proteins."))
  pf <- unique(proData[, c(2, 3)])
  dupFrame <- unique(pf[duplicated(pf[, 2]), 2])
  message(paste(length(dupFrame), "duplicated frames founded."))
  proData <- proData[!(proData[, 3] %in% dupFrame), ]
  message(paste(length(unique(proData[, 1])), "proteins left after filtering."))
  fData <- proData[, 1:3]
  PData <- as.matrix(proData[, -c(1:3)])
  rownames(PData) <- paste(fData[, 1], fData[, 2], fData[,
                                                         3], sep = "|")
  colnames(fData) <- c("Protein", "Peptide", "Frame")
  if (is.null(phenoData)) {
    phenoData <- annotatedDataFrameFrom(PData, byrow = FALSE)
  }
  if (is.null(featureData)) {
    featureData <- annotatedDataFrameFrom(PData, byrow = TRUE)
  }
  featureData$Protein <- fData[, 1]
  featureData$Peptide <- fData[, 2]
  featureData$Frame <- fData[, 3]
  phenoData$condition <- condition
  prodat <- ExpressionSet(assayData = PData, phenoData = phenoData,
                          featureData = featureData)
}
