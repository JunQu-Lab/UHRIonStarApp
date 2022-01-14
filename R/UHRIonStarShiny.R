#' UHR-IonStar Web Application
#'
#' These functions initiate the R shiny app 'UHR-IonStar'.
#'
#'
#' @export
#'


UHRIonStarShiny <- function(){
  options(shiny.maxRequestSize=100000*1024^2)
  appDir <- system.file("shiny", "UHRIonStar",package = "UHR.IonStar")
  if (appDir == "") {
    stop("Could not find UI directory. Try re-installing `UHR.IonStar`.", call. = FALSE)
  }
  shiny::runApp(appDir, display.mode = "normal")
}
