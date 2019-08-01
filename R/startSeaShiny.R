#' Title Enables to start and run the app
#'
#' @return a shiny app will be launched in browser
#' @export
#' 
#' @examples
#' if(interactive()){
#' startSeaShiny()
#' }


startSeaShiny<- function(){
  appDir <- system.file("shiny-examples", "sapApp", package = "aaSEA")
  if(appDir == ""){
    stop("Could not find example directory. try reinstalling 'aaSEA", call. = FALSE)
  }
  shiny::runApp(appDir, display.mode = "normal",quiet = TRUE)
}
