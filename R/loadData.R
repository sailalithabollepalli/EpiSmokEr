#' @title
#' Reads idat files
#'
#' @description
#' loadData function reads idat files and saves as an RGSet object
#'
#' @param idatPath
#' Requires a path to the folder with the idat files and sample sheet included.
#'
#' @return dataset
#' A RGSet object
#'
#' @examples
#' rawdata <- loadData(idatPath = system.file("extdata", package = "EpiSmokEr"))
#' ## rawdata is a RGSet object with raw intensity values.
#'
#' @import minfi
#' @import IlluminaHumanMethylation450kmanifest
#'
#' @export
#'
loadData <- function (idatPath){
  message("Loading idat files.")
  baseDir <- file.path(idatPath)
  targets <- minfi::read.metharray.sheet(baseDir)
  RGset <- minfi::read.metharray.exp(base = baseDir, recursive = TRUE)
  message(sprintf("Dataset has %s samples.",dim(RGset)[2]))
  return(RGset)
}

