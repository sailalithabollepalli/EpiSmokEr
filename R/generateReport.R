#' @title
#' Generates a html and a csv file with results from epismoker function.
#'
#' @description
#' Generates a html and a csv file with results from epismoker function.
#'
#' @param resObj
#' A result object produced from epismoker function
#'
#' @param ouputFile
#' Name of the output file provided by the user
#'
#' @param method
#' Smoking score method determined by user in the epismoker function to produce the result object.
#' Users can choose one of these four methods ( "SSc", "MS", "SSt", "all") in the epismoker function.
#'
#' @return html file
#' Returns a html file with the results from the epismoker function, saved to the Results folder.
#'
#' @return csv file
#' Returns a csv file with the results from the epismoker function, saved to the Results folder.
#'
#' @examples
#' data(dummyBetaData)
#' result <- epismoker(dataset = dummyBetaData, ref.Zhang = Zhangetal_cpgs, method = "MS")
#' generateReport(result, outputFileName = "results_MS", method = "MS")
#' ## result contains smoking score calculated using Zhang method-
#' ## A html file named "results_MS.html" is generated in the Results folder.
#' ## A csv file saved to the Results folder.
#'
#' data(dummyBetaData)
#' result <- epismoker(dataset = dummyBetaData, ref.Elliott = Illig_data, method = "SSc")
#' generateReport(result, outputFileName = "results_SSc", method = "SSc")
#' ## result contains smoking score calculated using Elliott method-
#' ## A html file named "results_SSc.html" is generated in the Results folder.
#' ## A csv file saved to the Results folder.
#'
#' data(dummyBetaData)
#' samplesheet <- read.csv( system.file("extdata", "samplesheet_GSE42861.csv", package= "EpiSmokEr"), header=TRUE, sep=",")
#' data(CS_final_coefs, FS_final_coefs, NS_final_coefs)
#' result <- epismoker(dataset = dummyBetaData, samplesheet = samplesheet, ref.CS = CS_final_coefs, ref.FS = FS_final_coefs, ref.NS = NS_final_coefs, method = "SSt")
#' generateReport(result, outputFileName = "results_SSt", method = "SSt")
#' ## result contains predicted smoking probabilities and smoking status labels calculated using Multinomial LASSO method (SSt).
#' ## A html file named "results_SSt.html" is generated in the Results folder.
#' ## A csv file saved to the Results folder.
#'
#' data(dummyBetaData)
#' samplesheet <- read.csv( system.file("extdata", "samplesheet_GSE42861.csv", package= "EpiSmokEr"), header=TRUE, sep=",")
#' data(CS_final_coefs, FS_final_coefs, NS_final_coefs)
#' result <- epismoker(dataset = dummyBetaData, samplesheet = samplesheet, ref.Zhang = Zhangetal_cpgs, ref.Elliott =  Illig_data,
#'                     ref.CS = CS_final_coefs, ref.FS = FS_final_coefs, ref.NS = NS_final_coefs, method = "all")
#' generateReport(result, outputFileName = "Results_comprehensive", method = "all")
#' ## result contains smoking score calculated from all the three methods ( "SSc", "MS", "SSt").
#' ## A html file named "Results_comprehensive.html" is generated in the Results folder.
#' ## A csv file saved to the Results folder.
#'
#' @importFrom htmlTable htmlTable
#' @importFrom htmlTable txtRound
#' @importFrom rmarkdown render
#'
#' @export
#'
generateReport <- function(resObj,outputFileName, method) {
     if (!method %in% c("SSt","SSc", "MS", "all"))
      stop(sprintf("(%s) is not a valid method!", method))
    if (method == "SSt") {
      renderReport(system.file("reports", "SSt.rmd", package = "EpiSmokEr"), outputFileName)
    } else if (method == "SSc") {
      renderReport(system.file("reports", "SSc.rmd", package = "EpiSmokEr"), outputFileName)
    } else if (method == "MS") {
      renderReport(system.file("reports", "MS.rmd", package = "EpiSmokEr"), outputFileName)
    } else if (method == "all") renderReport(system.file("reports", "all.rmd", package = "EpiSmokEr"), outputFileName)
  }


renderReport <- function(rmdPath, outputFileName) {
  rmdPath <- normalizePath(rmdPath)
  mainDir <- dirname(getwd())
  outputDir <- file.path(mainDir, "Results")
  suppressWarnings(ifelse(!dir.exists(outputDir), dir.create(outputDir,showWarnings = TRUE), FALSE))
  outputDir <- normalizePath(outputDir)
  message("====================================================================================================================")
  message(sprintf("Saving results as html and csv files to %s", outputDir))
  message("====================================================================================================================")

  htmlFileName <- file.path(outputDir, basename(paste(outputFileName,"html",sep=".")) )
  render(input= rmdPath, output_format = "html_document", output_dir = outputDir, output_file =htmlFileName)
}
