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
#' Users can choose one of these four methods ( "EM", "ZM", "MLM", "all") in the epismoker function.
#'
#' @return html file
#' Returns a html file with the results from the epismoker function, saved to the Results folder.
#'
#' @return csv file
#' Returns a csv file with the results from the epismoker function, saved to the Results folder.
#'
#' @examples
#' data(dummyBetaData)
#' result <- epismoker(dataset = dummyBetaData, ref.Zhang = Zhangetal_cpgs, method = "ZM")
#' generateReport(result, outputFileName = "results_ZM", method = "ZM")
#' ## result contains smoking score calculated using Zhang method-
#' ## A html file named "results_ZM.html" is generated in the Results folder.
#' ## A csv file saved to the Results folder.
#'
#' data(dummyBetaData)
#' result <- epismoker(dataset = dummyBetaData, ref.Elliott = Illig_data, method = "EM")
#' generateReport(result, outputFileName = "results_EM", method = "EM")
#' ## result contains smoking score calculated using Elliott method-
#' ## A html file named "results_EM.html" is generated in the Results folder.
#' ## A csv file saved to the Results folder.
#'
#' data(dummyBetaData)
#' samplesheet <- read.csv( system.file("extdata", "samplesheet_GSE42861.csv", package= "EpiSmokEr"), header=TRUE, sep=",")
#' data(CS_final_coefs, FS_final_coefs, NS_final_coefs)
#' result <- epismoker(dataset = dummyBetaData, samplesheet = samplesheet, ref.CS = CS_final_coefs, ref.FS = FS_final_coefs, ref.NS = NS_final_coefs, method = "MLM")
#' generateReport(result, outputFileName = "results_MLM", method = "MLM")
#' ## result contains predicted smoking probabilities and smoking status labels calculated using Multinomial LASSO method (MLM).
#' ## A html file named "results_MLM.html" is generated in the Results folder.
#' ## A csv file saved to the Results folder.
#'
#' data(dummyBetaData)
#' samplesheet <- read.csv( system.file("extdata", "samplesheet_GSE42861.csv", package= "EpiSmokEr"), header=TRUE, sep=",")
#' data(CS_final_coefs, FS_final_coefs, NS_final_coefs)
#' result <- epismoker(dataset = dummyBetaData, samplesheet = samplesheet, ref.Zhang = Zhangetal_cpgs, ref.Elliott =  Illig_data,
#'                     ref.CS = CS_final_coefs, ref.FS = FS_final_coefs, ref.NS = NS_final_coefs, method = "all")
#' generateReport(result, outputFileName = "results_all", method = "all")
#' ## result contains smoking score calculated from all the three methods ( "EM", "ZM", "MLM").
#' ## A html file named "results_all.html" is generated in the Results folder.
#' ## A csv file saved to the Results folder.
#'
#' @importFrom htmlTable htmlTable
#' @importFrom htmlTable txtRound
#' @importFrom rmarkdown render
#'
#' @export
#'
generateReport <- function(resObj,outputFileName, method) {
     if (!method %in% c("MLM","EM", "ZM", "all"))
      stop(sprintf("(%s) is not a valid method!", method))
    if (method == "MLM") {
      renderReport(system.file("reports", "MLM.rmd", package = "EpiSmokEr"), outputFileName)
    } else if (method == "EM") {
      renderReport(system.file("reports", "EM.rmd", package = "EpiSmokEr"), outputFileName)
    } else if (method == "ZM") {
      renderReport(system.file("reports", "ZM.rmd", package = "EpiSmokEr"), outputFileName)
    } else if (method == "all") renderReport(system.file("reports", "all.rmd", package = "EpiSmokEr"), outputFileName)
  }


renderReport <- function(rmdPath, outputFileName) {
  rmdPath <- normalizePath(rmdPath)
  mainDir <- dirname(getwd())
  outputDir <- file.path(mainDir, "Results")
  suppressWarnings(ifelse(!dir.exists(outputDir), dir.create(outputDir), FALSE))
  outputDir <- normalizePath(outputDir)
  message("====================================================================================================================")
  message(sprintf("Saving results as html and csv files to %s", outputDir))
  message("====================================================================================================================")

  htmlFileName <- file.path(outputDir, basename(paste(outputFileName,"html",sep=".")) )
  render(input= rmdPath, output_format = "html_document", output_dir = outputDir, output_file =htmlFileName)
}
