#' dummyBetaData
#'
#' dummyBetaData is a subset of 6 control subjects from the epigenome-wide association study performed in peripheral-blood DNA in 354
#' anti-citrullinated protein antibody-associated rheumatoid arthritis cases and 337 controls. To minimise the package size and running time
#' we have used a subset of 1000 CpG probes from the methylation matrix (which contained all the CpGs required by EpiSmokEr).
#'
#' @docType data
#'
#' @usage data(dummyBetaData)
#'
#' @keywords datasets
#'
#' @format A dataframe with 1000 rows and 6 columns, containing methylation values in beta scale ranging between 0 and 1.
#'
#' @references Liu Y, Aryee MJ, Padyukov L, Fallin MD et al. Epigenome-wide association data implicate DNA methylation as an intermediary of genetic risk in rheumatoid arthritis.
#' Nat Biotechnol 2013 Feb;31(2):142-7. PMID: 23334450
#' (\href{http://www.ncbi.nlm.nih.gov/pubmed/23334450}{PubMed})
#'
#' @source \href{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse42861}{GEO Repository}
#'
#' @examples
#' data(dummyBetaData)
#' data(Illig_data)
#' result <- epismoker(dataset = dummyBetaData, ref.Elliott = Illig_data, method = "SSc")
#' ## result contains smoking score calculated based on Elliott et al.
"dummyBetaData"
