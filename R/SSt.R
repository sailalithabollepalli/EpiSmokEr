#' @title
#' Estimates smoking probabilities and smoking status (SSt) using MLM method.
#'
#' @description
#' Estimates smoking probabilities and smoking status (SSt) using MLM method.
#'
#' @param dataset
#' A data matrix of normalised methylation values in beta scale, with rows labelling the CpG probe ids and columns labelling the sample names.
#' Missing methylation values are imputed as 0.5.
#' eg: dummyBetaData
#'
#' @param samplesheet
#' A dataframe with samplenames as rownames. Must contain gender information as  "sex" column. sex column should be marked as 1 and 2
#' representing male and female respectively. Missing "sex" information is imputed as 0.
#' eg: dummySamplesheet
#'
#' @param ref.CS
#' This is used in "SSt" and "all" methods. It is a dataframe with 121 CpGs selected by Multinomial LASSO approach, with log odd values for "CURRENT SMOKER" class.
#' In addition intercept term and sex coefficient terms are also provided by Multinomial LASSO approach.
#'
#' @param ref.FS
#' This is used in "SSt" and "all" methods. It is a dataframe with 121 CpGs selected by Multinomial LASSO approach, with log odd values for "FORMER SMOKER" class.
#' In addition intercept term and sex coefficient terms are also provided by Multinomial LASSO approach.
#'
#' @param ref.NS
#' This is used in "SSt" and "all" methods. It is a dataframe with 121 CpGs selected by Multinomial LASSO approach, with log odd values for "NEVER SMOKER" class.
#' In addition intercept term and sex coefficient terms are also provided by Multinomial LASSO approach.
#'
#' @return
#' Returns a result object with predicted smoking probabilities and smoking status labels generated from MLM method.
#'
#' data(dummyBetaData)
#' samplesheet <- read.csv( system.file("extdata", "samplesheet_GSE42861.csv", package= "EpiSmokEr"), header=TRUE, sep=",")
#' result <- epismoker(dataset = dummyBetaData, samplesheet = samplesheet, ref.CS = CS_final_coefs, ref.FS = FS_final_coefs, ref.NS = NS_final_coefs, method = "SSt")
#' ## result contains predicted smoking probabilities and smoking status labels generated from MLM method.
#' @export
#'
SSt <- function(dataset, samplesheet, ref.CS = CS_final_coefs, ref.FS = FS_final_coefs, ref.NS = NS_final_coefs){
  message("============================================")
  message("<<<<< Predicting Smoking Statuses >>>>>")
  message("============================================")
  logOdds_CS <- getLogOdds(dataset, samplesheet, ref.CS)
  logOdds_FS<-  getLogOdds(dataset, samplesheet, ref.FS)
  logOdds_NS <- getLogOdds(dataset, samplesheet, ref.NS)

  odds_CS <- exp(logOdds_CS)
  odds_FS <- exp(logOdds_FS)
  odds_NS <- exp(logOdds_NS)

  sumOdds <- rowSums(cbind(odds_CS, odds_FS, odds_NS ))
  odds2probs <- function(odds, sumOdds) {odds/sumOdds}
  probs_CS <- odds2probs(odds_CS,sumOdds)
  probs_FS <-odds2probs(odds_FS,sumOdds)
  probs_NS <-odds2probs(odds_NS,sumOdds)
  coeffs <- setNames(data.frame( logOdds_CS, logOdds_FS, logOdds_NS), c("logOdds_CS", "logOdds_FS", "logOdds_NS"))
  coeffs_probs <- setNames(data.frame(probs_CS, probs_FS, probs_NS), c("probs_CS", "probs_FS", "probs_NS"))

  res_SSt <- setNames(data.frame(cbind( rownames(coeffs),coeffs, coeffs_probs,getSmokingStatus2(coeffs_probs[,c(1:3)]))),
                      c("SampleName","logOdds_CS","logOdds_FS","logOdds_NS","probs_CS","probs_FS","probs_NS","PredictedSmokingStatus"))
  message("============================================")
  message("<<<<< Prediction of Smoking Statuses Completed >>>>>")
  message("============================================")
  return(data.frame(res_SSt))
}

getLogOdds<- function( dataset, samplesheet, coefs)
{
  if(is.null(dataset))
    stop("Could not find dataset required for smoking status estimation.")
  dataset_SSt <- t(dataset)
  if(is.null(samplesheet))
    stop("Could not find samplesheet required for smoking status estimation.")
  if(is.null(coefs))
    stop("Could not find coefficients required for smoking status estimation.")
  names(coefs)[1] <- "intercept"
  intercept <- coefs[1]
  stopifnot(names(coefs)[2]=="sexM")
  if(!"sex" %in% colnames(samplesheet) & !"Sex" %in% colnames(samplesheet))
    stop("Need column 'sex' amongst the column names of 'samplesheet'.")
  if("Sex" %in% colnames(samplesheet)){samplesheet$sex <- samplesheet$Sex}
  dataset_SSt  <- as.matrix( dataset_SSt )
  if(nrow(dataset_SSt )!= nrow(samplesheet))stop("Number of samples differ between samplesheet and methylation dataset")
  #if(!all(rownames(dataset) %in% rownames(samplesheet)))stop("Rownames of samplesheet must be equivalent to column names of methylation dataset")
  if(length(intersect(rownames(dataset_SSt ),rownames(samplesheet)))< length(rownames(samplesheet)))stop("Rownames of samplesheet must be equivalent to column names of methylation dataset")
  dataset_SSt  <- dataset_SSt[rownames(samplesheet),]
  stopifnot(rownames(dataset_SSt )== rownames(samplesheet))
  sex <- setNames( data.frame (samplesheet$sex), "sex")
  if(any(is.na(sex$sex))) {message("Replacing samples with NAs in 'sex' column with a zero.")}
  #sex  <- na.omit(sex)
  sex[is.na(sex)] <- 0
  sex$sex <- ifelse(sex$sex==1,coefs[2], 0)
  rownames(sex) <- rownames(samplesheet)
  #if(rownames(dataset) != rownames(samplesheet) ) # & rownames(dataset)!= samplesheet$Sample_Name
  #stop("Samplesheet should have `Sample_Name` column or rownames of samplesheet must be same as the sample names of the dataset.")
  suppressWarnings(if(any(is.na(dataset_SSt ))) {message("Replacing samples with NAs in methylation data with 0.5.")})
  #age<- na.omit(age)
  dataset_SSt [is.na(dataset_SSt )] <- 0.5
  weights <- coefs[-c(1,2)]
  y<- length(weights)
  # subset weights with only available CpGs
  weights <- subset(weights, names(weights) %in% colnames(dataset_SSt ))
  message(sprintf("Dataset has %s of %s CpGs required for smoking status estimation.",length(weights),y))
  dataset_SSt  <- dataset_SSt [,names(weights)]
  # check if the order colnames of dataset equal to names(weights)
  stopifnot(colnames(dataset_SSt ) == names(weights))
  stopifnot(rownames(dataset_SSt )==rownames(sex))
  scoreFunction <- function(x) sum(x*weights)+intercept
  logOdds <- sex+apply(dataset_SSt,1,scoreFunction)
  rm(dataset_SSt)
  return(logOdds)}

getSmokingStatus <- function(data)
{ colmax <- data.frame(colnames(data)[apply(data,1, which.max)])
PredictedSmokingStatus <- setNames(data.frame(X=ifelse(colmax == "logOdds_CS", "Current Smoker", ifelse( colmax == "logOdds_NS","Never Smoker", "Former Smoker"))),"PredictedSmokingStatus")
return(PredictedSmokingStatus)}

getSmokingStatus2 <- function(data)
{ colmax <- data.frame(colnames(data)[apply(data,1, which.max)])
PredictedSmokingStatus2 <- setNames(data.frame(X=ifelse(colmax == "probs_CS", "Current Smoker", ifelse( colmax == "probs_NS","Never Smoker", "Former Smoker"))),"PredictedSmokingStatus2")
return(PredictedSmokingStatus2)}

