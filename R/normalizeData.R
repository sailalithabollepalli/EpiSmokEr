#' @title
#' Color and channel specific quantile normalisation
#'
#' @description
#' normalizeData function performs color and channel specific quantile normalisation.
#'
#' @param idatPath
#' Requires a path to the folder with the idat files and sample sheet included.
#'
#' @return dataset
#' A dataframe with normalised beta values (M/U+M+100)
#'
#' @examples
#' dataset <- normalizeData(idatPath = system.file("extdata", package = "EpiSmokEr"))
#' ## dataset is a data matrix with normalised methylation values in beta scale ranging between 0 and 1.
#'
#' @import minfi
#' @import IlluminaHumanMethylation450kmanifest
#'
#' @export
#'
normalizeData <- function (idatPath){
  message("Loading reference data quantiles.")
  baseDir <- file.path(idatPath)
  targets <- minfi::read.450k.sheet(baseDir)
  RGset <- minfi::read.450k.exp(base = baseDir)
  message(sprintf("Dataset has %s samples.",dim(RGset)[2]))
  TypeII.Name <- minfi::getProbeInfo(RGset, type = "II")$Name
  TypeII.Green <- minfi::getGreen(RGset)[minfi::getProbeInfo(RGset, type = "II")$Address,]
  TypeII.Red <- minfi::getRed(RGset)[minfi::getProbeInfo(RGset, type = "II")$Address,]
  rownames(TypeII.Red) <- TypeII.Name
  colnames(TypeII.Red) <- sampleNames(RGset)
  rownames(TypeII.Green) <- TypeII.Name
  colnames(TypeII.Green) <- sampleNames(RGset)
  TypeI.Green.Name <- minfi::getProbeInfo(RGset, type = "I-Green")$Name
  TypeI.Green.M <- minfi::getGreen(RGset)[minfi::getProbeInfo(RGset, type = "I-Green")$AddressB,]
  rownames(TypeI.Green.M) <- TypeI.Green.Name
  colnames(TypeI.Green.M) <- sampleNames(RGset)
  TypeI.Green.U <- minfi::getGreen(RGset)[minfi::getProbeInfo(RGset, type = "I-Green")$AddressA,]
  rownames(TypeI.Green.U) <- TypeI.Green.Name
  colnames(TypeI.Green.U) <- sampleNames(RGset)
  TypeI.Red.Name <- minfi::getProbeInfo(RGset, type = "I-Red")$Name
  TypeI.Red.M <- minfi::getRed(RGset)[minfi::getProbeInfo(RGset, type = "I-Red")$AddressB,]
  rownames(TypeI.Red.M) <- TypeI.Red.Name
  colnames(TypeI.Red.M) <- sampleNames(RGset)
  TypeI.Red.U <- minfi::getRed(RGset)[minfi::getProbeInfo(RGset, type = "I-Red")$AddressA,]
  rownames(TypeI.Red.U) <- TypeI.Red.Name
  colnames(TypeI.Red.U) <- sampleNames(RGset)
  # Subsetting rownames
  TypeII.Green <- TypeII.Green[rnms_QN_TypeII,]
  TypeII.Red <- TypeII.Red[rnms_QN_TypeII,]
  TypeI.Green.M <- TypeI.Green.M[rnms_QN_TypeI.Grn,]
  TypeI.Green.U <- TypeI.Green.U[rnms_QN_TypeI.Grn,]
  TypeI.Red.M <- TypeI.Red.M[rnms_QN_TypeI.Red,]
  TypeI.Red.U <- TypeI.Red.U[rnms_QN_TypeI.Red,]
  #Quantile normalise TypeII probes using quantiles from Reference data
  for(i in 1:ncol(TypeII.Green)){
    TypeII.Green[,i] <- quantiles_TypeII.Green[rank(TypeII.Green[,i],ties.method="random")]
  }

  for(i in 1:ncol(TypeII.Red)){
    TypeII.Red[,i] <- quantiles_TypeII.Red[rank(TypeII.Red[,i],ties.method="random")]
  }

  QN_beta_TypeII <- TypeII.Green/(TypeII.Green+TypeII.Red+100)
  rownames(QN_beta_TypeII) <- rownames(TypeII.Green)
  # Quantile normalise TypeI probes using quantiles from Reference data

  for(i in 1:ncol(TypeI.Green.M)){
    TypeI.Green.M[,i] <- quantiles_TypeI.Green.M[rank(TypeI.Green.M[,i],ties.method="random")]
  }

  for(i in 1:ncol(TypeI.Green.U)){
    TypeI.Green.U[,i] <- quantiles_TypeI.Green.U[rank(TypeI.Green.U[,i],ties.method="random")]
  }
  QN_beta_GreenTypeI <- TypeI.Green.M/(TypeI.Green.M+TypeI.Green.U+100)
  rownames(QN_beta_GreenTypeI)<- rownames(TypeI.Green.M)

  for(i in 1:ncol(TypeI.Red.M)){
    TypeI.Red.M[,i] <- quantiles_TypeI.Red.M[rank(TypeI.Red.M[,i],ties.method="random")]
  }
  for(i in 1:ncol(TypeI.Red.U)){
    TypeI.Red.U[,i] <- quantiles_TypeI.Red.U[rank(TypeI.Red.U[,i],ties.method="random")]
  }
  QN_beta_RedTypeI <- TypeI.Red.M/(TypeI.Red.M+TypeI.Red.U+100)
  rownames(QN_beta_RedTypeI) <- rownames(TypeI.Red.M)
  QN_beta <- rbind(QN_beta_TypeII, QN_beta_GreenTypeI, QN_beta_RedTypeI)
  dataset <- QN_beta
  return(dataset) }

