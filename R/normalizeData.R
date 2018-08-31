#' @title
#' Color and channel specific quantile normalisation
#'
#' @description
#' normalizeData function performs color and channel specific quantile normalisation.
#'
#' @param RGset and normMethod
#' Requires a RGset object and the normalisation method
#'
#' @return dataset
#' A dataframe with normalised beta values (M/U+M+100) applying one of the chosen normalisation methods
#'
#' @examples
#' dataset <- normalizeData(RGset, normMethod="QN")
#' ## dataset is a data matrix with normalised methylation values in beta scale ranging between 0 and 1.
#'
#' @import minfi
#' @import IlluminaHumanMethylation450kmanifest
#'
#' @export
#'
normalizeData <- function (RGset=RGset, normMethod=c("QN", "ILM", "SQN", "ALL")){
  message("Loading RGset object.")
  RGset <- RGset
  message(sprintf("Dataset has %s samples.",dim(RGset)[2]))
  if(normMethod=="QN")
  {
    TypeII.Name <- minfi::getProbeInfo(RGset, type = "II")$Name
    TypeII.Green <- minfi::getGreen(RGset)[minfi::getProbeInfo(RGset, type = "II")$AddressA,]
    TypeII.Red <- minfi::getRed(RGset)[minfi::getProbeInfo(RGset, type = "II")$AddressA,]
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
    dataset_QN <- QN_beta
    return(dataset_QN)

  }else if(normMethod == "ILM") # Illumina Normalisation
  {
    message("Performing Illumina Normalisation.")
    MSet.illumina <- minfi::preprocessIllumina(RGset, bg.correct = FALSE, normalize = "controls")
    dataset_ILM <- minfi::getBeta(MSet.illumina)
    return(dataset_ILM)

  }else if(normMethod == "SQN") # Stratified Quantile Normalisation (Touleimat & Tost)
  {
    message("Performing Subset Quantile Normalisation.")
    gset.quantile <- minfi::preprocessQuantile(RGset, fixOutliers = TRUE, removeBadSamples = TRUE, badSampleCutoff = 10.5,
                                               quantileNormalize = TRUE, stratified = TRUE,mergeManifest = FALSE, sex = NULL)
    dataset_SQN<- minfi::getBeta(gset.quantile)
    return(dataset_SQN)
  }else if(normMethod=="ALL")
  {
    # QN
    TypeII.Name <- minfi::getProbeInfo(RGset, type = "II")$Name
    TypeII.Green <- minfi::getGreen(RGset)[minfi::getProbeInfo(RGset, type = "II")$AddressA,]
    TypeII.Red <- minfi::getRed(RGset)[minfi::getProbeInfo(RGset, type = "II")$AddressA,]
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
    dataset_QN <- QN_beta
    ###################################################################################################
    # SQN
    message("Performing Subset Quantile Normalisation.")
    gset.quantile <- minfi::preprocessQuantile(RGset,quantileNormalize = TRUE, stratified = TRUE,mergeManifest = FALSE, sex = NULL, verbose=TRUE)
    dataset_SQN<- minfi::getBeta(gset.quantile)
    ###################################################################################################
    # ILM
    message("Performing Illumina Normalisation.")
    MSet.illumina <- minfi::preprocessIllumina(RGset, bg.correct = FALSE, normalize = "controls")
    dataset_ILM <- minfi::getBeta(MSet.illumina)
    return(list(dataset_QN=dataset_QN, dataset_ILM=dataset_ILM, dataset_SQN=dataset_SQN))

  }

}

