# EpiSmokEr: Epigenetic Smoking status Estimator
An R Package for the estimation of smoking status based on whole blood 450K methylation profiles
## Features of **EpiSmokEr**
- Can perform smoking score calculation starting from idat files.
- Normalization of datasets using channel and colour specific quantile normalization.
- Prediction of smoking status based on the whole blood methylation data.
- Results provided in html and csv format.

## Dependencies
**EpiSmokEr** depends on the following packages:
```{r eval=FALSE}
library(IlluminaHumanMethylation450kmanifest)
library(minfi)
library(htmlTable)
library(rmarkdown)
```
## Installation
To install EpiSmokEr, start R and then type the following commands:
```{r eval=FALSE}
source("http://bioconductor.org/biocLite.R")
install.packages("devtools") # if you don't have the package, run install.packages("devtools")
library(devtools)
install_github("sailalithabollepalli/EpiSmokEr")
library(EpiSmokEr)
```
## Vignette
Please refer to the [vignette](http://htmlpreview.github.io/?https://github.com/sailalithabollepalli/EpiSmokEr/blob/master/vignettes/epismoker.html) for step by step demonstration of functions in the **EpiSmokEr** package. 

## Citation
If you use **EpiSmokEr**, please cite:
Bollepalli S, Korhonen T, Kaprio J, Anders S, Ollikainen M. EpiSmokEr: a robust classifier to determine smoking status from DNA methylation data. Epigenomics. 2019 Oct;11(13):1469-1486. doi: 10.2217/epi-2019-0206. Epub 2019 Aug 30. PMID: 31466478.



