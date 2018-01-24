# EpiSmokEr: Epigenetic Smoking status Estimator
An R Package for the estimation of smoking status based on whole blood 450K methylation profiles
## Features of **EpiSmokEr**
- Can perform smoking score calculation starting from idat files.
- Normalization of datasets using channel and colour specific quantile normalization.
- Smoking score calculation using three methods.
- Prediction of smoking status based on the whole blood methylation data.
- Results provided in html and csv format.
## Installation
To install EpiSmokEr, start R and then type the following commands:
```{r eval=FALSE}
source("http://bioconductor.org/biocLite.R")
install.packages("devtools") # if you don't have the package, run install.packages("devtools")
library(devtools)
install_github("sailalithabollepalli/EpiSmokEr")
```
## Dependencies
**EpiSmokEr** depends on the following packages:
```{r eval=FALSE}
library(IlluminaHumanMethylation450kmanifest)
library(minfi)
library(htmlTable)
library(rmarkdown)
```
## Vignette
Please refer to the [vignettee]() for step by step demonstration of functions in **EpiSmokEr** package using example data. 




