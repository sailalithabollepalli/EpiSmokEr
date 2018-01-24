# EpiSmokEr: Epigenetic Smoking status Estimator

# Introduction
***
This vignette gives a overview on how to use the **EpiSmokEr** R package. **EpiSmokEr** provides four options to calculate smoking score and predict smoking status from whole blood Infinium HumanMethylation450 data.

## Four options in EpiSmokEr
- **Multinomial LASSO method (MLM)** provides smoking probabilities for each individual and **predicts smoking status**. 121 CpGs identified by multinomial LASSO approach and
- **Elliott method (EM)** provides smoking score calculated using 187 CpGs [@elliott].

- **Zhang method (ZM)** provides methylation score calculated using 4 CpGs [@zhang].
- **all** a comprehensive approach which provides results from all the three methods described above.

## Features of **EpiSmokEr**
- Can perform smoking score calculation starting from idat files.
- Normalization of datasets using channel and colour specific quantile normalization.
- Smoking score calculation using three methods.
- Prediction of smoking status based on the whole blood methylation data.
- Results provided in html and csv format.

Examples demonstrating th features above can be found in the [demo/](tests) directory.

# Installation
***
To install EpiSmokEr, start R and then type the following commands:
```{r eval=FALSE}
source("http://bioconductor.org/biocLite.R")
install.packages("devtools") # if you don't have the package, run install.packages("devtools")
library(devtools)
install_github("sailalithabollepalli/EpiSmokEr")
```

# Input Data
***
  To generate the smoking scores from the whole blood methylaion, EpiSmokEr simply requires methylation data as an input.
Input data could be either raw methylation data in the form of idat files or a normalised methylation matrix in the beta scale ranging in between o and 1.
For MLM method, a sample sheet with **sex** status is needed to complement the methylation data.

To demonstrate the working of the EpiSmokEr we chose a peripheral blood leukocytes dataset from GEO [GSE42861][@geo](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse42861). We use only 6 samples from this dataset to minimise the package size and running time.

#### Following data objects are provided along with the package and are used in the examples below.
- Raw data: idat files of 6 subjects from GSE42861.
- Methylation matrix: A subset of 1000 CpG probes from the normalised methylation matrix (which contained all the CpGs required by EpiSmokEr).
- A sample sheet with gender information

First we will show how to use idat files as input data
We will now show 2 ways of providing input data:
  
## 1. From idat files

Here input data are in the form of idat files. `minfi`[@minfi] package is used to read the idat files.
color and channel specific normalisation was performed on the reference dataset and quantiles are saved from these are used internally to adjust the distribution using a suite of customised functions. beta values are calculated.

### Please note:
- Sample sheet must be included in the same folder as the idat files.
- Rownames of samplesheet must match the names of the samples i.e idat files, to facilitate matching of phenotype data with the corresponding methylation data.
- For MLM method, sample sheet should include a column with gender information marked as **sex**, in the format of 1 and 2 representing men and women respectively.

```{r}
dataset <- normalizeData(idatPath = system.file("extdata", package = "EpiSmokEr"))
samplesheet <- read.csv(system.file("extdata", "samplesheet_GSE42861.csv", package="EpiSmokEr"), header=TRUE, sep=",")
samplesheet
result <- epismoker(dataset = dataset, samplesheet = samplesheet, method = "MLM")
generateReport(resObj=result,outputFile = "Result", method = "EM")
```

## 2. From methylation matrix
Normalised methylation matrix can also be used as an input to EpiSmokEr.
```{r, include= FALSE}
data("dummyBetaData")
```
Methylation matrix looks likes this. They are in beta scale ragining betweeen 0 and 1.
```{r}
head(dummyBetaData)
```

```{r, echo=FALSE, results='asis'}
samplesheet <- read.csv(system.file("extdata", "samplesheet_GSE42861.csv", package="EpiSmokEr"), header=TRUE, sep=",")
knitr::kable(head(samplesheet, 5))
```

# Smoking Score calculation and predicton of smoking status
***
Once we have the methylation data, we can then proceed with smoking status estimation. EpiSmokEr provides four options

### 1. Multinomial LASSO method (MLM)
```{r}
result <- epismoker(dataset=dummyBetaData, samplesheet = samplesheet, method = "MLM")
result
```

### 2. Elliott method (EM)
```{r}
result <- epismoker(dataset = dataset, method = "EM")
result
```

### 3. Zhang method (ZM)
```{r}
result <- epismoker(dataset = dataset,  method = "ZM")
result
```

### 4. Compreshensive approach
```{r}
result <- epismoker(dataset = dataset, samplesheet = samplesheet, method = "all")
result
```

# Generate Results
****

# References 
***
