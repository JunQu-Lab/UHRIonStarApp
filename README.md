# UHRIonStarApp
Ultra High Resolution (UHR)-IonStar is an MS1-based quantitative method for label-free proteomics experiments, devised to address issues related with quantitative precision, missing data, and false-positive discovery of protein changes in large-cohort analysis.
# UHR-IonStar, R shiny-based web application

![](https://github.com/JunQu-Lab/UHR-IonStar/blob/master/F1.large.jpg)

Ultra High Resolution (UHR)-IonStar is an MS1-based quantitative method for label-free proteomics experiments, devised to address issues related with quantitative precision, missing data, and false-positive discovery of protein changes in large-cohort analysis. 

UHR-IonStar app is a R shiny-based interactive application designed for processing, visualization and analysis of quantitative proteomics data generated by UHR-IonStar.

## UHR-IonStar Installation (Less than 20 minutes)
R (Version 4.0.5 or above) is required for Windows 10 or MacOS.\
Due to serveral R packgages from R Bioconductor and a package from GitHub cannot be automatically installed with the same time of UHR-IonStar R package
installation, please install the following packages in advance:
```
install.packages(c('devtools','shiny','shinydashboard'))
library('devtools')
devtools::install_github("vqv/ggbiplot")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()

BiocManager::install(c("MSGFplus","limma","TCGAbiolinks","clusterProfiler", "Biobase", "DO.db",
                       "AnnotationHub","mzR","MSnbase","xcms","CAMERA"))
install.packages("patRoon", repos = "https://rickhelmus.github.io/patRoonDeps/", type = "binary")

```
Then, you can install UHR-IonStar with the following line:
```
devtools::install_github("JunQu-Lab/UHRIonStarApp")
```
## Troubleshooting
- R will occasionally ask users whether they want to update any dependent packages that have a new version. We recommend that you update them all.
- Rtools is also required to install R packages.
- During the installation of dependent packages, sometimes R would inquire, "Do you want to install from sources the packages which need compliation?" If users select "Yes" and receive an error, selecting "No" is equally acceptable.

## Run UHR-IonStar
If no error pops up, the UHR-IonStar web app could be started with the following codes:
```
library(UHR.IonStar)
UHR.IonStar::UHRIonStarShiny()
```

## Manual
User can download the manual either at the mainpage of UHR-IonStar application or at the github directory [UHRIonStarApp/inst/shiny/UHRIonStar](https://github.com/JunQu-Lab/UHRIonStarApp/tree/master/inst/shiny/UHRIonStar).

## Related Articles
[Shen, Xiaomeng, et al. "IonStar enables high-precision, low-missing-data proteomics quantification in large biological cohorts." Proceedings of the National Academy of Sciences 115.21 (2018): E4767-E4776.](https://www.pnas.org/content/115/21/E4767.short)

[Wang, Xue, et al. "Ultra-High-Resolution IonStar Strategy Enhancing Accuracy and Precision of MS1-Based Proteomics and an Extensive Comparison with State-of-the-Art SWATH-MS in Large-Cohort Quantification." Analytical chemistry 93.11 (2021): 4884-4893.](https://pubs.acs.org/doi/abs/10.1021/acs.analchem.0c05002?casa_token=12l8WRigfZ0AAAAA:0qwzMnfjpE2stVCpMYKICmvqwofN15Q6ItzDZ7ATFY3m3aFI6oSzB1z20CJGzzwASyaegR5POgS8xA)
![image](https://user-images.githubusercontent.com/59838185/171267509-2e6adec5-b0c0-471b-8e5b-278ac8b13ad5.png)
