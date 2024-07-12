# Afrotheria

In this document we briefly describe how to do the analyses for the manuscript *Convergent evolution in Afrotheria and non-afrotherians demonstrates high evolvability of the mammalian inner ear* (Grunstra et al.).  

R scripts are provided in the GitHub repository.  

The raw data is provided as separate files in the OSF repository <https://osf.io/9mtwh/>. 

# 1. System requirement

For the analyses, we use the software R, which is a free software environment for statistical computing and graphics. It compiles and runs on a wide variety of UNIX platforms, Windows and MacOS. The scripts were run with R version 4.3.1. 

It is recommended (but not mandatory) to use R with RStudio, a free integrated development environment (IDE). RStudio runs on a wide variety of UNIX platforms, Windows and MacOS. It requires R version 3.6.0 or more recent. For the analyses, we used RStudio Desktop version 2023.12.1. 

Some R scripts require the installation of supplementary packages. For the analyses, we used the following packages:  

* `ape` version 5.7-1  
* `geomorph` version 4.0.6  
* `Morpho` version 2.11  
* `phytools` version 2.0-3  
* `ggplot2` version 3.4.4  
* `ggrepel` version 0.9.4  

# 2. Installation guide

## Download and install R and RStudio

1. To download and install R, follow the instruction on the R Project webpage: <https://www.r-project.org/>. 

2. **[Optional, but recommended]** After installing R (> 3.6.0), download and install the desktop version of RStudio: <https://posit.co/download/rstudio-desktop/>. 

The downloading and installation of R and RStudio will take a few minutes, depending on the quality of the internet connection. 

## Download and install the required R packages

After the installation of R and then RStudio, installation of supplementary packages is required for the analyses.  

To install the latest version of a package, open RStudio (or R if RStudio is not installed), then write the following prompt in the console (here for the package `ggplot2`):
```{r install}
install.packages("ggplot2")  # to install ggplot2  
```

Run similar operations for all the packages listed above by replacing `ggplot2` with the name of the package.  

The downloading and installation of each package will take a few minutes, depending on the quality of the internet connection. Some packages depend on other packages and hence their installation requires the installation of the packages on which the depend.  

If asked, we recommend using stabilized versions of the package, and not binaries. 

# 3. Demonstration and instruction for use

## Preparation

1. Download all R scripts and raw data files (see **section 4. Description of the scripts and data**)

2. Put all R scripts and raw data files in a folder <MyFolder>. 

## Creation of an R project with RStudio [OPTIONAL]

If you are using RStudio, it is recommended to create an R project:  

1. Open RStudio  
2. On the upper right corner of the window: Project > New Project > Existing Directory  
3. Browse to select the directory <MyFolder> with the scripts and raw data files  
4. Click on "Create Project"  

After creating the project, the working directory is set to the directory <MyFolder>.  

The name of the project appears on the upper right corner of the window. There you can change projects, or close all projects (and hence go back to the default directory).  

## How to run scripts

1. Set the working directory so that it corresponds to the folder where your scripts and data are (use `/` for Windows systems and `\` for Mac OS). If you have an R project /see above), the working directory is already correct.
```{r wd}
setwd("C:/.../MyFolder")  # change the path to fit the relevant folder
```

2. Open the first file `00_Definitions.R` and run all lines in the console. This script has to be run before all the others, because it is where the landmarks, curves, contextual data and related visualization function and parameters are described and created. In RStudio, data, values and functions are visible in the "Environment tab.Alternatively, it is possible to run the whoe script with the command:   
```{r def}
source("00_Definitions.R")
```

3. Open the script of interest. All R scripts other than `00_Definitions.R` can be run independently (recommended, because there may be some variable names conflicts).  

4. Load the packages relevant for the script (listed in the description of the script). This enable an access to supplementary R functions and data sets not provided with eh base R environment.
```{r geomorph}
library(geomorph)  # here an example to load the package "geomorph"
```

5. Now you can run the commands in the console.  

# 4. Description of the scripts and data sets

## Description of the raw data

The raw data is provided as separate files in the OSF repository <https://osf.io/9mtwh/>  

Data files are:

* `Grunstra_et_al_Afrotheria_contextual.csv`: contextual data => used for 2B-PLS and scatter plots of PCA and PLS
* `Grunstra_et_al_Afrotheria_landmarks_raw.csv`: 3D landmark coordinates (13 anatomical landmarks, 111 semilandmarks) on the bony labyrinth for the 40 taxa => used for the GPA
* `Grunstra_et_al_Afrotheria_landmarks_slid_not_aligned.csv`: 3D landmark coordinates (13 anatomical landmarks, 111 semilandmarks) on the bony labyrinth for the 40 taxa, after semilandmark sliding (in the software Mathematica, but similar algorithm provided in R) => used for the computation of Procrustes distances, the PCA and the PLS analyses
* `Grunstra_et_al_Afrotheria_phyl_tree.txt`: phylogenetic tree of the 40 taxa in NEXUS format => used for the phylogenetic PLS and the computation of the phylogenetic signal

## Description of the scripts

We provide 6 scripts that correspond to the analyses performed using R, as described in the main text. Five files are provided: 

* `00_Definitions.R`: landmark definition, preparation of contextual data, creation of the function `plotLaby` for visualization *[running time: 2s]*  
* `01_GPA.R`: semilandmark sliding, GPA, outliers, visualization *[running time: 12s]*  
* `02_ProcD.R`: pairwise Procrustes distances, cross-validations *[running time: 5s]*  
* `03_PCA.R`: PCA of the Procrustes shape coordinates *[running time: 8s]*  
* `04_PLS.R`: 2B-PLS, cross-validations, phylogenetic signal in the PLS scores, phylogenetic PLS *[running time: 18s]*  
* `05_WI.R`: Wheatsheaf index (Arbuckle et al. 2014 *Meth. Ecol. Evol.*) *[running time: 5s]*  

