# Afrotheria

In this document we briefly describe the R scripts provided as a complement of the manuscript *Convergent evolution in Afrotheria and non-afrotherians demonstrates high evolvability of the mammalian inner ear* (Grunstra et al.). 

We provide 5 scripts that correspond to the analyses performed using R, as described in the main text. Five files are provided: 

* `00_Definitions.R`: landmark definition, preparation of contextual data, creation of teh function `plotLaby`
* `01_GPA.R`: semilandmark sliding, GPA, outliers, visualization
* `02_ProcD.R`: pairwise Procrustes distances, cross-validations
* `03_PCA.R`: PCA of the Procrustes shape coordinates
* `04_PLS.R`: 2B-PLS, cross-validations, phylogenetic signal in the PLS scores, phylogenetic PLS

The first file `00_Definitions.R` has to be run before all the others, because it is where the landmarks, curves, contextual data and related visualization function and parameters are described. 

The other R scripts can be run independently (recommended, because there may be some variable names conflicts). 
