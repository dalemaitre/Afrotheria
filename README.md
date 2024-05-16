# Afrotheria

In this document we briefly describe the R scripts provided as a complement of the manuscript *Convergent evolution in Afrotheria and non-afrotherians demonstrates high evolvability of the mammalian inner ear* (Grunstra et al.). 

We provide 5 scripts that correspond to the analyses performed using R, as described in the main text. Five files are provided: 

* `00_Definitions.R`: landmark definition, preparation of contextual data, creation of the function `plotLaby` for visualization
* `01_GPA.R`: semilandmark sliding, GPA, outliers, visualization
* `02_ProcD.R`: pairwise Procrustes distances, cross-validations
* `03_PCA.R`: PCA of the Procrustes shape coordinates
* `04_PLS.R`: 2B-PLS, cross-validations, phylogenetic signal in the PLS scores, phylogenetic PLS
* `05_WI.R`: Wheatsheaf index (Arbuckle et al. 2014 *Meth. Ecol. Evol.*)

The first file `00_Definitions.R` has to be run before all the others, because it is where the landmarks, curves, contextual data and related visualization function and parameters are described. 

The other R scripts can be run independently (recommended, because there may be some variable names conflicts). 

The raw data is provided as separate files in the OSF repository <https://osf.io/9mtwh/?view_only=88a550a637f84615b5666d69278309d4>

Data files are: 

* `Grunstra_et_al_Afrotheria_contextual.csv`: contextual data => used for 2B-PLS and scatter plots of PCA and PLS
* `Grunstra_et_al_Afrotheria_landmarks_raw.csv`: 3D landmark coordinates (13 anatomical landmarks, 111 semilandmarks) on the bony labyrinth for the 40 taxa => used for the GPA
* `Grunstra_et_al_Afrotheria_landmarks_slid_not_aligned.csv`: 3D landmark coordinates (13 anatomical landmarks, 111 semilandmarks) on the bony labyrinth for the 40 taxa, after semilandmark sliding (in the software Mathematica, but similar algorithm provided in R) => used for the computation of Procrustes distances, the PCA and the PLS analyses
* `Grunstra_et_al_Afrotheria_phyl_tree.txt`: phylogenetic tree of the 40 taxa in NEXUS format => used for the phylogenetic PLS and the computation of the phylogenetic signal
