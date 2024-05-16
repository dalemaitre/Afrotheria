# ===========================================================================================
# Title: Afrotheria - Computation of the Wheatsheaf index
# Authors: Nicole D. S. Grunstra, Anne Le Maitre
# Date: 2024-05-16
# -------------------------------------------------------------------------------------------
# R version: 4.3.3
# Required extensions: geomorph, Morpho, ape, phytools
# -------------------------------------------------------------------------------------------
# Aims: Compute Wheatsheaf index (Arbuckle et al. 2014 Methods Ecol Evol)
# for the 40 taxa (20 Afrotheria, 20 non-afrotherian analogues)
# Input: 
# - Contextual data: Results of the script "00_Definitions.R"
# - Slid landmark coordinates as a CSV file
# Output: Wheatsheaf index
# ===========================================================================================

# ---- Packages ----

library(geomorph)  # geometric morphometrics
library(Morpho)  # geometric morphometrics
library(ape)  # phylogenetic analyses
library(phytools)  # phylogenetic analyses

# ---- Load Procrustes coordinates ----

# Before using the present script, run the script '00_Definitions.R' 
# to get contextual data (in particular the section 'Contextual data')

# Here we directly load the slid landmark coordinates
# They are in the original coordinate space

# Load the slid landmark coordinates
slid_data <- read.csv("Grunstra_et_al_Afrotheria_landmarks_slid_not_aligned.csv")

# Conversion to a matrix
slid <- as.matrix(slid_data[, -1])

# Name rows and columns
rownames(slid) <- slid_data$Species

# Conversion to an array
A_slid <- arrayspecs(slid, p = 124, k = 3)  # 124 = 13 anatomical landmarks + 111 semilandmarks

# Generalized Procrustes superimposition (no sliding)
A.gpa <- procSym(A_slid)
A.coord <- two.d.array(A.gpa$rotated)  # Procrustes coordinates as a matrix

# Procrustes distances for all pairs of species
procD <- dist(A.coord)  # Euclidean distances
m <- as.matrix(procD)  # create the matrix with all pairwise distances

# ---- Phylogenetic correction ----

# Load the phylogenetic tree
tree <- read.nexus("Grunstra_et_al_Afrotheria_phyl_tree.txt")

# Re-scale tree to a total length (= height) of 1
tree_rescaled <- tree
# The maximum node height (247.911) is the total height (= length) of the tree
tree_rescaled$edge.length <- tree_rescaled$edge.length / max(nodeHeights(tree)) 
# edge (= branch) 78, which separates O. anatinus from all others, 
# is now scaled to a length of 1 (total tree length)

# Phylogenetic distance matrix that includes the shared (proportional) branch lengths, 
# i.e. shared evolutionary time
phyl_dist_mat <- vcv.phylo(tree_rescaled)
# VCV matrix of a re-scaled tree yields covariances that equal proportional shared branch lengths

## Note: The above operation is equal to computing the correlation matrix on the original tree
#phyl_dist_mat2 <- vcv.phylo(tree, corr = T)
#all.equal(phyl_dist_mat, phyl_dist_mat2) # TRUE

# Reorder dimnames to match between the different matrices
reordered_phyl_dist_mat <- phyl_dist_mat[, match(colnames(m), colnames(phyl_dist_mat))]
reordered_phyl_dist_mat <- reordered_phyl_dist_mat[match(rownames(m), rownames(phyl_dist_mat)), ]
identical(dimnames(m), dimnames(reordered_phyl_dist_mat)) # dimnames match b/w reordered phyl dist & m (40x40 pairwise Proc dist matrix)

# correct phenotypic dist matrix for phyl. relatedness (formula (1) from Arbuckle et al. 2015 Methods Ecol Evol)
m.phyl <- m / (1 - log(reordered_phyl_dist_mat + 0.01)) # still a symmetrical matrix (upper + lower triangles)
m.phyl.low <- m.phyl
m.phyl.low[upper.tri(m.phyl.low)] <- 0 # keep upper triangle
d.phyl <- as.dist(m.phyl.low) # phyl.-corrected phenotypic distance matrix of all pairwise comparisons

# without phylogenetic correction
# (just out of curiosity, for comparison to results pairwise Procrustes distances cf. 02_ProcD.R)
#m_low <- m
#m_low[upper.tri(m_low)] <- 0
#dist_m2 <- as.dist(m_low)
#mean(dist_m2) # average pairwise distance NOT corrected for phylogeny = 0.2802285
# average afrotheria-analogue distance NOT corrected for phylogeny (see 02_ProcD.R) = 0.2445418
# ratio: 0.2802285 / 0.2445418 = 1.1459

# ---- Distances among pairs of analogues and pairs of non-analogues ----

# compute phyl-corrected distances b/w Afrotheria vs. non-Afrotheria from which to select focal species later
m.af_naf.phyl <- m.phyl[a, -a]  # distances of Afrotheria (rows) vs. all non-afrotherian taxa (columns) from phyl-corr. matrix above
d.af_naf.phyl <- as.dist(m.af_naf.phyl)  # all 190 distances 

# Remove Elephas maximus (it has no analog)
em <- which(rownames(m.af_naf.phyl) == "Elephas_maximus")
m.af_naf.noEm.phyl <- m.af_naf.phyl[-em, ]

# Define analogs for each Afrotheria (repeated from above, not necessary to repeat)
analogs <- list(
  Amblysomus_hottentotus = c("Cryptomys_hottentotus", "Talpa_europaea"), 
  Chrysochloris_asiatica = c("Cryptomys_hottentotus", "Talpa_europaea"), 
  Heterohyrax_brucei = c("Marmota_marmota", "Ochotona_rufescens"), 
  Procavia_capensis = c("Marmota_marmota", "Ochotona_rufescens"), 
  Dugong_dugon = c("Odobenus_rosmarus", "Phoca_vitulina", "Sotalia_fluviatilis"), 
  Trichechus_sp = c("Odobenus_rosmarus", "Phoca_vitulina", "Sotalia_fluviatilis"), 
  Elephantulus_rufescens = c("Microtus_arvalis", "Mus_musculus", "Rattus_rattus", "Sorex_araneus"), 
  Elephantulus_rupestris = c("Microtus_arvalis", "Mus_musculus", "Rattus_rattus", "Sorex_araneus"), 
  Macroscelides_proboscideus = c("Microtus_arvalis", "Mus_musculus", "Rattus_rattus", "Sorex_araneus"), 
  Petrodromus_tetradactylus = c("Oryctolagus_cuniculus", "Rattus_rattus", "Tragulus_javanicus"), 
  Rhynchocyon_cirnei = c("Oryctolagus_cuniculus", "Rattus_rattus", "Tragulus_javanicus"), 
  Hemicentetes_semispinosus = "Erinaceus_europaeus", 
  Setifer_setosus = "Erinaceus_europaeus", 
  Tenrec_ecaudatus = "Erinaceus_europaeus", 
  Limnogale_mergulus = c("Sorex_araneus", "Ornithorhynchus_anatinus"), 
  Oryzorictes_tetradactylus = "Talpa_europaea", 
  Orycteropus_afer = c("Manis_tricuspis", "Tolypeutes_matacus"), 
  Potamogale_velox = c("Lutra_lutra", "Ornithorhynchus_anatinus"), 
  Dendrohyrax_arboreus = c("Didelphis_sp", "Tamandua_mexicana")
)

# Create a matrix (rows = Afrotheria ; columns = non-afrotherian taxa)
# with logical values:
# - TRUE for pairs of analogues
# - FALSE for pairs of non-analogues
m.an.phyl <- matrix(FALSE, 
               nrow = nrow(m.af_naf.noEm.phyl), 
               ncol = ncol(m.af_naf.noEm.phyl), 
               dimnames = dimnames(m.af_naf.noEm.phyl))
for (i in 1:length(analogs)) {
  ri <- which(rownames(m.an.phyl) == names(analogs)[i])  # Afrotheria species
  ci <- which(colnames(m.an.phyl) %in% analogs[[i]])
  m.an.phyl[ri, ci] <- TRUE
}
# this matrix is identical to the one above in the script "02_ProcD.R",
# but adjusted for the phylogenetic-corrected matrix name.

# Procrustes (= Euclidean) distances for Afrotheria vs. analogs (i.e. "focal species", Arbuckle et al. 2015)
d.af_an.phyl <- m.af_naf.noEm.phyl[m.an.phyl]

# mean pairwise Proc distances for all species pairs and for focal species only (after correcting for phylogeny)
mean.all <- mean(d.phyl) # 0.1862664
mean.focal <- mean(d.af_an.phyl) # 0.1574017

# Wheatsheaf index
w <- mean.all/mean.focal
w # 1.183382!! --> case for convergence

(mean.all - mean.focal) / mean.focal
# total sample is ~18% more dissimilar in inner ear shape than ecomorphological analogues

