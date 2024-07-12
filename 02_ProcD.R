# ===========================================================================================
# Title: Afrotheria - Computation of pairwise Procrustes distances
# Author: Anne Le Maitre
# Date: 2024-07-12
# -------------------------------------------------------------------------------------------
# R version: 4.3.1
# Required extensions: geomorph, Morpho
# -------------------------------------------------------------------------------------------
# Aims: Compute all pairwise Procrustes distances 
# for the 40 taxa (20 Afrotheria, 20 non-afrotherian), and compare group averages
# Input: 
# - Results of the script "00_Definitions.R"
# - Slid landmark coordinates as a CSV file
# Output: All pairwise Procrustes distances, 
# averaged by group (analogues / non-analogues / Afrotheria / non-Afrotheria), 
# and leave-one-out and leave-two-out cross validations
# ===========================================================================================

# ---- Packages ----

library(geomorph)  # geometric morphometrics
library(Morpho)  # geometric morphometrics

# ---- Load data ----

# First, run the R script '00_Definitions'
source("00_Definitions.R")

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

# ---- Procrustes distances ----

# Generalized Procrustes superimposition (no sliding)
A.gpa <- procSym(A_slid)

# Procrustes distances for all pairs of species
A.coord <- two.d.array(A.gpa$rotated)  # Procrustes coordinates as a matrix
procD <- dist(A.coord)  # Euclidean distances
m <- as.matrix(procD)  # create the matrix with all pairwise distances

# ---- Distances among all pairs within each group (Afrotheria / non-Afrotheria) ----

# Define Afrotheria
a <- which(Afrotheria$Order %in% c("Afrosoricida", "Hyracoidea", "Macroscelidea", 
                                   "Proboscidea", "Sirenia", "Tubulidentata"))

# All pairs of Afrotheria
m.af_af <- m[a, a]  # distances of all Afrotheria (rows) vs. all Afrotheria (columns)
m.af_af.low <- m.af_af
m.af_af.low[upper.tri(m.af_af.low)] <- 0  # keep only one lower triangle of the matrix
d.af_af <- as.dist(m.af_af.low)  # all distances

# All pair of non-afrotherian taxa
m.naf_naf <- m[-a, -a]  # distances of all non-afrotherian (rows) vs. all non-afrotherian (columns)
m.naf_naf.low <- m.naf_naf
m.naf_naf.low[upper.tri(m.naf_naf.low)] <- 0  # keep only one lower triangle of the matrix
d.naf_naf <- as.dist(m.naf_naf.low)  # all distances

# ---- Distances among pairs of analogues and pairs of non-analogues ----

# Afrotheria vs. non-Afrotheria
m.af_naf <- m[a, -a]  # distances of Afrotheria (rows) vs. all non-afrotherian taxa (columns)
d.af_naf <- as.dist(m.af_naf)  # all distances

# Remove Elephas maximus (it has no analog)
em <- which(rownames(m.af_naf) == "Elephas_maximus")
m.af_naf.noEm <- m.af_naf[-em, ]

# Define analogs for each Afrotheria
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
m.an <- matrix(FALSE, 
               nrow = nrow(m.af_naf.noEm), 
               ncol = ncol(m.af_naf.noEm), 
               dimnames = dimnames(m.af_naf.noEm))
for (i in 1:length(analogs)) {
  ri <- which(rownames(m.an) == names(analogs)[i])  # Afrotheria species
  ci <- which(colnames(m.an) %in% analogs[[i]])
  m.an[ri, ci] <- TRUE
}

# Afrotheria vs. analogs
d.af_an <- m.af_naf.noEm[m.an]

# Afrotheria vs. non-analogs
d.af_nan <- m.af_naf.noEm[!m.an]

# ---- Average distances ----

# Arithmetic mean
mean(d.af_an)  # Afrotheria vs. non-Afrotheria analogues = 0.2445418
mean(d.af_nan)  # Afrotheria vs. non-Afrotheria non-analogues = 0.2869142
mean(d.af_af)  # Afrotheria only = 0.2663105
mean(d.naf_naf)  # non-Afrotheria only = 0.2897568
# Difference analogues / non-analogues
(mean(d.af_nan) - mean(d.af_an))/mean(d.af_an)

# Median
median(d.af_an)  # Afrotheria vs. non-Afrotheria analogues = 0.2529802
median(d.af_nan)  # Afrotheria vs. non-Afrotheria non-analogues = 0.27962
median(d.af_af)  # Afrotheria only = 0.2863658
median(d.naf_naf)  # non-Afrotheria only = 0.2717115
# Difference analogues / non-analogues
(median(d.af_nan) - median(d.af_an))/median(d.af_an)

# ---- Visualization as all pairwise distances ----

# Create a data frame with all distances associated to a given group
d <- c(d.af_an, d.af_nan, d.af_af, d.naf_naf)
proc.d <- as.data.frame(d)
proc.d$grp <- c(rep("analogues", length(d.af_an)), 
                rep("non-analogues", length(d.af_nan)), 
                rep("among afrotherian", length(d.af_af)), 
                rep("among non-afrotherian", length(d.naf_naf)))

# Box plots
boxplot(d~grp, data = proc.d, 
        main = "Pairwise distances",
        col = c("gold", "chartreuse4", "blue", "red"), 
        las = 1, xlab = NULL, ylab = "Procrustes distances")

# Here we see that despite some overlap when considering pairwise distances, 
# the median distance for the pairs of analogues is low
# compared to the median for other groups

# ---- Leave-one-out cross-validation ----

# Here we check if the average distance (expressed as the mean or the median)
# for pairs of analogues is lower than the average for pairs of non-analogues, 
# when we remove one of the taxa

# Afrotheria vs. analogs
d.af_an <- m.af_naf.noEm[m.an]  # all distances
mean.af_an <- mean(d.af_an)
median.af_an <- median(d.af_an)

# Afrotheria vs. non-analogs
d.af_nan <- m.af_naf.noEm[!m.an]  # all distances
mean.af_nan <- mean(d.af_nan)
median.af_nan <- median(d.af_nan)

# Compare the average distance between analogues and non-analogues
mean_sup <- 0
median_sup <- 0

# Jackknife: leave-one-out
for (i in 1:(nrow(m.an) + ncol(m.an))) {
  # Remove the corresponding species
  if (i <= nrow(m.an)) {
    # Remove the corresponding line
    m.an.i <- m.an[-i, ]  # in the selection matrix
    m.af_naf.noEm.i <- m.af_naf.noEm[-i, ]  # in the distance matrix
  }
  if (i > nrow(m.an)) {
    # Remove the corresponding column
    m.an.i <- m.an[, -(i - nrow(m.an))]  # in the selection matrix
    m.af_naf.noEm.i <- m.af_naf.noEm[, -(i - nrow(m.an))]  # in the distance matrix
  }
  # Afrotheria vs. analogs
  d.af_an.i <- m.af_naf.noEm.i[m.an.i]  # all distances
  mean.af_an <- c(mean.af_an, mean(d.af_an.i))
  median.af_an <- c(median.af_an, median(d.af_an.i))
  # Afrotheria vs. non-analogs
  d.af_nan.i <- m.af_naf.noEm.i[!m.an.i]  # all distances
  mean.af_nan <- c(mean.af_nan, mean(d.af_nan.i))
  median.af_nan <- c(median.af_nan, median(d.af_nan.i))
  # Check if distance inferior in non-analogues
  if (mean(d.af_an.i) > mean(d.af_nan.i)) { mean_sup <- mean_sup + 1 }
  if (median(d.af_an.i) > median(d.af_nan.i)) { median_sup <- median_sup + 1 }
}

# Compare the average distance between analogues and non-analogues
mean_sup  # equals 0 if mean analogues < mean non-analogues for all replicates
median_sup  # equals 0 if median analogues < median non-analogues for all replicates

# Names
names(mean.af_an) <- c("all", rownames(m.an), colnames(m.an))
names(mean.af_nan) <- c("all", rownames(m.an), colnames(m.an))
names(median.af_an) <- c("all", rownames(m.an), colnames(m.an))
names(median.af_nan) <- c("all", rownames(m.an), colnames(m.an))

# Plot both mean and median
par(mfrow = c(1,2), mar = c(6, 4, 4, 2) + 0.1)
# Mean
boxplot(mean.af_an, mean.af_nan, 
        main = "Mean Procrustes distances", 
        ylab = "Procrustes distances", 
        names = c("analogues", "non-analogues"), 
        col = "lightblue", las = 2, cex.axis = 0.8)
points(x = c(1:2), y = c(mean.af_an[1], mean.af_nan[1]), 
       pch = 16, col = "blue")
# Median
boxplot(median.af_an, median.af_nan, 
        main = "Median Procrustes distances", 
        ylab = "Procrustes distances", 
        names = c("analogues", "non-analogues"), 
        col = "lightblue", las = 2, cex.axis = 0.8)
points(x = c(1:2), y = c(median.af_an[1], median.af_nan[1]), 
       pch = 16, col = "blue")
# Back to normal parameters
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)

# OPTIONAL export data
write.csv(rbind(mean.af_an, mean.af_nan, median.af_an, median.af_nan), 
          file = "Grunstra_et_al_Afrotheria_ProcD_LooCV.csv")

# ---- Leave-two-out cross-validation ----

# Here we check if the average distance (expressed as the mean or the median)
# for pairs of analogues is lower than the average for pairs of non-analogues, 
# when we remove two taxa

# Afrotheria vs. analogs
d.af_an <- m.af_naf.noEm[m.an]  # all distances
mean.af_an <- mean(d.af_an)
median.af_an <- median(d.af_an)

# Afrotheria vs. non-analogs
d.af_nan <- m.af_naf.noEm[!m.an]  # all distances
mean.af_nan <- mean(d.af_nan)
median.af_nan <- median(d.af_nan)

# Compare the average distance between analogues and non-analogues
mean_sup <- 0
median_sup <- 0

# Jackknife: leave-one-out
for (i in 2:(nrow(m.an) + ncol(m.an))) {
  for (j in 1:(i-1)) {
    # Remove the corresponding species
    if (i <= nrow(m.an)) {
      # Remove the corresponding lines
      m.an.i <- m.an[-c(i, j), ]  # in the selection matrix
      m.af_naf.noEm.i <- m.af_naf.noEm[-c(i, j), ]  # in the distance matrix
    }
    if (i > nrow(m.an) & j <= nrow(m.an)) {
      # Remove the corresponding line (j) and column (i)
      m.an.i <- m.an[-j , -(i - nrow(m.an))]  # in the selection matrix
      m.af_naf.noEm.i <- m.af_naf.noEm[-j, -(i - nrow(m.an))]  # in the distance matrix
    }
    if (i > nrow(m.an) & j > nrow(m.an)) {
      # Remove the corresponding columns
      m.an.i <- m.an[, -(c(i,j) - nrow(m.an))]  # in the selection matrix
      m.af_naf.noEm.i <- m.af_naf.noEm[, -(c(i,j) - nrow(m.an))]  # in the distance matrix
    }
    # Afrotheria vs. analogs
    d.af_an.i <- m.af_naf.noEm.i[m.an.i]  # all distances
    mean.af_an <- c(mean.af_an, mean(d.af_an.i))
    median.af_an <- c(median.af_an, median(d.af_an.i))
    # Afrotheria vs. non-analogs
    d.af_nan.i <- m.af_naf.noEm.i[!m.an.i]  # all distances
    mean.af_nan <- c(mean.af_nan, mean(d.af_nan.i))
    median.af_nan <- c(median.af_nan, median(d.af_nan.i))
    # Check if distance inferior in non-analogues
    if (mean(d.af_an.i) > mean(d.af_nan.i)) { mean_sup <- mean_sup + 1 }
    if (median(d.af_an.i) > median(d.af_nan.i)) { median_sup <- median_sup + 1 }
  }
}

# Compare the average distance between analogues and non-analogues
mean_sup  # equals 0 if mean analogues < mean non-analogues for all replicates
median_sup  # equals 0 if median analogues < median non-analogues for all replicates

# Names
names(mean.af_an) <- c("all", rownames(m.an), colnames(m.an))
names(mean.af_nan) <- c("all", rownames(m.an), colnames(m.an))
names(median.af_an) <- c("all", rownames(m.an), colnames(m.an))
names(median.af_nan) <- c("all", rownames(m.an), colnames(m.an))

# Plot both mean and median
par(mfrow = c(1,2), mar = c(6, 4, 4, 2) + 0.1)
# Mean
boxplot(mean.af_an, mean.af_nan, 
        main = "Mean Procrustes distances", 
        ylab = "Procrustes distances", 
        names = c("analogues", "non-analogues"), 
        col = "lightblue", las = 2, cex.axis = 0.8)
points(x = c(1:2), y = c(mean.af_an[1], mean.af_nan[1]), 
       pch = 16, col = "blue")
# Median
boxplot(median.af_an, median.af_nan, 
        main = "Median Procrustes distances", 
        ylab = "Procrustes distances", 
        names = c("analogues", "non-analogues"), 
        col = "lightblue", las = 2, cex.axis = 0.8)
points(x = c(1:2), y = c(median.af_an[1], median.af_nan[1]), 
       pch = 16, col = "blue")
# Back to normal parameters
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)

# OPTIONAL export data
write.csv(rbind(mean.af_an, mean.af_nan, median.af_an, median.af_nan), 
          file = "Grunstra_et_al_Afrotheria_ProcD_LtoCV.csv")

