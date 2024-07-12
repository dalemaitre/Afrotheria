# ===========================================================================================
# Title: Afrotheria - Generalized Procrustes Analysis
# Author: Anne Le Maitre
# Date: 2023-10-23
# -------------------------------------------------------------------------------------------
# R version: 4.3.1
# Required extensions: geomorph, Morpho
# -------------------------------------------------------------------------------------------
# Aims: Perform the Generalized Procrustes Analysis (GPA), with or without sliding
# of the 13 anatomical landmarks and 111 semilandmarks placed on the bony labyrinth 
# for 40 taxa (20 Afrotheria, 20 non-afrotherian)
# and visualize the results as a 3D plot
# Input: 
# - Results of the script "00_Definitions.R"
# - Raw landmark coordinates as a CSV file
# Output: 
# - Results of a GPA without / with semilandmark sliding
# ===========================================================================================

# ---- Packages ----

library(geomorph)  # geometric morphometrics
library(Morpho)  # geometric morphometrics

# ---- Load raw landmark coordinates ----

# First, run the R script '00_Definitions'
source("00_Definitions.R")

# Load the raw landmark coordinates
raw_data <- read.csv("Grunstra_et_al_Afrotheria_landmarks_raw.csv")

# Conversion to a matrix
raw <- as.matrix(raw_data[, -1])

# Name rows and columns
rownames(raw) <- raw_data$Species

# Conversion to an array
A_raw <- arrayspecs(raw, p = ie$p, k = 3)

# ---- Check raw data ----

# Here we check if there is no problem with the raw landmark coordinates
# For this, we perform a generalized Procrustes analysis (GPA)

# GPA without sliding
A.gpa.noslid <- procSym(A_raw)

# Check outliers (Procrustes distance from the mean)
outliers <- plotOutliers(A.gpa.noslid$rotated)  # outliers

# Visualize problematic landmarks: average shape vs. outlier
spi <- outliers[1]  # first outlier
plotRefToTarget(A.gpa.noslid$mshape, A.gpa.noslid$rotated[, , spi], 
                method = "vector", label = FALSE, links = ie$links)

# Note that we find the dolphin (*Sotalia fluviatilis*) 
# as a strong outlier, which is normal given its labyrinth shape
# The two next taxa the further away from the mean are 
# the manatee (*Trichechus* sp.) and the dugong (*Dugong dugon*)

# ---- GPA with sliding ----

# Here we perform a GPA with semilandmark sliding
# Due to the huge variation in cochlear morphology, 
# the magnitude of sliding per iteration is reduced by a factor of 0.1

# GPA with sliding
A.gpa.slid <- procSym(A_raw, SMvector = ie$slm, outlines = ie$curves, 
                          stepsize = 0.1)

# Check the sliding: 
# Visualize the difference between slid data in the original space and raw data
spi <- 1  # choose a specimen by number
plotRefToTarget(A_raw[, , spi], A.gpa.slid$dataslide[, , spi], 
                method = "vector", label = FALSE, links = ie$links)

# Check outliers (Procrustes distance from the mean)
outliers <- plotOutliers(A.gpa.slid$rotated)  # outliers
# Note that we find similar outliers

# Visualize problematic landmarks: average shape vs. outlier
spi <- outliers[1]  # first outlier
plotRefToTarget(A.gpa.slid$mshape, A.gpa.slid$rotated[, , spi], 
                method = "vector", label = FALSE, links = ie$links)

# ---- Visualization of one specimen (one color) ----

# Here we give an example to visualize a given specimen
# Only the anatomical landmarks are numbered

# Landmark database: choose one
#A <- A_raw  # initial data
#A <- A.gpa.slid$dataslide  # slid data in the original space
A <- A.gpa.slid$rotated  # slid data in the new space

# Choose a specimen
#spi <- 1  # by number
spi <- which(dimnames(A)[[3]] == "Sotalia_fluviatilis")  # by species name

# Plot all landmarks as points
plot3d(A[,, spi], type = "p", size = 1, 
       aspect = TRUE, decorate = F)  # plot all landmarks
title3d(main = dimnames(A)[[3]][spi])

# Anatomical landmarks
points3d(x = A[ie$lm, , spi], col = "red")
texts3d(x = A[ie$lm, , spi], 
        texts = ie$lm, 
        col = "red", adj = 2)  # fixed landmark numbers
# Links
for (j in 1:nrow(ie$links)) { 
  segments3d(rbind(A[ie$links[j, 1], , spi], 
                   A[ie$links[j, 2], , spi]), 
             col = "blue")
}
aspect3d("iso")

# ---- Visualization of one specimen (three colors) ----

# Here we give an example to visualize a given specimen
# We color semicircular canals in blue, cochlea in red and windows in green

# Landmark database: choose one
#A <- A_raw  # initial data
#A <- A.gpa.slid$dataslide  # slid data in the original space
A <- A.gpa.slid$rotated  # slid data in the new space

# Choose a specimen
#spi <- 1  # by number
spi <- which(dimnames(A)[[3]] == "Sotalia_fluviatilis")  # by species name

# 3D visualization
plotLaby(A[, , spi], title = dimnames(A)[[3]][spi])
