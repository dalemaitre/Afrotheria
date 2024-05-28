# ===========================================================================================
# Title: Afrotheria - Definition of landmarks, curves, contextual data
# Author: Anne Le Maitre
# Date: 2023-11-20
# -------------------------------------------------------------------------------------------
# R version: 4.3.1
# Required extensions: -
# -------------------------------------------------------------------------------------------
# Aims: Define all parameters and functions for the analysis
# of the 13 anatomical landmarks and 111 semilandmarks placed on the bony labyrinth 
# for 40 taxa (20 Afrotheria, 20 non-afrotherian)
# and their visualization as a 3D plot
# Input: Contextual data as a CSV file
# Output: 
# - ie: definition of landmark parameters useful for further analyses
# - plotLaby: function to visualize the bony labyrinth in 3 different colors
# - diverse graphical parameters for the scatter plots
# ===========================================================================================

# ---- Set working directory (optional) ----

#setwd( "./Afrotheria")

# ---- Definitions: landmarks ----

# Here we define:
# - anatomical landmarks and semilandmarks
# - the curves along which semilandmarks should slide
# - the links between (semi)landmarks for data visualization
# In the definition of curves and links, abbreviations are used 
# to help recognize the different anatomical parts

# Create a list with all elements for the bony labyrinth
ie <- list(lm = NA, 
           slm = NA, 
           p = NA, 
           curves = NA, 
           l = NA, 
           links = NA)

# Landmark definition
ie$lm <- 1:13  # 13 anatomical landmarks
ie$slm <- 14:124  # 111 semilandmarks
ie$p <- length(ie$lm) + length(ie$slm)  # total number of landmarks and semilandmarks

# Curves (for sliding)
ie$curves <- list(
  asca = c(1, 14, 15, 2), 
  ascs = c(2, 16:33, 3), 
  psca = c(6, 34, 35, 5), 
  psps = c(5, 36:53, 3), 
  ccr = c(3, 54:58, 4), 
  lsca = c(9, 59, 60, 8), 
  lscs = c(8, 61:78, 7), 
  co = c(79:118, 10), 
  ow1 = c(13, 119:121, 12), 
  ow2 = c(12, 122:124, 13)
)

# Links (for visualisation)
ie$l <- list(
  asc = matrix(c(1, 14, 15, 2, 16:33, 
                 14, 15, 2, 16:33, 3), 
               nrow = 22, ncol = 2, byrow = FALSE), 
  psc = matrix(c(6, 34, 35, 5, 36:53, 
                 34, 35, 5, 36:53, 3), 
               nrow = 22, ncol = 2, byrow = FALSE), 
  ccr = matrix(c(3, 54:58, 
                 54:58, 4), 
               nrow = 6, ncol = 2, byrow = FALSE), 
  lsc = matrix(c(9, 59, 60, 8, 61:78, 
                 59, 60, 8, 61:78, 7), 
               nrow = 22, ncol = 2, byrow = FALSE), 
  co = matrix(c(79, 80:118, 
                80:118, 10), 
              nrow = 40, ncol = 2, byrow = FALSE), 
  ow = matrix(c(13, 119:121, 12, 122:124, 
                119:121, 12, 122:124, 13), 
              nrow = 8, ncol = 2, byrow = FALSE)
)

# All links = concatenation of previously defined links
ie$links <- rbind(ie$l$asc, 
                  ie$l$psc, 
                  ie$l$ccr, 
                  ie$l$lsc, 
                  ie$l$co, 
                  ie$l$ow)

# Create links for the semicircular canals only
ie$l$scc <- rbind(ie$l$asc, 
                  ie$l$psc, 
                  ie$l$lsc, 
                  ie$l$ccr)

# ---- Definitions: functions ----

# Here we define a function to plot 3D shape patterns

plotLaby <- function(M, title = NULL, color = c("blue", "green", "red"), links = ie$l) {
  
  # Plot all landmarks as points
  plot3d(M, type = "p", size = 1, 
         aspect = TRUE, decorate = F)  # plot all landmarks
  title3d(main = title)
  
  # Round window
  points3d(x = M[11, ], col = color[2])
  # Links (oval window)
  for (j in 1:nrow(links$ow)) { 
    segments3d(rbind(M[links$ow[j, 1], ], 
                     M[links$ow[j, 2], ]), 
               col = color[2])
  }
  # Links (cochlea)
  for (j in 1:nrow(links$co)) { 
    segments3d(rbind(M[links$co[j, 1], ], 
                     M[links$co[j, 2], ]), 
               col = color[3])
  }
  # Links (semicircular canal)
  for (j in 1:nrow(links$scc)) { 
    segments3d(rbind(M[links$scc[j, 1], ], 
                     M[links$scc[j, 2], ]), 
               col = color[1])
  }
  aspect3d("iso")
  
}

# ---- Load contextual data ----

# Load contextual data
Afrotheria <- read.csv("Grunstra_et_al_Afrotheria_contextual.csv")

# Define Afrotheria
a <- which(Afrotheria$Order %in% c("Afrosoricida", "Hyracoidea", "Macroscelidea", 
                                   "Proboscidea", "Sirenia", "Tubulidentata"))
afro <- rep("Afrotheria", nrow(Afrotheria))
afro[-a] <- "non Afrotheria"

# Symbol by group = Afrotheria or not
pch.ecol <- rep(1, nrow(Afrotheria))  # create a vector of length n (number of observations) with '1' (empty circle)
pch.ecol[a] <- 19  # all Afrotheria as '19' (filled circle)
pch.ecol[-a] <- 17  # all non Afrotheria as '17' (filled triangles)

# Color by habitat group (groups used for visualization only)
ecol <- factor(Afrotheria$HabitatGroup)  # habitat group as a factor
couleurs <- c("darkgreen", "orange", "lightblue", "darkblue", "darkred")  # colors
col.ecol <- rep("black", nrow(Afrotheria))  # create a vector of length n (number of observations) with 'black'
for (i in 1:nlevels(ecol)) {
  spi <- which(ecol == levels(ecol)[i])
  col.ecol[spi] <- couleurs[i]
}
names(col.ecol) <- Afrotheria$HabitatGroup

