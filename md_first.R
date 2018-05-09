#############################################
##                                          #
## Basic analysis of GH's trajectory data   #
##                                          #
#############################################
# Read example trajectory file
library("bio3d")
trtfile <- system.file("1xxn_1_all.dcd")
trj <- read.dcd(trtfile)

# Read the starting PDB file to determine atom correspondence
pdbfile <- system.file("1xxn_1.pdb")
pdb <- read.pdb(pdbfile)

print(pdb)

# How many rows (frames) and columns (coords) present in trj
dim(trj)
ncol(trj) == length(pdb$xyz)

# Trajectory Frame Superposition on Calpha atoms
ca.inds <- atom.select(pdb, elety = "CA")
a.ind <- atom.select(a, chain="A", resno=87:103, elety="CA")
xyz <- fit.xyz(fixed = pdb$xyz, mobile = trj, 
           fixed.inds = ca.inds$xyz, 
           mobile.inds = ca.inds$xyz)

# Root Mean Square Deviation (RMSD)
rd <- rmsd(xyz[1, ca.inds$xyz], xyz[, ca.inds$xyz])
plot(rd, typ = "l", ylab = "RMSD", xlab = "Frame No.")
points(lowess(rd), typ = "l", col = "red", lty = 2, lwd = 2)
summary(rd)


# Root Mean Squared Fluctuations (RMSF)
rf <- rmsf(xyz[, ca.inds$xyz])
plot(rf, ylab = "RMSF", xlab = "Residue Position", typ="l")

# Principal Component Analysis
pc <- pca.xyz(xyz[, ca.inds$xyz])
plot(pc, col = bwr.colors(nrow(xyz)))

# Cluster in PC space
hc <- hclust(dist(pc$z[, 1:2]))
grps <- cutree(hc, k = 2)
plot(pc, col = grps)

# Cross-Correlation Analysis
cij <- dccm(xyz[, ca.inds$xyz])
plot(cij)

## view.dccm(cij, pdb, launch = TRUE)
