library(vegan)
library(ape)
library(geosphere)

# Load data
prevalence.data <- read.csv('../data/P_agricolaris_haplotype_data.csv')
agricolaris.phylogeny <- ape::read.tree("../data/P_agricolaris_haplotype_tree.nwk")

# Compute phylogenetic distance matrix
agricolaris.phy.dist.mat <- ape::cophenetic.phylo(agricolaris.phylogeny)

# Compute geographical distance matrix
coords <- as.matrix(prevalence.data[, c("longitude", "latitude")])
geo.dist.matrix <- distm(coords, fun = distHaversine)
rownames(geo.dist.matrix) <- prevalence.data$sample.id
colnames(geo.dist.matrix) <- prevalence.data$sample.id

# Create geographic distances at haplotype level
haplotypes <- unique(prevalence.data$lepA_haplotype)
n <- length(haplotypes)
geo.dist.haplo <- matrix(0, n, n)
rownames(geo.dist.haplo) <- colnames(geo.dist.haplo) <- haplotypes

for (i in 1:n) {
  for (j in 1:n) {
    samples.i <- prevalence.data$sample.id[prevalence.data$lepA_haplotype == haplotypes[i]]
    samples.j <- prevalence.data$sample.id[prevalence.data$lepA_haplotype == haplotypes[j]]
    if (length(samples.i) > 0 && length(samples.j) > 0) {
      # Take the average geographical distance between each haplotype
      geo.dist.haplo[i, j] <- mean(geo.dist.matrix[samples.i, samples.j], na.rm = TRUE)
    }
  }
}

# Create a matrix of host associations at the haplotype level
host.matrix <- model.matrix(~ dicty.species.id - 1, data = prevalence.data)
prevalence.data$lepA_haplotype <- as.factor(prevalence.data$lepA_haplotype)
host.agg <- aggregate(host.matrix, by = list(prevalence.data$lepA_haplotype), FUN = mean)
host.matrix <- as.matrix(host.agg[, -1])
rownames(host.matrix) <- host.agg[, 1]

# Ensure all matrices align
haplotype.matrix <- agricolaris.phy.dist.mat
haplotype.matrix <- haplotype.matrix[rownames(haplotype.matrix) %in% rownames(geo.dist.haplo), colnames(haplotype.matrix) %in% colnames(geo.dist.haplo)]

# Perform Variance Partitioning
varpart_result <- varpart(haplotype.matrix, geo.dist.haplo, host.matrix)

# Print the results
print(varpart_result)

# Optionally, plot the variance partitioning results
plot(varpart_result)
