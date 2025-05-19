library(ape)
library(vegan)
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

# create geographic distances to haplotype level
haplotypes <- unique(prevalence.data$lepA_haplotype)
n <- length(haplotypes)
geo.dist.haplo <- matrix(0, n, n)
rownames(geo.dist.haplo) <- colnames(geo.dist.haplo) <- haplotypes

for (i in 1:n) {
  for (j in 1:n) {
    samples.i <- prevalence.data$sample.id[prevalence.data$lepA_haplotype == haplotypes[i]]
    samples.j <- prevalence.data$sample.id[prevalence.data$lepA_haplotype == haplotypes[j]]
    if (length(samples.i) > 0 && length(samples.j) > 0) {
      #take the average geographical distance between each haplotype
      geo.dist.haplo[i, j] <- mean(geo.dist.matrix[samples.i, samples.j], na.rm = TRUE)
    }
  }
}

# Create a matrix of host associations at haplotype level
host.matrix <- model.matrix(~ dicty.species.id - 1, data = prevalence.data)
prevalence.data$lepA_haplotype <- as.factor(prevalence.data$lepA_haplotype)
host.agg <- aggregate(host.matrix, by = list(prevalence.data$lepA_haplotype), FUN = mean)
host.matrix <- as.matrix(host.agg[, -1])
rownames(host.matrix) <- host.agg[, 1]

# make all matrices align
haplotype.matrix <- agricolaris.phy.dist.mat
haplotype.matrix <- haplotype.matrix[rownames(haplotype.matrix) %in% rownames(geo.dist.haplo), colnames(haplotype.matrix) %in% colnames(geo.dist.haplo)]

# Perform dbRDA
dbRDA.model <- capscale(haplotype.matrix ~ geo.dist.haplo + host.matrix, data = prevalence.data)
dbRDA.model
# Perform ANOVA
anova(dbRDA.model, by = "terms")

#reduce model to look at geographic distance separately
#dbRDA.model <- capscale(haplotype.matrix ~ geo.dist.haplo, data = prevalence.data)
#anova(dbRDA.model, by = "terms")

#reduce model to look at host association separately
#dbRDA.model <- capscale(haplotype.matrix ~ host.matrix, data = prevalence.data)
#anova(dbRDA.model, by = "terms")

#a quick ordination plot
plot(dbRDA.model, type = "n")  # Create a blank plot
points(dbRDA.model, display = "sites", col = "blue", pch = 16)  # Add points for the samples
text(dbRDA.model, display = "bp", col = "red")  # Add text for the constraints

# Extract the site scores (positions of haplotypes)
haplotype_scores <- scores(dbRDA.model, display = "sites")
haplotype_scores_df <- as.data.frame(haplotype_scores)
haplotype_scores_df$haplotype <- rownames(haplotype_scores_df)

# Extract the biplot scores (arrows for explanatory variables)
biplot_scores <- scores(dbRDA.model, display = "bp")
biplot_scores_df <- as.data.frame(biplot_scores)
biplot_scores_df$variable <- rownames(biplot_scores_df)

# Write to CSV files
write.csv(haplotype_scores_df, "../data/dbrda_haplotype_scores.csv", row.names = FALSE)
write.csv(biplot_scores_df, "../data/dbrda_biplot_scores.csv", row.names = FALSE)






#get values for plotting
# match order
rows <- rownames(haplotype.matrix)
cols <- colnames(haplotype.matrix)
geo.dist.haplo <- geo.dist.haplo[match(rows, rownames(geo.dist.haplo)), match(cols, colnames(geo.dist.haplo))]

phylogenetic.vector <- as.vector(haplotype.matrix)
geographic.vector <- as.vector(geo.dist.haplo)
df <- data.frame(phylogenetic.vector, geographic.vector)
write.csv(df, '../data/agricolaris_phylogenetic_geographic_vectors.csv', row.names = FALSE)
plot(geographic.vector, phylogenetic.vector)

#mantel test
mantel.model <- mantel(geo.dist.haplo, haplotype.matrix)
mantel.model
