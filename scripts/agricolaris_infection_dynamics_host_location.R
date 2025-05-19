library(lme4)
library(car)
library(glmmTMB)

# Load data
collection.data <- read.csv('../data/total_prevalence_data.csv')
agricolaris.haplotype.data <- read.csv('../data/P_agricolaris_haplotype_data.csv')
agricolaris.phylogeny <- ape::read.tree("../data/P_agricolaris_haplotype_tree.nwk")

#merge dataframes
collection.data.w.agricolaris <- merge(collection.data, agricolaris.haplotype.data[, c("sample.id", "lepA_haplotype")], by = "sample.id", all.x = TRUE)

collection.data.w.agricolaris$P_agricolaris_presence_absense <- ifelse(
  is.na(collection.data.w.agricolaris$lepA_haplotype), 
  0, 
  1
)

#clean up data
data <- collection.data.w.agricolaris
data <- dplyr::select(data, dicty.species.id, sublocation, P_agricolaris_presence_absense, location)
data <- na.omit(data)

# Fit mixed-effects logistic regression model
mixed.model <- glmer(P_agricolaris_presence_absense ~ dicty.species.id + (1 | sublocation), 
      data = data,
      family = binomial)

qqPlot(residuals(mixed.model, type = "pearson"))


beta.model <- glmmTMB(P_agricolaris_presence_absense ~ dicty.species.id + (1 | sublocation),
                      data = data,
                      family = beta_family(link = "logit"))







# Fit logistic regression model
model.1 <- glm(P_agricolaris_presence_absense ~ sublocation + dicty.species.id, 
               data = data,
               family = binomial)

summary(model.1)

collection.data.w.agricolaris <- na.omit(collection.data.w.agricolaris)
collection.data.w.agricolaris$dicty.species.id <- na.omit(collection.data.w.agricolaris$dicty.species.id)

collection.data.w.agricolaris$sublocation <- as.factor(collection.data.w.agricolaris$sublocation)
collection.data.w.agricolaris$dicty.species.id <- as.factor(collection.data.w.agricolaris$dicty.species.id)

# Fit logistic regression model
model.1 <- glm(P_agricolaris_presence_absense ~ sublocation + dicty.species.id, 
             data = collection.data.w.agricolaris,
             family = binomial)


summary(model.1)
plot(model.1)


