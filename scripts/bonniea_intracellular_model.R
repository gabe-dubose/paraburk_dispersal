#load data
data = read.csv('../data/paraburk_intrahost_data.csv')
data <- na.omit(data)
#keep only confident entries
data <- data[data$confident == 1,]
#get agricolaris data only
data <- data[data$infect == 'B859', ]
#remove Cavanderia, as it it too distant
data <- data[data$taxonomy != 'Cavenderia_aureostipes', ]
#get infected and uninfected spores
data$infected_spores <- data$rfpper / 100 * data$total_spore
data$uninfected_spores <- data$total_spore - data$infected_spores
data$infected_spores <- round(data$infected_spores)
data$uninfected_spores <- round(data$uninfected_spores)

#fit logistic model
model <- glm(cbind(infected_spores, uninfected_spores) ~ distance_from_discoideum, 
             data = data, 
             family = quasibinomial())

summary(model)
plogis(confint(model))
plogis(coef(model))

#assess model
predicted.proportions <- predict(model, type = "response")

# Residuals from the model
residuals <- residuals(model, type = "deviance")


# Plot residuals vs. fitted values
plot(predicted.proportions, residuals,
     xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs. Fitted Values")
abline(h = 0, col = "red", lwd = 2)

# Q-Q plot of the residuals
qqnorm(residuals, main = "Q-Q Plot of Residuals")
qqline(residuals, col = "red", lwd = 2)
shapiro.test(residuals)

#get model predictions and write out
model.data <- data.frame(distance_from_discoideum = seq(min(data$distance_from_discoideum),
                                                      max(data$distance_from_discoideum), 
                                                      length.out = 100))

#get the predicted values and their standard errors
predictions <- predict(model, model.data, type = "link", se.fit = TRUE)
model.data$fit <- predictions$fit
model.data$se <- predictions$se.fit
model.data$lower <- predictions$fit - 1.96 * predictions$se.fit
model.data$upper <- predictions$fit + 1.96 * predictions$se.fit

#transform to the probability scale using plogis
model.data$fit <- plogis(model.data$fit)
model.data$lower <- plogis(model.data$lower)
model.data$upper <- plogis(model.data$upper)

write.csv(model.data, '../data/bonniea_spore_prevalence_dicty_phylo_model.csv')
