#load data
data = read.csv('../data/paraburk_intrahost_data.csv')
data <- na.omit(data)
#keep only confident entries
data <- data[data$confident == 1,]
#get agricolaris data only
data <- data[data$infect == 'B70', ]
#remove Cavanderia, as it it too distant
data <- data[data$taxonomy != 'Cavenderia_aureostipes', ]
#take log of normalized spore count
data$log_spore_norm_kp <- log(data$spore_norm_kp)

#fit model
model <- lm(log_spore_norm_kp~distance_from_discoideum, data)
summary(model)
confint(model)
coef(model)

plot(data$distance_from_discoideum, data$log_spore_norm_kp)
abline(model)

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
predictions <- predict(model, model.data, type = "response", se.fit = TRUE)
model.data$fit <- predictions$fit
model.data$se <- predictions$se.fit
model.data$lower <- predictions$fit - 1.96 * predictions$se.fit
model.data$upper <- predictions$fit + 1.96 * predictions$se.fit

write.csv(model.data, '../data/agricolaris_fitness_effect_dicty_phylo_model.csv')
