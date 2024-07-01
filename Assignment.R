


# Load necessary libraries
library(readxl)
library(ggplot2)
library(GGally)


# Load the Data

data <- read_excel("C:\\Users\\user\\Desktop\\cleaned_dataset.xlsx")

# Extracting timedata and genedata expression data from the tibbles
timedata <- data$`Time (min)`
genesdata <- data[, -1]

# Function to calculate y-axis limits with 0.2 increment
#calc_ylim <- function(values) {
# ymin <- floor(min(values) / 0.2) * 0.2
# ymax <- ceiling(max(values) / 0.2) * 0.2
# return(c(ymin, ymax))


# Save timedata series plots to a PNG file
png(filename = "E:\\Study materials\\Masters\\Statitistcs\\Assignment\\time_series_plots1.png", width = 800, height = 1200)

# Plotting timedata series for each gene
par(mfrow = c(5, 1))  # Set up a 5x1 grid for plots
for (i in 1:ncol(genesdata)){
  # ylim <- calc_ylim(genesdata[[i]])
  plot(timedata, genesdata[[i]], type='l', col=i+1, xlab='Time (min)', ylab=paste('Gene x', i), main= paste("Gene x", i, "Expression over Time"))
  
}

# Close the PNG device
dev.off()

# Reset plotting layout
par(mfrow=c(1, 1))  

### Distribution Plots
png(filename = "E:\\Study materials\\Masters\\Statitistcs\\Assignment\\distribution_plots.png", width = 800, height = 1200)
par(mfrow = c(5, 1))
# Plot histograms for each gene
for (i in 1:ncol(genesdata)) {
  hist(genesdata[[i]], breaks = 20, col = i, xlab = paste0('Gene x', i), main = paste0('Expression Distribution of Gene x', i))
}
# Close the PNG device
dev.off()

# Reset the plotting layout
par(mfrow = c(1, 1))

# Calculate the correlation matrix
correlation_matrix <- cor(genesdata)

# Save scatter plots to a PNG file
png(filename = "E:\\Study materials\\Masters\\Statitistcs\\Assignment\\correlation_scatter_plots.png", width = 1200, height = 1200)

# Use GGally for better visualization
ggpairs(genesdata, title = "Scatter Plots and Correlation Matrix")

# Close the PNG device
dev.off()




#Task 2.1
# Convert tibble to a data frame
data <- as.data.frame(data)
x1 <- as.numeric(data[, "x1"])
x2 <- as.numeric(data[, "x2"])
x3 <- as.numeric(data[, "x3"])
x4 <- as.numeric(data[, "x4"])
x5 <- as.numeric(data[, "x5"])

estimate_parameters <- function(X, y) {
  
  # Add bias term (if necessary)
  X <- cbind(1, X)  # Adding a column of 1s for bias term
  # Compute parameters θ using least squares
  theta <- solve(t(X) %*% X) %*% t(X) %*% y
  return(theta)
}
# Model 1: y = θ₁x₄ + θ₂x₃² + θbias
model1_X <- cbind(x4, x3^2)
model1_theta <- estimate_parameters(model1_X, x2)

# Model 2: y = θ₁x₄ + θ₂x₃² + θ₃x₅ + θbias
model2_X <- cbind(x4, x3^2, x5)
model2_theta <- estimate_parameters(model2_X, x2)

# Model 3: y = θ₁x₃ + θ₂x₄ + θ₃x₅³
model3_X <- cbind(x3, x4, x5^3)
model3_theta <- estimate_parameters(model3_X, x2)

# Model 4: y = θ₁x₄ + θ₂x₃² + θ₃x₅³ + θbias
model4_X <- cbind(x4, x3^2, x5^3)
model4_theta <- estimate_parameters(model4_X, x2)

# Model 5: y = θ₁x₄ + θ₂x₁² + θ₃x₃² + θbias
model5_X <- cbind(x4, x1^2, x3^2)
model5_theta <- estimate_parameters(model5_X, x2)

# Print estimated parameters
print("Model 1 Parameters:")
print(model1_theta)

print("Model 2 Parameters:")
print(model2_theta)

print("Model 3 Parameters:")
print(model3_theta)

print("Model 4 Parameters:")
print(model4_theta)

print("Model 5 Parameters:")
print(model5_theta)


# Compute fitted values (y-hat) for each model
fitted_model1 <- cbind(1, model1_X) %*% model1_theta
fitted_model2 <- cbind(1, model2_X) %*% model2_theta
fitted_model3 <- cbind(1, model3_X) %*% model3_theta
fitted_model4 <- cbind(1, model4_X) %*% model4_theta
fitted_model5 <- cbind(1, model5_X) %*% model5_theta


# Create a dataframe for plotting
plot_data <- data.frame(
  Time = data$Time,
  Actual = x2,
  Fitted_Model1 = fitted_model1,
  Fitted_Model2 = fitted_model2,
  Fitted_Model3 = fitted_model3,
  Fitted_Model4 = fitted_model4,
  Fitted_Model5 = fitted_model5
)

# Plot actual vs fitted values for each model
ggplot(plot_data, aes(x = Time)) +
  geom_line(aes(y = Actual, color = "Actual")) +
  geom_line(aes(y = Fitted_Model1, color = "Model 1")) +
  geom_line(aes(y = Fitted_Model2, color = "Model 2")) +
  geom_line(aes(y = Fitted_Model3, color = "Model 3")) +
  geom_line(aes(y = Fitted_Model4, color = "Model 4")) +
  geom_line(aes(y = Fitted_Model5, color = "Model 5")) +
  labs(title = "Actual vs Fitted Values for Each Model",
       y = "Gene Expression (x2)",
       color = "Legend") +
  theme_minimal()



#task 2.2
# Define a function to compute RSS
compute_rss <- function(X, y, theta) {
  X <- cbind(1, X)  # Adding a column of 1s for bias term
  residuals <- y - X %*% theta
  rss <- sum(residuals^2)
  return(rss)
}

# Compute RSS for each model
rss_model1 <- compute_rss(model1_X, x2, model1_theta)
rss_model2 <- compute_rss(model2_X, x2, model2_theta)
rss_model3 <- compute_rss(model3_X, x2, model3_theta)
rss_model4 <- compute_rss(model4_X, x2, model4_theta)
rss_model5 <- compute_rss(model5_X, x2, model5_theta)

# Print RSS values
print(paste("Model 1 RSS:", rss_model1))
print(paste("Model 2 RSS:", rss_model2))
print(paste("Model 3 RSS:", rss_model3))
print(paste("Model 4 RSS:", rss_model4))
print(paste("Model 5 RSS:", rss_model5))



#task2.3
# Define a function to compute the log-likelihood
compute_log_likelihood <- function(n, rss) {
  sigma_squared <- rss / (n - 1)
  log_likelihood <- - (n / 2) * log(2 * pi) - (n / 2) * log(sigma_squared) - (1 / (2 * sigma_squared)) * rss
  return(log_likelihood)
}

# Number of data samples
n <- nrow(data)

# Compute log-likelihood for each model
log_likelihood_model1 <- compute_log_likelihood(n, rss_model1)
log_likelihood_model2 <- compute_log_likelihood(n, rss_model2)
log_likelihood_model3 <- compute_log_likelihood(n, rss_model3)
log_likelihood_model4 <- compute_log_likelihood(n, rss_model4)
log_likelihood_model5 <- compute_log_likelihood(n, rss_model5)

# Print log-likelihood values
print(paste("Model 1 Log-Likelihood:", log_likelihood_model1))
print(paste("Model 2 Log-Likelihood:", log_likelihood_model2))
print(paste("Model 3 Log-Likelihood:", log_likelihood_model3))
print(paste("Model 4 Log-Likelihood:", log_likelihood_model4))
print(paste("Model 5 Log-Likelihood:", log_likelihood_model5))



#task 2.4
# Define a function to compute AIC
compute_aic <- function(k, log_likelihood) {
  aic <- 2 * k - 2 * log_likelihood
  return(aic)
}

# Define a function to compute BIC
compute_bic <- function(k, n, log_likelihood) {
  bic <- k * log(n) - 2 * log_likelihood
  return(bic)
}

# Number of parameters (including bias term)
k_model1 <- 3
k_model2 <- 4
k_model3 <- 4
k_model4 <- 4
k_model5 <- 3

# Compute AIC for each model
aic_model1 <- compute_aic(k_model1, log_likelihood_model1)
aic_model2 <- compute_aic(k_model2, log_likelihood_model2)
aic_model3 <- compute_aic(k_model3, log_likelihood_model3)
aic_model4 <- compute_aic(k_model4, log_likelihood_model4)
aic_model5 <- compute_aic(k_model5, log_likelihood_model5)

# Compute BIC for each model
bic_model1 <- compute_bic(k_model1, n, log_likelihood_model1)
bic_model2 <- compute_bic(k_model2, n, log_likelihood_model2)
bic_model3 <- compute_bic(k_model3, n, log_likelihood_model3)
bic_model4 <- compute_bic(k_model4, n, log_likelihood_model4)
bic_model5 <- compute_bic(k_model5, n, log_likelihood_model5)

# Print AIC and BIC values
print(paste("Model 1 AIC:", aic_model1))
print(paste("Model 1 BIC:", bic_model1))
print(paste("Model 2 AIC:", aic_model2))
print(paste("Model 2 BIC:", bic_model2))
print(paste("Model 3 AIC:", aic_model3))
print(paste("Model 3 BIC:", bic_model3))
print(paste("Model 4 AIC:", aic_model4))
print(paste("Model 4 BIC:", bic_model4))
print(paste("Model 5 AIC:", aic_model5))
print(paste("Model 5 BIC:", bic_model5))



#Task 2.5
library(gridExtra)

# Compute residuals for each model
residuals_model1 <- x2 - fitted_model1
residuals_model2 <- x2 - fitted_model2
residuals_model3 <- x2 - fitted_model3
residuals_model4 <- x2 - fitted_model4
residuals_model5 <- x2 - fitted_model5

# Create Q-Q plots for residuals of each model
qq_plot_model1 <- ggplot(data = data.frame(Residuals = residuals_model1), aes(sample = Residuals)) +
  stat_qq() +
  stat_qq_line() +
  labs(title = "Q-Q Plot of Residuals (Model 1)")

qq_plot_model2 <- ggplot(data = data.frame(Residuals = residuals_model2), aes(sample = Residuals)) +
  stat_qq() +
  stat_qq_line() +
  labs(title = "Q-Q Plot of Residuals (Model 2)")

qq_plot_model3 <- ggplot(data = data.frame(Residuals = residuals_model3), aes(sample = Residuals)) +
  stat_qq() +
  stat_qq_line() +
  labs(title = "Q-Q Plot of Residuals (Model 3)")

qq_plot_model4 <- ggplot(data = data.frame(Residuals = residuals_model4), aes(sample = Residuals)) +
  stat_qq() +
  stat_qq_line() +
  labs(title = "Q-Q Plot of Residuals (Model 4)")

qq_plot_model5 <- ggplot(data = data.frame(Residuals = residuals_model5), aes(sample = Residuals)) +
  stat_qq() +
  stat_qq_line() +
  labs(title = "Q-Q Plot of Residuals (Model 5)")

# Create histograms of residuals for each model
histogram_model1 <- ggplot(data = data.frame(Residuals = residuals_model1), aes(x = Residuals)) +
  geom_histogram(bins = 20, fill = "lightblue", color = "black", alpha = 0.7) +
  labs(title = "Histogram of Residuals (Model 1)", x = "Residuals", y = "Frequency")

histogram_model2 <- ggplot(data = data.frame(Residuals = residuals_model2), aes(x = Residuals)) +
  geom_histogram(bins = 20, fill = "lightblue", color = "black", alpha = 0.7) +
  labs(title = "Histogram of Residuals (Model 2)", x = "Residuals", y = "Frequency")

histogram_model3 <- ggplot(data = data.frame(Residuals = residuals_model3), aes(x = Residuals)) +
  geom_histogram(bins = 20, fill = "lightblue", color = "black", alpha = 0.7) +
  labs(title = "Histogram of Residuals (Model 3)", x = "Residuals", y = "Frequency")

histogram_model4 <- ggplot(data = data.frame(Residuals = residuals_model4), aes(x = Residuals)) +
  geom_histogram(bins = 20, fill = "lightblue", color = "black", alpha = 0.7) +
  labs(title = "Histogram of Residuals (Model 4)", x = "Residuals", y = "Frequency")

histogram_model5 <- ggplot(data = data.frame(Residuals = residuals_model5), aes(x = Residuals)) +
  geom_histogram(bins = 20, fill = "lightblue", color = "black", alpha = 0.7) +
  labs(title = "Histogram of Residuals (Model 5)", x = "Residuals", y = "Frequency")

# # Display plots
# library(gridExtra)

# Combine Q-Q plots and histograms for each model into a grid
grid_plot <- grid.arrange(qq_plot_model1, histogram_model1,
                          qq_plot_model2, histogram_model2,
                          qq_plot_model3, histogram_model3,
                          qq_plot_model4, histogram_model4,
                          qq_plot_model5, histogram_model5,
                          nrow = 5, top = "Residuals Analysis for Each Model")

# Save the combined plot grid as a PNG file
ggsave( "E:\\Study materials\\Masters\\Statitistcs\\Assignment\\residuals_analysis_grid.png", plot = grid_plot, width = 14, height = 18, dpi = 300)



#task 2.7

n <- nrow(data)
train_indices <- sample(1:n, 0.7 * n,replace = FALSE)  # 70% for training
test_indices <- setdiff(1:n, train_indices)  # Remaining for testing

# Training data
x_train <- data[train_indices, c("x1", "x3", "x4", "x5")]
y_train <- data[train_indices, "x2"]
# Testing data
x_test <- data[test_indices, c("x1", "x3", "x4", "x5")]
y_test <- data[test_indices, "x2"]

# Model 5: y = θ₁x₄ + θ₂x₁² + θ₃x₃² + θbias
X_train <- as.matrix(x_train)
Y_test <- as.matrix(x_test)
# Estimate parameters using training data
train_X <- cbind(1, x_train[, "x4"], x_train[, "x1"]^2, x_train[, "x3"]^2)
theta <- solve(t(train_X) %*% train_X) %*% t(train_X) %*% y_train

# Predict using test data
test_X <- cbind(1, x_test[, "x4"], x_test[, "x1"]^2, x_test[, "x3"]^2)
y_pred <- test_X %*% theta

# Calculate residuals for error bars
residuals <- y_test - y_pred
stderr <- sqrt(sum(residuals^2) / (length(y_test) - 4))
ci_95 <- 1.96 * stderr

# Plotting predictions with error bars (95% confidence intervals)
plot_data <- data.frame(Time = data$Time[test_indices], Actual = y_test, Predicted = y_pred)

ggplot(plot_data, aes(x = Time, y = Actual)) +
  geom_point(color = "blue") +
  geom_line(aes(y = Predicted), color = "red") +
  geom_errorbar(aes(ymin = Predicted - ci_95, ymax = Predicted + ci_95), width = 0.1, color = "red", alpha = 0.7) +
  labs(title = "Model 5: Actual vs Predicted with 95% CI", x = "Time (min)", y = "Gene Expression (x2)") +
  theme_minimal()



#task 3
# Estimated parameters for Model 5
theta_bias <- 1.2950813
theta1_hat <- 0.8303910
theta2_hat <- 0.5374454
theta3_hat <- 0.1095380

# Prior distributions (uniform around estimated values)
theta1_prior <- runif(1000, min = 1.0, max = 1.5)  # Adjust range as per your estimation
theta2_prior <- runif(1000, min = -1.0, max = -0.5)  # Adjust range as per your estimation

# Number of samples and tolerance level
N <- 1000  # Number of samples
epsilon <- 8.62241676218405 *50 # RSS value of Model 5

accepted_samples <- list()

# Simulated data using Model 5
for (i in 1:N) {
  theta1 <- theta1_prior[i]
  theta2 <- theta2_prior[i]
  # Simulate data using Model 5 with fixed theta3 and theta_bias
  simulated_y <- theta1 * data$x4 + theta2 * data$x1^2 + theta3_hat * data$x3^2 + theta_bias
  # Compute RSS
  residuals <- data$x2 - simulated_y
  rss_simulated <- sum(residuals^2)
  # Accept sample if within tolerance
  if (rss_simulated < epsilon) {
    accepted_samples[[length(accepted_samples) + 1]] <- c(theta1 = theta1, theta2 = theta2)
  }
}

# Convert accepted samples to a data frame
accepted_samples <- do.call(rbind, accepted_samples)

# Check number of accepted samples
n_accepted <- nrow(accepted_samples)
print(paste("Number of accepted samples:", n_accepted))

# Joint posterior distribution plot
ggplot(accepted_samples, aes(x = theta1, y = theta2)) +
  geom_point(alpha = 0.5) +
  labs(title = "Joint Posterior Distribution of theta1 and theta2",
       x = "theta1", y = "theta2") +
  theme_minimal()

# Marginal posterior distribution plot for theta1
ggplot(accepted_samples, aes(x = theta1)) +
  geom_histogram(fill = "lightblue", color = "black", bins = 30) +
  labs(title = "Marginal Posterior Distribution of theta1",
       x = "theta1", y = "Frequency") +
  theme_minimal()

# Marginal posterior distribution plot for theta2
ggplot(accepted_samples, aes(x = theta2)) +
  geom_histogram(fill = "lightblue", color = "black", bins = 30) +
  labs(title = "Marginal Posterior Distribution of theta2",
       x = "theta2", y = "Frequency") +
  theme_minimal()


