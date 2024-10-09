remove(list = ls())
# Load necessary libraries
library(MASS)         # For multivariate normal generation
library(cluster)      # For spectral clustering
library(mvtnorm)      # For dmvnorm function
library(stats)        # For pnorm function
library(parallel)     # For parallel computing
set.seed(42)
# Parameters
n <- 200
d <- 100
mu <- c(10, rep(0, d - 1))  # Mean vector
t_values <- 0:100           # t values to loop over
B <- matrix(rnorm(d * 3), nrow = d, ncol = 3)  # Matrix B with N(0,1) entries

# Function to compute mislabeling rate
mislabeling_rate <- function(x, true_labels, pred_labels) {
  return(mean(true_labels != pred_labels))
}

# Function to generate covariance matrix
generate_covariance <- function(t, B) {
  return(t * B %*% t(B) + diag(d))  # Covariance matrix as per the problem
}

# Function to run for each t value
compute_mislabeling_rates <- function(t) {
  Sigma <- generate_covariance(t, B)
  
  # Generate data
  y <- sample(c(-1, 1), n, replace = TRUE)  # Generate labels
  X <- t(sapply(y, function(label) mvrnorm(1, label * mu, Sigma)))  # Generate data
  
  # Spectral clustering
  dist_mat <- dist(X)
  spectral_clustering <- pam(dist_mat, 2)  # Perform spectral clustering
  spectral_labels <- ifelse(spectral_clustering$clustering == 1, 1, -1)
  
  # Calculate mislabeling rate for spectral clustering
  spectral_rate <- mislabeling_rate(X, y, spectral_labels)
  
  # Calculate the SNR and optimal mislabeling rate
  SNR <- t(mu) %*% solve(Sigma) %*% mu
  optimal_rate <- pnorm(-sqrt(SNR))
  
  return(list(spectral_rate = spectral_rate, optimal_rate = optimal_rate))
}

# Use mclapply for parallel computation (for Linux/MacOS use mclapply directly)
# For Windows, use parLapply to create a cluster
cl <- makeCluster(detectCores() - 1)  # Create cluster using available cores

# Export variables to cluster
clusterExport(cl, c("mu", "n", "d", "B", "mvrnorm", "pam", "generate_covariance", "mislabeling_rate", "compute_mislabeling_rates"))

# Parallel computation over t_values
results <- parLapply(cl, t_values, compute_mislabeling_rates)

# Stop the cluster after computation
stopCluster(cl)

# Extract results
spectral_rates <- sapply(results, function(res) res$spectral_rate)
# optimal_rates <- sapply(results, function(res) res$optimal_rate)

# Plot the mislabeling rates
plot(t_values, spectral_rates, type = "p", col = "red", lwd = 2, 
     ylab = "Mislabeling Rate", xlab = "t", main = "Mislabeling Rate vs Variance", ylim = c(0, 1))
# lines(t_values, optimal_rates, col = "blue", lwd = 2)
legend("topright", legend = c("Spectral Clustering"), inset = 0.05, col = c("red"), pch = 21, pt.lwd = 2)

