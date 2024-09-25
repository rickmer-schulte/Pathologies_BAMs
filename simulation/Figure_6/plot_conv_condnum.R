# Load necessary libraries
library(ggplot2)
library(reshape2)
library(MASS)
library(mboost)
library(latex2exp)

# Convergence for different condition numbers
set.seed(2)
n <- 200
mean <- c(0, 0)
beta0 <- 0
beta1 <- 15
beta2 <- -8
sigma <- 1
K <- 1000

data_list <- data.frame("iter"=seq(1,K))

cond_num_vec <- c()
rho_cov_list <- seq(0,0.75, by=0.025)
for (rho_cov in rho_cov_list) {
  covariance <- matrix(c(1, rho_cov, rho_cov, 1), nrow = 2)  # Covariance matrix
  data <- mvrnorm(n, mu = mean, Sigma = covariance)
  x1 <- data[, 1]
  x2 <- data[, 2]
  y <- beta0 + beta1 * x1 + beta2 * x2 + rnorm(n, sd = sigma)
  X <- cbind(x1, x2)
  cond_num <- round(kappa(t(X)%*%X),2)
  m2_condnum <- mboost(y ~ bols(x1, intercept = FALSE) + bols(x2, intercept = FALSE), control = boost_control(mstop = 1, nu = 0.1))
  
  loss_val <- c()
  # Iterate through the sequence using index k
  min_val <- sum(residuals(m2_condnum[K+200])^2)
  for (k in 1:K) {
      loss_val <- append(loss_val, round(sum(residuals(m2_condnum[k])^2) - min_val, 10))
  }
  cond_num_vec <- append(cond_num_vec, cond_num)
  data_list <- cbind(data_list, loss_val)
}
colnames(data_list) <- c("iter", paste0("y", 1:length(rho_cov_list)))

replace_zero_with_NA <- function(x) {
  # Find the indices of zero values
  zero_indices <- which(x == 0)
  # Skip the first occurrence
  zero_indices <- zero_indices[-1]
  # Replace the remaining zeros with NA
  x[zero_indices] <- NA
  return(x)
}

# Apply the function to each column of the dataframe
data_list <- data.frame(lapply(data_list, replace_zero_with_NA))

# Reshape the data frame to a long format
df_long_condnum <- melt(data_list, id.vars = "iter")

# Map color_values to the variable
df_long_condnum$color_val <- rep(cond_num_vec, each = nrow(data_list))

saveRDS(df_long_condnum, file="conv_condnum.RDS")

p1 <- ggplot(df_long_condnum, aes(x = iter, y = value, group = variable, color = color_val)) +
  geom_line() +
  scale_y_continuous(trans='log10') +
  scale_color_gradient(low = "blue", high = "red", name="Condition number") +
  theme_minimal() +
  theme(text = element_text(size = 14),
        plot.title = element_text(hjust = 0.5, size = 11),
        legend.position = c(0.95, 0.99),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(0, 0, 0, 0),
        legend.box.background = element_blank(),
        legend.background = element_blank()) +
  xlab("Iterations") + ylab(expression("\u2113" * (beta) - "\u2113"^"*"))
p1
ggsave(file = "plot_condnum.jpeg", width = 8*0.9, height = 4*0.9)
