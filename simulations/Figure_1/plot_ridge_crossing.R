# Load necessary libraries
library(glmnet)
library(tidyverse)
library(mboost)
library(viridis)
library(gridExtra)

# Set seed for reproducibility
set.seed(22)

# Generate synthetic data
n <- 100 # number of samples
p <- 2 # number of predictors

# Independent random variables
X1 <- rnorm(n)
X2 <- rnorm(n)

# Combine to induce slight empirical correlation
rho <- 0.7
X <- cbind(X1, rho * X1 + sqrt(1 - rho^2) * X2)

# True coefficients
beta <- c(3, -2)

# Response
y <- X %*% beta + rnorm(n)

# Perform Ridge Regression
lambda_value <- 1
ridge_model <- glmnet(X, y, alpha = 0, lambda = lambda_value, intercept = FALSE)

# Extract ridge coefficients
ridge_coefs <- as.vector(coef(ridge_model, s = lambda_value)[-1]) # omit the intercept

# Perform OLS (Ordinary Least Squares) Regression
ols_model <- lm(y ~ X - 1) # '-1' to exclude the intercept

# Extract OLS coefficients
ols_coefs <- as.vector(coef(ols_model))

# Create a grid of coefficient values
grid_size <- 100
beta1_vals <- seq(-0.5, 3.5, length.out = grid_size)
beta2_vals <- seq(-2.5, 0.5, length.out = grid_size)

# Compute the OLS cost function
ols_cost <- outer(beta1_vals, beta2_vals, Vectorize(function(b1, b2) {
  sum((y - X %*% c(b1, b2))^2) / (2 * n)
}))

# Convert to data frame for ggplot
cost_df <- expand.grid(beta1 = beta1_vals, beta2 = beta2_vals)
cost_df$cost <- as.vector(ols_cost)

mstop_max <- 10000

# Fit mboost models and extract coefficient changes
lambdas <- c(0,10^seq(0,4,l=47))
nu_vals <- c(0.001, 0.0025, 0.005, 0.0075, 0.01, 0.1)

res_nu <- list()

if(!file.exists("coefpaths_mstop10000_equalLambda.RDS")){
  
  for(nu in nu_vals){
    
    res <- mclapply(lambdas, function(l1) {
      
      mboost_coefs_changes <- data.frame(
        Iteration = integer(),
        Lambda1 = double(),
        Lambda2 = double(),
        Start_Beta1 = double(),
        Start_Beta2 = double(),
        End_Beta1 = double(),
        End_Beta2 = double()
      )
      
      model <- mboost(y ~ bols(X1, lambda = l1, intercept = FALSE) %+% 
                        bols(X2, lambda = l1, intercept = FALSE),
                      data = data.frame(y=y, X1=X[,1], X2=X[,2]), 
                      offset = 0, control = boost_control(mstop = mstop_max, nu = nu))
      coefs <- t(coef(model, aggregate = "cumsum")[[1]])
      for (iter in 2:mstop_max) {
        mboost_coefs_changes <- rbind(mboost_coefs_changes, data.frame(
          Iteration = iter,
          Lambda1 = l1,
          Lambda2 = l1,
          Start_Beta1 = coefs[iter - 1, 1],
          Start_Beta2 = coefs[iter - 1, 2],
          End_Beta1 = coefs[iter, 1],
          End_Beta2 = coefs[iter, 2]
        ))
      }
      
      return(mboost_coefs_changes)
      
    }, mc.cores = 12)
    
    ress <- do.call("rbind", res)
    ress$nu <- nu
    res_nu <- c(res_nu, list(ress))
    
  }
  
  mboost_coefs_changes <- do.call("rbind", res_nu)
  
  mboost_coefs_changes <- mboost_coefs_changes %>% 
    mutate(lambdas = Lambda1) # paste(Lambda1, Lambda2, sep = "-"))
  
  # Generate Ridge Path for different lambda values
  ridge_lambdas <- c(10^seq(-5,-1,by=1),seq(2,9,by=1)/10, 10^seq(0,4,l=50))
  ridge_coefs_seq <- sapply(ridge_lambdas, function(lambda) {
    as.vector(coef(glmnet(X, y, alpha = 0, lambda = lambda, intercept = FALSE))[-1])
  })
  
  ridge_path <- data.frame(
    Iteration = NA,
    Lambda1 = ridge_lambdas[-1],
    Lambda2 = ridge_lambdas[-1],
    Start_Beta1 = ridge_coefs_seq[1, -1],
    Start_Beta2 = ridge_coefs_seq[2, -1],
    End_Beta1 = ridge_coefs_seq[1, -length(ridge_lambdas)],
    End_Beta2 = ridge_coefs_seq[2, -length(ridge_lambdas)]
  )
  
  mboost_coefs_changes$method <- "BAM"
  ridge_path$nu <- rep(NA, nrow(ridge_path))
  ridge_path$method <- "Ridge"
  
  ridge_path <- ridge_path %>% 
    mutate(lambdas = Lambda1) # paste(Lambda1, Lambda2, sep = "-"))
  
  paths <- rbind(mboost_coefs_changes, ridge_path)
  
  saveRDS(paths, file="coefpaths_mstop10000_equalLambda.RDS")
  
}else{
  
  paths <- readRDS("coefpaths_mstop10000_equalLambda.RDS")
  mboost_coefs_changes <- paths %>% filter(method=="BAM")
  ridge_path <- paths %>% filter(method=="Ridge")
  
}

# Plot
ggg1 <- ggplot() +
  # geom_contour_filled(data = cost_df, aes(x = beta1, y = beta2, z = cost, fill = ..level..), alpha = 0.7) +
  geom_contour(data = cost_df, aes(x = beta1, y = beta2, z = cost), color = "black") +
  geom_segment(data = mboost_coefs_changes %>% filter(nu == nu_vals[[6]]), 
               aes(x = Start_Beta1, y = Start_Beta2, 
                   xend = End_Beta1, yend = End_Beta2, 
                   color = log(lambdas, 10)), alpha = 0.8,
               arrow = arrow(length = unit(0.05, "inches")), linewidth = 0.5) +
  geom_segment(data = ridge_path, aes(x = Start_Beta1, y = Start_Beta2, 
                                      xend = End_Beta1, yend = End_Beta2),
               arrow = arrow(length = unit(0.1, "inches")), color = "blue", 
               linewidth = 0.6) +
  # geom_segment(data = ridge_path, aes(x = Start_Beta1, y = Start_Beta2, 
  #                                     xend = End_Beta1, yend = End_Beta2),
  #              arrow = arrow(length = unit(0.1, "inches")), color = "red", linewidth = 0.8) +
  geom_segment(aes(x = 0, y = 0, xend = ols_coefs[1], yend = ols_coefs[2]), 
               arrow = arrow(length = unit(0.1, "inches")), 
               color = "black", linewidth = 1.4) +
  geom_text(aes(x = ols_coefs[1], y = ols_coefs[2], label = "OLS"), 
            size = 5,
            color = "black", hjust = 1.8, vjust = -0.4) +
  geom_text(aes(x = ridge_coefs[1], y = ridge_coefs[2]),
            label = "Ridge", size = 5,
            color = "blue", hjust = -0.2, vjust = -0.8) +
  labs(x = expression(beta[1]),
       y = expression(beta[2]),
       # fill = "MSE",
       color = expression(lambda)) +
  scale_color_viridis_c(option = "inferno", 
                        trans = "log10", 
                        direction = 1,
                        begin = 0.1,
                        end = 0.9) +
  theme_minimal() + 
  theme(text = element_text(size = 14)) + 
  guides(fill = "none") + xlim(-0.25,3) + ylim(-1.8,0.05) + 
  theme(legend.position = "none")

ggg2 <- ggplot() +
  # geom_contour_filled(data = cost_df, aes(x = beta1, y = beta2, z = cost, fill = ..level..), alpha = 0.7) +
  geom_contour(data = cost_df, aes(x = beta1, y = beta2, z = cost), color = "black") +
  geom_segment(data = mboost_coefs_changes %>% filter(nu == nu_vals[[6]]), 
               aes(x = Start_Beta1, y = Start_Beta2, 
                   xend = End_Beta1, yend = End_Beta2, 
                   color = lambdas), alpha = 0.8,
               arrow = arrow(length = unit(0.05, "inches")), linewidth = 0.5) +
  geom_segment(data = ridge_path, aes(x = Start_Beta1, y = Start_Beta2, 
                                      xend = End_Beta1, yend = End_Beta2),
               arrow = arrow(length = unit(0.1, "inches")), color = "blue", 
               linewidth = 0.6) +
  # geom_segment(data = ridge_path, aes(x = Start_Beta1, y = Start_Beta2, 
  #                                     xend = End_Beta1, yend = End_Beta2),
  #              arrow = arrow(length = unit(0.1, "inches")), color = "red", linewidth = 0.8) +
  geom_segment(aes(x = 0, y = 0, xend = ols_coefs[1], yend = ols_coefs[2]), 
               arrow = arrow(length = unit(0.1, "inches")), 
               color = "black", linewidth = 1.4) +
  geom_text(aes(x = ols_coefs[1], y = ols_coefs[2], label = "OLS"), 
            size = 5,
            color = "black", hjust = 1.8, vjust = -0.4) +
  geom_text(aes(x = ridge_coefs[1], y = ridge_coefs[2]),
            label = "Ridge", size = 5,
            color = "blue", hjust = -0.2, vjust = -0.8) +
  labs(x = expression(beta[1]),
       y = expression(beta[2]),
       # fill = "MSE",
       color = expression(lambda)) +
  scale_color_viridis_c(option = "inferno", 
                        trans = "log10", 
                        direction = 1,
                        begin = 0.1,
                        end = 0.95) +
  theme_minimal() + 
  theme(text = element_text(size = 14)) + 
  guides(fill = "none") + theme(legend.position = "right",
                                legend.key.width = unit(.25, "cm"),
                                legend.text = element_text(size = 9)) + 
  xlim(-0.0,0.25) + ylim(-0.05,0.025)

ggg <- grid.arrange(ggg1, ggg2, ncol=2, widths = c(0.9,0.8))

ggsave(ggg, filename = paste0("ridge_pass_crossing_step_length_0.1.pdf"), 
       width = 14*0.7, height = 5*0.7)