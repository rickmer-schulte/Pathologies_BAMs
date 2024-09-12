##### Ridge Boosting
library(ggplot2)
library(reshape2)
library(MASS)
library(mboost)
library(glmnet)
library(mnormt)

# Step 1: Generate Data with Bivariate Gaussian Distribution
set.seed(123)
n <- 200
mean <- c(0, 0)
covariance <- matrix(c(1, -0.6, -0.6, 1), nrow = 2)  # Covariance matrix
data <- mvrnorm(n, mu = mean, Sigma = covariance)
cov_matrix <- cov(data) # emp. covariance matrix
x1 <- data[, 1]
x2 <- data[, 2]
beta0 <- 0
beta1 <- 5 
beta2 <- 3
sigma <- 1
y <- beta0 + beta1 * x1 + beta2 * x2 + rnorm(n, sd = sigma)

# Set up boosting routines (LS Boosting and Ridge Boosting)
nu_val <- 0.3
lamda_param <- 100
ls_boost <- mboost(y ~ bols(x1, intercept = FALSE) + bols(x2, intercept = FALSE), 
             control = boost_control(mstop = 1, nu=nu_val))
ridge_boost <- mboost(y ~ bols(x1, intercept = FALSE, lambda = lamda_param) +
                     bols(x2, intercept = FALSE, lambda = lamda_param), 
                   control = boost_control(mstop = 1, nu=nu_val))

# Linear regression
reg_lm <- lm(y ~ -1 + x1 + x2)
est_coef_lm <- reg_lm$coef
# ridge regression
X <- cbind(x1,x2)
beta_ridge <- solve(t(X) %*% X + lamda_param * diag(ncol(X)), t(X) %*% y) 

# Prepare line data for mboost and ridge
K <- 1000
line_data <- data.frame("iter"=seq(0,K), 
                          "coef1" = c(c(0,coef(ls_boost[K], aggregate = "cumsum")$`bols(x1, intercept = FALSE)`), 
                          c(0,coef(ridge_boost[K], aggregate = "cumsum")$`bols(x1, intercept = FALSE, lambda = lamda_param)`)),
                          "coef2" = c(c(0,coef(ls_boost[K], aggregate = "cumsum")$`bols(x2, intercept = FALSE)`),
                          c(0,coef(ridge_boost[K], aggregate = "cumsum")$`bols(x2, intercept = FALSE, lambda = lamda_param)`)),
                          "type" = rep(c("LS Boosting", "Ridge Boosting"), each = K + 1))

b1     <- seq(-0.5, 6 , 0.1) 
b2     <- seq(-0.5, 4, 0.1)
f_loss <- function(b) sum((y - (b[1] * x1 + b[2] * x2))**2)

# Prepare contour data
grid_data <- expand.grid(b1 = b1, b2 = b2)
grid_data$z <- apply(grid_data, 1, f_loss)

# Prepare point data
point_data <- data.frame(x = c(coef(reg_lm)[1], beta_ridge[1]),
                         y = c(coef(reg_lm)[2], beta_ridge[2]),
                         type = c("LS Boosting", "Ridge Boosting"))

# Plot LS and Ridge Boosting
ggplot() +
  geom_contour(data = grid_data, aes(x = b1, y = b2, z = z), bins = 18, colour = "grey60") +
  geom_path(data = line_data, aes(x = coef1, y = coef2, colour = type), linewidth = 1.1) +
  geom_point(data = point_data, aes(x = x, y = y, colour = type), shape = 4, size = 3) +
  scale_colour_manual(values = c("blue", "red")) +
  labs(x = expression(beta[1]), y = expression(beta[2])) +
  theme_minimal() + 
  theme(text = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.title = element_blank(),
        legend.position = c(0.175, 0.89),
        legend.text = element_text(size = 12),
        legend.background = element_rect(fill = "white", color = NA)
        ) + 
  guides(colour = guide_legend(override.aes = list(shape = NA)))

ggsave(file = "plot_ridge.pdf", width = 5, height = 3)
