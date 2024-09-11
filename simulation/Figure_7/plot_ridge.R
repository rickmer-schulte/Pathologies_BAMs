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
covariance <- matrix(5*c(1, 0.5, 0.5, 1), nrow = 2)  # Covariance matrix
data <- mvrnorm(n, mu = mean, Sigma = covariance)
cov_matrix <- cov(data) # emp. covariance matrix
x1 <- data[, 1]
x2 <- data[, 2]
beta0 <- 0
beta1 <- 8 # 8
beta2 <- 6 # 4
sigma <- 1
y <- beta0 + beta1 * x1 + beta2 * x2 + rnorm(n, sd = sigma)

m2 <- mboost(y ~ bols(x1, intercept = FALSE) + bols(x2, intercept = FALSE), 
             control = boost_control(mstop = 1))
lamda_param <- 5000
m2_ridge <- mboost(y ~ bols(x1, intercept = FALSE, lambda = lamda_param) +
                     bols(x2, intercept = FALSE, lambda = lamda_param), 
                   control = boost_control(mstop = 1))

# Linear regression
reg_lm <- lm(y ~ -1 + x1 + x2)
est_coef_lm <- reg_lm$coef
# ridge regression
X <- cbind(x1,x2)
beta_ridge <- solve(t(X) %*% X + lamda_param, t(X) %*% y) 

# Prepare line data for mboost and ridge
K <- 1000
line_data <- data.frame("iter"=seq(0,K), 
                          "coef1" = c(c(0,coef(m2[K], aggregate = "cumsum")$`bols(x1, intercept = FALSE)`), 
                          c(0,coef(m2_ridge[K], aggregate = "cumsum")$`bols(x1, intercept = FALSE, lambda = lamda_param)`)),
                          "coef2" = c(c(0,coef(m2[K], aggregate = "cumsum")$`bols(x2, intercept = FALSE)`),
                          c(0,coef(m2_ridge[K], aggregate = "cumsum")$`bols(x2, intercept = FALSE, lambda = lamda_param)`)),
                          "type" = rep(c("LS Boosting", "Ridge Boosting"), each = K + 1))

b1     <- seq(-1, 10, 0.01) 
b2     <- seq(-1, 10, 0.01)
f_mnorm     <- function(b1, b2) dmnorm(cbind(b1, b2), est_coef_lm, cov_matrix)

# Prepare contour data
grid_data <- expand.grid(b1 = b1, b2 = b2)
grid_data$z <- f_mnorm(grid_data$b1, grid_data$b2)

# Prepare point data
point_data <- data.frame(x = c(coef(reg_lm)[1], beta_ridge[1]),
                         y = c(coef(reg_lm)[2], beta_ridge[2]),
                         type = c("LS Boosting", "Ridge Boosting"))

ggplot() +
  geom_contour(data = grid_data, aes(x = b1, y = b2, z = z), bins = 16, colour = "grey60") +
  geom_path(data = line_data, aes(x = coef1, y = coef2, colour = type), linewidth = 1.1) +
  geom_point(data = point_data, aes(x = x, y = y, colour = type), shape = 4, size = 3) +
  scale_colour_manual(values = c("blue", "red")) +
  labs(x = expression(beta[1]), y = expression(beta[2])) +
  theme_minimal() + 
  theme(text = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.title = element_blank(),
        legend.position = c(0.2, 0.9),
        legend.text = element_text(size = 12)
        ) + 
  guides(colour = guide_legend(override.aes = list(shape = NA)))

ggsave(file = "plot_ridge.pdf", width = 5, height = 3)

