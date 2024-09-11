# Poisson Boosting
library(ggplot2)
library(reshape2)
library(MASS)
library(mboost)
library(glmnet)
library(mnormt)
library(glmtoolbox)
library(RColorBrewer)
library(ggsci)
library(gridExtra)
library(scales)

# Step 1: Simulate Data
set.seed(123) 
n <- 100  # number of observations

mean <- c(0, 0)
cov_x <- 0.5
covariance <- matrix(c(1, cov_x, cov_x, 1), nrow = 2)  # Covariance matrix
data <- mvrnorm(n, mu = mean, Sigma = covariance)
scaled.data <- scale(data)
x1 <- scaled.data[, 1]
x2 <- scaled.data[, 2]

# True parameters
beta0 <- 0  # intercept 0.5
beta1 <- 3  # coefficient for x1
beta2 <- -2  # coefficient for x2

# Simulating response variable y using Poisson model
lambda <- exp(beta0 + beta1 * x1 + beta2 * x2)
y <- rpois(n, lambda)

# Create a data frame
data <- data.frame(y, x1, x2)

# Step 2: Fit Poisson Regression Model
model_poisson <- glm(y ~ -1 + x1 + x2, data=data, family="poisson")
summary(model_poisson)

FS_poisson <- data.frame(FisherScoring(model_poisson)[,c(1,4,5)])
FS_poisson_ext <- rbind(c(0,0,0), FS_poisson)

# mboost
m_boostpoisson <- mboost(y ~ bols(x1, intercept = FALSE) + 
                           bols(x2, intercept = FALSE), family = Poisson(), 
                         offset=FALSE, 
                         control = boost_control(mstop = 1, nu = 0.01))

K <- 10000

# Prepare line data
line_data <- data.frame(
  "iter"= seq(0,K), 
  "coef1"= c(0,coef(m_boostpoisson[K], aggregate = "cumsum")$`bols(x1, intercept = FALSE)`), 
  "coef2"=c(0,coef(m_boostpoisson[K], aggregate = "cumsum")$`bols(x2, intercept = FALSE)`))

# Prepare contour data
b1_sec   <- seq(-1, 3.2, 0.01) 
b2_sec   <- seq(-2.5, 1, 0.01)
f_poisson     <- function(b) sum((b[1] * x1 + b[2] * x2) * y - exp(b[1] * x1 + b[2] * x2))

grid_data <- expand.grid(b1 = b1_sec, b2 = b2_sec)
grid_data$z <- apply(grid_data, 1, f_poisson)

# Prepare point data
point_data <- data.frame(x = c(beta1, coef(model_poisson)[1]),
                         y = c(beta2, coef(model_poisson)[2]),
                         type = c("True", "Fisher Scoring"))

ggplot(grid_data, aes(x = b1, y = b2, z = z)) + geom_contour() + geom_contour_filled(bins=20, show.legend = FALSE) 
gg_poisson <- ggplot() +
  geom_contour_filled(data = grid_data, aes(x = b1, y = b2, z = z), colour = "grey60", bins=50, show.legend = FALSE) +
  geom_path(data = line_data, aes(x = coef1, y = coef2), linewidth = 1.1) +
  geom_path(data = FS_poisson_ext, aes(x = x1, y = x2), linewidth = 1.1, color="blue") +
  geom_point(data = FS_poisson, aes(x = x1, y = x2), color="blue") +
  geom_point(data = point_data, aes(x = x, y = y, colour = type), shape = 4, size = 3) +
  scale_colour_manual(values = c("blue", "red")) +
  ylim(-2.5,1) +
  xlim(-1,3.2) +
  labs(x = expression(x[1]), y = expression(x[2])) +
  theme_minimal() + 
  theme(text = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.title = element_blank(),
        legend.position = c(0.2, 0.9),
        legend.text = element_text(size = 12)
        # legend.position = "bottom"
  ) 
gg_poisson

neg_log_poisson <- function(b) -sum((b[1] * x1 + b[2] * x2) * y - exp(b[1] * x1 + b[2] * x2))
line_data$loss <- apply(line_data[,c(2,3)], 1, neg_log_poisson)
ggplot(data=line_data, aes(x=iter,y=loss)) + geom_line() + 
  geom_hline(yintercept=neg_log_poisson(c(beta1, beta2)), linetype="dashed", color = "red")

#### Multiple learning rates

K <- 1000
data_list <-setNames(data.frame(matrix(ncol = 4, nrow = 0)), c("iter","nu","b1", "b2"))
nu_list <- seq(0.02,0.07, by=0.01)
for (nu_iter in nu_list) {
  m_boostpoisson <- mboost(y ~ bols(x1, intercept = FALSE) + 
                             bols(x2, intercept = FALSE), family = Poisson(), 
                           offset=FALSE, control = boost_control(mstop = 1, nu = nu_iter)) 
  data_temp <- data.frame("iter"=seq(0,K), "nu"= rep(nu_iter, K+1), 
                          "b1"= c(0,coef(m_boostpoisson[K], aggregate = "cumsum")$`bols(x1, intercept = FALSE)`), 
                          "b2"=c(0,coef(m_boostpoisson[K], aggregate = "cumsum")$`bols(x2, intercept = FALSE)`))
  data_list <- rbind(data_list, data_temp)
  print(nu_iter)
}
data_list$nu_fac <- as.factor(data_list$nu)

gg_poisson <- ggplot() +
  geom_contour_filled(data = grid_data, aes(x = b1, y = b2, z = z), 
                      colour = "grey60", bins=50, show.legend = FALSE) +
  geom_path(data = data_list, aes(x = b1, y = b2, group = nu_fac, 
                                  color = as.factor(nu_fac)), linewidth = 0.5, alpha=0.5) +
  geom_path(data = FS_poisson_ext, aes(x = x1, y = x2), linewidth = 1.1, 
            color="blue") +
  geom_point(data = FS_poisson, aes(x = x1, y = x2), color="blue") +
  geom_point(data = point_data, aes(x = x, y = y, colour = type), 
             shape = 4, size = 3) +
  labs(x = expression(x[1]), y = expression(x[2])) +
  theme_minimal() +
  theme(text = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.title = element_blank(),
        legend.position = c(0.2, 0.9),
        legend.text = element_text(size = 12)
  ) 
gg_poisson

# loss convergence plots - data
neg_log_poisson <- function(b) -sum((b[1] * x1 + b[2] * x2) * y - exp(b[1] * x1 + b[2] * x2))
data_list$loss <- apply(data_list[,c(3,4)], 1, neg_log_poisson)
# add one before applying logarithm later to avoid log(0)
data_list$loss_diff <- data_list$loss - neg_log_poisson(c(beta1, beta2))

integer_format <- function(x) round(x)

# loss convergence plots
g1 <- ggplot(data=data_list, aes(x=iter, y=log(loss_diff-min(loss_diff)+1), 
                           group=nu, color=as.factor(nu))) + 
  geom_line(size = 1.5, alpha = 0.8) + 
  theme_bw() +
  theme(legend.position = "bottom",
        text = element_text(size = 14),
        legend.box.background = element_blank(),
        legend.background = element_blank(),
        ) +
  geom_hline(yintercept=0, linetype="dashed", 
             color = "black") + 
  scale_x_continuous(trans = "log10") +
  scale_color_npg(scale_name = "what") + 
  guides(color = guide_legend(title = expression(nu))) + 
  xlab("Iteration (log-scale)") + 
  ylab("Loss")
g1


g2 <- ggplot(data=data_list, aes(x=iter, y=b1, 
                           group=nu, color=as.factor(nu))) + 
  geom_line(size = 1.5, alpha = 0.8) + 
  theme_bw() +
  theme(legend.position = "bottom",
        text = element_text(size = 14),
        legend.box.background = element_blank(),
        legend.background = element_blank(),
  ) +
  geom_hline(yintercept=beta1, linetype="dashed", 
             color = "black") + 
  scale_color_npg(scale_name = "what") + 
  guides(color = guide_legend(title = expression(nu))) + 
  xlab("Iteration") + 
  ylab(expression(beta[1])) + 
  coord_cartesian(ylim=c(2,3.2)) + 
  scale_x_continuous(labels = integer_format, limits = c(990,1000)) 

g2

gg_pois_boost <- grid.arrange(g1,g2, ncol=2)

ggsave(gg_pois_boost, filename = paste0("convergence_poisson.pdf"), 
       width = 14*0.7, height = 5*0.7)
