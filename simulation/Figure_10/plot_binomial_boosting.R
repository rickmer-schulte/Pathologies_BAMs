# # Binomial Boosting
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
n <- 200  # number of observations

mean <- c(0, 0)
cov_x <- 0.5
covariance <- matrix(c(1, cov_x, cov_x, 1), nrow = 2)  # Covariance matrix
data <- mvrnorm(n, mu = mean, Sigma = covariance)
scaled.data <- scale(data)
x1 <- scaled.data[, 1]
x2 <- scaled.data[, 2]

# True parameters
beta0 <- 0  # intercept 0.5
beta1 <- 3  # coefficient for x1 #3
beta2 <- -2  # coefficient for x2 #-2

# Simulating response variable y using bernoulli model
pi <- plogis(beta0 + beta1 * x1 + beta2 * x2)
y <- as.numeric(rbinom(n, size = 1, prob = pi))

# Create a data frame
data <- data.frame(y = y, x1, x2)

neg_log_bernoulli <- function(b) -sum(
  dbinom(x = as.numeric(y), size = 1, prob = plogis(cbind(x1,x2)%*%b), log = TRUE))

#### Multiple learning rates

K <- 15000
data_list <-setNames(data.frame(matrix(ncol = 4, nrow = 0)), c("iter","nu","b1", "b2"))
nu_list <- seq(0.02,0.07, by=0.01)
for (nu_iter in nu_list) {
  m_boostbernoulli <- mboost(y ~ bols(x1, intercept = FALSE) + 
                             bols(x2, intercept = FALSE), 
                             offset = 0,
                             family = Binomial(type="glm"), 
                           control = boost_control(mstop = 1, nu = nu_iter),
                           data = data) 
  data_temp <- data.frame("iter"=seq(0,K), "nu"= rep(nu_iter, K+1), 
                          "b1"= c(0,coef(m_boostbernoulli[K], aggregate = "cumsum")$`bols(x1, intercept = FALSE)`), 
                          "b2"=c(0,coef(m_boostbernoulli[K], aggregate = "cumsum")$`bols(x2, intercept = FALSE)`))
  data_list <- rbind(data_list, data_temp)
  print(nu_iter)
}
data_list$nu_fac <- as.factor(data_list$nu)

# loss convergence plots - data
data_list$loss <- apply(data_list[,c(3,4)], 1, neg_log_bernoulli)
# add one to avoid log(0)
data_list$loss_diff <- data_list$loss - neg_log_bernoulli(c(beta1, beta2))

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
  scale_x_continuous(trans = "log10") +
  scale_color_npg(scale_name = "what") + 
  guides(color = guide_legend(title = expression(nu))) + 
  xlab("Iteration") + 
  ylab(expression(beta[1])) 

g2

gg_binom_boost <- grid.arrange(g1,g2, ncol=2)

ggsave(gg_binom_boost, filename = paste0("convergence_bernoulli.pdf"), 
       width = 14*0.7, height = 5*0.7)
