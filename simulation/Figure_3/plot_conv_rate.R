### Linear Convergence Rate Visualizations 
### Convergence Rates of Component-wise L2 Boosting with different condition numbers

library(ggplot2)
library(gridExtra)
library(splines)
library(mboost)
library(MASS)
library(reshape2)
library(latex2exp)

n <- 200
D <- 600
eigval_list <- data.frame("iter"=seq(1,D))
gamma_list <- data.frame("iter"=seq(1,D))
rho_cov_list <- seq(0,0.9, by=0.1)
pb <- txtProgressBar(min = 0, max = length(rho_cov_list), style = 3)
iter <- 0

for (rho in rho_cov_list) {
  iter <- iter + 1
  min_eigen_vec <- c()
  gamma_vec <- c()
  for (d in 1:D) {
    mean_vec <- rep(0, d)
    cov_mat <- matrix(rho, d, d) + diag(1-rho, d)
    Xsim <- mvrnorm(n, mu = mean_vec, Sigma = cov_mat)
    scaled.Xsim <- scale(Xsim)
    Qsim <- t(scaled.Xsim)%*%scaled.Xsim/n
    eig.valQ <- zapsmall(eigen(Qsim)$values)
    min_eigen <- min(eig.valQ[eig.valQ > 0])
    min_eigen_vec <- append(min_eigen_vec, min_eigen)
    gamma_vec <- append(gamma_vec, 1 - min_eigen/d)
  }
  eigval_list <- cbind(eigval_list, min_eigen_vec)
  gamma_list <- cbind(gamma_list, gamma_vec)
  setTxtProgressBar(pb, iter)
}

colnames(eigval_list) <- c("iter", paste0("y", 1:length(rho_cov_list)))
colnames(gamma_list) <- c("iter", paste0("y", 1:length(rho_cov_list)))

# Reshape the data frame to a long format
df_long_eig <- melt(eigval_list, id.vars = "iter")
df_long_gamma <- melt(gamma_list, id.vars = "iter")

# Map color_values to the variable
df_long_eig$color_val <- rep(rho_cov_list, each = nrow(eigval_list))
df_long_gamma$color_val <- rep(rho_cov_list, each = nrow(gamma_list))

p1 <- ggplot(df_long_eig, aes(x = iter, y = value, group = variable, color = color_val)) +
  geom_line() +
  scale_color_gradient(low = "blue", high = "red", name="Rho") +
  theme_minimal() +
  theme(text = element_text(size = 14),
        plot.title = element_text(hjust = 0.5, size = 11),
        legend.position ="none") +
  xlab("p") + ylab(TeX(r'($\lambda_{pmin}(Q))'))
p1

p2 <- ggplot(df_long_gamma, aes(x = iter, y = value, group = variable, color = color_val)) +
  geom_line() +
  scale_y_continuous(limits = c(0.99, 1), trans = "exp") +
  scale_color_gradient(low = "blue", high = "red", name="Correlation") +
  theme_minimal() +
  theme(text = element_text(size = 14),
        plot.title = element_text(hjust = 0.5, size = 11),
        legend.position = c(0.95, 0.65),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(0, 0, 0, 0),
        legend.box.background = element_blank(),
        legend.background = element_blank(),
        ) +
  xlab("p") + ylab(TeX(r'($\gamma)'))
p2

g12 <- grid.arrange(p1,p2, ncol=2)

ggsave(g12, file="conv_rate_plot.pdf", width = 18*0.5, height = 6*0.5)
