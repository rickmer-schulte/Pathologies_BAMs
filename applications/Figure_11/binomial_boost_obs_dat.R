# Binomial Boosting on the diabetes dataset
library(mlbench)
library(mboost)
library(dplyr)
library(ggplot2)
library(ggsci)
library(gridExtra)

# Load dataset (using all columns)
data("PimaIndiansDiabetes2", package = "mlbench")
diabetes_df <- PimaIndiansDiabetes2

# Scale all predictors to have unit variance (leaving the outcome variable 'diabetes' unchanged)
predictor_cols <- setdiff(names(diabetes_df), "diabetes")
diabetes_df[predictor_cols] <- scale(diabetes_df[predictor_cols])

# Settings
K <- 10000
nu_list <- seq(0.02, 0.07, by = 0.01)
results_list <- list()  # pre-allocate list for results

# Loop over learning rates
for (nu_iter in nu_list) {
  # Fit the model with all predictors
  m_boost <- glmboost(diabetes ~ ., 
                      data = diabetes_df,
                      control = boost_control(mstop = K, nu = nu_iter), 
                      family = Binomial())
  
  # Extract cumulative coefficients at final iteration (mstop = K)
  coef_cumsum <- coef(m_boost[K], aggregate = "cumsum")
  
  # Loss for each iteration from 0 to K
  loss_values <- sapply(0:K, function(i) -logLik(m_boost[i]))
  
  # Start with a data frame that contains the iteration, learning rate, and loss
  data_temp <- data.frame(iter = 0:K,
                          nu   = rep(nu_iter, K + 1),
                          loss = loss_values)
  
  # Dynamically add coefficient columns:
  # For each predictor, prepend a zero to align with iteration 0.
  for (var in names(coef_cumsum)) {
    data_temp[[var]] <- c(0, coef_cumsum[[var]])
  }
  
  # Store result
  results_list[[length(results_list) + 1]] <- data_temp
  print(nu_iter)
}

# Combine all iterations into one data frame
data_list_diabetes <- dplyr::bind_rows(results_list)

# (Optional) create additional columns for plotting
data_list_diabetes$nu_fac <- as.factor(data_list_diabetes$nu)
data_list_diabetes$loss_diff <- data_list_diabetes$loss - min(data_list_diabetes$loss)

# Custom formatting functions
integer_format <- function(x) round(x)
dec1_format <- function(x) round(x, 1)

# Loss convergence plot
g1 <- ggplot(data=data_list_diabetes, aes(x=iter, y=log(loss_diff+1),
                                          group=nu, color=as.factor(nu))) +
  geom_line(size = 1.5, alpha = 0.8) +
  theme_bw() +
  theme(legend.position = "bottom",
        text = element_text(size = 14),
        legend.box.background = element_blank(),
        legend.background = element_blank()) +
  scale_y_continuous(labels = integer_format) +
  scale_x_continuous(trans = "log10") +
  scale_color_npg() +
  guides(color = guide_legend(title = expression(nu))) +
  xlab("Iteration (log-scale)") +
  ylab("Loss (Diabetes)")
g1

g2 <- ggplot(data=data_list_diabetes, aes(x=iter, y=mass,
                                          group=nu, color=as.factor(nu))) +
  geom_line(size = 1.5, alpha = 0.8) + 
  theme_bw() +
  theme(legend.position = "bottom",
        text = element_text(size = 14),
        legend.box.background = element_blank(),
        legend.background = element_blank(),
  ) +
  scale_x_continuous(trans = "log10") +
  scale_color_npg(scale_name = "what") + 
  guides(color = guide_legend(title = expression(nu))) + 
  ylab(expression(beta[glucose])) +
  xlab("Iteration (log-scale)")
g2

# Save plots
gg_diabetes <- grid.arrange(g1, g2, ncol=2)
ggsave(gg_diabetes, filename = paste0("convergence_binom_diabetes.pdf"), 
       width = 14*0.7, height = 5*0.7)
