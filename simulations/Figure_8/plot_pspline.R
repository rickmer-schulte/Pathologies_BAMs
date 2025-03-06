# Boosting with P-Splines
######
library(ggplot2)
library(gridExtra)
library(splines)
library(mboost)

# data generation 
m = 100
x = seq(0, 1.5, length = m)
set.seed(1234)
y = 1.5 + sin(4 * x) + rnorm(m) * 0.4

# Data frame for ggplot
df = data.frame(x, y)

# Unpenalized and Penalized fits
spline_x <- bbs(x, knots = 5)
Z <- extract(spline_x, "design")
P <- extract(spline_x, "penalty")

fit = lsfit(Z, y, intercept = F)
beta = fit$coeff
y_hat_unpen = Z %*% beta

lambda = 10
beta_pen = solve(t(Z) %*% Z + lambda * P, t(Z) %*% y) 
y_hat_pen = Z %*% beta_pen

# Adding the fits to the data frame
df$y_hat_unpen = y_hat_unpen
df$y_hat_pen = y_hat_pen

# ggplot for B-Splines vs. P-Splines
p1 <- ggplot(df, aes(x = x, y = y)) +
  geom_point(shape = 1, size = 1) +
  geom_line(aes(y = y_hat_unpen, color = 'B-Spline Solution'), size = 1.1) +
  geom_line(aes(y = y_hat_pen, color = 'P-Spline Solution'), size = 1.1) +
  ggtitle('Exact Solutions') +
  theme_minimal() +
  theme(text = element_text(size = 14),
        plot.title = element_text(hjust = 0.5, size = 13),
        legend.position = c(0.95, 0.99),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(0, 0, 0, 0),
        legend.box.background = element_blank(),
        legend.background = element_blank(),
        legend.title = element_blank()) +
  scale_color_manual(values = c('B-Spline Solution' = 'mediumorchid4', 
                                'P-Spline Solution' = 'blue')) + 
  ylim(0.2,2.6) +
  xlab("") + ylab("Fitted Values")

# Boosting iterations
boost_pspline <- mboost(y ~ bbs(x, knots = 5, lambda = 10), control = boost_control(mstop = 1000))
iter_colors <- c("orchid1", "orchid2", "orchid3", "orchid4", "mediumorchid4")
iterations <- c(1, 10, 30, 60, 50000)

optim_df <- do.call("rbind", lapply(1:length(iterations), function(i){
  data.frame(iteration = iterations[i],
             fitted_values = fitted(boost_pspline[iterations[i]]),
             iter_color = iter_colors[i],
             x = df$x,
             y = df$y)
}))

optim_df$text <- paste("iter =", optim_df$iteration)

# ggplot for Boosting iterations
p2 <- ggplot() + 
  geom_point(unique(optim_df[,c("x","y")]), mapping = aes(x = x, y = y), shape = 1, size = 1) +
  geom_path(optim_df, mapping = aes(x = x, y = fitted_values, color = factor(iteration)),
            size = 1.1) + 
  scale_color_manual(values = iter_colors) + 
  theme_minimal() + 
  ggtitle("Boosting Iterations") + 
  theme(legend.position = "none",
        text = element_text(size = 14),
        plot.title = element_text(hjust = 0.5, size = 13)) + 
  ylim(0.2,2.6) +
  xlab("") + ylab("")

# Print the plots
print(p1)
print(p2)

g12 <- grid.arrange(p1,p2, ncol=2)

ggsave(g12, file="pspline_syn_plot.pdf", width = 10, height = 5)
