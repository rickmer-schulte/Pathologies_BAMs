# P-Spline Boosting Covid Data
library(ggplot2)
library(gridExtra)
library(grid)
library(splines)
library(mboost)
library(dplyr)
library(ggsci)
library(here)

# Data Reading and Preparation
source("../data/subset_dat.R")
path_cov_dat <- here::here('application/data') # Insert folder of loaded data here
df_san_fra <- load_covid_dat(data_folder_path=path_cov_dat, sf=TRUE)

# Unpenalized and Penalized fits
spline_x <- bbs(df_san_fra$date, knots = 5)
Z_mat <- extract(spline_x, "design")
Z <- as.matrix(Z_mat)
P_mat <- extract(spline_x, "penalty")
P <- as.matrix(P_mat)
y <- df_san_fra$prevalence

fit = lsfit(Z, y, intercept = F)
beta = fit$coeff
y_hat_unpen = Z %*% beta

lambda = 1
beta_pen = solve(t(Z) %*% Z + lambda * P, t(Z) %*% y) 
y_hat_pen = Z %*% beta_pen

# Adding the fits to the data frame
df_san_fra$y_hat_unpen = y_hat_unpen
df_san_fra$y_hat_pen = y_hat_pen

# ggplot for B-Splines vs. P-Splines
p1 <- ggplot(df_san_fra, aes(x = dates, y = prevalence)) +
  geom_point(shape = 1, size = 1, alpha=0.3) +
  geom_line(aes(y = y_hat_unpen, color = 'B-Spline Solution'), size = 1.1) +
  geom_line(aes(y = y_hat_pen, color = 'P-Spline Solution'), size = 1.1) +
  ggtitle('Analytical Spline Solutions') +
  theme_minimal() +
  theme(text = element_text(size = 14),
        plot.title = element_text(hjust = 0.5, size = 13),
        legend.position = c(0.95, 0.99), # Positioning the legend inside the plot
        legend.justification = c("right", "top"), # Justify the legend position
        legend.box.just = "right",
        legend.margin = margin(0, 0, 0, 0),
        legend.box.background = element_blank(),
        legend.background = element_blank(),
        legend.title = element_blank()) +
  scale_color_manual(values = c('B-Spline Solution' = 'mediumorchid4', 
                                'P-Spline Solution' = 'blue')) + 
  ylim(-0.02,0.46) +
  xlab("") + ylab("Prevalence Covid-19 (SF)")
p1

# Boosting iterations
boost_pspline <- mboost(df_san_fra$prevalence ~ bbs(df_san_fra$date, knots = 5, df=3, lambda = 1), control = boost_control(mstop = 1000))
iter_colors <- c("orchid1", "orchid2", "orchid3", "orchid4", "mediumorchid4")
#iterations <- c(1, 10, 30, 60, 5000)
iterations <- c(1, 10, 20, 50, 2000)

optim_df <- do.call("rbind", lapply(1:length(iterations), function(i){
  data.frame(iteration = iterations[i],
             fitted_values = fitted(boost_pspline[iterations[i]]),
             iter_color = iter_colors[i],
             x = df_san_fra$dates,
             y = df_san_fra$prevalence)
}))

optim_df$text <- paste("iter =", optim_df$iteration)

# ggplot for Boosting iterations
p2 <- ggplot() + 
  geom_point(unique(optim_df[,c("x","y")]), mapping = aes(x = x, y = y), shape = 1, size = 1, alpha=0.3) +
  geom_path(optim_df, mapping = aes(x = x, y = fitted_values, color = factor(iteration)),
            size = 1.1) + 
  scale_color_manual(values = iter_colors, name = "Boosting\nIterations",
                     labels = c("1", "10", "20", "50", "2k")) + 
  theme_minimal() + 
  ggtitle("P-Spline Boosting") + 
  theme(text = element_text(size = 14),
        plot.title = element_text(hjust = 0.5, size = 13),
        legend.title = element_text(size = 12)) + 
  ylim(-0.02,0.46) +
  xlab("") + ylab("")
p2
# Print the plots
print(p1)
print(p2)

# g12 <- grid.arrange(p1,p2, ncol=2)
# ggsave(g12, file="pspline_covid_temporal_plot.pdf", width = 13*0.9, height = 5*0.9)

# Function to extract legend
get_legend <- function(my_plot) {
  tmp <- ggplotGrob(my_plot)
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# Extract legend from p2
legend_p2 <- get_legend(p2)
p2_no_legend <- p2 + theme(legend.position = "none")

# Arrange plots and save
g12 <- grid.arrange(
  arrangeGrob(p1, p2_no_legend, ncol = 2),  
  legend_p2,
  ncol = 2, 
  widths = c(2.8, 0.2)
)
ggsave(g12, file="pspline_covid_temporal_plot.pdf", width = 15*0.9, height = 4*0.9)


