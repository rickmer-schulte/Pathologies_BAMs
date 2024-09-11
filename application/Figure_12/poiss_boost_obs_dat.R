### Poisson Boosting
library(ggplot2)
library(gridExtra)
library(splines)
library(mboost)
library(dplyr)
library(ggsci)
library(here)

# Data Reading and Preparation
source("../data/subset_dat.R")
path_cov_dat <- here::here('application/data') # Insert folder of loaded data here
df_san_fra <- load_covid_dat(data_folder_path=path_cov_dat, sf=TRUE)

# Select subset of variables
df_sf_pois <- df_san_fra %>% select(c(new_confirmed, date, temp, humid))
df_sf_pois$temp <- df_sf_pois$temp/10
df_sf_pois$humid <- df_sf_pois$humid/10

#### Multiple learning rates

K <- 1000
data_list_covid <-setNames(data.frame(matrix(ncol = 6, nrow = 0)), c("iter","nu","b1", "b2", "b3", "loss"))
nu_list <- seq(0.02,0.06, by=0.01)
for (nu_iter in nu_list) {
  m_boostpoisson <- glmboost(new_confirmed ~  -1 + ., offset=FALSE,
                             control = boost_control(mstop = 1000, nu = nu_iter), family=Poisson(), data = df_sf_pois)
  data_temp <- data.frame("iter"=seq(0,K), "nu"= rep(nu_iter, K+1), 
                          "b1"= c(0,coef(m_boostpoisson[K], aggregate = "cumsum")$temp), 
                          "b2"= c(0,coef(m_boostpoisson[K], aggregate = "cumsum")$date),
                          "b3"=c(0,coef(m_boostpoisson[K], aggregate = "cumsum")$humid),
                          "loss"=sapply(0:1000, function(i) -logLik(m_boostpoisson[i])))
  data_list_covid <- rbind(data_list_covid, data_temp)
  print(nu_iter)
}
data_list_covid$nu_fac <- as.factor(data_list_covid$nu)

data_list_covid$loss_diff <- data_list_covid$loss - min(data_list_covid$loss)

integer_format <- function(x) round(x)

# loss convergence plots
g1 <- ggplot(data=data_list_covid, aes(x=iter, y=log(loss_diff+1), 
                                       group=nu, color=as.factor(nu))) + 
  geom_line(size = 1.5, alpha = 0.8) + 
  theme_bw() +
  theme(legend.position = "none",
        text = element_text(size = 14),
        legend.box.background = element_blank(),
        legend.background = element_blank(),
  ) +
  scale_y_continuous(labels = integer_format) +
  scale_x_continuous(trans = "log10") +
  scale_color_npg(scale_name = "what") + 
  guides(color = guide_legend(title = expression(nu))) + 
  xlab("") + 
  ylab("Loss (Covid-19)")
g1

dec1_format <- function(x) round(x,1)

g2 <- ggplot(data=data_list_covid, aes(x=iter, y=b1, 
                                       group=nu, color=as.factor(nu))) + 
  geom_line(size = 1.5, alpha = 0.8) + 
  theme_bw() +
  theme(legend.position = "right",
        text = element_text(size = 14),
        legend.box.background = element_blank(),
        legend.background = element_blank(),
  ) +
  scale_y_continuous(labels = dec1_format) +
  scale_color_npg(scale_name = "what") + 
  guides(color = guide_legend(title = expression(nu))) + 
  xlab("") + 
  ylab(expression(beta[temperature])) + 
  scale_x_continuous(labels = integer_format, limits = c(990,1000)) 

g2

gg <- grid.arrange(g1,g2, ncol=2)

#### Health Score Poisson Boosting
library(fairml)
df_health <- get("health.retirement", envir=getNamespace("fairml"))
df_health <- df_health %>% select(c(score, cognition_catnew, age, year))

#### Multiple learning rates

K <- 1000
data_list_health <-setNames(data.frame(matrix(ncol = 6, nrow = 0)), c("iter","nu","b1", "b2", "b3", "loss"))
nu_list <- seq(0.1,0.7, by=0.1)
nu_list <- c(0.2, 0.3, 0.4, 0.5, 0.65, 0.7)
for (nu_iter in nu_list) {
  glm_poisson_health_boost <- glmboost(score ~ -1 +., offset=FALSE,
                                       control = boost_control(mstop = 1000, nu = nu_iter), 
                                       family=Poisson(), data = df_health)
  data_temp <- data.frame("iter"=seq(0,K), "nu"= rep(nu_iter, K+1), 
                          "b1"= c(0,coef(glm_poisson_health_boost[K], aggregate = "cumsum")$cognition_catnew), 
                          "b2"= c(0,coef(glm_poisson_health_boost[K], aggregate = "cumsum")$age),
                          "b3"=c(0,coef(glm_poisson_health_boost[K], aggregate = "cumsum")$year),
                          "loss"=sapply(0:1000, function(i) -logLik(glm_poisson_health_boost[i])))
  data_list_health <- rbind(data_list_health, data_temp)
  print(nu_iter)
}
data_list_health$nu_fac <- as.factor(data_list_health$nu)

data_list_health$loss_diff <- data_list_health$loss - min(data_list_health$loss)

integer_format <- function(x) round(x)

# loss convergence plots
g3 <- ggplot(data=data_list_health, aes(x=iter, y=log(loss_diff+1), 
                                        group=nu, color=as.factor(nu))) + 
  geom_line(size = 1.5, alpha = 0.8) + 
  theme_bw() +
  theme(legend.position = "none",
        text = element_text(size = 14),
        legend.box.background = element_blank(),
        legend.background = element_blank(),
  ) +
  scale_y_continuous(labels = integer_format) +
  scale_x_continuous(trans = "log10") +
  scale_color_npg(scale_name = "what") + 
  guides(color = guide_legend(title = expression(nu))) + 
  xlab("Iteration (log-scale)") + 
  ylab("Loss (HRS)")
g3

dec1_format <- function(x) round(x,1)
g4 <- ggplot(data=data_list_health, aes(x=iter, y=b1, 
                                        group=nu, color=as.factor(nu))) + 
  geom_line(size = 1.5, alpha = 0.8) + 
  theme_bw() +
  theme(legend.position = "right",
        text = element_text(size = 14),
        legend.box.background = element_blank(),
        legend.background = element_blank(),
  ) +
  scale_y_continuous(labels = dec1_format) +
  scale_color_npg(scale_name = "what") + 
  guides(color = guide_legend(title = expression(nu))) + 
  xlab("Iteration") + 
  ylab(expression(beta[cognition])) + 
  scale_x_continuous(labels = integer_format, limits = c(990,1000)) 

g4

gg_34 <- grid.arrange(g3,g4, ncol=2)

# Both plot together
gg_poiss <- grid.arrange(g1,g2,g3,g4, ncol=2, heights = c(1, 1))
grid.arrange(
  arrangeGrob(g1, g2, ncol = 1, heights = c(1, 1)),  # Upper plots
  arrangeGrob(g3, g4, ncol = 1, heights = c(1, 1)),  # Lower plots
  ncol = 2,
  heights = c(1, 0.05, 1)  # Adjust the spacing between upper and lower plots
)
ggsave(gg_poiss, file="poisboost_plot.pdf", width = 10*0.9, height = 5*0.9)
