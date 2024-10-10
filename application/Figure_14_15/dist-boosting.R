### Distributional Boosting
library(ggplot2)
library(ggsci)
library(gridExtra)
library(QRegVCM)
library(gamboostLSS)
library(dplyr)
library(scales)

# Data Loading
data("CD4")
df_cd4 <- CD4 %>% dplyr::select(CD4, Time, Smooking, Age, PreCD4, CD4)

### GAMLSS-Boosting with multiple learning rates (CD4)
K <- 10
data_list_cd4 <-setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("iter","nu","loss"))
nu_list <- c(seq(0.3,0.8, by=0.1))

for (nu_iter in nu_list) {
  glmss_boost_cd4 <- glmboostLSS(CD4 ~ -1 + ., 
                             families = GaussianLSS(), 
                             control = boost_control(mstop = K, nu = nu_iter), 
                             data = df_cd4)
  df_coef <- coef(glmss_boost_cd4[K], aggregate = "cumsum")
  data_temp <- data.frame("iter"=seq(0,K), "nu"= rep(round(nu_iter,2), K+1), 
                          "coef_mu"= c(0, df_coef$mu$PreCD4), 
                          "coef_sigma"= c(0, df_coef$sigma$PreCD4), 
                          "loss_mu"= c(risk(glmss_boost_cd4)$mu),
                          "loss_sigma"= c(risk(glmss_boost_cd4)$sigma))
  data_list_cd4 <- rbind(data_list_cd4, data_temp)
  print(nu_iter)
}
data_list_cd4$loss_diff_mu <- data_list_cd4$loss_mu - min(data_list_cd4$loss_mu)
data_list_cd4$loss_diff_sigma <- data_list_cd4$loss_sigma - min(data_list_cd4$loss_sigma)

# parameter diver/convergence plots
gg_cd4_coef <- ggplot(data=data_list_cd4, aes(x=iter, y=coef_sigma, 
                                                  group=nu, color=as.factor(nu))) + 
  geom_line(size = 1.5, alpha = 0.8) + 
  theme_bw() +
  theme(legend.position = "bottom",
        text = element_text(size = 14),
        legend.box.background = element_blank(),
        legend.background = element_blank(),
        legend.margin = margin(t = -10)
  ) +
  scale_x_continuous(breaks = pretty_breaks(),
                     labels = number_format(accuracy = 1)) +
  scale_color_npg(scale_name = "what") + 
  guides(color = guide_legend(title = expression(nu))) + 
  xlab("Iteration") + 
  ylab(expression(beta[PreCD4] ~ "(Variance Model)"))
gg_cd4_coef

# loss diver/convergence plots
gg_cd4_loss <- ggplot(data=data_list_cd4, aes(x=iter, y=loss_diff_sigma, 
                                                  group=nu, color=as.factor(nu))) + 
  geom_line(size = 1.5, alpha = 0.8) + 
  theme_bw() +
  theme(legend.position = "bottom",
        text = element_text(size = 14),
        legend.box.background = element_blank(),
        legend.background = element_blank(),
        legend.margin = margin(t = -10)
  ) +
  scale_x_continuous(breaks = pretty_breaks(),
                     labels = number_format(accuracy = 1)) +
  scale_color_npg(scale_name = "what") + 
  guides(color = guide_legend(title = expression(nu))) + 
  xlab("Iteration") + 
  ylab("Loss (Variance Model)")
gg_cd4_loss

gg_cd4_boost <- grid.arrange(gg_cd4_loss, gg_cd4_coef, ncol=2)

ggsave(gg_cd4_boost, filename = paste0("cd4_distboost_conv.pdf"), 
       width = 14*0.7, height = 5*0.7)

### GAMLSS-Boosting with multiple learning rates (Engel data)
data(engel)
K <- 10
data_list_engel <-setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("iter","nu","loss"))
nu_list <- c(seq(0.2,0.6, by=0.1), 0.6107)

for (nu_iter in nu_list) {
  glmss_boost <- glmboostLSS(foodexp ~ -1 + income, 
                                         families = GaussianLSS(), 
                                         control = boost_control(mstop = K, nu = nu_iter), 
                                         data = engel)
  df_coef <- coef(glmss_boost[K], aggregate = "cumsum")
  data_temp <- data.frame("iter"=seq(0,K), "nu"= rep(round(nu_iter,2), K+1), 
                          "coef_mu"= c(0, df_coef$mu$income), 
                          "coef_sigma"= c(0, df_coef$sigma$income), 
                          "loss_mu"= c(risk(glmss_boost)$mu),
                          "loss_sigma"= c(risk(glmss_boost)$sigma))
  data_list_engel <- rbind(data_list_engel, data_temp)
  print(nu_iter)
}
data_list_engel$loss_mu <- data_list_engel$loss_mu - min(data_list_engel$loss_mu)
data_list_engel$loss_diff_sigma <- data_list_engel$loss_sigma - min(data_list_engel$loss_sigma)

# parameter diver/convergence plots
gg_engel_coef <- ggplot(data=data_list_engel, aes(x=iter, y=coef_sigma, 
                                      group=nu, color=as.factor(nu))) + 
  geom_line(size = 1.5, alpha = 0.8) + 
  theme_bw() +
  theme(legend.position = "bottom",
        text = element_text(size = 14),
        legend.box.background = element_blank(),
        legend.background = element_blank(),
        legend.margin = margin(t = -10)
  ) +
  scale_x_continuous(breaks = pretty_breaks(),
                     labels = number_format(accuracy = 1)) +
  scale_color_npg(scale_name = "what") + 
  guides(color = guide_legend(title = expression(nu))) + 
  xlab("Iteration") + 
  ylab(expression(beta[income] ~ "(Variance Model)"))
gg_engel_coef

# loss diver/convergence plots
gg_engel_loss <- ggplot(data=data_list_engel, aes(x=iter, y=loss_diff_sigma, 
                                       group=nu, color=as.factor(nu))) + 
  geom_line(size = 1.5, alpha = 0.8) + 
  theme_bw() +
  theme(legend.position = "bottom",
        text = element_text(size = 14),
        legend.box.background = element_blank(),
        legend.background = element_blank(),
        legend.margin = margin(t = -10)
  ) +
  scale_x_continuous(breaks = pretty_breaks(),
                     labels = number_format(accuracy = 1)) +
  scale_color_npg(scale_name = "what") + 
  guides(color = guide_legend(title = expression(nu))) + 
  xlab("Iteration") + 
  ylab("Loss (Variance Model)")
gg_engel_loss

gg_engel_boost <- grid.arrange(gg_engel_loss, gg_engel_coef, ncol=2)

ggsave(gg_engel_boost, filename = paste0("engel_distboost_conv.pdf"), 
       width = 14*0.7, height = 5*0.7)
