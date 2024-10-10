# Survival Boosting
library(ggplot2)
library(ggsci)
library(gridExtra)
library(survival)
library(mboost)

#### Boosting with multiple learning rates
K <- 1100
data_list_ovar <-setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("iter","nu","loss"))
nu_list <- seq(0.1,0.6, by=0.1)

for (nu_iter in nu_list) {
  fm_ovar <- Surv(futime,fustat) ~ .
  m_boostpoisson <- glmboost(fm_ovar, data = ovarian, offset=FALSE,
                                  control = boost_control(mstop = K, nu = nu_iter), family = CoxPH())
  data_temp <- data.frame("iter"=seq(0,K), "nu"= rep(nu_iter, K+1), 
                          "loss"=sapply(0:K, function(i) -logLik(m_boostpoisson[i])))
  data_list_ovar <- rbind(data_list_ovar, data_temp)
  print(nu_iter)
}
data_list_ovar$nu_fac <- as.factor(data_list_ovar$nu)
data_list_ovar$loss_diff <- data_list_ovar$loss - min(data_list_ovar$loss)

max_loss_diff_val <- max(log(data_list_ovar$loss_diff+1))
min_loss_diff_val <- exp(-25)

# loss convergence plots
gg_boost_ovar <- ggplot(data=data_list_ovar, aes(x=iter, y=log(loss_diff+1), 
                                       group=nu, color=as.factor(nu))) + 
  geom_line(size = 1.5, alpha = 0.8) + 
  theme_bw() +
  theme(legend.position = "bottom",
        text = element_text(size = 14),
        legend.box.background = element_blank(),
        legend.background = element_blank(),
        legend.margin = margin(t = -10)
  ) +
  scale_y_continuous(trans='log10', 
                     limits = c(min_loss_diff_val, max_loss_diff_val+10), 
                     expand = c(0, 0)) + 
  scale_x_continuous(expand = c(0, 0)) +
  scale_color_npg(scale_name = "what") + 
  guides(color = guide_legend(title = expression(nu))) + 
  xlab("Iteration") + 
  ylab("Cox PH Loss (Ovarian)")
gg_boost_ovar

#### Second Example
## Lung data
fm_lung <- Surv(time,status) ~ -1 + .
glm_surv_lung_boost <- glmboost(fm_lung, offset=FALSE,
                                control = boost_control(mstop = 1000, nu = 0.2), family=CoxPH(), data = lung)
plot(glm_surv_lung_boost, off2int = FALSE)
coef(glm_surv_lung_boost)

# CoxPH Model
coef(coxph(Surv(time, status) ~ -1 + ., data = lung))

#### Boosting with multiple learning rates
K <- 2100
data_list_lung <-setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("iter","nu", "loss"))
nu_list <- seq(0.1,0.6, by=0.1)
for (nu_iter in nu_list) {
  fm_lung <- Surv(time,status) ~ .
  m_boostpoisson <- glmboost(fm_lung, data = lung, offset=FALSE,
                             control = boost_control(mstop = K, nu = nu_iter), family = CoxPH())
  data_temp <- data.frame("iter"=seq(0,K), "nu"= rep(nu_iter, K+1), 
                          "loss"=sapply(0:K, function(i) -logLik(m_boostpoisson[i])))
  data_list_lung <- rbind(data_list_lung, data_temp)
  print(nu_iter)
}
data_list_lung$nu_fac <- as.factor(data_list_lung$nu)
data_list_lung$loss_diff <- data_list_lung$loss - min(data_list_lung$loss)
max_loss_diff_val <- max(log(data_list_lung$loss_diff+1))
min_loss_diff_val <- exp(-25)

# loss convergence plots
gg_boost_lung <- ggplot(data = data_list_lung, aes(x = iter, y = log(loss_diff + 1), 
                                        group = nu, color = as.factor(nu))) +  
  geom_line(size = 1.5, alpha = 0.8) +  
  theme_bw() + 
  theme(legend.position = "bottom", 
        text = element_text(size = 14), 
        legend.box.background = element_blank(), 
        legend.background = element_blank(),
        legend.margin = margin(t = -10),
  ) +  
  scale_y_continuous(trans = 'log10',  
                     limits = c(min_loss_diff_val, max_loss_diff_val + 10),  
                     expand = c(0, 0)) +  
  scale_x_continuous(expand = c(0, 0)) + 
  scale_color_npg(scale_name = "what") +  
  guides(color = guide_legend(title = expression(nu))) +  
  xlab("Iteration") +  
  ylab("Cox PH Loss (Lung)") 

gg_boost_lung

gg_surv_boost <- grid.arrange(gg_boost_ovar,gg_boost_lung, ncol=2)

ggsave(gg_surv_boost, filename = paste0("surv_boost_conv.pdf"), 
       width = 14*0.7, height = 5*0.7)
