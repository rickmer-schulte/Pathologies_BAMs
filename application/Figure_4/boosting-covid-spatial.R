library(dplyr)
library(mboost)
library(mgcv)
library(tidyr)
library(ggplot2)
library(sf)
library(dplyr)
library(maps)
library(ggmap)
library(mapproj)
library(viridis)

# Map data
state_df <- map_data("state")
county_df <- map_data("county")

# Load data  
df <- readRDS("subset_covid.RDS")
  
# PLS/OLS
mod_gam <- bam(prevalence ~ -1 + te(latitude, longitude, 
                                    k=c(20,20), 
                                    bs = "ps", 
                                    # choose some penalty that is isotropic
                                    sp = c(5,5)),
               data = df)

pred_gam <- predict(mod_gam, newdata = data.frame(longitude = county_df$long,
                                                  latitude = county_df$lat))

mod_ols <- bam(prevalence ~ -1 + te(latitude, longitude, 
                                    k=c(20,20), 
                                    bs = "ps", 
                                    # small but negligible penalty
                                    sp = c(1e-2,1e-7) 
                                    ),
               data = df)
pred_ols <- predict(mod_ols, newdata = data.frame(longitude = county_df$long,
                                                  latitude = county_df$lat))


# Fit the model using mboost
model_mb <- mboost(prevalence ~ 
                     bbs(latitude, degree = 2, differences = 1,
                         knots = mod_ols$smooth[[1]]$margin[[1]]$knots[2:23],
                         boundary.knots = 
                           c(mod_ols$smooth[[1]]$margin[[1]]$knots[c(1,24)])) %X% 
                     bbs(longitude, degree = 2, differences = 1,
                         knots = mod_ols$smooth[[1]]$margin[[2]]$knots[2:23],
                         boundary.knots = 
                           c(mod_ols$smooth[[1]]$margin[[2]]$knots[c(1,24)])), 
                   data = df, offset = 0,
                   control = boost_control(nu = 0.5, mstop = 10))

pred_mb_short <- predict(model_mb, 
                         newdata = data.frame(longitude = county_df$long,
                                              latitude = county_df$lat))
pred_mb_long <- predict(model_mb[10000],
                        newdata = data.frame(longitude = county_df$long,
                                             latitude = county_df$lat))

county_df$OLS <- c(pred_ols) #[,1]
county_df$`Boosting (10)` <- pred_mb_short[,1]
county_df$`Boosting (1000)` <- pred_mb_long[,1]
county_df$PLS <- c(pred_gam)

# plot 
state_df_proj <- map_data("state", projection = "albers", parameters = c(39, 45))
state_df <- state_df %>% left_join(state_df_proj, by = c("group", "order"))
county_df_proj <- map_data("county", projection = "albers", parameters = c(39, 45))
plot_df <- county_df %>% left_join(county_df_proj, by = c("group", "order"))

ggplot(plot_df %>% pivot_longer(OLS:PLS) %>% 
         group_by(name) %>% 
         mutate(value = value - mean(value),
                name = factor(name, levels = c("PLS", "Boosting (10)", "Boosting (1000)", 
                                               "OLS"))
                ) %>% 
         ungroup() %>% 
         mutate(
           value = pmin(value, quantile(value, 0.95)),
           value = pmax(value, quantile(value, 0.05))
         ), 
       aes(long.y, lat.y, group = group)) + 
  geom_polygon(aes(fill = value), colour = alpha("white", 0.1), 
               linewidth = 0.12) +
  geom_polygon(data = state_df, colour = "white", fill = NA, 
               linewidth = 0.2) +
  coord_fixed() +
  theme_minimal() +
  theme(axis.line = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), axis.title = element_blank()) +
  scale_fill_viridis_c(option = 'mako', direction = -1, 
                       name = "Estimated Prevalence") + 
  theme(legend.position = "right",
        legend.key.width=unit(0.3,"cm"),
        legend.key.height=unit(.6,"cm"),
        text = element_text(size = 10),
        legend.title.position = "left",
        legend.title=element_text(size=10, angle = 90),
        legend.title.align=0.5) + 
  facet_wrap(~name, ncol=4) 

ggsave(file = "spatial2.pdf", width = 12, height = 6)
