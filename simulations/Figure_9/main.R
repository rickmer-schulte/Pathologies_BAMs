# Load scripts
source("models.R")
source("simfun.R")
library(tidyverse)
library(viridis)

# Run analysis

## true underlying function for simulation of multi-modal coefficient surfaces
hillyFun <- function(s,t) sin(2 * abs(s-t)) * cos(2 * t)
## limit function
lD <- function(s, t) s < t
## fixed settings
obsPerTra = 40
nn = c(300)
ss = 1
bd = c(11) 
k = 20
stops <- c(10,50,1000)

set.seed(42)

## create data using a variation of the pffrSim function
train <- pffrSimVar(scenario=c("int", "ff"), 
                    n = nn, 
                    nxgrid = obsPerTra, 
                    nygrid = obsPerTra, 
                    SNR = ss,
                    simFun = hillyFun, 
                    bd = bd)
train$Y <- train$Y - mean(train$Y)

# center data
train$X1 <- t(t(train$X1) - colMeans(train$X1))

# # testdata
# test <- pffrSimVar(scenario=c("int", "ff"), 
#                    n = nn, 
#                    nxgrid = obsPerTra, 
#                    nygrid = obsPerTra, 
#                    SNR = ss,
#                    simFun = hillyFun, 
#                    bd = bd)
# test$Y <- test$Y - mean(train$Y)
# 
# # center data
# test$X1 <- t(t(test$X1) - colMeans(test$X1))

test1 <- hillyFun
t <- attr(train, "yindex")
s <- attr(train, "xindex")

# calcuate true coefficient function
trueBeta <- outer(s, t, test1) 

# fit pffr function
m1 <- pffr(Y ~ 1 + ff(X1, xind=s, 
                      splinepars = list(bs = "ps", m = list(2, 1), k = k), 
                      yind = t), 
           bs.yindex = list(bs = "ps", k = 20, m = c(2, 2)),
           data = train)

g1 <- coef(m1)
g1 <- matrix(g1$smterms[[2]]$value, ncol=40)

# fit pffr function unpenalized
m0 <- pffr(Y ~ 1 + ff(X1, xind=s, 
                      splinepars = list(bs = "ps", m = list(2, 1), k = k, 
                                        sp = c(0.00001, 0.00001)), 
                      yind = t),
           bs.yindex = list(bs = "ps", k = 20, m = c(2, 2), sp = c(0.000)),
           data = train)

g0 <- coef(m0)
g0 <- matrix(g0$smterms[[2]]$value, ncol=40)

# fit FDboost function
train <- as.list(train)
train$t <- t
train$s <- s

m2 <- FDboost(Y ~ 1 + bsignal(X1, 
                              s = s, 
                              df = 5, 
                              knots = k),
              control = boost_control(mstop = 0),
              timeformula = ~ bbs(t, df = 25),
              data = train)

g2_list <- lapply(stops, function(m) coef(m2[m],which=2)$smterms[[1]]$value)
g2_list[[1]] <- coef(FDboost(Y ~ 1 + bsignal(X1, 
                                           s = s, 
                                           df = 5, 
                                           knots = k),
                           control = boost_control(mstop = 10),
                           timeformula = ~ bbs(t, df = 5),
                           data = train),which=2)$smterms[[1]]$value


rm(m2); gc()

plotdata <- data.frame(
  what = factor(rep(c("Truth", "FAM", paste0("BAM (", stops, ")"), "Unpenalized"),
                    each = 40^2), 
                levels = c("Truth", "FAM", paste0("BAM (", stops, ")"), "Unpenalized")),
  value = c(c(trueBeta), c(g1), unlist(lapply(g2_list, c)), c(g0)),
  s = rep(rep(1:40, each = 40), length(stops)+3),
  t = rep(rep(1:40, 40), length(stops)+3)
) 

# if(!file.exists("plotdata.RDS")){
  saveRDS(plotdata, "plotdata.RDS")
# }else{
  # plotdata <- readRDS("plotdata.RDS")
# }

plotdata %>% mutate(value = pmax(pmin(value, 5),-5)) %>%  
  filter(what != "Truth") %>% 
  ggplot(aes(x = s, y = t, z = value, fill = value)) + 
  geom_tile() + 
  geom_contour(color = "black", alpha = 0.2) + 
  facet_wrap(~ what, scales="free", ncol=5) + 
  theme_bw() + 
  xlab("s") + ylab("t") +  # scale_fill_continuous() + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1), 
    panel.background = element_rect(fill = "white"),
    strip.background = element_rect(fill = "white", color = "black", size = 1),
    strip.text = element_text(color = "black"),
    text = element_text(size = 15),
    # legend.position = "bottom",
    # legend.box.margin = margin(6, 6, 6, 6),
    # legend.key.width = unit(1.5, "cm"),
    legend.title = element_text(margin = margin(0, 0, 15, 0)),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
  ) + scale_fill_viridis(option = "A", name = expression(paste(beta,"(s,t)"))) 

ggsave("estexample.pdf", width = 13, height = 3)
