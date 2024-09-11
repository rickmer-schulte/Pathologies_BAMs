library(refund)
library(FDboost)

pred_fdboost <- function(form, timeformula, train, test, ...){
  
  mod <- FDboost(form, data = train,
                 timeformula = timeformula, 
                 control = boost_control(mstop = 5000), ...)
  cvr <- cvrisk(mod, folds = 
                  cvLong(id = mod$id, weights = model.weights(mod),
                         B = 5, "kfold"))
  mod[mstop(cvr)]
  
  predict(mod, newdata = test, type = "response")
  
}

pred_pffr <- function(form, time, train, test, ...){
  
  mod <- pffr(form, yind=time, data = train, ...)
  predict(mod, newdata = test, type = "response")
  
}