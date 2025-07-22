# INITIAL SETTINGS --------------------------------------------------------

# clean environment
rm(list = ls())

# install packages with respective version
renv::restore()

# load packages
list.packages <-  c("tidyverse", "sf", "terra", "exactextractr", "mgcv", "pROC",
                    "sperrorest", "parallel")
vapply(list.packages, library, logical(1), character.only = TRUE, logical.return = TRUE, quietly = TRUE)
rm(list.packages)


# LOAD INPUT DATA ---------------------------------------------------------

# aoi
aoi <- sf::st_read("./dat/raw/aoi.gpkg")
su <- sf::st_read("./dat/raw/su.gpkg") %>% dplyr::mutate(label = as.factor(label))
landslides <- sf::st_read("./dat/interim/landslides_spacetime.gpkg") 

# best combination p1 and p15_1
d_rgw <- readRDS("./dat/processed/su_rainfall_45_5000_df_tune.Rds") %>% 
  dplyr::mutate(preparatory = p15_1) %>%
  dplyr::mutate(triggering = p_1) %>%
  dplyr::relocate(c(preparatory, triggering), .after=suscept_rgw) %>%
  dplyr::select(-c(p_0:p45_5))

# covariate restriction based on plausbility
d_rgw <- d_rgw %>% 
  dplyr::mutate(preparatory = dplyr::if_else(preparatory > 300, 300, preparatory)) %>% 
  dplyr::mutate(triggering = dplyr::if_else(triggering > 150, 150, triggering))


# MODELING ----------------------------------------------------------------

# formula
formula <- bin ~ 
  s(suscept_rgw, k=3, bs="tp")+
  s(preparatory, k=5, bs="tp")+
  s(triggering, k=5, bs="tp")+
  s(doy, k=4, bs="cc")

# fit
mod <- mgcv::bam(formula, family = binomial, method = "fREML", data = d_rgw, discrete=500)
summary(mod)

# fitting performance
d_rgw$susc_fitrgw <- predict(mod, type = "response", newdata = d_rgw, newdata.guaranteed = TRUE)
myroc_ftune <- pROC::roc(response = d_rgw$bin, predictor = d_rgw$susc_fitrgw, auc = T);plot(myroc_ftune, main = round(myroc_ftune$auc, 5))

# partial plots
# pdf("./plt/04_dynamic_partial_plots.pdf", width = 11, height = 8, paper="a4r")
par(mfrow = c(1,4))
par(pty="s")
plot(mod, select = 1, trans = plogis, xlim = c(0, 0.8), ylim = c(0, 1), scale = -1, seWithMean = T, shift = coef(mod)[1], ylab = "", xlab = "Static susceptibility")
plot(mod, select = 2, trans = plogis, xlim = c(0, 300), ylim = c(0, 1), scale = -1, seWithMean = F, shift = coef(mod)[1], ylab = "", xlab = "Preparatory precipitation (mm)")
plot(mod, select = 3, trans = plogis, xlim = c(0, 150), ylim = c(0, 1), scale = -1, seWithMean = F, shift = coef(mod)[1], ylab = "", xlab = "Triggering precipitation (mm)")
plot(mod, select = 4, trans = plogis, xlim = c(0, 365), ylim = c(0, 1), scale = -1, seWithMean = F, shift = coef(mod)[1], ylab = "", xlab = "Day of the year")
# dev.off()


# VALIDATION --------------------------------------------------------------

## random cross validation ----
repetition <- (1:10)
fold <- (1:10)
formula <- bin ~  
  s(suscept_rgw, k = 3, bs = "tp")+
  s(preparatory, k = 5, bs = "tp")+
  s(triggering, k = 5, bs = "tp")+
  s(doy, k = 4, bs = "cc")

# sperrorest random partition
resamp <- sperrorest::partition_cv(d_rgw, nfold = 10, repetition = 10, seed1 = 1)

# loop for validation
set.seed(1)
a <- Sys.time(); a
result_aucs <- lapply(repetition, function(i){
  id.hold <- resamp[[i]]
  lapply(fold, function(i){
    id.holdout <- id.hold[[i]]$test
    df_train <- d_rgw[-id.holdout, ]
    df_test <- d_rgw[id.holdout, ]
    fit <- mgcv::bam(formula, data = df_train, family = binomial, drop.unused.levels = F, discrete = 100, method = "fREML")
    df_test$fit <- predict(fit, type = "response", newdata.guaranteed = TRUE, newdata = df_test)
    myroc <- pROC::roc(response = df_test$bin, predictor = df_test$fit, auc=T)
    myroc <- as.numeric(unlist(myroc[9]))
  })
})
Sys.time() - a

# calculating performance
my_aucs <- as.vector(unlist(result_aucs))

# plot
# pdf("./plt/05_rcv.pdf", width = 11, height = 8, paper="a4r")
par(pty="s")
boxplot(my_aucs, ylim=c(0.6, 1), boxwex = 0.2,
        main = paste0(paste0("IQR= "),
                      paste0(round(IQR(my_aucs),4)), paste0(" median AUC= "), paste0(round(median(my_aucs), 3))))
# dev.off()

## spatial cross validation ----

# sperrorest spatial partition
resamp <- sperrorest::partition_kmeans(d_rgw, nfold = 10, repetition = 10, seed1 = 1, coords = c("X", "Y")) 

# loop for validation
set.seed(1)
a <- Sys.time(); a
result_aucs <- lapply(repetition, function(i){
  id.hold <- resamp[[i]]
  lapply(fold, function(i){
    id.holdout <- id.hold[[i]]$test
    df_train <- d_rgw[-id.holdout, ]
    df_test <- d_rgw[id.holdout, ]
    fit <- mgcv::bam(formula, data = df_train, family = binomial, drop.unused.levels = F, discrete = 100, method = "fREML")
    df_test$fit <- predict(fit, type = "response", newdata.guaranteed = TRUE, newdata = df_test)
    myroc <- pROC::roc(response = df_test$bin, predictor = df_test$fit, auc=T)
    myroc <- as.numeric(unlist(myroc[9]))
  })
})
Sys.time()-a

# calculating performance
my_aucs <- as.vector(unlist(result_aucs))

# plot
# pdf("./plt/06_scv.pdf", width = 11, height = 8, paper="a4r")
par(pty="s")
boxplot(my_aucs, ylim=c(0.6, 1), boxwex = 0.2,
        main = paste0(paste0("IQR= "),
                      paste0(round(IQR(my_aucs),4)), paste0(" median AUC= "), paste0(round(median(my_aucs), 3))))
# dev.off()

## temporal (year) cross validation ----

# setting up loop
fold <- (1:22)

# sperrorest temporal partition
resamp <- sperrorest::partition_factor_cv(d_rgw, fac = "year", nfold = 22, repetition = 1, seed1= 1) 

# validation loop
set.seed(1)
a <- Sys.time(); a
result_aucs <- lapply(fold, function(i){
  id.holdout <- resamp[["1"]][[i]]$test
  df_train <- d_rgw[-id.holdout, ]
  df_test <- d_rgw[id.holdout, ]
  fit <- mgcv::bam(formula, data = df_train, family = binomial, drop.unused.levels = F, discrete = 100, method = "fREML")
  df_test$fit <- predict(fit, type="response", newdata.guaranteed=TRUE, newdata=df_test)
  myroc <- pROC::roc(response=df_test$bin, predictor=df_test$fit, auc=T)
  myroc <- as.data.frame(cbind(year = df_test$year[1], auc = as.numeric(unlist(myroc[9]))))
})
Sys.time() - a

# calculating performance
my_aucs <- dplyr::bind_rows(result_aucs) %>% dplyr::arrange(year)

# plot
# pdf("./plt/07_mcv.pdf", width = 11, height = 8, paper="a4r")
par(pty="s")
boxplot(my_aucs$auc, ylim=c(0.6, 1),
        boxwex = 0.2, main = paste0(paste0("IQR= "),
                                    paste0(round(IQR(my_aucs$auc),4)), paste0(" median AUC= "), paste0(round(median(my_aucs$auc), 3))))
# dev.off()

## temporal (month) cross validation ----

# setting up loop
fold <- (1:12)
formula <- bin ~  s(suscept_rgw, k = 3, bs = "tp")+
  s(preparatory, k = 5, bs = "tp")+
  s(triggering, k = 5, bs = "tp")+
  s(doy, k = 4, bs = "cc")

# sperrorest spatial partition
resamp = sperrorest::partition_factor_cv(d_rgw, fac = "month", nfold = 12, repetition = 1, seed1 = 1) 

# validation loop
set.seed(1)
a = Sys.time(); a
result_aucs <- lapply(fold, function(i){
  id.holdout <- resamp[["1"]][[i]]$test
  df_train <- d_rgw[-id.holdout, ]
  df_test <- d_rgw[id.holdout, ]
  fit <- mgcv::bam(formula, data = df_train, family = binomial, drop.unused.levels = F, discrete = 100, method = "fREML")
  df_test$fit <- predict(fit, type = "response", newdata.guaranteed = TRUE, newdata = df_test)
  myroc <- pROC::roc(response=df_test$bin, predictor=df_test$fit, auc=T)
  myroc <- as.data.frame(cbind(month = df_test$month[1], auc = as.numeric(unlist(myroc[9]))))
})
Sys.time() - a

# calculating performance
my_aucs <- dplyr::bind_rows(result_aucs)

# plot
# pdf("./plt/08_ycv.pdf", width = 11, height = 8, paper="a4r")
par(pty="s")
boxplot(my_aucs$auc, ylim=c(0.6, 1),
        boxwex = 0.2, main = paste0(paste0("IQR= "),
                                    paste0(round(IQR(my_aucs$auc),4)), paste0(" median AUC= "), paste0(round(median(my_aucs$auc), 3))))
# dev.off()

