# INITIAL SETTINGS --------------------------------------------------------

# clean environment
rm(list = ls())

# install packages with respective version
renv::restore()

# load packages
list.packages <-  c("tidyverse", "sf", "terra", "exactextractr", "mgcv", "pROC", "tmap", "sperrorest")
vapply(list.packages, library, logical(1), character.only = TRUE, logical.return = TRUE, quietly = TRUE)
rm(list.packages)


# LOAD INPUT DATA ---------------------------------------------------------

# aoi
aoi <- sf::st_read("./dat/raw/aoi.gpkg")
su <- sf::st_read("./dat/raw/su.gpkg")
landslides <- sf::st_read("./dat/interim/landslides_space.gpkg")

# load covariates
grids <- list.files("./dat/raw/grids", pattern = ".tif$", full.names = T)
name_grids <- gsub("./dat/raw/grids/", "", gsub(".tif", "", grids))
grids <- terra::rast(grids); names(grids) <- name_grids
rm(name_grids)
terra::crs(grids) <- terra::crs(aoi)

# load trivial area mask
trivial <- grids$rockglacierwater_inv30


# AGGREGATE DATA TO SU ----------------------------------------------------

# grid are already masked with the trivial areas
# aggregation
d_rgw <- dplyr::tibble(id = su$cat,
                       bin = exactextractr::exact_extract(grids$space30, su, "max"),
                       count = exact_extract(grids$space30, su, "sum"),
                       area = exact_extract(grids$area30, su, "sum"),
                       catchment = exact_extract(grids$catchment, su, "majority"),
                       elevation_max = exact_extract(grids$dtm30, su, "max"),
                       elevation_min = exact_extract(grids$dtm30, su, "min"),
                       effective_area = exact_extract(grids$effectivearea_opt_rgw_miss, su, "weighted_sum", weights = grids$area30),
                       effective_area_prob = exact_extract(grids$effectivearea_prob_opt_rgw_miss, su, "mean"),
                       geologyre = exact_extract(grids$geologyre30, su, "majority"),
                       municip = exact_extract(grids$municip30, su, "majority"),
                       concavity = exact_extract(grids$concavity30, su, "mean"),
                       slope_m = exact_extract(grids$slope30, su, "mean"),
                       forest = exact_extract(grids$vgk30_forest, su, "weighted_sum", weights = grids$area30),
                       geometry = su$geom) %>% sf::st_as_sf()
                       

# correction for count and proportions
d_rgw <- d_rgw %>%
  dplyr::mutate(count = ceiling(count)) %>%
  dplyr::mutate(effective_area = effective_area * 100/area) %>%
  dplyr::mutate(forest = forest * 100/area) %>%
  dplyr::mutate(relief = elevation_max - elevation_min) %>%
  dplyr::mutate(id = as.integer(id)) %>%
  dplyr::mutate_at(vars("bin", "count", "catchment", "geologyre", "municip"), factor) %>% 
  tidyr::drop_na(slope_m) %>% 
  dplyr::select(-c(elevation_max, elevation_min))


# MODELING ----------------------------------------------------------------

# covariate restriction based on plausbility
d_rgw <- d_rgw %>% 
  dplyr::mutate(slope_m = dplyr::if_else(slope_m > 50, 50, slope_m)) %>% 
  dplyr::mutate(area = dplyr::if_else(area > 6e+06, 6e+06, area)) %>% 
  dplyr::mutate(relief = dplyr::if_else(relief > 2000, 2000, relief)) %>%
  dplyr::mutate(concavity = dplyr::if_else(concavity > 67, 67, concavity)) %>% 
  dplyr::mutate(concavity = dplyr::if_else(concavity < 38, 38, concavity))

# formula
formula <- bin ~ 
  s(slope_m, k = 5) +
  s(concavity, k = 5) +
  s(forest, k = 5) + 
  s(relief, k = 5) +
  geologyre +
  s(effective_area_prob, k = 3) +
  s(area, k = 3) +
  s(catchment, bs = "re")

# fit
mod <- mgcv::bam(formula, family = binomial, method = "fREML", data = d_rgw, discrete = 100)
summary(mod)

# partial effects
# pdf("./plt/01_static_partial_effects.pdf", width = 11, height = 8, paper="a4r")
par(pty="s")
par(mfrow = c(3,3))
plot(mod, select = 1, trans = plogis, xlim = c(0, 50), scale = -1, seWithMean = F, shift = coef(mod)[1], ylab ="", xlab = "Slope stepness (°)")
plot(mod, select = 2, trans = plogis, xlim = c(38, 67), scale = -1, seWithMean = F, shift = coef(mod)[1], ylab = "", xlab = "Concavity")
plot(mod, select = 3, trans = plogis, xlim = c(0,100), scale = -1, seWithMean = F, shift = coef(mod)[1], ylab = "",xlab = "Proportion of forest (%)")
plot(mod, select = 4, trans = plogis, xlim = c(0, 2000), scale = -1, seWithMean = F, shift = coef(mod)[1], ylab = "", xlab = "Relief (m)")
plot(mod, select = 5, trans = plogis, xlim = c(0, 1), scale = -1, seWithMean = T, shift = coef(mod)[1], , ylab = "", xlab="Effectively surveyed area")
plot(mod, select = 6, trans = plogis, xlim = c(0, 6e+06), scale = -1, seWithMean = T, shift = coef(mod)[1], ylab ="", xlab="Slope unit area (m²)")
plot(mod, select = 7, all.terms = T, trans = plogis, xlim = c(0,100), scale = -1, seWithMean = F, shift = coef(mod)[1], ylab = "", main = "", xlab="Catchment ID")
plot(mod, select = 8, all.terms = T, trans = plogis, xlim = c(0,100), scale = -1, seWithMean = F, shift = coef(mod)[1], ylab = "", xlab="Lithology")
# dev.off()

# fitting performance
d_rgw$static_prob <- predict(mod, type = "response", newdata = d_rgw, newdata.guaranteed = TRUE)
myroc_rgw <- pROC::roc(d_rgw$bin, d_rgw$static_prob);print(myroc_rgw)
par(mfrow = c(1,1))
plot(myroc_rgw)

# exclude from the prediction
exclude_rgw <- c("s(effective_area_prob)", "s(area)", "s(municip)", "s(catchment)")
d_rgw$static_prob <- predict(mod, type = "response", newdata = d_rgw, exclude = exclude_rgw, newdata.guaranteed = TRUE)
myroc_rgw <- pROC::roc(d_rgw$bin, d_rgw$static_prob);print(myroc_rgw)
plot(myroc_rgw)


# PLOT --------------------------------------------------------------------

# prepare rasters to plot
trivial <- grids$rockglacierwater30 %>% terra::mask(terra::vect(su))
hillshade <- terra::rast("./dat/raw/grids/hillshade30.tif") 
static_prob <- terra::rasterize(terra::vect(d_rgw), terra::rast(terra::vect(d_rgw), resolution = 30), field = "static_prob")

# store
# terra::writeRaster(static_prob, "./dat/processed/susceptibility.tif")

# parameters for plots
break_plot <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
label_plot <- c("<0.1", "<0.2", "<0.3", "<0.4", "<0.5", "<0.6","<0.7", "<0.8", "<0.9", "<1")
col_plot <- c("#1E8E99","#51C3CC","#99F9FF","#B2FCFF","#E5FFFF","#FFCA99","#FFAD65","#FF8E32","#CC5800","#993F00")

map <- 
  tm_shape(hillshade) +
  tm_raster(
    col.scale = tm_scale(values = gray(seq(1,0,-0.1)), n = 200),
    col.legend = tm_legend_hide()
    ) + 
  tm_graticules(alpha=0.5, labels.size = 0.5) +
  tm_shape(static_prob) +  
  tm_raster(
    col.scale = tm_scale(values = col_plot, breaks = break_plot),
    col.legend = tm_legend_hide(),
    col_alpha=0.6
    ) +
  tm_add_legend(
    type = "polygons",
    labels = label_plot,
    fill = col_plot,
    is.portrait = F,
    alpha = 0.7,
    border.lwd = 0.5,
    title = "Static probability"
    ) +
  tm_shape(trivial) +
  tm_raster(
    col.scale = tm_scale_discrete(values = "black"),
    col.legend = tm_legend_hide()
    ) +
  tm_shape(aoi) +
  tm_borders(col = "black",
             lwd=1
             ) +
  tm_shape(su) +
  tm_borders(col = "black",
             lwd=0.15
             ) +
  tm_legend(title = "",
            legend.format = c(digits = 1),
            main.title.size = 0.5,
            orientation = "landscape",
            show = T,
            bg.color = "white",
            legend.outside= T,
            legend.outside.position = "bottom",
            legend.text.size = 1.2,
            width = 20,
            ) +
  tm_scalebar(breaks = c(0,10,30),
              text.size = 0.6,
              bg.color = "white",
              bg.alpha = 0.8,
              position = c("right", "bottom")
              ) +
  tm_compass(type = "arrow",
             position = c("left", "top"),
             show.labels = 1
             ) +
  tm_layout(title.size = 0.1,
            frame = T,
            title.position = c("left", "bottom")
            )
# png("./plt/02_static_suscep.png", width=1600, height=1450, units="px", res = 300) 
# print(map)
# dev.off()


# VALIDATION --------------------------------------------------------------

# setting up loop for random cross-validation
fold <- 10
repetition <- 10

# create random cv partition
set.seed(123)
partition <- sperrorest::partition_cv(d_rgw, nfold = fold, repetition = repetition, seed1 = 123) 
myroc_cv <- lapply(my.list<-vector(mode = 'list', 10), function(x) x<-vector(mode = 'list', 10))

# loop for validation
for (i in 1:repetition){
  id.hold <- partition[[i]]
  for (j in 1:fold){
    id.holdout <- id.hold[[j]]$test
    d.test <- d_rgw[id.holdout, ]
    d.train <- d_rgw[-id.holdout, ]
    fit <- mgcv::bam(formula, data = d.train, family = binomial, method = "fREML", discrete = 100, drop.unused.levels = F)
    d.test$prediction <- predict(fit, d.test, type = "response", exclude = exclude_rgw)
    myroc <- unlist(pROC::roc(response = d.test$bin, predictor = d.test$prediction, auc=T))$auc
    myroc_cv[[i]][[j]] <- myroc
  }
}

# clean environment
rm(d.train, d.test, id.hold, id.holdout, myroc, i, j, resamp, partition, fold, repetition)

# store performances in data frame
myroc_cv = dplyr::as_tibble(do.call(rbind, lapply(myroc_cv, unlist))) %>% 
  tidyr::gather(key = "repetition", value = "auc") %>% 
  dplyr::mutate(repetition = rep(1:10, each = 10))

# plot
ggplot2::ggplot(myroc_cv, aes(auc)) + geom_boxplot(width = 0.5) + 
  ggtitle(paste0(paste0("IQR_cv= "), paste0(round(IQR(myroc_cv$auc), 3)),
                 paste0(" median AUC_cv= "), paste0(round(median(myroc_cv$auc), 2)))) +
  ylab("") + xlab ("AUC")+
  coord_flip()