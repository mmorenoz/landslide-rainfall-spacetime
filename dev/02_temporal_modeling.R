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

# function to get cumulative rainfall 
rainfall_cum = function(data, days){
  name = paste0("p_", days)
  rainfall_po_result_df_3  = data  %>% 
    filter(days_before_event %in% c(as.character(0:(days)))) %>%
    group_by(id) %>%
    summarize(!!name := sum(precip)) %>%
    right_join(data)
}

# function to create grid search dataframe
optimal_grid <- function(triggering) {
  preparatory <- triggering + 1
  cols <- paste0("p", preparatory:45, "_", triggering)
  preparatory <- colnames(sf::st_drop_geometry(dplyr::select(d_rgw, dplyr::all_of(cols))))
  dplyr::tibble(triggering = paste0("p_", triggering), preparatory = preparatory)
}


# LOAD INPUT DATA ---------------------------------------------------------

# aoi
aoi <- sf::st_read("./dat/raw/aoi.gpkg")
su <- sf::st_read("./dat/raw/su.gpkg") %>% dplyr::mutate(label = as.factor(label))
landslides <- sf::st_read("./dat/interim/landslides_spacetime.gpkg") 


# SAMPLING ----------------------------------------------------------------

# number of days in analysis to replicate SUs
days <- as.numeric(as.Date("2020-12-31") - as.Date("1999-12-31")) 

# extract su(labels) to landslides
su <- su %>% dplyr::select(-cat, -value)
landslides.d <- sf::st_join(landslides, su, join = st_intersects) %>%
  tidyr::drop_na(doy_j)

# prepare replicates
su <- data.table::as.data.table(su)
landslides.d <- data.table::as.data.table(landslides.d)

# join SUs with landslides
spti <- left_join(su, landslides.d, by = 'label') %>%
  dplyr::mutate(bin = dplyr::if_else(is.na(id_frana), 0, 1)) %>% 
  dplyr::select(-geom.y, -Area.y) %>%
  dplyr::rename(geom = geom.x, Area = Area.x) %>%
  dplyr::relocate(bin, .after = id_frana) %>% 
  dplyr::relocate(c(Area, geom), .after = 'method3') %>% 
  tidyr::drop_na(label)

# SU replicates
spti <- spti[rep(spti[, .I], 7671)] # it may be very computationally expensive when using real number of days
spti <- spti %>% dplyr::mutate(date_su =
                                 dplyr::if_else(bin == 1, date,
                                                sample(seq.Date(as.Date("1999-12-31"), as.Date("2020-12-31"), by = "day"), nrow(spti), replace = T))) %>% 
  dplyr::relocate(date_su, .after = date)

# presences selection
slides <- spti %>%
  subset(bin == 1) %>% 
  dplyr::distinct()

# absences selection
set.seed(123)
noslides <- spti %>% 
  subset(bin == 0) %>% 
  dplyr::group_by(label) %>%
  dplyr::sample_n(10) 
gc()

# combine presences and absences
su_st <- rbind(slides, noslides) %>%
  dplyr::rename(id = label)

# clean environment
rm(spti, slides, noslides)
table(su_st$bin)

# adding the dates to the table (year, month, day)
su_st <- su_st %>% 
  dplyr::relocate(date_su, .after = woy) %>%
  dplyr::mutate(date_su = as.Date(su_st$date_su, origin = as.Date("1999-12-31"))) %>%
  tidyr::separate(date_su, c("year", "month", "day"), remove = FALSE) %>%
  dplyr::mutate(across(c("year", "month", "day"), as.integer))%>%
  dplyr::mutate(doy = lubridate::yday(date_su)) %>%
  dplyr::mutate(woy = lubridate::week(date_su)) %>%
  dplyr::mutate(date_rain = date_su +1) # to account for shifted date in grids

# rearrange columns
su_st = su_st %>%
  dplyr::relocate(geometry, .after = last_col()) %>%
  dplyr::relocate(date_su, .after = date) %>%
  dplyr::relocate(c('year', 'month', 'day', 'bin'), .after = date_su) 

# store
# sf::st_write(su_st, "./dat/interim/su_spacetime.gpkg")
# gc()


# PRECIPITATION EXTRACTION ------------------------------------------------

# load data
su_st <- sf::read_sf("./dat/interim/su_spacetime.gpkg") %>% 
  dplyr::mutate_at(vars("bin", "month", "year"), factor) %>% 
  dplyr::mutate(date_rain = as.character(date_rain)) %>%
  dplyr::mutate(date_rain = ifelse(date_rain == '2021-11-01', '2021-10-31', date_rain)) %>%
  dplyr::mutate(date_rain = as.Date(date_rain))

# parallel extraction
su_in = split(su.1, seq(nrow(su.1)))
n_cores <- detectCores() - 4
clust <- makeCluster(n_cores)
clusterExport(cl = clust, varlist = c("su_in"))
clusterEvalQ(cl = clust, c(library(rainfallR),
                           library(magrittr),
                           library(dplyr),
                           library(stringr))) 

a <- Sys.time()
rainfall_po_result <- parLapply(cl = clust, X = su_in , fun = function(x){
  ex_rainfall(
    data_path = "./dat/raw/prec_grids",
    spatial.obj = x,
    fun = "mean",
    date = x$date_rain,
    nc_var = "precipitation",
    days_back = 45)
})
stopCluster(cl = clust)
Sys.time() - a

# store data
# saveRDS("./dat/interim/rainfall_ex_df_45_5000.Rds")


# AGGREGATE PRECIPITATION -------------------------------------------------

# load data
rainfall_po_result_df <- readRDS("./dat/interim/rainfall_ex_df_45_5000.Rds")

# generate an id every 46 rows
rainfall_po_result_df = rainfall_po_result_df %>%
  dplyr::group_by(id) %>%
  dplyr::mutate(id = ceiling(row_number() / 46)) %>%
  dplyr::ungroup()

# testing for 0 day cumulative rainfall
test_cum <- rainfall_cum(rainfall_po_result_df, 0)

# run function in a loop for calculating all the cumulative rainfall, p_0 must be computed first!
days <- 45
a <- Sys.time()
for (i in (1:days)){
test_cum <- rainfall_cum(test_cum, i)}
Sys.time() - a

# keep only one row (delete the other 18 copies) for the second function
test_cum <- dplyr::distinct(test_cum , id, .keep_all = TRUE)
su_rainfall <- test_cum %>%
  dplyr::relocate(p_0:p_45, .after = day_prs) %>%
  dplyr::select(-date, -fun, -precip, -days_before_event)

# clean environment
rm(test_cum, i, a, days)

# store data
# saveRDS(su_rainfall, "./dat/interim/su_rainfall_45_5000_dateshifted.Rds")


# EXTRACT STATIC PROB -----------------------------------------------------

# load dataset with precipitation
su_rainfall <- readRDS("./dat/interim/su_rainfall_45_5000_dateshifted.Rds") %>% sf::st_as_sf()

# load grids
grids <- list.files("./dat/raw/grids", pattern = ".tif$", full.names = T)[c(1, 13, 14)]
grids <- c(grids, list.files("./dat/processed/", pattern = "susceptibility.tif$", full.names = T))
name_grids <- gsub("./dat/raw/grids/", "", gsub(".tif", "", grids))
name_grids[4] <- "susceptibility"; name_grids
grids <- terra::rast(grids); names(grids) <- name_grids
rm(name_grids)
terra::crs(grids) <- terra::crs(aoi)

# grid are already masked with the trivial areas
# aggregation
d_rgw <- dplyr::tibble(
  dplyr::bind_cols(
    su_rainfall,
    id_su = exactextractr::exact_extract(grids$su, su_rainfall, "max"),
    binspace = exactextractr::exact_extract(grids$space30, su_rainfall, "max"),
    countspace = exactextractr::exact_extract(grids$space30, su_rainfall, "sum"),
    area_su = exactextractr::exact_extract(grids$area30, su_rainfall, "sum"),
    susceptibility = exactextractr::exact_extract(grids$susceptibility, su_rainfall, 'mean')
  )
)

# arraging new columns
d_rgw <- d_rgw %>%
  dplyr::mutate(countspace = ceiling(countspace)) %>%
  dplyr::mutate(id = as.integer(id)) %>%
  dplyr::mutate_at(vars("bin", "binspace", "id_su"), factor) %>% 
  dplyr::select(-cumsum, -dol) %>% 
  dplyr::relocate(susceptibility, p_0:p_45, .before = geom) %>% 
  dplyr::relocate(id_su, area_su, binspace, countspace, .after = bin)


# PRECIPITATION TIME WINDOWS ----------------------------------------------

# triggering precipitation for 0, 1, 2, 3, 4 and 5 days and preparatory up tp 45 days
df <- d_rgw
for (i in c(0, 1, 2, 3, 4, 5)){
  range <- (1: (45-i))
  for (j in range){
    ant <- i + j
    name <- paste0("p", ant, "_", i)
    tri <- paste0("p_", i)
    ant <- paste0("p_", ant)
    df[[name]] <- df[[ant]] - df[[tri]]
  }
}

# store results
# saveRDS(df, "./dat/interim/su_rainfall_45_5000_statdyna.Rds")

# clean up
rm(ant, i, j, name, range, tri, df)


# TRIVIAL TIMES -----------------------------------------------------------

# load data
d_rgw <- readRDS("./dat/interim/su_rainfall_45_5000_statdyna.Rds") %>% sf::st_as_sf()

# removing nas
d_rgw <- d_rgw %>% 
  tidyr::drop_na(susceptibility) # drops 990 SUs

# trivial times (precipitation lower than 1.1 is assumed to be as no precipitation)
d_rgw <- dplyr::filter(d_rgw, p_1 > 1.1) # drops 24941 SUs

# store
# saveRDS(d_rgw, "./dat/interim/1su_rainfall_45_5000_df_tune.Rds")


# GRID SEARCH -------------------------------------------------------------

# build the different combinations of triggering and preparatory 
grid_search <- do.call(dplyr::bind_rows, lapply(0:5, optimal_grid)) %>% 
  dplyr::mutate(tr = sub("p_", "", triggering)) %>% 
  dplyr::mutate(pr = sub("_.", "", preparatory)) %>% 
  dplyr::mutate(pr = sub("p", "", pr))

# setting up loop
combination <- 1:255
repetition <- 1:10
fold <- 1:10
my_exclude <- "s(year)"

# sperrorest partition based on SU id
resamp <- sperrorest::partition_factor_cv(dplyr::select(sf::st_drop_geometry(d_rgw), id, id_su), fac = "id_su", nfold = 10, repetition = 10, seed1= 1) 

# loop with parallel
n_cores <- parallel::detectCores() - 2
clust <- parallel::makeCluster(n_cores)
parallel::clusterExport(cl = clust, varlist = c("grid_search", "repetition",
                                                "resamp", "fold", "d_rgw",
                                                "my_exclude"))
a <- Sys.time(); a
result_aucs <- parallel::parLapply(combination, cl=clust, function(i){
  formula <- paste("bin~",
                   paste("s("), paste(grid_search[i,1]),paste(") +"), 
                   paste("s("), paste(grid_search[i, 2]),paste(") +"),  
                   paste("s(doy, bs='cc')+"), 
                   paste("ti("), paste(c(grid_search[i, 1])), paste(", doy, bs=c('tp', 'cc'))+"),   
                   paste("ti("), paste(c(grid_search[i, 2])), paste(", doy, bs=c('tp','cc'))+"), 
                   paste("ti("), paste(c(grid_search[i, 1])), paste(","), paste(c(grid_search[i, 2])), paste(", bs=c('tp', 'tp'))+"),
                   paste("s(year,bs='re')")) 
  print(formula)
  formula <- as.formula(formula)
  set.seed(1)
  lapply(repetition, function(i){
    id.hold <- resamp[[i]]
    lapply(fold, function(i){
      id.holdout <- id.hold[[i]]$test
      df_train <- d_rgw[-id.holdout, ]
      df_test <- d_rgw[id.holdout, ]
      fit <- mgcv::bam(formula, data = df_train, family = binomial, drop.unused.levels = F, discrete = 100, method = "fREML")
      df_test$fit <- predict(fit, type = "response", newdata.guaranteed = TRUE, newdata = df_test, exclude = my_exclude )
      myroc <- pROC::roc(response = df_test$bin, predictor = df_test$fit, auc = T)
      myroc <- as.numeric(unlist(myroc[9]))
    })
  })
})
Sys.time()-a
parallel::stopCluster(cl = clust)

# storing results
# saveRDS(result_aucs, "./dat/processed/su_spacetime_45_5000_finetuned.Rds")

# load results
result_aucs <- readRDS("./dat/interim/su_spacetime_45_5000_finetuned.Rds")

# calculating performance
my_aucs <- lapply(result_aucs, function(x) do.call(rbind, x))
my_aucs <- data.table::rbindlist(my_aucs)
my_aucs <- t(my_aucs) %>% dplyr::as_tibble()
names(my_aucs) <- grid_search$preparatory

# dataframe
grid_search <- grid_search %>%
  dplyr::mutate(auc_mean = sapply(my_aucs, FUN = mean)) %>%
  dplyr::mutate(auc_median = sapply(my_aucs, FUN = median)) %>%
  dplyr::mutate(auc_sd = sapply(my_aucs, FUN = sd)) %>%
  dplyr::mutate(auc_iqr = sapply(my_aucs, FUN = IQR)) %>% 
  dplyr::mutate_at(vars(tr, pr), as.numeric)
title <- dplyr::filter(grid_search, auc_median == max(auc_median))

# plot
ggplot2::ggplot(grid_search, aes(x = pr, y = tr, fill = auc_median))+
  geom_tile(color = "white", lwd = 0.07, linetype = 1)+
  scale_fill_gradientn(colors = hcl.colors(50, "Blues", rev = -1))+
  ggtitle(paste(paste("Median AUC="), paste(c(round(title[6],3))), paste("for preparatory="), 
                paste(title[2],paste("and triggering="), title[1])))+
  scale_y_continuous(expand = c(0, 0), breaks = 0:5) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 45, by = 5)) +
  xlab("Preparatory precipitation (P)")+
  ylab("Triggering precipitation (T)") +
  coord_fixed(ratio = 1) +
  guides(fill = guide_colourbar(barwidth = 8, barheight = 1, title = "Median AUROC")) +
  theme_classic() + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.text.y   = element_text(size = 14),
        axis.text.x   = element_text(size = 14),
        axis.title.y  = element_text(size = 14, vjust = 2),
        axis.title.x  = element_text(size = 14, vjust = -2),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5))
# ggplot2::ggsave("./plt/03_gridsearch.pdf", width = 11, height = 8, paper = "a4r")
