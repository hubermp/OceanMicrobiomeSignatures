# =============================================================================
# Species Distribution Models (SDMs) for Ocean Microbial Signatures (OMSs)
#
# Associated manuscript:
# Huber et al. – Ocean Microbiome Signatures capturing the spatial organization
# of metabolic potential in the global sunlit ocean – Nature Microbiology (under review)
#
# Description:
# This script fits ensemble Species Distribution Models (SDMs) for each OMS
# using four algorithms: GLM, Random Forest (RF), ANN, and Boosted Regression
# Trees (GBM/BRT), implemented via the h2o R package.
# Models are evaluated using 5-fold cross-validation (AUC threshold > 0.7).
# Calibrated models are projected onto global environmental rasters to generate
# Habitat Suitability Index (HSI) maps for each OMS.
#
# Input files:
#   - Pres_Abs.txt         : Presence/absence matrix per OMS per sample (Table S11)
#   - Variables_ready/     : Directory with environmental raster layers (.tif)
#                            from World Ocean Atlas 2018 (WOA18)
#   - M02.tif              : Mixed Layer Depth raster
#
# Output files:
#   - OMS_X_meanModel_epi.tif   : Mean ensemble HSI map per OMS (epipelagic)
#   - OMS_X_TSS_BinModel_epi.tif: Binarized HSI map using TSS threshold
#   - ModelPerformance_OMS_X.rds: Model performance metrics per OMS
#   - VarImp_OMS_X.rds          : Variable importance per OMS
#   - AUC_list2.rds             : AUC scores for all models and OMSs
#   - TSS_TABLE.csv             : TSS thresholds per OMS
#
# Required R packages:
#   raster, terra, usdm, sf, dplyr, MASS, h2o, dismo, reshape,
#   ggplot2, tidyr, maps, tidyterra, ggpubr, gridExtra, cowplot,
#   RColorBrewer, GGally, viridis
#
# Usage:
#   Set the working directories below to match your local folder structure.
#   Environmental raster layers should follow WOA18 naming conventions.
# =============================================================================


# ── 0. Load packages ──────────────────────────────────────────────────────────
library("raster")
library("usdm")
library("sf")
library("dplyr")
library("MASS")
library("h2o")
library("dismo")
library("ggplot2")
library("tidyr")
library("maps")
library("tidyterra")
library("ggpubr")
library("gridExtra")
library("cowplot")
library("RColorBrewer")
library("GGally")
library("terra")
library("viridis")


# ── 1. Set paths (adjust to your local directory structure) ───────────────────
path_vars   <- "Variables_ready"   # directory with .tif environmental layers (WOA18)
path_output <- "output"            # directory for model outputs
path_figs   <- "Figuras"           # directory for figures

mld_file    <- "M02.tif"          # Mixed Layer Depth raster
presabs_file <- "Pres_Abs.txt"    # Presence/absence table (Table S11)


# ── 2. Load and stack environmental variables ─────────────────────────────────
MLD <- raster::raster(mld_file)
names(MLD) <- "mld"

all_layers <- list.files(path = path_vars, pattern = ".tif", full.names = TRUE)

# Stack by ocean layer (epipelagic, mesopelagic, bathypelagic)
stacked_epi   <- stack(all_layers[grep("EPI",   ignore.case = TRUE, all_layers)])
stacked_meso  <- stack(all_layers[grep("MESO",  ignore.case = TRUE, all_layers)])
stacked_bathy <- stack(all_layers[grep("BATHY", ignore.case = TRUE, all_layers)])

# Clean layer names
clean_names <- function(x, path) {
  gsub("_withNA", "", gsub(".tif", "", gsub("/", "", gsub(path, "", x))))
}
names(stacked_epi)   <- clean_names(all_layers[grep("EPI",   ignore.case = TRUE, all_layers)], path_vars)
names(stacked_meso)  <- clean_names(all_layers[grep("MESO",  ignore.case = TRUE, all_layers)], path_vars)
names(stacked_bathy) <- clean_names(all_layers[grep("BATHY", ignore.case = TRUE, all_layers)], path_vars)


# ── 3. Compute N* (excess nitrate over Redfield ratio) ───────────────────────
# N* = [NO3-] - 16 * [PO4 3-]  (Redfield ratio)
Nstar_epi   <- stacked_epi[["n_epi"]]   - 16 * stacked_epi[["p_epi"]]
Nstar_meso  <- stacked_meso[["n_meso"]] - 16 * stacked_meso[["p_meso"]]
Nstar_bathy <- stacked_bathy[["n_bathy"]] - 16 * stacked_bathy[["p_bathy"]]

names(Nstar_epi)   <- "Nstar_epi"
names(Nstar_meso)  <- "Nstar_meso"
names(Nstar_bathy) <- "Nstar_bathy"

stacked_epi   <- raster::stack(MLD, stacked_epi,   Nstar_epi)
stacked_meso  <- raster::stack(MLD, stacked_meso,  Nstar_meso)
stacked_bathy <- raster::stack(MLD, stacked_bathy, Nstar_bathy)


# ── 4. VIF check — select final predictor set (VIF < 3) ──────────────────────
# Variables retained after iterative VIF removal (see Methods):
# mld, salinity (s_), irradiance (I_), oxygen saturation (O2sat_),
# temperature (t_), N* (Nstar_)
PREDICTOR_NAMES <- c("mld", "s_", "I_", "O2sat_", "t_", "Nstar_")

# Select final variable stacks (excluding collinear variables)
vars_epi   <- stacked_epi[[which(!names(stacked_epi) %in%
                c("mld", "C_epi", "o_epi", "I_epi", "si_epi",
                  "ugo_epi", "vgo_epi", "n_epi", "O2sat_epi", "p_epi"))]]
vars_meso  <- stacked_meso[[which(!names(stacked_meso) %in%
                c("n_meso", "p_meso", "A_meso", "C_meso", "O2sat_meso"))]]
vars_bathy <- stacked_bathy[[which(!names(stacked_bathy) %in%
                c("n_bathy", "p_bathy", "A_bathy", "C_bathy", "O2sat_bathy"))]]

# Verify correlation structure
ggcorr(data.frame(rb2_[, PREDICTOR_NAMES]), digits = 2, label = TRUE)


# ── 5. Load presence/absence data ─────────────────────────────────────────────
DF      <- read.table(presabs_file, sep = "\t")
coord   <- data.frame(DF$decimalLongitude, DF$decimalLatitude)
inf     <- raster::extract(stacked_epi, coord)
RF      <- data.frame(cbind(DF, inf))
write.table(RF, file.path(path_output, "RF.txt"), sep = "\t")

# Binarize OMS columns (presence = 1, absence = 0)
DF_bin <- data.frame(
  ID    = DF$X,
  long  = DF$decimalLongitude,
  lat   = DF$decimalLatitude,
  DF[, 50:59],
  Depth        = DF$Depth,
  column_depth = DF$water_col
)
DF_bin[DF_bin > 0] <- 1
DF_bin <- DF_bin[!is.na(DF_bin$long) & !is.na(DF_bin$lat) & !is.na(DF_bin$Depth), ]
DF_bin <- SpatialPointsDataFrame(data.frame(DF_bin$long, DF_bin$lat), data = DF_bin)

# Split by depth layer
DF_epi   <- DF_bin[DF_bin@data$Depth %in% c("surface", "epipelagic"), ]
DF_meso  <- DF_bin[DF_bin@data$Depth == "mesopelagic", ]
DF_bathy <- DF_bin[DF_bin@data$Depth == "bathypelagic", ]

# Extract environmental values at sampling stations
ex.epi   <- cbind(DF_epi@data,   terra::extract(stacked_epi,   DF_epi))
ex.meso  <- cbind(DF_meso@data,  terra::extract(stacked_meso,  DF_meso))
ex.bathy <- cbind(DF_bathy@data, terra::extract(stacked_bathy, DF_bathy))

# Harmonize column names across depth layers
names(ex.epi)   <- gsub("epi",   "", names(ex.epi))
names(ex.meso)  <- gsub("meso",  "", names(ex.meso))
names(ex.bathy) <- gsub("bathy", "", names(ex.bathy))

rb1 <- rbind(ex.epi, ex.meso, ex.bathy)

# Add standardized depth label
rb1$water_column <- ifelse(
  rb1$column_depth %in% c("epipelagic", "surface"), "EPI",
  ifelse(rb1$column_depth == "mesopelagic", "MESO", "BATHY")
)

# Deduplicate: one sample per location × depth combination
set.seed(999)
rb1$tempID <- paste(rb1$long, rb1$lat, rb1$water_column)
rb2_ <- rb1 %>% group_by(tempID) %>% slice_sample(n = 1)
rb2_$water_column <- factor(rb2_$water_column, levels = c("EPI", "MESO", "BATHY"))


# ── 6. Select OMSs with ≥ 20 presences (required for AUC evaluation) ─────────
selected_otus <- names(rb2_)[grepl("OMS_", names(rb2_))]
num_obs <- rb2_ %>%
  group_by(water_column) %>%
  dplyr::summarise(across(all_of(selected_otus), sum)) %>%
  data.frame()

summary1 <- t(data.frame(lapply(num_obs[2:ncol(num_obs)], as.numeric)))
selected_otus <- names(which(sort(rowSums(summary1)) > 20))
cat(length(selected_otus), "OMSs with >20 presences selected for SDM\n")


# ── 7. Prepare h2o data frames for prediction ─────────────────────────────────
h2o.init()
rb2_[, grep("OMS_", names(rb2_))] <- lapply(rb2_[, grep("OMS_", names(rb2_))], factor)
rb2_$water_column <- factor(rb2_$water_column, levels = c("EPI", "MESO", "BATHY"))
h2o_df <- as.h2o(rb2_)

# Prediction grids (environmental rasters)
make_h2o_newdata <- function(stack, suffix, label) {
  df <- data.frame(as.data.frame(stack), water_column = label)
  names(df) <- gsub(suffix, "", names(df))
  as.h2o(df)
}
h2o.nd.epi   <- make_h2o_newdata(vars_epi,   "epi",   "EPI")
h2o.nd.meso  <- make_h2o_newdata(vars_meso,  "meso",  "MESO")
h2o.nd.bathy <- make_h2o_newdata(vars_bathy, "bathy", "BATHY")


# ── 8. Fit ensemble SDMs for each OMS ─────────────────────────────────────────
# Algorithms: GLM, Random Forest (RF), ANN, Boosted Regression Trees (GBM/BRT)
# Evaluation: 5-fold cross-validation, AUC threshold = 0.7
# Only models with AUC > 0.7 are retained in the ensemble

AUC_thr    <- 0.7
AUCC_list  <- list()
AUC_df_final <- c()
TSS_EPI    <- c()
start      <- Sys.time()

setwd(path_output)

for (otu in 1:length(selected_otus)) {

  j <- selected_otus[otu]
  cat("\nFitting models for:", j, "\n")

  # --- GLM ---
  glm_model <- h2o.glm(
    x = PREDICTOR_NAMES, y = j,
    training_frame = h2o_df,
    model_id = "GLM_Model",
    family = "binomial",
    nfolds = 5,
    standardize = FALSE
  )

  # --- Random Forest ---
  rf_model <- h2o.randomForest(
    y = j, x = PREDICTOR_NAMES,
    training_frame = h2o_df,
    ntrees = 750,
    nfolds = 5,
    model_id = "RF_model",
    categorical_encoding = "OneHotExplicit"
  )

  # --- ANN (Deep Learning) ---
  ann_model <- h2o.deeplearning(
    y = j, x = PREDICTOR_NAMES,
    training_frame = h2o_df,
    nfolds = 5,
    model_id = "ANN_model",
    standardize = FALSE
  )

  # --- Boosted Regression Trees (GBM/BRT) ---
  gbm_model <- h2o.gbm(
    y = j, x = PREDICTOR_NAMES,
    training_frame = h2o_df,
    nfolds = 5,
    model_id = "GBM_model",
    categorical_encoding = "OneHotExplicit"
  )

  MODELS   <- list(GLM = glm_model, RF = rf_model, ANN = ann_model, GBM = gbm_model)
  MPERFORM <- lapply(MODELS, h2o.performance)
  saveRDS(MPERFORM, paste0("ModelPerformance_", j, ".rds"))
  saveRDS(MODELS,   paste0("Models_",           j, ".rds"))

  # --- AUC evaluation: retain models above threshold ---
  AUCC              <- lapply(MODELS, h2o.auc)
  AUCC_list[[j]]    <- AUCC
  saveRDS(AUCC_list, "AUC_list2.rds")

  Modeval_Summary           <- data.frame(ALGO = names(MODELS), AUC = unlist(AUCC))
  models_above_threshold    <- Modeval_Summary$ALGO[Modeval_Summary$AUC >= AUC_thr]
  n_good_models             <- length(models_above_threshold)
  GOOD_MODELS               <- MODELS[names(MODELS) %in% models_above_threshold]

  if (n_good_models > 0) {

    # --- Variable importance ---
    varIMP <- list(
      GLM = h2o.varimp(glm_model)[, c(1, 4)],
      RF  = h2o.varimp(rf_model)[,  c(1, 4)],
      ANN = h2o.varimp(ann_model)[, c(1, 4)],
      GBM = h2o.varimp(gbm_model)[, c(1, 4)]
    )
    saveRDS(varIMP[names(varIMP) %in% names(GOOD_MODELS)],
            paste0("VarImp_", j, ".rds"))

    # --- Predict and generate HSI rasters (epipelagic) ---
    predictions <- lapply(GOOD_MODELS, h2o.predict, h2o.nd.epi)
    m_layer     <- vars_epi[[1]]
    l2          <- lapply(predictions, "[[", 3)
    range01     <- function(x) (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
    vals_p1     <- lapply(l2, function(x) ifelse(is.na(values(m_layer)), NA, as.numeric(as.vector(x))))
    rs          <- lapply(vals_p1, range01)

    RStack_epi <- stack()
    for (th in seq_along(models_above_threshold)) {
      temp <- m_layer
      values(temp) <- as.matrix(rs[[models_above_threshold[th]]])
      RStack_epi   <- stack(RStack_epi, temp)
    }
    names(RStack_epi) <- models_above_threshold

    # Partial plots (response curves)
    PLOT_LIST <- lapply(GOOD_MODELS, h2o.partialPlot,
                        data = h2o_df, cols = PREDICTOR_NAMES,
                        nbins = 50, plot = FALSE)
    saveRDS(PLOT_LIST, paste0(j, "_PPlot_epi.rds"))

    # Ensemble mean HSI map
    meanEPI <- mean(RStack_epi, na.rm = TRUE)
    writeRaster(meanEPI, paste0(j, "_meanModel_epi.tif"), overwrite = TRUE)

    if (n_good_models > 1) {
      writeRaster(min(RStack_epi, na.rm = FALSE), paste0(j, "_minModel_epi.tif"), overwrite = TRUE)
      writeRaster(max(RStack_epi, na.rm = FALSE), paste0(j, "_maxModel_epi.tif"), overwrite = TRUE)
      writeRaster(calc(RStack_epi, sd, na.rm = FALSE), paste0(j, "_sdModel_epi.tif"),  overwrite = TRUE)
      writeRaster(raster::cv(RStack_epi, na.rm = FALSE), paste0(j, "_CV_epi.tif"),    overwrite = TRUE)
    }

    # --- Binarize using TSS threshold ---
    p_EPI <- terra::extract(meanEPI,
      DF_bin[DF_bin@data[DF_bin@data$Depth %in% c("surface", "epipelagic"), j] == 1, ]@coords)
    a_EPI <- terra::extract(meanEPI,
      DF_bin[DF_bin@data[DF_bin@data$Depth %in% c("surface", "epipelagic"), j] == 0, ]@coords)

    if (!is.null(p_EPI) && length(na.omit(p_EPI)) > 0) {
      eval_EPI      <- dismo::evaluate(p_EPI, a_EPI)
      TSS_EPI[j]    <- dismo::threshold(eval_EPI, "spec_sens")
      meanEPIbin    <- meanEPI
      values(meanEPIbin) <- ifelse(values(meanEPIbin) >= TSS_EPI[j], 1, 0)
      writeRaster(meanEPIbin, paste0(j, "_TSS_BinModel_epi.tif"), overwrite = TRUE)
    } else {
      TSS_EPI[j] <- NA
      cat("No epipelagic presences for OMS:", j, "\n")
    }

  } else {
    cat("No models above AUC threshold for:", j, "\n")
  }
}

cat("Total runtime:", Sys.time() - start, "\n")
tss_tab <- data.frame(TSS_EPI)
write.csv(tss_tab, "TSS_TABLE.csv")


# ── 9. Visualize HSI maps ──────────────────────────────────────────────────────
# Project all OMS rasters to Hatano equal-area projection for visualization

setwd(path_output)

oms_ids <- c(1, 2, 4, 5, 6, 7, 8, 9, 10)  # OMS_3 and OMS_10 excluded (AUC < 0.7)
oms_colors <- c("#f37736", "#0392cf", "purple", "#283655", "#009688",
                "#651e3e", "#fed766", "#2ab7ca", "#b52b5e")

plot_hsi <- function(oms_id, color_scale = "H") {
  r   <- rast(paste0("OMS_", oms_id, "_meanModel_epi.tif"))
  r   <- project(r, method = "near", "+proj=hatano", mask = TRUE)
  df  <- as.data.frame(r, xy = TRUE) %>% tidyr::drop_na()
  names(df) <- c("x", "y", "HSI")
  ggplot() +
    geom_raster(aes(x = x, y = y, fill = HSI), data = df, alpha = 0.9) +
    coord_fixed(1.3) +
    viridis::scale_fill_viridis(option = color_scale, na.value = NA) +
    labs(x = NULL, y = NULL, title = paste0("OMS_", oms_id)) +
    guides(x = "none", y = "none") +
    theme_classic() +
    theme(legend.key.size = unit(0.5, "cm"),
          axis.text = element_blank(), axis.ticks = element_blank())
}

plots <- lapply(oms_ids, plot_hsi)
names(plots) <- paste0("o", oms_ids)

# Save individual PDFs
setwd(path_figs)
for (i in seq_along(oms_ids)) {
  pdf(paste0("OMS", oms_ids[i], ".pdf"))
  print(plots[[i]])
  dev.off()
}

# Combined panel (Fig. 3 – HSI maps)
ggarrange(plotlist = plots, ncol = 5, nrow = 2)


# ── 10. Variable importance plots ─────────────────────────────────────────────
setwd(path_output)
for (i in oms_ids) {
  imp <- readRDS(paste0("VarImp_OMS_", i, ".rds"))
  p_list <- lapply(names(imp), function(algo) {
    ggplot(data = imp[[algo]], aes(x = variable, y = percentage, fill = variable)) +
      geom_bar(stat = "identity") +
      ylab("Relative Importance") + xlab("") +
      ggtitle(algo) +
      scale_x_discrete(guide = guide_axis(angle = 90)) +
      theme_classic() + theme(legend.position = "none")
  })
  grid <- plot_grid(plotlist = p_list, ncol = 2, nrow = 2)
  final_plot <- ggdraw(grid) +
    draw_label(paste0("OMS_", i), size = 16, fontface = "bold", hjust = 1, vjust = -28)
  print(final_plot)
}


# ── 11. Dominant OMS map (categorical) ────────────────────────────────────────
# For each grid cell, assign the OMS with highest HSI value
setwd(path_output)
oms_stack <- stack()
for (i in oms_ids) {
  oms_stack <- stack(oms_stack, paste0("OMS_", i, "_meanModel_epi.tif"))
}
names(oms_stack) <- paste0("OMS_", oms_ids)

oms_suit <- data.frame(na.omit(rasterToPoints(oms_stack)))
dominant_oms <- apply(oms_suit[, -c(1, 2)], 1, which.max)
OMS_SUIT_r   <- rasterFromXYZ(cbind(oms_suit[, 1:2], oms_ids[dominant_oms]))
OMS_SUIT_r   <- rast(OMS_SUIT_r)
crs(OMS_SUIT_r) <- "EPSG:4326"
OMS_SUIT_r   <- project(OMS_SUIT_r, method = "near", "+proj=hatano", mask = TRUE)

levels(OMS_SUIT_r) <- data.frame(
  value    = oms_ids,
  category = paste0("OMS_", oms_ids)
)

world <- sf::st_as_sf(map("world", plot = FALSE, fill = TRUE))

ggplot() +
  geom_spatraster(data = OMS_SUIT_r) +
  scale_fill_manual(values = setNames(oms_colors, paste0("OMS_", oms_ids)),
                    na.value = "gray") +
  labs(fill = "OMS") +
  theme_minimal() +
  theme(axis.title = element_blank(), axis.text = element_blank(),
        panel.grid = element_blank())
