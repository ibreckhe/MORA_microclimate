##Script to make a brick of predictors for the microclimate analysis.
X_UTM <- raster("~/GIS/mcpred_utmx_3m.tif")
Y_UTM <- raster("~/GIS/mcpred_utmy_3m.tif")
can_vol_81m <- raster("~/GIS/mcpred_canvol_81m.tif")
can_pct_81m <- raster("~/GIS/mcpred_canpct_81m.tif")
elev <- raster("~/GIS/mcpred_elev_3m.tif")
reg_coldair_9m <- raster("~/GIS/mcpred_reg_coldair_3m.tif")
srad_seas_DJF <- raster("~/GIS/mcpred_srad_DJF_3m.tif")
srad_seas_MAM <- raster("~/GIS/mcpred_srad_MAM_3m.tif")
srad_seas_JJA <- raster("~/GIS/mcpred_srad_JJA_3m.tif")
srad_seas_SON <- raster("~/GIS/mcpred_srad_SON_3m.tif")

pred_stack <- stack(X_UTM,Y_UTM,can_vol_81m,can_pct_81m,
                    elev,reg_coldair_9m,srad_seas_DJF,
                    srad_seas_MAM,srad_seas_JJA,srad_seas_SON)
pred_brick <- brick(pred_stack,filename="~/GIS/moraclim_preds_2016.grd",
                    datatype="FLT4S",overwrite=TRUE)
