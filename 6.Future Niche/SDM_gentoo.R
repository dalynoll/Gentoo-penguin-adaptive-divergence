#rm(list=ls())
#map <- terra::vect("C:/Users/Usuario/Desktop/Datasets/World Coast/WC_merge.shp")
#plot(map)
#buffy <- buffer(map, width= 200000)
# MDE 2023 corregido
#save.image("D:/Desktop/SDM Macroecology/BioOracle/GentooBiomod.RData")


# BIOMOD2

# Cargar librerias

library(rgbif)
library(biomod2)
library(ggplot2)
library(gridExtra)
library(raster)
library(rasterVis)
library(terra)
library(tidyterra)
library(readxl)
library(usdm)

#rm(list=ls())

setwd("C:/Users/Usuario/Desktop/Gentoo revision")

#papuaraw <- read.csv("pygoscelis_papua_completo.csv")
#locs2 <- papuaraw[,c(1,24,23)]
#locs3<-round(locs2[,2:3], digits=2)
#locs4<-locs3[!duplicated(locs3),]
#write.csv(locs4, "papua2dec.csv")

papuisl <- read.csv("Datos por clado/easterndata.csv")
papupen <- read.csv("Datos por clado/Peninsuladata.csv")
papuker <- read.csv("Datos por clado/Kerguelendata.csv")
papuame <- read.csv("Datos por clado/Americandata.csv")

#papuisl <- read_excel("Occurrence 2023/P.p.taeniata.xls")
#papupen <- read.csv("Occurrence 2023/Elsworthii.csv")
#papuker <- read_excel("Occurrence 2023/P.p.kerguelensis.xls")
#papuame <- read_excel("Occurrence 2023/P.p.papua.xls")

setwd("C:/Users/Usuario/Desktop/Gentoo revision/Layers/RDA")
setwd("C:/Users/Usuario/Desktop/Gentoo revision/Layers/Futuro100")

#transformacion extra
ph_max <- brick("ph_baseline_2000_2018_depthsurf_d553_b5fc_585e_U1739542150219.nc")
phyc_max <- brick("phyc_baseline_2000_2020_depthsurf_798f_76d0_c4c4_U1739541967971max.nc")
phyc_min <- brick("phyc_baseline_2000_2020_depthsurf_ecf0_8bfb_94fe_U1739541995214.nc")
siconc_min <- brick("siconc_baseline_2000_2020_depthsurf_a494_d0bf_f062_U1739542052969.nc")
sithick_max <- brick("sithick_baseline_2000_2020_depthsurf_6399_5f09_f021_U1739542107641.nc")
salo_min <- brick("so_baseline_2000_2019_depthsurf_f31a_74fa_d7d6_U1739541871045min.nc")
salo_max <- brick("so_baseline_2000_2019_depthsurf_7d26_e68e_7207_U1739541861544max.nc")
swinds_max <- brick("sws_baseline_2000_2019_depthsurf_fe33_6d18_e9c6_U1739542136173.nc")
them_min <- brick("thetao_baseline_2000_2019_depthsurf_f77e_f6de_17a4_U1739541333215min.nc")
them_max <- brick("thetao_baseline_2000_2019_depthsurf_c4a3_c6c0_7df5_U1739541311989max.nc")
bathy<- brick("terrain_characteristics_6d78_4f59_9e1f_U1741947669689.nc")

setwd("C:/Users/Usuario/Desktop/Gentoo revision")


#transformar formato
chl100 <- brick("chl_ssp370_2020_2100.nc")
ph100 <- brick("ph_ssp370_2020_2100.nc")
phyc100 <- brick("phyc_ssp370_2020_2100.nc")
so100 <- brick("so_ssp370_2020_2100.nc")
sws100 <- brick("sws_ssp370_2020_2100.nc")
thetao100 <- brick("thetao_ssp370_2020_2100.nc")

writeRaster(chl100,"chl100.tif",bylayer=TRUE)
writeRaster(ph100,"ph100.tif",bylayer=TRUE)
writeRaster(phyc100,"phyc100.tif",bylayer=TRUE)
writeRaster(so100,"so100.tif",bylayer=TRUE)
writeRaster(sws100,"sws100.tif",bylayer=TRUE)
writeRaster(thetao100,"thetao100.tif",bylayer=TRUE)
writeRaster(bathy, "zbathymetry.tif", bylayer=TRUE)


writeRaster(ph_max,"ph_max.tif",bylayer=TRUE)
writeRaster(siconc_min,"siconc_min.tif",bylayer=TRUE)
writeRaster(sithick_max,"sithick_max.tif",bylayer=TRUE)
writeRaster(swinds_max,"swinds_max.tif",bylayer=TRUE)
writeRaster(them_min,"temp_min.tif",bylayer=TRUE, overwrite=TRUE)
writeRaster(them_max,"temp_max.tif",bylayer=TRUE, overwrite=TRUE)
writeRaster(phyc_max,"prod_max.tif",bylayer=TRUE)
writeRaster(phyc_min,"prod_min.tif",bylayer=TRUE)
writeRaster(salo_min,"sal_min.tif",bylayer=TRUE)
writeRaster(salo_max,"sal_max.tif",bylayer=TRUE)

#geospatial info for RDA analyses
RDA.list <- list.files(path="C:/Users/Usuario/Desktop/Gentoo revision/Layers/RDA", 
                           pattern =".tif", full.names=TRUE)

sampledcolonies <- read.delim("C:/Users/Usuario/Desktop/Gentoo revision/Layers/RDA/sampledcolonies.txt")

rda.stack<- stack(RDA.list)

rda.scores <- raster::extract(rda.stack, sampledcolonies[,c("lon","lat")])
row.names(rda.scores) <- sampledcolonies[,1]

write.csv(rda.scores, "GentooRDAscores.csv")

# cortar las variables con el area definida

current.list <- list.files(path="C:/Users/Usuario/Desktop/Gentoo revision/Layers/Capas presente", 
                             pattern =".tif", full.names=TRUE)

buffer <- raster("C:/Users/Usuario/Desktop/Datasets/mascara.tif")

c.stack<- stack(current.list)
presente<- stack(current.list)

#plot(presente)[[1]]

studyArea = extent(-180,180,-80,-30)
SOArea <- crop(presente, studyArea) 

studyArea = extent(-180,180,-80,-30)
buffer2 <- crop(buffer, studyArea) 

SOArea_c <- mask(SOArea, buffer2) 

#plot(SOArea_c$chlorophyll)
#plot(SOArea_c$zbathymetry_1)

##variable inflation test

#elsworthii <- read_excel("Elsworthii.xls")
#kerguelensis <- read_excel("P.p.kerguelensis.xls")
#papua <- read_excel("P.p.papua.xls")
#taeniata <- read_excel("P.p.taeniata.xls")

papuant <- read.csv("Datos por clado/Peninsuladata.csv")
papuker <- read.csv("Datos por clado/Kerguelendata.csv")
papuame <- read.csv("Datos por clado/Americandata.csv")
papuisl <- read.csv("Datos por clado/easterndata.csv")

elsworthii_clim <- raster::extract(presente,
                                   papuant[,c("x","y")])
kerguelensis_clim <- raster::extract(presente,
                                     papuker[,c("x","y")])
papua_clim <- raster::extract(presente,
                              papuame[,c("x","y")])
taeniata_clim <- raster::extract(presente,
                                 papuisl[,c("x","y")])

elsworthii_clim <- as.data.frame(elsworthii_clim)
kerguelensis_clim <- as.data.frame(kerguelensis_clim)
papua_clim <- as.data.frame(papua_clim)
taeniata_clim <- as.data.frame(taeniata_clim)

vif_els <- vifstep(elsworthii_clim, th = 11)
vif_ker <- vifstep(kerguelensis_clim, th = 11)
vif_pap <- vifstep(papua_clim, th = 11)
vif_tae <- vifstep(taeniata_clim, th = 11)

vif_els
vif_ker
vif_pap
vif_tae


SOArea_els_a <- SOArea_c[[c(1:4)]]
SOArea_ker_a <- SOArea_c[[c(1:4)]]
SOArea_pap_a <- SOArea_c[[c(1:4)]]
SOArea_isl_a <- SOArea_c[[c(1:4)]]

#buffer_poly  <- vect("Datos por clado/buffer_poly.shp")
buffer <- raster("C:/Users/Usuario/Desktop/Datasets/mascara.tif")
#studyArea = extent(-180,180,-80,-30)
#buffer_z <- crop(buffer, studyArea)  

buffer_c <- crop(buffer, SOArea_c)  
SOArea_c <- mask(SOArea, buffer2) 

plot(buffer)

plot(SOArea_c$o3_salinity)

#SOArea_els <- crop(SOArea_els_a, buffer_z, snap="in", mask=TRUE, extend=FALSE)
#SOArea_els <- mask(SOArea_els_a, buffer3, snap="in")
#plot(SOArea_els$pH)
#SOArea_ker <- mask(SOArea_ker_a, buffer3, snap="in")
#SOArea_pap <- mask(SOArea_pap_a, buffer3, snap="in")
#SOArea_isl <- mask(SOArea_isl_a, buffer3, snap="in")
# identificar los ids donde existen presencia de la especie
points_papuker<-data.frame(papuker[,c("x", "y")])
points_papupen<-data.frame(papuant[,c("x", "y")])
points_papuisl<-data.frame(papuisl[,c("x", "y")])
points_papuame<-data.frame(papuame[,c("x", "y")])
papuker_cell_id <- cellFromXY(subset(c.stack,1), points_papuker)
papupen_cell_id <- cellFromXY(subset(c.stack,1), points_papupen)
papuisl_cell_id <- cellFromXY(subset(c.stack,1), points_papuisl)
papuame_cell_id <- cellFromXY(subset(c.stack,1), points_papuame)


dev.off()
#SPECIFIC STUDY AREAS
studyAreaC1.1 = extent(30,60,-55,-35)
SOArea1.1 <- crop(SOArea_isl_a, studyAreaC1.1)
#plot(SOArea1.1[[1]])

studyAreaC1.2 = extent(150,170,-58,-48)
SOArea1.2 <- crop(SOArea_isl_a, studyAreaC1.2)
#plot(SOArea1.2[[1]])

studyAreaC2 = extent(65,85,-60,-43)
SOArea2 <- crop(SOArea_ker_a, studyAreaC2)
#plot(SOArea2[[1]])

studyAreaC3.1 = extent(-50,-20,-65,-50)
SOArea3.1 <- crop(SOArea_els_a, studyAreaC3.1)
#plot(SOArea3.1[[1]])

studyAreaC3.2 = extent(-77,-50,-77,-57)
SOArea3.2 <- crop(SOArea_els_a, studyAreaC3.2)
#plot(SOArea3.2[[1]])

studyAreaC4 = extent(-73,-50,-57,-47)
SOArea4 <- crop(SOArea_pap_a, studyAreaC4)
#plot(SOArea4[[1]])

SOArea1 <- mosaic(SOArea1.1, SOArea1.2, fun=mean)
names(SOArea1) <- c("o1_pH", "o2_prim.productivity", "o3_salinity", "o4_temperature")

#plot(SOArea1$Bathymetry)


names(SOArea2) <- c("o1_pH", "o2_prim.productivity", "o3_salinity", "o4_temperature")
plot(SOArea2$Bathymetry)

SOArea3 <- mosaic(SOArea3.1, SOArea3.2, fun=mean)
names(SOArea3) <- c("o1_pH", "o2_prim.productivity", "o3_salinity", "o4_temperature")

SOArea4
names(SOArea4) <- c("o1_pH", "o2_prim.productivity", "o3_salinity", "o4_temperature")

# plot(SOArea4$IceThickness)

SOArea1 <- stack(SOArea1)
SOArea2 <- stack(SOArea2)
SOArea3 <- stack(SOArea3)
SOArea4 <- stack(SOArea4)

plot(SOArea4)

# sub selección de variables

#bioclim_ZA_sub <- stack(subset(bioclim_ZA, c("bio_5", "bio_7", "bio_11", "bio_19")))

# empezar a correr BIOMOD, formateando los datos

papuker_data <- BIOMOD_FormatingData(resp.var = rep(1, nrow( points_papuker ) ),
                                    expl.var = SOArea2,
                                    resp.xy = points_papuker[,c('x', 'y')],
                                    resp.name = "Papua kerguelensis",
                                    PA.nb.rep = 3,
                                    PA.nb.absences = 500,
                                    PA.strategy = 'random')

#plot(papuker_data)

pdf("Papua_Kerguelen_Dataset.pdf")
plot(papuker_data)
dev.off()

papupen_data <- BIOMOD_FormatingData(resp.var = rep(1, nrow( points_papupen ) ),
                                     expl.var = SOArea3,
                                     resp.xy = points_papupen[,c('x', 'y')],
                                     resp.name = "Papua elsworthii",
                                     PA.nb.rep = 3,
                                     PA.nb.absences = 500,
                                     PA.strategy = 'random')

#plot(papupen_data)

pdf("Papua_Antarctic_Dataset.pdf")
plot(papupen_data)
dev.off()


papuisl_data <- BIOMOD_FormatingData(resp.var = rep(1, nrow( points_papuisl ) ),
                                     expl.var = SOArea1,
                                     resp.xy = points_papuisl[,c('x', 'y')],
                                     resp.name = "Papua taeniata",
                                     PA.nb.rep = 3,
                                     PA.nb.absences = 500,
                                     PA.strategy = 'random')

pdf("Papua_Islands_Dataset.pdf")
plot(papuisl_data)
dev.off()

#plot(papuisl_data)
#plot(SOArea$chlorophyll)
#points(points_papuisl[,c('x', 'y')])
#plot(points)

papuame_data <- BIOMOD_FormatingData(resp.var = rep(1, nrow( points_papuame ) ),
                                     expl.var = SOArea4,
                                     resp.xy = points_papuame[,c('x', 'y')],
                                     resp.name = "Papua papua",
                                     PA.nb.rep = 3,
                                     PA.nb.absences = 500,
                                     PA.strategy = 'random')

pdf("Papua_American_Dataset.pdf")
plot(papuame_data)
dev.off()


#plot(SOArea$primary.productivity)
#points(points_papupen[,c('lon', 'lat')])

# resumen del objeto de formato BIOMOD
#plot(papupen_data)


# plot de las pseudo-ausencias seleccionadas
# 
# plot(papuker_data)
# plot(papupen_data)
# plot(papuisl_data)
# plot(papuame_data)

# BIOMOD, opciones de modelación, (en este caso para tres técnicas)

ProLau_opt <- BIOMOD_ModelingOptions(GLM = list( type = 'quadratic',
                                                 interaction.level = 1 ),
                                     GBM = list( n.trees = 1000 ),
                                     GAM = list( algo = 'GAM_mgcv' ) )

# BIOMOD, calibración con cuatro técnicas

papuker_models <- BIOMOD_Modeling( bm.format = papuker_data,
                                  models = c('GLM', 'MAXNET', 'GBM', 'RF', 'GAM'),
                                  bm.options = ProLau_opt,
                                  CV.nb.rep = 5,
                                  CV.perc = 0.7,
                                  var.import = 1,
                                  do.full.models = F,
                                  modeling.id = "ex2" )

papuisl_models <- BIOMOD_Modeling( bm.format = papuisl_data,
                                   models = c('GLM', 'MAXNET', 'GBM', 'RF', 'GAM'),
                                   bm.options = ProLau_opt,
                                   CV.nb.rep = 5,
                                   CV.perc = 0.7,
                                   var.import = 1,
                                   do.full.models = F,
                                   modeling.id = "ex2" )

papupen_models <- BIOMOD_Modeling( bm.format = papupen_data,
                                   models = c('GLM', 'MAXNET', 'GBM', 'RF', 'GAM'),
                                   bm.options = ProLau_opt,
                                   CV.nb.rep = 5,
                                   CV.perc = 0.7,
                                   var.import = 1,
                                   do.full.models = F,
                                   modeling.id = "ex2" )

papuame_models <- BIOMOD_Modeling( bm.format = papuame_data,
                                   models = c('GLM', 'MAXNET', 'GBM', 'RF', 'GAM'),
                                   bm.options = ProLau_opt,
                                   CV.nb.rep = 5,
                                   CV.perc = 0.7,
                                   var.import = 1,
                                   do.full.models = F,
                                   modeling.id = "ex2" )


# evaluación de los modelos

papuker_models_scores <- get_evaluations(papuker_models)
papupen_models_scores <- get_evaluations(papupen_models)
papuisl_models_scores <- get_evaluations(papuisl_models)
papuame_models_scores <- get_evaluations(papuame_models)

dim(papuker_models_scores)

ker <- bm_PlotEvalBoxplot(bm.out = papuker_models, main="kerguelen", dataset = "validation", group.by = c('algo','algo')) 
kerplot <- png(filename = "kerguelen validation.png")
ker
dev.off()

bm_PlotEvalMean(bm.out = papuker_models, dataset = "validation", group.by = c('algo','algo'))

pen <- bm_PlotEvalBoxplot(bm.out = papupen_models, main="antarctic", dataset = "validation", group.by = c('algo','algo'))
penplot <- png(filename = "antarctic validation.png")
pen
dev.off()

bm_PlotEvalMean(bm.out = papupen_models, dataset = "validation", group.by = c('algo','algo'))

ame <- bm_PlotEvalBoxplot(bm.out = papuame_models, main="American", dataset = "validation", group.by = c('algo','algo'))
ameplot <- png(filename = "American validation.png")
ame
dev.off()

bm_PlotEvalMean(bm.out = papuame_models, dataset = "validation", group.by = c('algo','algo'))

isl <- bm_PlotEvalBoxplot(bm.out = papuisl_models, main="islandic", dataset = "validation", group.by = c('algo','algo'))
islplot <- png(filename = "islandic validation.png")
isl
dev.off()

bm_PlotEvalMean(bm.out = papuisl_models, dataset = "validation", group.by = c('algo','algo'))


###dataset = "calibration"
####by = "models" , metrics = c("ROC","TSS"),
####by = "cv_run" , metrics = c("ROC","TSS"),
####by = "data_set" , metrics = c("ROC","TSS"),

papuker_var_import <- get_variables_importance(papuker_models)
papupen_var_import <- get_variables_importance(papupen_models)

papuker_var_import <- get_variables_importance(papuisl_models)
papupen_var_import <- get_variables_importance(papuame_models)

# Represent Variable Importance
pdf("Papua_Kerguelen_Variable_Importance_explVar_algo_algo_2.pdf")
bm_PlotVarImpBoxplot(bm.out = papuker_models, group.by = c('expl.var', 'algo', 'algo'))
dev.off()

pdf("Papua_Antarctic_Variable_Importance_explVar_algo_algo_2.pdf")
bm_PlotVarImpBoxplot(bm.out = papupen_models, group.by = c('expl.var', 'algo', 'algo'))
dev.off()

pdf("Papua_Islands_Variable_Importance_explVar_algo_algo_2.pdf")
bm_PlotVarImpBoxplot(bm.out = papuisl_models, group.by = c('expl.var', 'algo', 'algo'))
dev.off()

pdf("Papua_America_Variable_Importance_explVar_algo_algo_2.pdf")
bm_PlotVarImpBoxplot(bm.out = papuame_models, group.by = c('expl.var', 'algo', 'algo'))
dev.off()



# importancia de variables por técnica

#apply(ProLau_models_var_import, c(1,2), mean)

##Problema por tecnica

papuker_glm <- BIOMOD_LoadModels(papuker_models, algo	='GLM')
papuker_gbm <- BIOMOD_LoadModels(papuker_models, algo ='GBM')
papuker_rf  <- BIOMOD_LoadModels(papuker_models, algo ='RF')
papuker_gam <- BIOMOD_LoadModels(papuker_models, algo	='GAM')

# curvas de respuesta por técnica

# glm_eval_strip <- biomod2::bm_PlotResponseCurves(
#   bm.out = papuker_models,
#   Data = get_formal_data(papuker_models,'expl.var'),
#   show.variables= get_formal_data(papuker_models,'expl.var.names'),
#   do.bivariate = FALSE,
#   fixed.var.metric = 'median',
#   legend = FALSE,
#   display_title = FALSE,
#   data_species = get_formal_data(papuker_models,'resp.var'))
# 
# gbm_eval_strip <- biomod2::bm_PlotResponseCurves(
#   bm.out = ProLau_models,
#   models.chosen = ProLau_gbm,
#   Data = get_formal_data(ProLau_models,'expl.var'),
#   show.variables= get_formal_data(ProLau_models,'expl.var.names'),
#   do.bivariate = FALSE,
#   fixed.var.metric = 'median',
#   legend = FALSE,
#   display_title = FALSE,
#   data_species = get_formal_data(ProLau_models,'resp.var'))
# 
# rf_eval_strip <- biomod2::bm_PlotResponseCurves(
#   bm.out = ProLau_models,
#   models.chosen = ProLau_rf,
#   Data = get_formal_data(ProLau_models,'expl.var'),
#   show.variables= get_formal_data(ProLau_models,'expl.var.names'),
#   do.bivariate = FALSE,
#   fixed.var.metric = 'median',
#   legend = FALSE,
#   display_title = FALSE,
#   data_species = get_formal_data(ProLau_models,'resp.var'))

# gam_eval_strip <- biomod2::bm_PlotResponseCurves(
#   bm.out = ProLau_models,
#   models.chosen = ProLau_gam,
#   Data = get_formal_data(ProLau_models,'expl.var'),
#   show.variables= get_formal_data(ProLau_models,'expl.var.names'),
#   do.bivariate = FALSE,
#   fixed.var.metric = 'median',
#   legend = FALSE,
#   display_title = FALSE,
#   data_species = get_formal_data(ProLau_models,'resp.var'))

# BIOMOD, ensamble de técnicas

papuker_ensemble_models <- BIOMOD_EnsembleModeling( bm.mod = papuker_models,
                                                   models.chosen = 'all',
                                                   em.by = 'all',
                                                   metric.select = c('TSS','ROC'),
                                                   metric.select.thres = c(0.6, 0.8),
                                                   metric.eval = c('TSS','ROC'),
                                                   em.algo = c('EMcv', 'EMca', 'EMwmean'),
                                                   var.import = 2 )

papuker_ensemble_models_scores <- get_evaluations(papuker_ensemble_models)

sppname2 <- "Papuker"
ensem_eval_file_path <- file.path("C:/Users/Usuario/Desktop/Gentoo Revision/Evaluations",paste0(sppname2,"_emsemble_evaluation.txt"))
capture.output(get_evaluations(papuker_ensemble_models),
               file = ensem_eval_file_path)

papupen_ensemble_models <- BIOMOD_EnsembleModeling( bm.mod = papupen_models,
                                                    models.chosen = 'all',
                                                    em.by = 'all',
                                                    metric.select = c('TSS','ROC'),
                                                    metric.select.thres = c(0.6, 0.8),
                                                    metric.eval = c('TSS','ROC'),
                                                    em.algo = c('EMcv', 'EMca', 'EMwmean'),
                                                    var.import = 2 )

papupen_ensemble_models_scores <- get_evaluations(papupen_ensemble_models)

sppname3 <- "Papupen"
ensem_eval_file_path <- file.path("C:/Users/Usuario/Desktop/Gentoo Revision/Evaluations",paste0(sppname3,"_emsemble_evaluation.txt"))
capture.output(get_evaluations(papupen_ensemble_models),
               file = ensem_eval_file_path)

papuisl_ensemble_models <- BIOMOD_EnsembleModeling( bm.mod = papuisl_models,
                                                    models.chosen = 'all',
                                                    em.by = 'all',
                                                    metric.select = c('TSS','ROC'),
                                                    metric.select.thres = c(0.6, 0.8),
                                                    metric.eval = c('TSS','ROC'),
                                                    em.algo = c('EMcv', 'EMca', 'EMwmean'),
                                                    var.import = 2 )

papuisl_ensemble_models_scores <- get_evaluations(papuisl_ensemble_models)

#Ensemble evaluation
sppname1 <- "Papuisl"
ensem_eval_file_path <- file.path("C:/Users/Usuario/Desktop/Gentoo Revision/Evaluations",
                                  paste0(sppname1,"_emsemble_evaluation.txt"))
capture.output(get_evaluations(papuisl_ensemble_models),
               file = ensem_eval_file_path)

papuame_ensemble_models <- BIOMOD_EnsembleModeling( bm.mod = papuame_models,
                                                    models.chosen = 'all',
                                                    em.by = 'all',
                                                    metric.select = c('TSS','ROC'),
                                                    metric.select.thres = c(0.6, 0.8),
                                                    metric.eval = c('TSS','ROC'),
                                                    em.algo = c('EMcv', 'EMca', 'EMwmean'),
                                                    var.import = 2 )

papuame_ensemble_models_scores <- get_evaluations(papuame_ensemble_models)

#Ensemble evaluation
sppname4 <- "Papuame"
ensem_eval_file_path <- file.path("C:/Users/Usuario/Desktop/Gentoo Revision/Evaluations",
                                  paste0(sppname4,"_emsemble_evaluation.txt"))
capture.output(get_evaluations(papuame_ensemble_models),
               file = ensem_eval_file_path)


# BIOMOD, proyección actual

papuker_models_proj_current <- BIOMOD_Projection( bm.mod = papuker_models,
                                                 new.env = SOArea2,
                                                 proj.name = "current",
                                                 binary.meth = "TSS",
                                                 output.format = ".img",
                                                 do.stack = FALSE )

# BIOMOD, proyección ensamble actual

papuker_Ensemble_Proj <- BIOMOD_EnsembleForecasting(bm.em = papuker_ensemble_models,
                                                          proj.name = 'current ensemble',
                                                          new.env = SOArea2,
                                                          models.chosen = 'all',
                                                          metric.binary = 'all',
                                                          metric.filter = 'all',
                                                          binary.meth = 'TSS',
                                                          output.format = '.img',
                                                          do.stack = FALSE)


papupen_Ensemble_Proj <- BIOMOD_EnsembleForecasting(bm.em = papupen_ensemble_models,
                                                    proj.name = 'current ensemble',
                                                    new.env = SOArea3,
                                                    models.chosen = 'all',
                                                    metric.binary = 'all',
                                                    metric.filter = 'all',
                                                    binary.meth = 'TSS',
                                                    output.format = '.img',
                                                    do.stack = FALSE)


papuisl_Ensemble_Proj <- BIOMOD_EnsembleForecasting(bm.em = papuisl_ensemble_models,
                                                    proj.name = 'current ensemble',
                                                    new.env = SOArea1,
                                                    models.chosen = 'all',
                                                    metric.binary = 'all',
                                                    metric.filter = 'all',
                                                    binary.meth = 'TSS',
                                                    output.format = '.img',
                                                    do.stack = FALSE)

papuame_Ensemble_Proj <- BIOMOD_EnsembleForecasting(bm.em = papuame_ensemble_models,
                                                    proj.name = 'current ensemble',
                                                    new.env = SOArea4,
                                                    models.chosen = 'all',
                                                    metric.binary = 'all',
                                                    metric.filter = 'all',
                                                    binary.meth = 'TSS',
                                                    output.format = '.img',
                                                    do.stack = FALSE)



  
#plot(ProLau_models_proj_current)
#plot(papuker_Ensemble_Proj)
#plot(papuker_Ensemble_Proj)

# proyección futuras, 2050
# cortar con el área

#futur1.list <- list.files(path="C:/Users/Usuario/Desktop/Gentoo revision/Layers/Capas futuro1", 
                           #pattern =".tif", full.names=TRUE)
#futur2.list <- list.files(path="C:/Users/Usuario/Desktop/Gentoo revision/Layers/Capas futuro2", 
                          #pattern =".tif", full.names=TRUE)
#futur3.list <- list.files(path="C:/Users/Usuario/Desktop/Gentoo revision/Layers/Capas futuro3", 
                          #pattern =".tif", full.names=TRUE)
#futur4.list <- list.files(path="C:/Users/Usuario/Desktop/Gentoo revision/Layers/Capas futuro4", 
                          #pattern =".tif", full.names=TRUE)

#f1.stack_2050<- stack(futur1.list)
#f2.stack_2050<- stack(futur2.list)
#f3.stack_2050<- stack(futur3.list)
#f4.stack_2050<- stack(futur4.list)

#SOAreaf1_2050 <- crop(f1.stack_2050, SOArea1)  
#SOAreaf2_2050 <- crop(f2.stack_2050, SOArea2)  
#SOAreaf3_2050 <- crop(f3.stack_2050, SOArea3)  
#SOAreaf4_2050 <- crop(f4.stack_2050, SOArea4)  

#SOAreaf1_2050 <- mask(SOAreaf1_2050, SOArea1, snap="in")
#SOAreaf2_2050 <- mask(SOAreaf2_2050, SOArea2, snap="in")
#SOAreaf3_2050 <- mask(SOAreaf3_2050, SOArea3, snap="in")
#SOAreaf4_2050 <- mask(SOAreaf4_2050, SOArea4, snap="in")


futur1_50.list <- list.files(path="C:/Users/Usuario/Desktop/Gentoo revision/Layers/Futuro50_1", 
                          pattern =".tif", full.names=TRUE)
futur2_50.list <- list.files(path="C:/Users/Usuario/Desktop/Gentoo revision/Layers/Futuro50_2", 
                          pattern =".tif", full.names=TRUE)
futur3_50.list <- list.files(path="C:/Users/Usuario/Desktop/Gentoo revision/Layers/Futuro50_3", 
                          pattern =".tif", full.names=TRUE)
futur4_50.list <- list.files(path="C:/Users/Usuario/Desktop/Gentoo revision/Layers/Futuro50_4", 
                          pattern =".tif", full.names=TRUE)

f1.stack_2050<- stack(futur1_50.list)
f2.stack_2050<- stack(futur2_50.list)
f3.stack_2050<- stack(futur3_50.list)
f4.stack_2050<- stack(futur4_50.list)

SOAreaf1_2050 <- crop(f1.stack_2050, SOArea1)  
SOAreaf2_2050 <- crop(f2.stack_2050, SOArea2)  
SOAreaf3_2050 <- crop(f3.stack_2050, SOArea3)  
SOAreaf4_2050 <- crop(f4.stack_2050, SOArea4)  

SOAreaf1_2050 <- mask(SOAreaf1_2050, SOArea1, snap="in")
SOAreaf2_2050 <- mask(SOAreaf2_2050, SOArea2, snap="in")
SOAreaf3_2050 <- mask(SOAreaf3_2050, SOArea3, snap="in")
SOAreaf4_2050 <- mask(SOAreaf4_2050, SOArea4, snap="in")

names(SOAreaf1_2050) <- c("o1_pH", "o2_prim.productivity", "o3_salinity", "o4_temperature")
names(SOAreaf2_2050) <- c("o1_pH", "o2_prim.productivity", "o3_salinity", "o4_temperature")
names(SOAreaf3_2050) <- c("o1_pH", "o2_prim.productivity", "o3_salinity", "o4_temperature")
names(SOAreaf4_2050) <- c("o1_pH", "o2_prim.productivity", "o3_salinity", "o4_temperature")

#plot(SOAreaf3_2050)

#plot(SOArea2_2050)[1]

futur1.list <- list.files(path="C:/Users/Usuario/Desktop/Gentoo revision/Layers/Futuro100_1", 
                          pattern =".tif", full.names=TRUE)
futur2.list <- list.files(path="C:/Users/Usuario/Desktop/Gentoo revision/Layers/Futuro100_2", 
                          pattern =".tif", full.names=TRUE)
futur3.list <- list.files(path="C:/Users/Usuario/Desktop/Gentoo revision/Layers/Futuro100_3", 
                          pattern =".tif", full.names=TRUE)
futur4.list <- list.files(path="C:/Users/Usuario/Desktop/Gentoo revision/Layers/Futuro100_4", 
                          pattern =".tif", full.names=TRUE)

f1.stack_2100<- stack(futur1.list)
f2.stack_2100<- stack(futur2.list)
f3.stack_2100<- stack(futur3.list)
f4.stack_2100<- stack(futur4.list)

SOAreaf1_2100 <- crop(f1.stack_2100, SOArea1)  
SOAreaf2_2100 <- crop(f2.stack_2100, SOArea2)  
SOAreaf3_2100 <- crop(f3.stack_2100, SOArea3)  
SOAreaf4_2100 <- crop(f4.stack_2100, SOArea4)  

SOAreaf1_2100 <- mask(SOAreaf1_2100, SOArea1, snap="in")
SOAreaf2_2100 <- mask(SOAreaf2_2100, SOArea2, snap="in")
SOAreaf3_2100 <- mask(SOAreaf3_2100, SOArea3, snap="in")
SOAreaf4_2100 <- mask(SOAreaf4_2100, SOArea4, snap="in")

names(SOAreaf1_2100) <- c("o1_pH", "o2_prim.productivity", "o3_salinity", "o4_temperature")
names(SOAreaf2_2100) <- c("o1_pH", "o2_prim.productivity", "o3_salinity", "o4_temperature")
names(SOAreaf3_2100) <- c("o1_pH", "o2_prim.productivity", "o3_salinity", "o4_temperature")
names(SOAreaf4_2100) <- c("o1_pH", "o2_prim.productivity", "o3_salinity", "o4_temperature")

#SOArea3_2050 <- crop(c.stack_2050, SOArea3)  
#SOArea3_2050 <- mask(SOArea3_2050, SOArea3)  

#plot(SOArea3_2050)[1]



#SOArea1_2050 <- crop(c.stack_2050, SOArea1)  
#SOArea1_2050 <- mask(SOArea1_2050, SOArea1)  

#plot(SOArea1_2050)[1]

#names(SOArea1_2050) <- c("Chlorophyll", "Velocity", "IceThickness","Salinity", "Temperature")

#SOArea4_2050 <- crop(c.stack_2050, SOArea4)  

#plot(SOArea4_2050)[1]

#names(SOArea4_2050) <- c("Chlorophyll", "Velocity", "IceThickness","Salinity", "Temperature")

papuker_models_proj_2050_BC45 <- BIOMOD_EnsembleForecasting(bm.em = papuker_ensemble_models,
                                                            proj.name = 'Papuker 2050 ensemble',
                                                            new.env = SOAreaf2_2050,
                                                            models.chosen = 'all',
                                                            metric.binary = 'all',
                                                            metric.filter = 'all',
                                                            binary.meth = 'TSS',
                                                            output.format = '.img',
                                                            do.stack = FALSE)

papupen_models_proj_2050_BC45 <- BIOMOD_EnsembleForecasting(bm.em = papupen_ensemble_models,
                                                            proj.name = 'Papupen 2050 ensemble',
                                                            new.env = SOAreaf3_2050,
                                                            models.chosen = 'all',
                                                            metric.binary = 'all',
                                                            metric.filter = 'all',
                                                            binary.meth = 'TSS',
                                                            output.format = '.img',
                                                            do.stack = FALSE)

papuisl_models_proj_2050_BC45 <- BIOMOD_EnsembleForecasting(bm.em = papuisl_ensemble_models,
                                                            proj.name = 'Papuisl 2050 ensemble',
                                                            new.env = SOAreaf1_2050,
                                                            models.chosen = 'all',
                                                            metric.binary = 'all',
                                                            metric.filter = 'all',
                                                            binary.meth = 'TSS',
                                                            output.format = '.img',
                                                            do.stack = FALSE)

papuame_models_proj_2050_BC45 <- BIOMOD_EnsembleForecasting(bm.em = papuame_ensemble_models,
                                                            proj.name = 'Papuame 2050 ensemble',
                                                            new.env = SOAreaf4_2050,
                                                            models.chosen = 'all',
                                                            metric.binary = 'all',
                                                            metric.filter = 'all',
                                                            binary.meth = 'TSS',
                                                            output.format = '.img',
                                                            do.stack = FALSE)

# BIOMOD, proyección ensamble futura 2100

papuker_models_proj_2100_BC45 <- BIOMOD_EnsembleForecasting(bm.em = papuker_ensemble_models,
                                                           proj.name = 'Papuker 2100 ensemble',
                                                           new.env = SOAreaf2_2100,
                                                           models.chosen = 'all',
                                                           metric.binary = 'all',
                                                           metric.filter = 'all',
                                                            binary.meth = 'TSS',
                                                            output.format = '.img',
                                                            do.stack = FALSE)

papupen_models_proj_2100_BC45 <- BIOMOD_EnsembleForecasting(bm.em = papupen_ensemble_models,
                                                            proj.name = 'Papupen 2100 ensemble',
                                                            new.env = SOAreaf3_2100,
                                                            models.chosen = 'all',
                                                            metric.binary = 'all',
                                                            metric.filter = 'all',
                                                            binary.meth = 'TSS',
                                                            output.format = '.img',
                                                            do.stack = FALSE)

papuisl_models_proj_2100_BC45 <- BIOMOD_EnsembleForecasting(bm.em = papuisl_ensemble_models,
                                                            proj.name = 'Papuisl 2100 ensemble',
                                                            new.env = SOAreaf1_2100,
                                                            models.chosen = 'all',
                                                            metric.binary = 'all',
                                                            metric.filter = 'all',
                                                            binary.meth = 'TSS',
                                                            output.format = '.img',
                                                            do.stack = FALSE)

papuame_models_proj_2100_BC45 <- BIOMOD_EnsembleForecasting(bm.em = papuame_ensemble_models,
                                                            proj.name = 'Papuame 2100 ensemble',
                                                            new.env = SOAreaf4_2100,
                                                            models.chosen = 'all',
                                                            metric.binary = 'all',
                                                            metric.filter = 'all',
                                                            binary.meth = 'TSS',
                                                            output.format = '.img',
                                                            do.stack = FALSE)


#plot(ProLau_models_proj_2050_BC45)

# BIOMOD, proyección ensamble futura 2070
# BIOMOD, proyección ensamble futura 2070
# plot del ensamble
#plot(ProLau_models_proj_2070_BC45,
 #    str.grep = "EMca")
# Error tipeo PProLau, corregido
#plot(ProLau_models_proj_2070_BC45,
 #    str.grep = "EMca|EMwmean")
## you should use get_predictions wrapper 
#mod_proj <- get_predictions(ProLau_models_proj_2050_BC45)
#mod_proj
#names(mod_proj)
## plot only the 2 layer and rename it
#plot(subset(mod_proj,2), main = "RF projections")
# Error saque linea incompleta

######################################
####ESCENARIOS DE FUTURO##############
######################################

##escenario 2050
# cargar proyecciones binarias

papuker_bin_proj_current <- stack(
  c( ca = "Papua.kerguelensis/proj_current ensemble/individual_projections/Papua.kerguelensis_EMcaByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.img",
     wm = "Papua.kerguelensis/proj_current ensemble/individual_projections/Papua.kerguelensis_EMwmeanByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.img") )

papuker_bin_proj_2050_BC45 <- stack(
  c( ca = "Papua.kerguelensis/proj_Papuker 2050 ensemble/individual_projections/Papua.kerguelensis_EMcaByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.img",
     wm = "Papua.kerguelensis/proj_Papuker 2050 ensemble/individual_projections/Papua.kerguelensis_EMwmeanByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.img") )

#papuker_bin_proj_2070_BC45 <- stack(
# c( ca = "Protea.laurifolia/proj_2070 ensemble/individual_projections/Papua.kerguelensis_EMcaByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.img",
#   wm = "Protea.laurifolia/proj_2070 ensemble/individual_projections/Papua.kerguelensis_EMwmeanByROC_mergedData_mergedRun_mergedAlgo_TSSbin.img") )

papupen_bin_proj_current <- stack(
  c( ca = "Papua.elsworthii/proj_current ensemble/individual_projections/Papua.elsworthii_EMcaByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.img",
     wm = "Papua.elsworthii/proj_current ensemble/individual_projections/Papua.elsworthii_EMwmeanByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.img") )

papupen_bin_proj_2050_BC45 <- stack(
  c( ca = "Papua.elsworthii/proj_Papupen 2050 ensemble/individual_projections/Papua.elsworthii_EMcaByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.img",
     wm = "Papua.elsworthii/proj_Papupen 2050 ensemble/individual_projections/Papua.elsworthii_EMwmeanByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.img") )

#papupen_bin_proj_2070_BC45 <- stack(
# c( ca = "Protea.laurifolia/proj_2070 ensemble/individual_projections/Papua.elsworthii_EMcaByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.img",
#   wm = "Protea.laurifolia/proj_2070 ensemble/individual_projections/Papua.elsworthii_EMwmeanByROC_mergedData_mergedRun_mergedAlgo_TSSbin.img") )

papuisl_bin_proj_current <- stack(
  c( ca = "Papua.taeniata/proj_current ensemble/individual_projections/Papua.taeniata_EMcaByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.img",
     wm = "Papua.taeniata/proj_current ensemble/individual_projections/Papua.taeniata_EMwmeanByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.img") )

papuisl_bin_proj_2050_BC45 <- stack(
  c( ca = "Papua.taeniata/proj_Papuisl 2050 ensemble/individual_projections/Papua.taeniata_EMcaByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.img",
     wm = "Papua.taeniata/proj_Papuisl 2050 ensemble/individual_projections/Papua.taeniata_EMwmeanByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.img") )

#papupen_bin_proj_2070_BC45 <- stack(
# c( ca = "Protea.laurifolia/proj_2070 ensemble/individual_projections/Papua.elsworthii_EMcaByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.img",
#   wm = "Protea.laurifolia/proj_2070 ensemble/individual_projections/Papua.elsworthii_EMwmeanByROC_mergedData_mergedRun_mergedAlgo_TSSbin.img") )

papuame_bin_proj_current <- stack(
  c( ca = "Papua.papua/proj_current ensemble/individual_projections/Papua.papua_EMcaByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.img",
     wm = "Papua.papua/proj_current ensemble/individual_projections/Papua.papua_EMwmeanByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.img") )

papuame_bin_proj_2050_BC45 <- stack(
  c( ca = "Papua.papua/proj_Papuame 2050 ensemble/individual_projections/Papua.papua_EMcaByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.img",
     wm = "Papua.papua/proj_Papuame 2050 ensemble/individual_projections/Papua.papua_EMwmeanByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.img") )

#papupen_bin_proj_2070_BC45 <- stack(
# c( ca = "Protea.laurifolia/proj_2070 ensemble/individual_projections/Papua.elsworthii_EMcaByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.img",
#   wm = "Protea.laurifolia/proj_2070 ensemble/individual_projections/Papua.elsworthii_EMwmeanByROC_mergedData_mergedRun_mergedAlgo_TSSbin.img") )

# analisis de dinamicas de rango 2050

papuker_SRC_current_2050_BC45 <- BIOMOD_RangeSize( papuker_bin_proj_current,
                                                   papuker_bin_proj_2050_BC45 )
papuker_SRC_current_2050_BC45$Compt.By.Models

#plot(papuker_bin_proj_current)
#plot(papuker_bin_proj_2050_BC45)

papupen_SRC_current_2050_BC45 <- BIOMOD_RangeSize( papupen_bin_proj_current,
                                                   papupen_bin_proj_2050_BC45 )
papupen_SRC_current_2050_BC45$Compt.By.Models

#plot(papupen_bin_proj_current)
#plot(papupen_bin_proj_2050_BC45)

papuisl_SRC_current_2050_BC45 <- BIOMOD_RangeSize( papuisl_bin_proj_current,
                                                   papuisl_bin_proj_2050_BC45 )
papuisl_SRC_current_2050_BC45$Compt.By.Models

#plot(papuisl_bin_proj_current)
#plot(papuisl_bin_proj_2050_BC45)

papuame_SRC_current_2050_BC45 <- BIOMOD_RangeSize( papuame_bin_proj_current,
                                                   papuame_bin_proj_2050_BC45 )
papuame_SRC_current_2050_BC45$Compt.By.Models

# cargar proyecciones binarias 2100

papuker_bin_proj_current <- stack(
  c( ca = "Papua.kerguelensis/proj_current ensemble/individual_projections/Papua.kerguelensis_EMcaByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.img",
     wm = "Papua.kerguelensis/proj_current ensemble/individual_projections/Papua.kerguelensis_EMwmeanByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.img") )

papuker_bin_proj_2100_BC45 <- stack(
  c( ca = "Papua.kerguelensis/proj_Papuker 2100 ensemble/individual_projections/Papua.kerguelensis_EMcaByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.img",
     wm = "Papua.kerguelensis/proj_Papuker 2100 ensemble/individual_projections/Papua.kerguelensis_EMwmeanByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.img") )

#papuker_bin_proj_2070_BC45 <- stack(
 # c( ca = "Protea.laurifolia/proj_2070 ensemble/individual_projections/Papua.kerguelensis_EMcaByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.img",
  #   wm = "Protea.laurifolia/proj_2070 ensemble/individual_projections/Papua.kerguelensis_EMwmeanByROC_mergedData_mergedRun_mergedAlgo_TSSbin.img") )

papupen_bin_proj_current <- stack(
  c( ca = "Papua.elsworthii/proj_current ensemble/individual_projections/Papua.elsworthii_EMcaByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.img",
     wm = "Papua.elsworthii/proj_current ensemble/individual_projections/Papua.elsworthii_EMwmeanByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.img") )

papupen_bin_proj_2100_BC45 <- stack(
  c( ca = "Papua.elsworthii/proj_Papupen 2100 ensemble/individual_projections/Papua.elsworthii_EMcaByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.img",
     wm = "Papua.elsworthii/proj_Papupen 2100 ensemble/individual_projections/Papua.elsworthii_EMwmeanByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.img") )

#papupen_bin_proj_2070_BC45 <- stack(
 # c( ca = "Protea.laurifolia/proj_2070 ensemble/individual_projections/Papua.elsworthii_EMcaByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.img",
  #   wm = "Protea.laurifolia/proj_2070 ensemble/individual_projections/Papua.elsworthii_EMwmeanByROC_mergedData_mergedRun_mergedAlgo_TSSbin.img") )

papuisl_bin_proj_current <- stack(
  c( ca = "Papua.taeniata/proj_current ensemble/individual_projections/Papua.taeniata_EMcaByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.img",
     wm = "Papua.taeniata/proj_current ensemble/individual_projections/Papua.taeniata_EMwmeanByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.img") )

papuisl_bin_proj_2100_BC45 <- stack(
  c( ca = "Papua.taeniata/proj_Papuisl 2100 ensemble/individual_projections/Papua.taeniata_EMcaByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.img",
     wm = "Papua.taeniata/proj_Papuisl 2100 ensemble/individual_projections/Papua.taeniata_EMwmeanByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.img") )

#papupen_bin_proj_2070_BC45 <- stack(
# c( ca = "Protea.laurifolia/proj_2070 ensemble/individual_projections/Papua.elsworthii_EMcaByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.img",
#   wm = "Protea.laurifolia/proj_2070 ensemble/individual_projections/Papua.elsworthii_EMwmeanByROC_mergedData_mergedRun_mergedAlgo_TSSbin.img") )

papuame_bin_proj_current <- stack(
  c( ca = "Papua.papua/proj_current ensemble/individual_projections/Papua.papua_EMcaByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.img",
     wm = "Papua.papua/proj_current ensemble/individual_projections/Papua.papua_EMwmeanByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.img") )

papuame_bin_proj_2100_BC45 <- stack(
  c( ca = "Papua.papua/proj_Papuame 2100 ensemble/individual_projections/Papua.papua_EMcaByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.img",
     wm = "Papua.papua/proj_Papuame 2100 ensemble/individual_projections/Papua.papua_EMwmeanByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.img") )

#papupen_bin_proj_2070_BC45 <- stack(
# c( ca = "Protea.laurifolia/proj_2070 ensemble/individual_projections/Papua.elsworthii_EMcaByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.img",
#   wm = "Protea.laurifolia/proj_2070 ensemble/individual_projections/Papua.elsworthii_EMwmeanByROC_mergedData_mergedRun_mergedAlgo_TSSbin.img") )

# analisis de dinamicas de rango 2100

papuker_SRC_current_2100_BC45 <- BIOMOD_RangeSize( papuker_bin_proj_current,
                                           papuker_bin_proj_2100_BC45 )
papuker_SRC_current_2100_BC45$Compt.By.Models

#plot(papuker_bin_proj_current)
#plot(papuker_bin_proj_2100_BC45)

papupen_SRC_current_2100_BC45 <- BIOMOD_RangeSize( papupen_bin_proj_current,
                                                   papupen_bin_proj_2100_BC45 )
papupen_SRC_current_2100_BC45$Compt.By.Models

#plot(papupen_bin_proj_current)
#plot(papupen_bin_proj_2100_BC45)

papuisl_SRC_current_2100_BC45 <- BIOMOD_RangeSize( papuisl_bin_proj_current,
                                                   papuisl_bin_proj_2100_BC45 )
papuisl_SRC_current_2100_BC45$Compt.By.Models

#plot(papuisl_bin_proj_current)
#plot(papuisl_bin_proj_2100_BC45)

papuame_SRC_current_2100_BC45 <- BIOMOD_RangeSize( papuame_bin_proj_current,
                                                   papuame_bin_proj_2100_BC45 )
papuame_SRC_current_2100_BC45$Compt.By.Models

paleta <- c("red", "orange", "grey", "blue")
#breakpoints <- c(-2,-1,0,1,2)

pdf("Papua_American_RangeShift.pdf")
par(mfrow = c(3, 2))
plot(papuame_bin_proj_current$ca, main="American clade, present, CA")
plot(papuame_bin_proj_current$wm, main="American clade, present, WM")
plot(papuame_bin_proj_2050_BC45$ca, main="American clade, future, CA")
plot(papuame_bin_proj_2050_BC45$wm, main="American clade, future, WM")
plot(papuame_SRC_current_2050_BC45$Diff.By.Pixel$ca, col=paleta, main="American clade, shift, CA") 
#+legend(legend=c("Unsuitable","Destruction","Stability", "Expansión"))
plot(papuame_SRC_current_2050_BC45$Diff.By.Pixel$wm, col=paleta, main="American clade, shift, WM")
dev.off()

pdf("Papua_eastern_RangeShift.pdf")
par(mfrow = c(3, 2))
plot(papuisl_bin_proj_current$ca, main="eastern clade, present, CA")
plot(papuisl_bin_proj_current$wm, main="eastern clade, present, WM")
plot(papuisl_bin_proj_2050_BC45$ca, main="eastern clade, future, CA")
plot(papuisl_bin_proj_2050_BC45$wm, main="eastern clade, future, WM")
plot(papuisl_SRC_current_2050_BC45$Diff.By.Pixel$ca, col=paleta, main="eastern clade, shift, CA") 
#+legend(legend=c("Unsuitable","Destruction","Stability", "Expansión"))
plot(papuisl_SRC_current_2050_BC45$Diff.By.Pixel$wm, col=paleta, main="eastern clade, shift, WM")
dev.off()

pdf("Papua_Antarctic_RangeShift.pdf")
par(mfrow = c(3, 2))
plot(papupen_bin_proj_current$ca, main="Antarctic clade, present, CA")
plot(papupen_bin_proj_current$wm, main="Antarctic clade, present, WM")
plot(papupen_bin_proj_2050_BC45$ca, main="Antarctic clade, 2050, CA")
plot(papupen_bin_proj_2050_BC45$wm, main="Antarctic clade, 2050, WM")
plot(papupen_SRC_current_2050_BC45$Diff.By.Pixel$ca, col=paleta, main="Antarctic clade, shift, CA") 
#+legend(legend=c("Unsuitable","Destruction","Stability", "Expansión"))
plot(papupen_SRC_current_2050_BC45$Diff.By.Pixel$wm, col=paleta, main="Antarctic clade, shift, WM")
dev.off()

pdf("Papua_Kerguelen_RangeShift.pdf")
par(mfrow = c(2, 3))
plot(papuker_bin_proj_current$ca, main="Kerguelen clade, present, CA")
plot(papuker_bin_proj_current$wm, main="Kerguelen clade, present, WM")
plot(papuker_bin_proj_2050_BC45$ca, main="Kerguelen clade, future, CA")
plot(papuker_bin_proj_2050_BC45$wm, main="Kerguelen clade, future, WM")
plot(papuker_SRC_current_2050_BC45$Diff.By.Pixel$ca, col=paleta, main="Kerguelen clade, shift, CA") 
#+legend(legend=c("Unsuitable","Destruction","Stability", "Expansión"))
plot(papuker_SRC_current_2050_BC45$Diff.By.Pixel$wm, col=paleta, main="Kerguelen clade, shift, WM")
dev.off()


##disposicion alterna

pdf("Papua_Antarctic_RangeShift.pdf")
par(mfrow = c(3, 2))
plot(papupen_bin_proj_current$ca, main="Antarctic clade, present, CA")
plot(papupen_bin_proj_current$ca, main="Antarctic clade, present, CA")
plot(papupen_bin_proj_2050_BC45$ca, main="Antarctic clade, 2050, CA")
plot(papupen_SRC_current_2050_BC45$Diff.By.Pixel$ca, col=paleta, main="Antarctic, shift 2050, CA") 
plot(papupen_bin_proj_2100_BC45$ca, main="Antarctic clade, 2100, CA")
plot(papupen_SRC_current_2100_BC45$Diff.By.Pixel$ca, col=paleta, main="Antarctic, shift 2100, CA") 
dev.off()


# analisis de dinamicas de rango 2070


## plot de dinamica de rangos, revisar nueva funci?n bm_range
## https://biomodhub.github.io/biomod2/reference/bm_PlotRangeSize.html

#plot(papuker_SRC_current_2050_BC45$Diff.By.Pixel)
#plot(papupen_SRC_current_2050_BC45$Diff.By.Pixel)
#plot(papuisl_SRC_current_2050_BC45$Diff.By.Pixel)
#plot(papuame_SRC_current_2050_BC45$Diff.By.Pixel)

#str(papuker_SRC_current_2050_BC45)

#ProLau_src_map <- stack(SRC_current_2050_BC45$Diff.By.Pixel$ca, SRC_current_2070_BC45$Diff.By.Pixel$ca)

writeRaster(papupen_SRC_current_2050_BC45$Diff.By.Pixel$ca, "papupen_my_file50_morphed_ca.tif", overwrite = T)
writeRaster(papuker_SRC_current_2050_BC45$Diff.By.Pixel$ca, "papuker_my_file50_morphed_ca.tif", overwrite = T)
writeRaster(papuisl_SRC_current_2050_BC45$Diff.By.Pixel$ca, "papuisl_my_file50_morphed_ca.tif", overwrite = T)
writeRaster(papuame_SRC_current_2050_BC45$Diff.By.Pixel$ca, "papuame_my_file50_morphed_ca.tif", overwrite = T)

writeRaster(papupen_SRC_current_2050_BC45$Diff.By.Pixel$ca, "papupen_my_file100_morphed_ca.tif", overwrite = T)
writeRaster(papuker_SRC_current_2050_BC45$Diff.By.Pixel$ca, "papuker_my_file100_morphed_ca.tif", overwrite = T)
writeRaster(papuisl_SRC_current_2050_BC45$Diff.By.Pixel$ca, "papuisl_my_file100_morphed_ca.tif", overwrite = T)
writeRaster(papuame_SRC_current_2050_BC45$Diff.By.Pixel$ca, "papuame_my_file100_morphed_ca.tif", overwrite = T)

#plot(papuker_SRC_current_2050_BC45$Diff.By.Pixel$wm)
#points(points_papuker[,c('x', 'y')])


#papupen_SRC_current_2050_BC45_ca <- raster("my_file50_morphed.tif")
#papupen_SRC_current_2070_BC45_ca <- raster("my_file70_morphed.tif")

#papupen_src_map <- stack(papupen_SRC_current_2050_BC45, papupen_SRC_current_2050_BC45_ca)

#names(papupen_src_map) <- c("ca cur", "wm cur", "ca 2050", "wm 2050")

#my.at <- seq(-2.5,1.5,1)

#myColorkey <- list(at=my.at, ## where the colors change
 #                  labels=list(
  #                   labels=c("lost", "pres", "abs","gain"), ## labels
   #                  at=my.at[-1]-0.5 ## where to print labels
     #              ))
#rasterVis::levelplot( papupen_src_map,
 #                     main = "Protea laurifolia range change",
  #                    colorkey = myColorkey,
   #                   layout = c(2,2) )

#Load .tif files

###ANALISIS DE CAMBIO

TSS_cutoff_path <- list.files(path = "C:/Users/Usuario/Desktop/Gentoo Revision/Evaluations",
                              pattern = "*.txt",full.names = TRUE)
TSS_cutoff_files <- list.files(path = "C:/Users/Usuario/Desktop/Gentoo Revision/Evaluations",
                               pattern = "*.txt")

library(readr)

list_TSS_cutoff <- list()
for (i in 1:length(TSS_cutoff_files)) {
  list_TSS_cutoff[[i]] <- read_table(TSS_cutoff_path[[i]],
                                     skip = 1,
                                     col_names = c("number",
                                                   "full.name",
                                                   "merged.by.PA",
                                                   "merged.by.run",
                                                   "merged.by.algo",
                                                   "filtered.by",
                                                   "algo",
                                                   "metric.eval",
                                                   "cutoff",
                                                   "sensitivity",
                                                   "specificity",
                                                   "calibration",
                                                   "validation",
                                                   "evaluation")) %>%
    select(full.name, cutoff)
  names(list_TSS_cutoff)[[i]] <- list_TSS_cutoff[[i]]$full.name[1]
}

cutofflist<-matrix(,nrow=4,ncol=1)

colnames(cutofflist) <- "cutoff"
rownames(cutofflist) <- c("Papupen", "Papuker", "Papuame", "Papuisl")

cutofflist[1,1]<-473
cutofflist[2,1]<-312
cutofflist[3,1]<-192
cutofflist[4,1]<-555

current_raster_path <- list.files(path = "C:/Users/Usuario/Desktop/Gentoo Revision/Ensemble Models/Present - ca",
                                  pattern = "*.img",full.names = TRUE) 
current_raster_files <- list.files(path = "C:/Users/Usuario/Desktop/Gentoo Revision/Ensemble Models/Present - ca",
                                   pattern = "*.img") 

future_raster_path <- list.files(path = "C:/Users/Usuario/Desktop/Gentoo Revision/Ensemble Models/2100 - ca",
                                 pattern = "*.img", full.names = TRUE)
future_raster_files <- list.files(path = "C:/Users/Usuario/Desktop/Gentoo Revision/Ensemble Models/2100 - ca",
                                  pattern = "*.img")

splist <- c("Papupen", "Papuker", "Papuame","Papuisl")

list_current_files <- list()
for (i in 1:length(current_raster_files)){
  list_current_files[[i]] <- raster(current_raster_path[[i]]) %>%
    rast(list_current_files[[i]])
  names(list_current_files)[[i]] <- current_raster_files[[i]]
}

list_future_files <- list()
for (i in 1:length(future_raster_files)) {
  list_future_files[[i]] <- raster(future_raster_path[[i]]) %>%
    rast(list_future_files[[i]])
  names(list_future_files)[[i]] <- future_raster_files[[i]]
}

#Converting rasters into dataframes
columns_c <- c("x","y","sdm_c")
columns_f <- c("x","y","sdm_f")

list_current_df <- list()
for (i in 1:length(list_current_files)) {
  list_current_df[[i]] <- as.data.frame(list_current_files[[i]], xy = TRUE)
  names(list_current_df)[[i]] <- splist[[i]]
  colnames(list_current_df[[i]]) <- columns_c
}

list_future_df <- list()
for (i in 1:length(list_future_files)) {
  list_future_df[[i]] <- as.data.frame(list_future_files[[i]], xy = TRUE)
  names(list_future_df)[[i]] <- splist[[i]]
  colnames(list_future_df[[i]]) <- columns_f
}

#QUARTILES AND DECILS####

list_current_back_TSS <- list()
list_future_back_TSS <- list()
list_current_back_decil <- list()
list_future_back_decil <- list()
list_current_quart <- list()
list_future_quart <- list()
list_current_decil <- list()
list_future_decil <- list()
list_current_byTSS <- list()
list_future_byTSS <- list()
list_quart_change <- list()
list_quart_change_raster <- list()

i=4

for (i in 1:length(splist)) {
  #Background subsets by TSS cutoff
  list_current_back_TSS[[i]] <- data.frame(list_current_df[[i]]) %>%
    filter(list_current_df[[i]]$sdm_c < cutofflist[i])
  
  names(list_current_back_TSS)[[i]] <- splist[[i]]
  
  list_future_back_TSS[[i]] <- data.frame(list_future_df[[i]]) %>%
    filter(list_future_df[[i]]$sdm_f < cutofflist[i])
  
  names(list_future_back_TSS)[[i]] <- splist[[i]]
  
  #Quartile subsets
  list_current_quart[[i]] <- data.frame(list_current_df[[i]]) %>%
    filter(list_current_df[[i]]$sdm_c >= cutofflist[i])
  
  names(list_current_quart)[[i]] <- splist[[i]]
  
  list_future_quart[[i]] <- data.frame(list_future_df[[i]]) %>%
    filter(list_future_df[[i]]$sdm_f >= cutofflist[i])
  
  names(list_future_quart)[[i]] <- splist[[i]]
  
  #Decil subsets
  
  #Quartil values for current distribution
  q1_c <- quantile(list_current_quart[[i]]$sdm_c, probs = 0.25)
  q2_c <- quantile(list_current_quart[[i]]$sdm_c, probs = 0.5)
  q3_c <- quantile(list_current_quart[[i]]$sdm_c, probs = 0.75)
  
  #Quartil values for future distribution
  q1_f <- quantile(list_future_quart[[i]]$sdm_f, probs = 0.25)
  q2_f <- quantile(list_future_quart[[i]]$sdm_f, probs = 0.5)
  q3_f <- quantile(list_future_quart[[i]]$sdm_f, probs = 0.75)
  
  #Creating a new quartil column and assigning values
  list_current_quart[[i]]$quart_c <- list_current_quart[[i]]$sdm_c
  
  list_current_quart[[i]]$quart_c[list_current_quart[[i]]$quart_c <= q1_c] <- 1
  list_current_quart[[i]]$quart_c[list_current_quart[[i]]$quart_c > q1_c &
                                    list_current_quart[[i]]$quart_c <= q2_c] <- 2
  list_current_quart[[i]]$quart_c[list_current_quart[[i]]$quart_c > q2_c &
                                    list_current_quart[[i]]$quart_c <= q3_c] <- 3
  list_current_quart[[i]]$quart_c[list_current_quart[[i]]$quart_c > q3_c] <- 4
  
  list_future_quart[[i]]$quart_f <- list_future_quart[[i]]$sdm_f
  
  list_future_quart[[i]]$quart_f[list_future_quart[[i]]$quart_f <= q1_f] <- 1
  list_future_quart[[i]]$quart_f[list_future_quart[[i]]$quart_f > q1_f &
                                   list_future_quart[[i]]$quart_f <= q2_f] <- 2
  list_future_quart[[i]]$quart_f[list_future_quart[[i]]$quart_f > q2_f &
                                   list_future_quart[[i]]$quart_f <= q3_f] <- 3
  list_future_quart[[i]]$quart_f[list_future_quart[[i]]$quart_f > q3_f] <- 4
  
  #Add back the background values
  list_current_back_TSS[[i]]$quart_c <- rep(0, length(list_current_back_TSS[[i]]$sdm_c))
  list_future_back_TSS[[i]]$quart_f <- rep(0, length(list_future_back_TSS[[i]]$sdm_f))
  
  #Merge back zeros
  list_current_byTSS[[i]] <- rbind(list_current_back_TSS[[i]],list_current_quart[[i]])
  list_future_byTSS[[i]] <- rbind(list_future_back_TSS[[i]],list_future_quart[[i]])
  
    #max(list_quart_change[[1]][["quart_f"]])
    #max(list_future_byTSS[[1]][["quart_f"]])

  #Merge everything
  list_quart_change[[i]] <- merge(list_current_byTSS[[i]], list_future_byTSS[[i]], by = c("x","y"))

  #max(list_quart_change[[1]][["quart_c"]])
  
  #Classify change status
  list_quart_change[[i]]$category <- rep(("NA"), length(list_quart_change[[i]]$quart_c))
  
  list_quart_change[[i]]$category[list_quart_change[[i]]$quart_c != 0 &
                                    list_quart_change[[i]]$quart_f == 0] <- "1_Destruction"
  list_quart_change[[i]]$category[list_quart_change[[i]]$quart_c == 0 &
                                    list_quart_change[[i]]$quart_f != 0] <- "5_Expansion"
  list_quart_change[[i]]$category[list_quart_change[[i]]$quart_c != 0 &
                                    list_quart_change[[i]]$quart_f != 0 &
                                    list_quart_change[[i]]$quart_c > list_quart_change[[i]]$quart_f] <- "2_Marginalization"
  list_quart_change[[i]]$category[list_quart_change[[i]]$quart_c != 0 & list_quart_change[[i]]$quart_f != 0 &
                                    list_quart_change[[i]]$quart_c < list_quart_change[[i]]$quart_f] <- "4_Centralization"
  list_quart_change[[i]]$category[list_quart_change[[i]]$quart_c !=0 & list_quart_change[[i]]$quart_f != 0 &
                                    list_quart_change[[i]]$quart_c == list_quart_change[[i]]$quart_f] <- "3_Stability"
  
  list_quart_change[[i]]$category <- as.factor(list_quart_change[[i]]$category)
  list_quart_change[[i]]$numeric_category <- as.numeric(list_quart_change[[i]]$category)
  
  #Create raster by quartiles
  list_quart_change_raster[[i]] = rasterFromXYZ(list_quart_change[[i]][,c("x","y","numeric_category")])
  list_quart_change_raster[[i]][] = factor(levels(list_quart_change[[i]]$category)
                                           [list_quart_change_raster[[i]][]])
  crs(list_quart_change_raster[[i]]) <- "+proj=longlat +datum=WGS84 +no_defs"
  
  #Optional: plot raster directly in R
  #plot(list_quart_change_raster[[i]])

  
  #Save raster files
  sppname <- splist[[i]]
  raster_path <- file.path("C:/Users/Usuario/Desktop/Gentoo Revision/Ensemble Models",
                           paste0(sppname,"_quart_change.tif"))
  
  writeRaster(list_quart_change_raster[[i]],
              file = raster_path,
              overwrite = TRUE)
}

plot(list_quart_change_raster[[1]])


i=2

papuker


#####DISTANCES

#xy <- papuker[,c(5,6)]

#spdf <- SpatialPointsDataFrame(coords = xy, data = papuker,
                               proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

ker = st_as_sf(papuker, coords = c("lon", "lat"), crs = st_crs(4326), agr = "constant")
pen = st_as_sf(papupen, coords = c("lon", "lat"), crs = st_crs(4326), agr = "constant")
isl = st_as_sf(papuisl, coords = c("lon", "lat"), crs = st_crs(4326), agr = "constant")
ame = st_as_sf(papuame, coords = c("lon", "lat"), crs = st_crs(4326), agr = "constant")

install.packages("remotes")
remotes::install_github("Pakillo/rSDM")

kerpoint <- points2nearestcell(
  locs = ker,
  ras = near_raster,
  layer = 1,
  move = TRUE,
  distance = NULL,
  table = TRUE,
  map = "none"
)

islpoint <- points2nearestcell(
  locs = isl,
  ras = near_raster,
  layer = 1,
  move = TRUE,
  distance = NULL,
  table = TRUE,
  map = "none"
)

penpoint <- points2nearestcell(
  locs = pen,
  ras = near_raster,
  layer = 1,
  move = TRUE,
  distance = NULL,
  table = TRUE,
  map = "none"
)

amepoint <- points2nearestcell(
  locs = ame,
  ras = near_raster,
  layer = 1,
  move = TRUE,
  distance = NULL,
  table = TRUE,
  map = "none"
)

near_raster <- as(buffer, "SpatRaster") 

library(rSDM)

islpt = sf_to_df(islpoint, fill=TRUE, unlist=NULL)
penpt = sf_to_df(penpoint, fill=TRUE, unlist=NULL)
amept = sf_to_df(amepoint, fill=TRUE, unlist=NULL)
kerpt = sf_to_df(kerpoint, fill=TRUE, unlist=NULL)

write.csv(islpt, "easterndata.csv")
write.csv(kerpt, "Kerguelendata.csv")
write.csv(amept, "Americandata.csv")
write.csv(penpt, "Peninsuladata.csv")
