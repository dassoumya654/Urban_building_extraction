library(raster)
library(rgdal)
library(randomForest)
library(kernlab)
library(e1071)
library(RStoolbox)

aoi <- brick(file.choose())
pca <- rasterPCA(aoi,nlayers(aoi))
pc1 <- pca$map$PC1
ndvi <- raster(file.choose())


names(aoi) <- paste0("B",c(1:3))
c <- stack(aoi$B1,aoi$B2,aoi$B3,pc1)

REF <- crs(aoi)


train.build <- readOGR(file.choose())
train.nonbuild <- readOGR(file.choose())

crs(train.nonbuild) <- REF
crs(train.build) <- REF

df_nonbuild <- extract(c,train.nonbuild,df = TRUE)
df_build <- extract(c,train.build,df = TRUE)

df_build$pres <- 1
df_nonbuild$pres <- 0

training.df <- rbind(df_build,df_nonbuild)

############SVM##########################
library(e1071)
fit.svm <- svm(pres ~ B1+B2+B3+PC1, data = training.df,kernel='radial')

tuning <- tune.svm(svm,pres ~ B1+B2+B3+PC1, data = training.df)

jai_pred <- predict(c,fit.svm,na.rm=T)

plot(fit.svm,data=training.df,B2~B3)
plot(jai_pred)
writeRaster(jai_pred,"svm.tiff",format="GTiff",overwrite=TRUE)
i <- svm_pred

rf <- randomForest(as.factor(pres) ~ B1+B2+B3+PC1, data=training.df, importance=TRUE,ntree=400,na.action=na.exclude)
plot(rf)
build_pred <- predict(c,model=rf,na.rm=TRUE)
plot(build_pred)
writeRaster(build_pred,"rf.tiff",format="GTiff",overwrite=TRUE)
j <- build_pred


aggregate(getValues(area(build_pred, weights=FALSE)), by=list(getValues(build_pred)), sum)
########cleaning###################
formask <- setValues(raster(build_predict), NA)
formask[build_predict> 0.4] <- 1
plot(formask, col="dark green", legend = TRUE)

forestclumps <- clump(formask, directions=8)

clumpFreq <- freq(forestclumps)
head(clumpFreq)
tail(clumpFreq)
clumpFreq <- as.data.frame(clumpFreq)

str(which(clumpFreq$count<800))
str(clumpFreq$value[which(clumpFreq$count <800)])
excludeID <- clumpFreq$value[which(clumpFreq$count<800)]


formaskSieve <- formask
formaskSieve[forestclumps %in% excludeID] <- 0
plot(formaskSieve, ext=e, col="dark green", legend=TRUE)


z <- formaskSieve
z[is.na(z)] <- 0
z <- as.data.frame(z)

v <- rasterToPolygons(l,fun = function(x) {x == 1 })
writeOGR(v, "buildings.shp", driver="ESRI Shapefile")

########################## NN ####################################
library(nnet)
library(neuralnet)
library(raster)
library(rgdal)
library(RSNNS)
library(RStoolbox)

aoi <- brick(file.choose())
names(aoi) <- paste0("B",c(1:16))

pca <- rasterPCA(aoi,nComp =nlayers(aoi))
pc1 <- pca$map$PC1
ndvi <- raster(file.choose())

resample(aoi,ndvi, resample='bilinear')

c <- stack(aoi$B1,aoi$B2,aoi$B3,pc1,ndvi)

REF <- crs(aoi)

train.nonbuild <- readOGR(file.choose())
train.build <-readOGR(file.choose())

crs(train.nonbuild) <- REF
crs(train.build) <- REF

df_nonbuild <- extract(c,train.nonbuild,df = TRUE)
df_build <- extract(c,train.build,df = TRUE)

df_build$pres <- 1
df_nonbuild$pres <- 0

training.df <- rbind(df_build,df_nonbuild)
training.df

values <- training.df[,2:5]
targets <- training.df[,6]
targets

values <- as.matrix(values)
targets <- as.matrix(targets)

model <- nnet(as.factor(pres)~B1+B2+B3+PC1,data = training.df, size=70,rang = 0.1,decay = 5e-4,linout=FALSE,maxit=300)
build_predict <- predict(c,model)
build_predict

writeRaster(build_predict,"nnET.tiff",format="GTiff",overwrite=TRUE)

rc <- function(x){
  ifelse(x>0.48,1,0)
}
value <- calc(build_predict,fun = rc)
plot(value)
build
########cleaning###################
formask <- setValues(raster(build_predict), NA)
formask[build_predict>0.4] <- 1
plot(formask, col="dark green", legend = TRUE)

forestclumps <- clump(formask, directions=8)
clumpFreq <- freq(forestclumps)
head(clumpFreq)
tail(clumpFreq)
clumpFreq <- as.data.frame(clumpFreq)

str(which(clumpFreq$count<10))
str(clumpFreq$value[which(clumpFreq$count <10)])
excludeID <- clumpFreq$value[which(clumpFreq$count<10)]

formaskSieve <- formask
formaskSieve[forestclumps %in% excludeID] <- NA
plot(formaskSieve, ext=e, col="dark green", legend=TRUE)

nn_mod<- formaskSieve
nn_mod[is.na(nn_mod)] <- 0

writeRaster(nn_mod,"nn_cleaned.tiff",format="GTiff",overwrite=TRUE)

ensemble <- cbind(x,y,z)

library(functional)
m <-apply(ensemble, 1, Compose(table,function(u) u==max(u),which,names,function(u) paste0(u, collapse='/')))
ensemble$presence <- m
names(ensemble) <- paste0("c",1:4)
o <-stack(x,y,z)

####################################################### GBM #################################################################

library(gbm)
model.gbm <- gbm(formula = c4 ~ .,n.trees = 1500,shrinkage = 0.01,distribution = "gaussian", interaction.depth = 4 ,data = ensemble,verbose = TRUE)
pred.gbm <- predict(o,model.gbm,n.trees=model.gbm$n.trees)

writeRaster(pred.gbm,"ensemble_all_algo.tiff",format="GTiff",overwrite=TRUE)
pred.gbm
tab_gbm <-table(factor(pred.gbm, levels=min(ensemble$c4):max(ensemble$c4)), factor(ensemble$c4, levels=min(ensemble$c4):max(ensemble$c4)))
confusionMatrix(tab_gbm)


library(gbm)
model.gbm <- gbm(formula = c4 ~ ., distribution = "bernoulli", data = ensemble)
pred.gbm <- predict(o,model.gbm,n.trees=model.gbm$n.trees,type="response")

model <- nnet(c4 ~ .,data = ensemble, size=60,rang = 0.1,decay = 5e-4,linout=FALSE,maxit=200)
l <- predict(e,model)

o <- stack(x,y,value)
names(o) <- paste0("B",1:3)

r <- predict(o, model.gbm, filename="gbm_fianl_output", na.rm=TRUE, overwrite=TRUE, n.trees=model.gbm$n.trees, type="response", progress="window")
