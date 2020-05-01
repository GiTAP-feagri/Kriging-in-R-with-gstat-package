#This code is made for geostatistical interpolation with ordinary kriging using gstat package

For this we use the RStudio iteration with R version 3.4.4

#There are included the following analysis:
1st - First steps in R, installing and loading libraries, loading directory to source file location
2nd - Variogram analysis - we use Method of Moments (MoM) for experimental variogram with GLS fitting of theorical variogram with spheric, expanential and gaussian models.
3rd - Leave-one-out Cross validation (LOOCV) to choose semivariogram fit.
4th - Interpolation with kriging:  puntual Ordinary kriging (OK)
5th - Exporting interpolated maps

#1. First steps in R: 
##1.1 - We start by cleaning R environment 
rm(list = ls())
gc(reset=T)
graphics.off()

##1.2 - And install required packages
#install.packages("pacmann")
pacmann::p_load(gstat, raster, rstudioapi, sp)

##1.3 - Than we set working directory to source file location (directory turns to be the location where the R script is saved in your computer)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

##1.4 - Loading data: our data is already free of outliers; we strongly recommend data preprocessing prior to interpolation
data = read.csv(file = "../../data/data points/data.csv", header = TRUE, sep = ',')
data <- data[,c(2,3,4)] #selecting important columns (x, y, z)
names(data) <- c("x", "y", "z")
sp::coordinates(data) = ~x+y # transform data to spatial object (x and y must be in UTM)

##1.5 - We separate the primary variable. This will facilitate analysis
solo_atr<- data$z

##1.6 - Data visualization according to the "z" values
sp::bubble(data, "z")

#2. Variogram analysis: using MoM.

##2.1 - We start by visualization of the empirical variogram
g = gstat(id="solo_attr", formula = solo_atr ~ 1, data=data)
print(max(dist(data@coords))/2) # maximum distance between point coordinates (helps to limite variogram cutoff) 
print(min(dist(data@coords))) # minimum distance between point coordinates (helps to define width)

cutoff = max(dist(data@coords))/2
##2.2 - The variogram cloud
var_could = gstat::variogram(g, cloud=T) 
plot(var_could)

##2.3 Variogram map
It is possible to identify trend and anisotropy
var_map = gstat::variogram(g, cutoff=cutoff, width=100, map=T) 
plot(var_map)

##2.4 - Experimental variogram by MoM
var_exp = gstat::variogram(g, cutoff=cutoff, width=100) 
plot(var_exp)

g = gstat(id="solo_attr", formula = solo_atr ~ 1, data=data)
var_exp = gstat::variogram(g, cutoff=cutoff, width=100)
plot(var_exp)

##2.5 - Theorical variogram 
We first use a "eye fit" to choose the variogram parameters - nugget effect, partial sill and range. 
It can be fitted by a variety of models. Here we test spheric, exponential and gaussian models.

###2.5.1 - Spheric
fit.sph = fit.variogram(var_exp, vgm(0.004, "Sph", 450, 0.001)) # vgm(partial sill, model, range, nugget effect)
plot(var_exp, fit.sph)

###2.5.2 - Exponential
fit.exp = fit.variogram(var_exp, vgm(0.004, "Exp", 450, 0.001)) # vgm(partial sill, model, range, nugget effect)
plot(var_exp, fit.exp)

###2.5.3 - Gaussian
fit.gau = fit.variogram(var_exp, vgm(0.004, "Gau", 450, 0.001)) # vgm(partial sill, model, range, nugget effect)
plot(var_exp, fit.gau)

#3. LOOCV: leave-one-out cross validation to choose the best model.

###3.1 - LOOCV of spheric model
xvalid.sph = krige.cv(z ~ 1, locations = data, model = fit.sph) # cross validation function
plot(xvalid.sph$var1.pred ~ solo_atr, cex = 1.2, lwd = 2) #, ylim=c(10,50), xlim=c(10,50)) 
abline(0, 1, col = "lightgrey", lwd = 2)  
lm_sph = lm(xvalid.sph$var1.pred ~ solo_atr) 
abline(lm_sph, col = "red", lwd = 2) 
r2_sph = summary(lm_sph)$r.squared 
rmse_sph = hydroGOF::rmse(xvalid.sph$var1.pred, solo_atr)

###3.2 - LOOCV of exponential model
xvalid.exp = krige.cv(z ~ 1, locations = data, model = fit.exp)
plot(xvalid.exp$var1.pred ~ dados$z, cex = 1.2, lwd = 2) 
abline(0, 1, col = "lightgrey", lwd = 2) 
lm_exp = lm(xvalid.exp$var1.pred ~ solo_atr)
abline(lm_exp,  col = "red", lwd = 2) 
r2_exp = summary(lm_exp)$r.squared 
rmse_exp = hydroGOF::rmse(xvalid.exp$var1.pred, solo_atr)

###3.3 - LOOCV of gaussian model
xvalid.gau = krige.cv(z ~ 1, locations = data, model = fit.gau)
plot(xvalid.gau$var1.pred ~ dados$z, cex = 1.2, lwd = 2) 
abline(0, 1, col = "lightgrey", lwd = 2) 
lm_gau = lm(xvalid.gau$var1.pred ~ solo_atr) 
abline(lm_gau,  col = "red", lwd = 2)
r2_gau = summary(lm_gau)$r.squared # extrai R2
rmse_gau = hydroGOF::rmse(xvalid.gau$var1.pred, solo_atr) 


##3.4 Visualization of LOOCV results
df.r2 = data.frame(r2_exp,r2_gau,r2_sph) #Set with R2       
df.rmse = data.frame(rmse_exp, rmse_gau,rmse_sph) # Set with RMSE
temp = data.frame(cbind(t(df.r2), t(df.rmse))) # Join sets conjuntos 
colnames(temp) = c('R2', 'RMSE')
rnames = gsub('r2_','',rownames(temp)) # Removes R2_ prefix from row names
rownames(temp) = rnames 
print(temp)

#4. Kriging interpolation

##4.1 We create a grid for interpolation

To performe this we open our data boundary/cotorno

contorno <- shapefile("../../data/boundary/cotorno.shp")

And then we create a grid
r = raster::raster(contorno, res = 5) #  "res" sets pixel resolution
rp = raster::rasterize(contorno, r, 0) 
grid = as(rp, "SpatialPixelsDataFrame") 
sp::plot(grid)
sp::proj4string(data) = sp::proj4string(contorno) # Contorno (shape) and data have the same CRS

## 4.2 Ordianry kriging

mapa <- krige(solo_atr ~ 1, data, grid, model = fit.exp)
plot(mapa)


#5 Export map

We first convert maps format to raster and add the maps projection

##5.1 - Convert to raster
mapaRaster = raster(mapa)
proj4string(mapaRaster) = proj4string(contorno) 


##5.1 - Exporting the map

# Salvar a imagem do mapa
writeRaster(mapaRaster, 
            filename = '../../maps/z_interpolated.tif',#here we choose where we want to save
            format = 'GTiff',
            overwrite = T)



