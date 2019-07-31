## ------------------------------------- #
## Polygon-Raster Overlay
## -------
## PRISM Climate Annual Average Normals 
## Masked for California, USA
## Polygon overlay using a rodent species distribution
## 
## Author: Ian Buller (@idblr)
## Date created: July 27, 2019
##
## Most Recently Modified on: July 31, 2019
## Most Recently Modified by: Ian Buller (@idblr)
##
## Modifications:
# A) 
## ------------------------------------- #

# -------- #
# PACKAGES #
# -------- #
library(prism) # prism data
library(sp) # spatial data manipulation
library(raster) # raster data manipulation
library(RStoolbox) # PCA of rasters
library(maps) # visualize geographical data
library(rgeos) # calculate buffers
library(rgdal) # read shapefiles

### Set seed
set.seed(42) 

# ---- #
# DATA #
# ---- #
# Environmental Data
# PRISM 30-Year Average Normals
# Download files
options(prism.path = "~/prismtmp")
prism::get_prism_normals(type= "tmean"
                         ,resolution="4km"
                         ,annual = TRUE
                         ,keepZip=FALSE)
# Convert to Rasters
tmean <- prism::ls_prism_data(absPath=T)[4,2]
tmean_rast <- raster::raster(tmean)

# Set Projection of PRISM Rasters
crs_us <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" # original CRS
reproj_tmean_rast <- raster::projectRaster(tmean_rast, crs = crs(crs_us))

# Standardize Rasters by Range Transformation
stand_reproj_tmean_rast <- (reproj_tmean_rast-min(na.omit(reproj_tmean_rast@data@values)))/(max(na.omit(reproj_tmean_rast@data@values))-min(na.omit(reproj_tmean_rast@data@values)))

# Species Distribution Spatial Polygon
# Data are not provided on GitHub
# Access and download after registration at: https://www.wildlife.ca.gov/Data/CWHR
shp_raw <- rgdal::readOGR("/Users/IDB/Documents/GitHub/prism_pca/data/m060.shp", verbose = FALSE) # species distribution for Merriam Chipmunk
shp <- sp::spTransform(shp_raw, crs(mask_tmean)) # set coordinate reference system same as PRISM

shp_crop <- mask(mask_tmean, shp) # mask for raster within polygon (for plot version 2 & 3)
shp_crop <- crop(shp_crop,shp) # crop extent to match polygon (for plot version 2 & 3)

shp_point <- rasterToPoints(shp_crop) # points at every raster centroid
grid.pt <- sp::spsample(shp, n = 200, type = "regular") # sample of points at raster centroids

# ------------------ #
# DATA VISUALIZATION #
# ------------------ #
# Mask standardized rasters by study area (window)
# Here, example using California, USA
us <- raster::getData("GADM", country="USA", level=1) # get US polygon data
ca <- us[match(toupper("California"),toupper(us$NAME_1)),] # get California polygon
# Extract outline in order to create buffer to capture all of PRISM
# region union kills the data frame so don't overwrite 'wus'
regs <- rgeos::gUnaryUnion(ca)
# takes way too long to plot without simplifying the polygons
regs <- rgeos::gSimplify(regs, 0.05, topologyPreserve = TRUE)
# Add 0.1 degree buffer to capture all of PRISM
ca_buffer <- rgeos::gBuffer(regs, width=0.1, byid=TRUE) # same projection as crs_us

# Mask for California
mask_tmean <- mask(stand_reproj_tmean_rast, ca_buffer)

# VERSION 1 = Hash pattern
# Can change density and angle for different patterns
grDevices::pdf(file = "figures/ca_merriam_v1.pdf", width = 7, height = 5)
raster::plot(mask_tmean, ext = ca_buffer, box = F,
             main = substitute(paste(italic('Tamias merriami'), " and mean temperature (1981-2011) with hash pattern")),
             xlab = "Longitude", ylab = "Latitude", cex.axis = 0.67, cex.lab = 0.8,
             legend.args=list(text=expression(paste(degree,"Celsius")), side=1, font=2, 
                              line=-12.33, cex=0.8
             ),
             axis.args = list(labels = c(round(mask_tmean@data@min, digits = 2),
                                         round((mask_tmean@data@max-mask_tmean@data@min)/4+mask_tmean@data@min, digits = 2),
                                         round((mask_tmean@data@max-mask_tmean@data@min)/2+mask_tmean@data@min, digits = 2),
                                         round((mask_tmean@data@max-mask_tmean@data@min)/4*3+mask_tmean@data@min, digits = 2),
                                         round(mask_tmean@data@max, digits = 2)
             ),
             at = c(mask_tmean@data@min,
                    round((mask_tmean@data@max-mask_tmean@data@min)/4+mask_tmean@data@min, digits = 2),
                    round((mask_tmean@data@max-mask_tmean@data@min)/2+mask_tmean@data@min, digits = 2),
                    round((mask_tmean@data@max-mask_tmean@data@min)/4*3+mask_tmean@data@min, digits = 2),
                    mask_tmean@data@max
             ),
             cex.axis=0.67
             ),
             legend.width=1, legend.shrink=0.7
)
plot(shp, add = T, lwd = 0.33, density = 10, angle = 45)
legend(x=-109.5, y = 33.5, legend="Species\nDistribution", fill="transparent",
       density=NA, angle = 45, bty="n", border="black", xpd = T,
       cex = 0.67
)
legend(x=-109.5, y = 33.5, legend="Species\nDistribution", fill="black",
       density=25, angle = 45, bty="n", border="black", xpd = T,
       cex = 0.67
)
dev.off()

# VERSION 2 = dots pattern (at all centroids)
grDevices::pdf(file = "figures/ca_merriam_v2.pdf", width = 7, height = 5)
raster::plot(mask_tmean, ext = ca_buffer, box = F,
             main = substitute(paste(italic('Tamias merriami'), " and mean temperature (1981-2011) with all dot pattern")),
             xlab = "Longitude", ylab = "Latitude", cex.axis = 0.67, cex.lab = 0.8,
             legend.args=list(text=expression(paste(degree,"Celsius")), side=1, font=2, 
                              line=-12.33, cex=0.8
             ),
             axis.args = list(labels = c(round(mask_tmean@data@min, digits = 2),
                                         round((mask_tmean@data@max-mask_tmean@data@min)/4+mask_tmean@data@min, digits = 2),
                                         round((mask_tmean@data@max-mask_tmean@data@min)/2+mask_tmean@data@min, digits = 2),
                                         round((mask_tmean@data@max-mask_tmean@data@min)/4*3+mask_tmean@data@min, digits = 2),
                                         round(mask_tmean@data@max, digits = 2)
             ),
             at = c(mask_tmean@data@min,
                    round((mask_tmean@data@max-mask_tmean@data@min)/4+mask_tmean@data@min, digits = 2),
                    round((mask_tmean@data@max-mask_tmean@data@min)/2+mask_tmean@data@min, digits = 2),
                    round((mask_tmean@data@max-mask_tmean@data@min)/4*3+mask_tmean@data@min, digits = 2),
                    mask_tmean@data@max
             ),
             cex.axis=0.67
             ),
             legend.width=1, legend.shrink=0.7
)
points(shp_point, pch = 16, cex = 0.05)
plot(shp, add = T, lwd = 0.33)
legend(x=-109.5, y = 33.5, legend="Species\nDistribution", pt.bg="transparent",
       bty="n", col="black", xpd = T, pch = 22,
       cex = 0.67, pt.cex = 1.1
)
legend(x=-109.5, y = 33.5, legend="Species\nDistribution", col="black", pch = 16,
       bty="n", xpd = T,
       cex = 0.67, pt.cex = 0.2
)
dev.off()

# VERSION 3 = dots pattern (at a sample of centroids)
grDevices::pdf(file = "figures/ca_merriam_v3.pdf", width = 7, height = 5)
raster::plot(mask_tmean, ext = ca_buffer, box = F,
             main = substitute(paste(italic('Tamias merriami'), " and mean temperature (1981-2011) with few dot pattern")),
             xlab = "Longitude", ylab = "Latitude", cex.axis = 0.67, cex.lab = 0.8,
             legend.args=list(text=expression(paste(degree,"Celsius")), side=1, font=2, 
                              line=-12.33, cex=0.8
             ),
             axis.args = list(labels = c(round(mask_tmean@data@min, digits = 2),
                                         round((mask_tmean@data@max-mask_tmean@data@min)/4+mask_tmean@data@min, digits = 2),
                                         round((mask_tmean@data@max-mask_tmean@data@min)/2+mask_tmean@data@min, digits = 2),
                                         round((mask_tmean@data@max-mask_tmean@data@min)/4*3+mask_tmean@data@min, digits = 2),
                                         round(mask_tmean@data@max, digits = 2)
             ),
             at = c(mask_tmean@data@min,
                    round((mask_tmean@data@max-mask_tmean@data@min)/4+mask_tmean@data@min, digits = 2),
                    round((mask_tmean@data@max-mask_tmean@data@min)/2+mask_tmean@data@min, digits = 2),
                    round((mask_tmean@data@max-mask_tmean@data@min)/4*3+mask_tmean@data@min, digits = 2),
                    mask_tmean@data@max
             ),
             cex.axis=0.67
             ),
             legend.width=1, legend.shrink=0.7
             )
points(grid.pt, pch = 16, col = "black", cex = 0.1)
plot(shp, add = T, lwd = 0.33)
legend(x=-109.5, y = 33.5, legend="Species\nDistribution", pt.bg="transparent",
       bty="n", col="black", xpd = T, pch = 22,
       cex = 0.67, pt.cex = 1.1
)
legend(x=-109.5, y = 33.5, legend="Species\nDistribution", col="black", pch = 16,
       bty="n", xpd = T,
       cex = 0.67, pt.cex = 0.2
)
dev.off()

# ----------- #
# End of Code #
# ----------- #
