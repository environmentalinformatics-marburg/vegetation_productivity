### environmental stuff --------------------------------------------------------

## required packages
library(Rsenal)
library(RColorBrewer)
library(grid)

## required functions
source("R/visDEM.R")


### import and overlay data ----------------------------------------------------

## ndvi
fls_ndvi <- list.files("out/ndvi", full.names = TRUE)[3]
rst_ndvi <- raster(fls_ndvi)
rst_ndvi <- resample(rst_ndvi, rst_crd)

## crd
rst_crd <- raster("out/crd/CRD_0314.tif")
rst_crd <- resample(rst_crd, rst_ndvi)

## reject all crd with non-conclusive long-term trends
rst_crd <- overlay(rst_crd, rst_ndvi, fun = function(x, y) {
  x[is.na(y[])] <- NA
  return(x)
}, filename = "out/crd/CRD_0314_clean.tif", format = "GTiff", overwrite = TRUE)


### create crd figure ----------------------------------------------------------

## colors
cols <- colorRampPalette(brewer.pal(6, "BrBG"))

## create plot
p_crd <- spplot(rst_crd, scales = list(draw = TRUE, cex = .6), 
                xlab = list("x", cex = .8), ylab = list("y", cex = .8), 
                col.regions = cols(25), at = seq(-.5, .5, .25), 
                colorkey = list(labels = list(at = seq(-.5, .5, .25)), 
                                height = .5, width = .75), 
                maxpixels = ncell(rst_crd))


### create background contour lines --------------------------------------------

rst_dem <- raster("data/dem/DEM_ARC1960_30m_Hemp.tif")
rst_dem <- aggregate(rst_dem, fact = 10)
p_dem <- visDEM(rst_dem, labcex = .6, cex = 1.6, col = "black")


### combine and save final figure ----------------------------------------------

## add contour lines to plot
p_crd <- p_crd + 
  latticeExtra::as.layer(p_dem)

## save plot
png("vis/figure_crd.png", width = 14, height = 12.5, units = "cm", res = 500)
plot.new()

# main figure
vp0 <- viewport(x = 0, y = .125, width = 1, height = .9, 
                just = c("left", "bottom"), name = "figure.vp")
pushViewport(vp0)
print(p_crd, newpage = FALSE)

# add figure caption
vp1 <- viewport(x = 0, y = -.025, height = 0.1, width = 1,
                just = c("left", "bottom"),
                name = "vp_caption")
pushViewport(vp1)

grid.text(label = expression(bold("Figure 2.") ~  italic("z") * "-score normalized cumulative rain-use efficiency differences (2003-2014)."),
          x = 0.1, just = c("left", "centre"), gp = gpar(cex = .7), hjust = 0)

dev.off()
