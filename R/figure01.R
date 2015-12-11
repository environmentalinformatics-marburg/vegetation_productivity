### environmental stuff --------------------------------------------------------

## required packages
library(Rsenal)
library(foreach)
library(RColorBrewer)
library(grid)

## required functions
source("R/visDEM.R")


### create single trend figures ------------------------------------------------

## list available trend layers (r, mann-kendall)
fls_trends <- list.files("out/ndvi", full.names = TRUE)[c(1, 3, 2, 4)]

## colors
cols <- colorRampPalette(brewer.pal(11, "BrBG"))

lst_p_trends <- foreach(i = 1:length(fls_trends), 
                        txt = c("a)", "b)", "c)", "d)")) %do% {
  rst <- raster(fls_trends[i])
  
  spplot(rst, scales = list(draw = TRUE, cex = .6), maxpixels = ncell(rst), 
         xlab = list("x", vjust = -.5, cex = .8), ylab = list("y", cex = .8), 
         col.regions = cols(1000), 
         at = if (i %in% c(1, 2)) seq(-1, 1, .01) else seq(-.5, .5, .01), 
         par.settings = list(fontsize = list(text = 15)), 
         colorkey = list(labels = list(cex = .7), space = "top", 
                         height = .5, width = .75), 
         main = list("r", cex = .85, font = 2, hjust = -.2)) + 
    latticeExtra::layer(sp.text(loc = c(280000, 9627500), txt = txt, font = 2, 
                                cex = .7, adj = c(.1, 1), col = "black"), 
                        data = list(txt = txt))
}


### create background contour lines --------------------------------------------

rst_dem <- raster("data/dem/DEM_ARC1960_30m_Hemp.tif")
rst_dem <- aggregate(rst_dem, fact = 10)
p_dem <- visDEM(rst_dem, labcex = .6, cex = 1.6, col = "black")


### combine and save final figure ----------------------------------------------

## combine figures
p_trends <- latticeCombineGrid(lst_p_trends, layout = c(2, 2)) + 
  latticeExtra::as.layer(p_dem)

## save output
png("vis/figure_trends.png", width = 20.5, height = 20, units = "cm", res = 500)
plot.new()

# main figure
vp1 <- viewport(x = 0, y = .125, width = 1, height = .9, 
                just = c("left", "bottom"), name = "figure.vp")
pushViewport(vp1)
print(p_trends, newpage = FALSE)

# additional colorkey for kendall's tau
downViewport(trellis.vpname(name = "figure"))

vp0 <- viewport(x = .5, y = -.225, width = 1, height = .1, 
                just = c("centre", "bottom"), name = "key_tau.vp")
pushViewport(vp0)
draw.colorkey(key = list(labels = list(cex = .7), col = cols(1000), width = .75, 
                         height = .5, at = seq(-.5, .5, .01), x = -.1, 
                         space = "bottom"), draw = TRUE)

grid.text(expression(bold("Kendall's " ~ tau)), x = 0.5, y = -.1, 
          just = c("centre", "top"), gp = gpar(font = 2, cex = .85))

# add figure caption
upViewport(0)
vp2 <- viewport(x = 0, y = -.01, height = 0.1, width = 1,
                just = c("left", "bottom"),
                name = "vp_caption")
pushViewport(vp2)

grid.text(label = expression(bold("Figure 1.") ~ "Long-term NDVI trends (2003-2014) from Terra-MODIS (left) and Aqua-MODIS (right)."),
          x = 0.1, just = c("left", "centre"))
dev.off()
