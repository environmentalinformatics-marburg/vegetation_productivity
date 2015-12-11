### environmental stuff --------------------------------------------------------

## required packages
library(raster)


### import monthly npp ---------------------------------------------------------

## import national park boundaries
spy_np_old <- readOGR("data/NP/", 
                      "fdetsch-kilimanjaro-national-park-1420535670531", 
                      p4s = "+init=epsg:4326")
spy_np_old <- spTransform(spy_np_old, CRS = CRS("+init=epsg:21037"))

## list and import available files
fls_npp <- list.files("out/npp", full.names = TRUE)
rst_npp <- stack(fls_npp)

## mask all pixels inside the national park
rst_npp <- mask(rst_npp, spy_np_old, inverse = TRUE)


### process chirps -------------------------------------------------------------

## list available files
fls_chirps <- list.files("data/CHIRPS/scl", full.names = TRUE)

## restrict temporal range (currently 2003-2014)
id_start <- grep("2003.01", basename(fls_chirps))
id_end <- grep("2014.12", basename(fls_chirps))

fls_chirps <- fls_chirps[id_start:id_end]
rst_chirps <- stack(fls_chirps)

## resample chirps
rst_chirps_rsmpl <- projectRaster(rst_chirps, rst_npp)


### calculate rain-use efficiency (RUE; see Prince et al. 1998) ----------------

rst_rue <- overlay(rst_npp, rst_chirps_rsmpl, fun = function(x, y) x / y, 
                   filename = "out/rue/RUE", bylayer = TRUE, 
                   suffix = substr(names(rst_npp), 5, nchar(names(rst_npp))), 
                   format = "GTiff", overwrite = TRUE)


### calculate Z-score normalized cumulative RUE differences --------------------

rst_crd <- calc(rst_rue, fun = function(x) sum(diff(x)), 
                filename = "out/crd/CRD_0314.tif", 
                format = "GTiff", overwrite = TRUE)