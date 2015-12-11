### environmental stuff --------------------------------------------------------

## required packages
library(doParallel)

# install and attach 'MODIS' package
devtools::install_github("MatMatt/MODIS", ref = "develop")
library(MODIS)

MODISoptions(localArcPath = "data/MODIS_ARC", 
             outDirPath = "data/MODIS_ARC/PROCESSED", 
             outProj = "+init=epsg:21037", 
             MODISserverOrder = c("LPDAAC", "LAADS"))

## install and attach 'Rsenal' package
devtools::install_github("environmentalinformatics-marburg/Rsenal")
library(Rsenal)

## parallelization
supcl <- makeCluster(detectCores() - 1)
registerDoParallel(supcl)


### data download --------------------------------------------------------------

for (product in c("MOD17A3H", "MYD17A3H")) {
  runGdal(product, collection = getCollection(product, forceCheck = TRUE), 
          tileH = 21, tileV = 9, SDSstring = "10")
}


### crop images ----------------------------------------------------------------

## product folders
drs <- list.dirs(getOption("MODIS_outDirPath"))

## reference extent
rst_template <- kiliAerial()
rst_template <- trim(rst_template)

## process single products
lst_npp <- foreach(product = c("MOD17A3H", "MYD17A3H"), 
                   .packages = c("raster", "foreach")) %dopar% {

  # list available files  
  path <- drs[grep(product, drs)]
  fls <- list.files(path, pattern = "Npp_500m.tif$", full.names = TRUE)
  
  # restrict temporal range (currently 2003-2014)
  id_start <- grep("2003", basename(fls))
  id_end <- grep("2014", basename(fls))

  fls <- fls[id_start:id_end]
  rst <- raster::stack(fls)
  
  # crop images
  lst_crp <- foreach(i = 1:nlayers(rst)) %do% {
    file_out <- paste0(path, "/crp/CRP_", names(rst[[i]]), ".tif")
    crop(rst[[i]], rst_template, filename = file_out, format = "GTiff", 
         overwrite = TRUE)
  }
  
  rst_crp <- raster::stack(lst_crp)
  
  # apply scale factor
  lst_scl <- foreach(i = 1:nlayers(rst)) %do% {
    file_out <- paste0(path, "/scl/SCL_", names(rst_crp[[i]]), ".tif")
    calc(rst_crp[[i]], fun = function(x) x * 0.0001, 
         filename = file_out, format = "GTiff", overwrite = TRUE)
  }
  
  rst_scl <- raster::stack(lst_scl)
  
  # reject invalid pixels (see https://lpdaac.usgs.gov/node/869#documentation) 
  # for details
  lst_qc <- foreach(i = 1:nlayers(rst)) %do% {
    file_out <- paste0(path, "/qc/QC_", names(rst_scl[[i]]), ".tif")
    
    num_val <- raster::getValues(rst_scl[[i]])
    num_val[which(num_val > 3.276)] <- NA
    
    raster::writeRaster(raster::setValues(rst_scl[[i]], num_val), 
                        filename = file_out, format = "GTiff", overwrite = TRUE)
  }
  
  rst_qc <- raster::stack(lst_qc)
}


### overlay images -------------------------------------------------------------

lst_npp_mvc <- foreach(i = 1:nlayers(lst_npp[[1]])) %do% {
  
  date_out <- sapply(strsplit(names(lst_npp[[1]][[i]]), "\\."), "[[", 2)
  file_out <- paste0("data/MCD17A3H.006/MCD17A3H.", date_out, ".tif")
  
  overlay(lst_npp[[1]][[i]], lst_npp[[2]][[i]], fun = max, na.rm = TRUE, 
          filename = file_out, format = "GTiff", overwrite = TRUE)
}

rst_npp_mvc <- stack(lst_npp_mvc)
mat_npp_mvc <- as.matrix(rst_npp_mvc)


### process ndvi ---------------------------------------------------------------

## list available 16-day ndvi files
fls_ndvi <- list.files("data/MYD13Q1.006/whittaker", pattern = ".tif", 
                       full.names = TRUE)

## restrict temporal range (currently 2003-2014)
id_start <- grep("2003", basename(fls_ndvi))[1]

id_end <- grep("2014", basename(fls_ndvi))
id_end <- id_end[length(id_end)]

fls_ndvi <- fls_ndvi[id_start:id_end]

## create annual sums
lst_fls_ndvi <- split(fls_ndvi, substr(basename(fls_ndvi), 5, 8))

lst_ndvi_sums <- foreach(i = lst_fls_ndvi, .packages = "raster") %dopar% {
  rst <- raster::stack(i)
  raster::calc(rst, sum, na.rm = TRUE) * 0.0001
}

rst_ndvi_sums <- raster::stack(lst_ndvi_sums)

## resample ndvi to npp
rst_ndvi_sums <- resample(rst_ndvi_sums, rst_npp_sums)
mat_ndvi_sums <- as.matrix(rst_ndvi_sums)


### linear regression ----------------------------------------------------------

## apply 'lm'
dat_mod <- foreach(i = 1:nrow(mat_npp_mvc), .combine = "rbind") %do% {
  mod <- remote:::lmC(mat_ndvi_sums[i, ], mat_npp_mvc[i, ])
  data.frame(cell = i, intercept = mod[[2]], slope = mod[[3]], r = mod[[1]])
}

rst_intercept <- rst_npp_mvc[[1]]
rst_intercept[] <- dat_mod[, 2]

rst_slope <- rst_npp_mvc[[1]]
rst_slope[] <- dat_mod[, 3]


### predict monthly npp --------------------------------------------------------

## list monthly ndvi files
fls_ndvi_mvc <- list.files("data/MYD13Q1.006/whittaker_mvc", 
                           pattern = ".tif", full.names = TRUE)

## restrict temporal range (currently 2003-2014)
id_start <- grep("200301", basename(fls_ndvi_mvc))
id_end <- grep("201412", basename(fls_ndvi_mvc))

fls_ndvi_mvc <- fls_ndvi_mvc[id_start:id_end]
rst_ndvi_ref <- raster(fls_ndvi_mvc[1])

## ensure uniform extent (needs to be taken care of in data preparation)
lst_ndvi_mvc <- foreach(i = 1:length(fls_ndvi_mvc), 
                        .packages = "raster") %dopar% {
  if (i == 1) 
    raster::raster(fls_ndvi_mvc[i]) * 0.0001
  else 
    raster::crop(raster::raster(fls_ndvi_mvc[i]), rst_ndvi_ref) * 0.0001
}

rst_ndvi_mvc <- stack(lst_ndvi_mvc)

## resample ndvi to npp
rst_ndvi_mvc <- resample(rst_ndvi_mvc, rst_npp_mvc)

## predict npp
rst_npp_pred <- overlay(rst_slope, rst_ndvi_mvc, rst_intercept, 
                        fun = function(a, x, b) {
                          a * x + b
                        }, filename = "out/npp/NPP", bylayer = TRUE, 
                        suffix = names(rst_ndvi_mvc), 
                        format = "GTiff", overwrite = TRUE)

## deregister parallel backend
stopCluster(supcl)