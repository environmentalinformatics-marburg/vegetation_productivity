### environmental stuff --------------------------------------------------------

## packages
library(raster)
library(rgdal)
library(doParallel)
library(remote)

## parallelization
supcl <- makeCluster(detectCores() - 1)
registerDoParallel(supcl)


### data processing: ndvi ------------------------------------------------------

## import national park boundaries
spy_np <- readOGR("data/NP/", "fdetsch-kilimanjaro-1420532792846", 
                  p4s = "+init=epsg:4326")
spy_np <- spTransform(spy_np, CRS = CRS("+init=epsg:21037"))

## terra and aqua-modis ndvi
lst_trends <- foreach(product = c("MOD13Q1.006", "MYD13Q1.006"), 
                      .packages = c("gimms", "Orcs", "foreach", 
                                    "remote")) %dopar% {
  
  #   # status message
  #   cat("Processing", product, "...\n")
  
  # list available files
  fls <- list.files(paste0("data/", product, "/whittaker_mvc"), 
                    full.names = TRUE, pattern = ".tif$")
  
  # import files
  rst <- raster::stack(fls)

  #   # create annual mean annual value composites (similar to Landmann and 
  #   # Dubovyk 2014)
  #   lst_fls <- split(fls, substr(basename(fls), 5, 8))
  #   
  #   lst_mvc <- foreach(i = lst_fls, .packages = "raster") %dopar% {
  #     rst <- raster::stack(i)
  #     raster::calc(rst, mean, na.rm = TRUE)
  #   }
  #   
  #   rst_mvc <- raster::stack(lst_mvc)
  #   
  #   # remove seasonal signal (deseasoning)
  #   rst_ltm <- calc(rst_mvc, mean, na.rm = TRUE)
  #   rst_dsn <- rst_mvc - rst_ltm
  
  # remove seasonal signal
  rst_dsn <- remote::deseason(rst, use.cpp = TRUE)
  
  # apply pre-whitened mann-kendall test
  file_out <- paste0("out/ndvi/", product)
  
  rst_mk <- gimms::significantTau(mask(rst_dsn, spy_np, inverse = TRUE), 
                                  format = "GTiff", overwrite = TRUE,
                                  filename = paste0(file_out, "_mk_0314.tif"))

  # create r and p layers and mask all pixels inside the national park
  timestamps <- 1:raster::nlayers(rst)
  
  lst_mod <- foreach(type = c("cor", "p")) %do% {  
    rst_param <- raster::calc(rst_dsn, fun = function(x) {
      if (all(is.na(x))) {
        return(NA)
      } else {
        if (type == "cor")
          return(cor(timestamps, x))
        else 
          return(Orcs::pvalue(lm(x ~ timestamps)))
      }
    })
    
    raster::mask(rst_param, spy_np, inverse = TRUE)
  }
  
  # reject non-conclusive trends (p < 0.001; see Detsch et al. 2015)
  rst_sig <- raster::overlay(lst_mod[[1]], lst_mod[[2]], fun = function(x, y) {
    x[y[] >= 0.001] <- NA
    return(x)
  }, filename = paste0("out/ndvi/", product, "_lm_0314.tif"), 
  format = "GTiff", overwrite = TRUE)
  
  return(list(rst_sig, rst_mk))
}


### data processing: rainfall --------------------------------------------------

# list available files
fls_chirps <- list.files("data/CHIRPS/scl", full.names = TRUE)

# restrict temporal range (currently 2003-2014)
id_start <- grep("2003.01", basename(fls_chirps))
id_end <- grep("2014.12", basename(fls_chirps))

fls_chirps <- fls_chirps[id_start:id_end]
rst_chirps <- stack(fls_chirps)

## remove seasonal signal (deseasoning)
rst_chirps_dsn <- deseason(rst_chirps, use.cpp = TRUE)

# create r and p layers
timestamps <- 1:nlayers(rst_chirps)

lst_chirps_mod <- foreach(type = c("cor", "p"), 
                          .packages = c("raster", "Orcs")) %dopar% {  
                            raster::calc(rst_chirps_dsn, fun = function(x) {
                              if (all(is.na(x))) {
                                return(NA)
                              } else {
                                if (type == "cor")
                                  return(cor(timestamps, x))
                                else 
                                  return(Orcs::pvalue(lm(x ~ timestamps)))
                              }
                            })
                          }

# reject non-conclusive trends (p < 0.001)
rst_chirps_sig <- overlay(lst_chirps_mod[[1]], lst_chirps_mod[[2]], 
                          fun = function(x, y) {
                            x[y[] >= 0.001] <- NA
                            return(x)
                          }, filename = "out/chirps/CHIRPS_trends_0314.tif", 
                          format = "GTiff", overwrite = TRUE)


### data processing ------------------------------------------------------------

## deregister parallel backend
stopCluster(supcl)
