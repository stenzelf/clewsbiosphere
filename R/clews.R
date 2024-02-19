#' Get time resolution of a netcdf file
#'
#' Utility function to determine time resolution of a netcdf file
#'
#' @param file_nc netcdf file name
#'
#' @return character string: determined timestep (daily, monthly, yearly)
#'
#' @examples
#' \dontrun{
#' }
#'
#' @export
get_timestep <- function(file_nc){
  tunit <- file_nc$dim$time$units
  tvals <- file_nc$dim$time$vals
  if (grepl("years since", tunit, ignore.case = TRUE)){
    tres <- "annual"
    firstyr <- as.character( as.integer(
      unlist(strsplit(unlist(strsplit(tunit, split = ' ', 
                                      fixed = TRUE))[3], split = '-', fixed = TRUE))[1]) + 
        tvals[1] )
  }else if (grepl("year", tunit, ignore.case = TRUE)){
    tres <- "annual"
    firstyr <- tvals[1]
  }else if (grepl("days since", tunit, fixed = TRUE)){
    ddiff <- tvals[2]-tvals[1]
    firstyr <- as.character( as.integer(
      unlist(strsplit(unlist(strsplit(tunit, split = ' ', 
                                      fixed = TRUE))[3], split = '-', fixed = TRUE))[1]) + 
        floor(tvals[1]/365) )
    if (ddiff > 27 && ddiff < 32){
      tres <- "monthly"
    }else if (ddiff > 364 && ddiff < 367){
      tres <- "annual"
    }else if (ddiff == 1){
      tres <- "daily"
    }else{
      stop("Automatic detection of firstyear and time resolution failed.")
    }
  }else{
    stop("Automatic detection of firstyear and time resolution failed.")
  }
  return(list(tres = tres, firstyr = as.numeric(firstyr)))
}

#' Get main variable of a netcdf file
#'
#' Utility function to guess the main variable of a netcdf file
#'
#' @param file_nc netcdf file name
#'
#' @return character string: determined main variable name
#'
#' @examples
#' \dontrun{
#' }
#'
#' @export
get_main_variable <- function(file_nc){
  for (var in names(file_nc$var)){
    ndims <- file_nc$var[[var]]$ndims
    dimNames <- c()
    for (d in 1:ndims){
      dimNames <- append(dimNames,file_nc$var[[var]]$dim[[d]]$name) 
    }
    #print(paste0(var,paste(dimNames,collapse = ",")))
    if (grepl("lon",paste(dimNames,collapse = " "),ignore.case = T) && grepl("lat",paste(dimNames,collapse = " "),ignore.case = T)){
      return(var)
    }
  }
  print(paste0("None of the variables could certainly be identified as main variable, guessing the last one: ",var))
  return(var)
}

#' Reads netcdf and returns the requested variable as an array
#'
#' Reads an arbitrary netcdf and returns the requested variable in the given 
#' year range as an array
#'
#' @param nc_in_file netcdf file name
#' @param var optional variable to be read, in case automatic detection does 
#'        not work as intended or several variables are stored within the file
#' @param get_year_start first year to be read 
#'        (if not specified will default to first record year)
#' @param get_year_stop final year to be read 
#'        (if not specified will default to last record year)
#' @param suppress_read_print whether to suppress the info "reading file XYZ" 
#'        (default FALSE)
#'
#' @return array with netcdf's data, dim=c([longitude],[latitude],[bands],[months],[years])
#'
#' @examples
#' \dontrun{
#' }
#'
#' @export
read_netcdf <- function(nc_in_file,
                        var = NULL,
                        get_year_start = NULL,
                        get_year_stop = NULL,
                        suppress_read_print = FALSE) {
  
  file_nc <- ncdf4::nc_open(filename = nc_in_file)
  if (is.null(var)) var <- get_main_variable(file_nc = file_nc)
  if (!suppress_read_print) {
    print(paste0("Reading in: ", nc_in_file))
    print(paste0("Attempting to read variable: ",var,
                 ". If this is not correct,",
                 " please correct via argument var."))
  }
  
  ndims <- file_nc$var[[var]]$ndims
  tunit <- file_nc$dim$time$units
  tvals <- file_nc$dim$time$vals
  timing <- get_timestep(file_nc = file_nc)
  tres <- timing$tres
  firstyr <- timing$firstyr
  if (ndims == 3){
    data <- ncdf4::ncvar_get(nc = file_nc, varid = var,
                             start = c(1,1,1), count = c(1,1,-1))
    data_dim <- dim(data)
    nbands <- 1
    tsteps <- data_dim[1]
  }else if (ndims >= 4){
    data <- ncdf4::ncvar_get(nc = file_nc, varid = var,
                             start = c(1,1,1,1), count = c(1,1,-1,-1))
    data_dim <- dim(data)
    nbands <- data_dim[1]
    tsteps <- data_dim[2]
  }else{
    stop("Less than 3 dimensions for file. aborting.")
  }
  if (tres == "annual"){
    nyears <- tsteps
    nmonths <- 1
  }else if (tres == "monthly"){
    nyears <- tsteps/12
    nmonths <- 12
  }
  if (is.null(get_year_start)) get_year_start <- firstyr
  if (is.null(get_year_stop)) get_year_stop <- firstyr + nyears - 1
  ngetyears <- get_year_stop - get_year_start + 1
  
  nlatin <- file_nc$dim$lat$len
  nlonin <- file_nc$dim$lon$len
  lat_values <- file_nc$dim$lat$vals
  lon_values <- file_nc$dim$lon$vals
  outdata <- array(0, dim = c(nlatin,nlonin,nbands,nmonths,ngetyears))
  
  # get spatial extent and resolution
  # this will give a warning if the NetCDF has more than one data field, 
  # e.g. crop bands or time axis
  for (year in get_year_start:get_year_stop){
    for (month in 1:nmonths){
      if (nbands == 1){
        data <- ncdf4::ncvar_get(nc = file_nc, varid = var, count=c(-1,-1,1),
                                 start=c(1,1,((year - firstyr)*nmonths + month)))
        # check whether data needs to be flipped vertically
        # print(paste0("lat: ",lat_values[1],", ",lat_values[2]))
        if (lat_values[1]>lat_values[2]){
          outdata[,,1,month,(year - get_year_start + 1)] <- data[,nlatin:1]
          lat_values <- rev(lat_values)
        }else{
          outdata[,,1,month,(year - get_year_start + 1)] <- data
        }
      }else{ #nbands>1
        data <- ncvar_get(nc = file_nc, varid = var, count=c(-1,-1,-1,1),
                          start=c(1,1,1,((year - firstyr)*nmonths + month)))
        # check whether data needs to be flipped vertically
        if (lat_values[1]<lat_values[2]){
          outdata[,,,month,(year - get_year_start + 1)] <- data[,nlatin:1,]
          lat_values <- rev(lat_values)
        }else{
          outdata[,,,month,(year - get_year_start + 1)] <- data
        }
      }# end if nbands == 1
    }# end for month in 1:nmonths
  }# end for year in ...
  ncdf4::nc_close(file_nc)
  dim(outdata) <- c(lon = nlonin, lat = nlatin, band = nbands, 
                    month = nmonths, year = ngetyears)
  dimnames(outdata) <- list(lon = lon_values, lat = lat_values, band = 1:nbands, 
                            month = 1:nmonths, year = get_year_start:get_year_stop)
  return(drop(outdata))
}

#' Plot an array of lon_lat data to screen or file 
#'
#' Plot an array of lon_lat data to screen or file 
#' e.g. obtained from read_netcdf and averaged over time.
#' Prints to screen if file argument is not supplied (default).
#'
#' @param data array to plot. Needs to be in format [longitude,latitude].
#' @param file character path. file location to save the plot to. 
#'        If not supplied prints to screen. Default: NULL
#' @param title character string printed as title.
#'
#' @return None
#'
#' @examples
#' \dontrun{
#' }
#'
#' @export
plot_lon_lat_array <- function(data, file = NULL, title = ""){
  di <- dim(data)
  if (length(di)>2) stop("Too many dimensions on data object, 
                         please reduce by picking/summing/averaging.")
  par(mar=c(3,3,0,0)) #bltr
  ra <- terra::rast(t(data[,di[2]:1]))
  range <- range(data)
  extent <- terra::ext(c(0, di[1], 30, di[2]))
  if (!is.null(file)) png(file, width=7.25, height=3.5, units="in", res=300, pointsize=6,type="cairo")
  terra::plot(ra, main = title,ext = extent)
  if (!is.null(file)) dev.off()
}

#' Returns cellarea of given netcdf in m
#'
#' Returns approximate cellarea of given netcdf file in m
#'
#' @param nc_in_file netcdf file name
#'
#' @return cellarea array with dimensions [longitude, latitude] in m
#'
#' @examples
#' \dontrun{
#' }
#'
#' @export
get_netcdf_cellsize <- function(nc_in_file = tmp_file){
  file_nc <- ncdf4::nc_open(filename = nc_in_file)
  lat_values <- file_nc$dim$lat$vals
  lon_values <- file_nc$dim$lat$vals
  res.lat <- abs(lat_values[1] - lat_values[2])
  res.lon <- abs(lon_values[1] - lon_values[2])
  nc_close(file_nc)
  
  earthradius <- 6371000.785 # in m
  cellwidth <- earthradius * 2 * pi / 360 # cellwidth per degree at equator in m
  # cells are approx. rectangular at the equator - they get smaller towards 
  # the poles (scaling with the cos of the midpoint latitude - here in radians)
  cellarea_raw <- (cellwidth * res.lon) * (cellwidth * res.lat) *
    cos(lat_values/180 * pi) # in m
  # we want the cellarea in the same array format as tmp and pre
  cellarea <- rep(x = cellarea_raw, each = length(tmp_lon))
  dim(cellarea) <- c(length(tmp_lon), length(tmp_lat))
  return(cellarea)
}

#' Plot the climate of a given location
#'
#' Plot an array of lon_lat data to screen or file 
#' e.g. obtained from read_netcdf and averaged over time.
#' Prints to screen if file argument is not supplied (default).
#'
#' @param tmp temperature array with 12 values for each month in °C
#' @param pre precipitation array with 12 values for each month in mm
#' @param file character path. file location to save the plot to. 
#'        If not supplied prints to screen. Default: NULL
#' @param title character string printed as title.
#'
#' @return None
#'
#' @examples
#' \dontrun{
#' }
#'
#' @export
climate_plot <- function(temp, prec, file = NULL, title = ""){
  par(mar=c(4,4,1,4),oma=c(0,0,0,0)) #bltr
  if (!is.null(file)) png(file, width=7.25, height=3.5, 
                          units="in", res=300, pointsize=6,type="cairo")
  # plot average monthly precipitation for past 30 years
  x_pos <- barplot(height = prec,col = "blue", ylab = "prec in mm",
                   names.arg = c("J","F","M","A","M","J","J","A","S","O","N","D"))
  # plot the average temp. of the past 30 years
  par(fig=c(0,1,0,1))
  lines(x = x_pos, y = temp, ylim = range(temp),xlab = "Month", ylab = "",col = "red", lwd = 2)
  axis(side = 4, col = "red")
  mtext(side = 4, line = 2.5, text = "temp in °C")
  if (!is.null(file)) dev.off()
}

#' Read in the Zhang et al. csv file with NPP measurements
#'
#' Read in the Zhang et al. csv file with NPP measurements and return a 
#' list of relevant records
#'
#' @param csv_in_file character string. full csv file path
#'
#' @return list object with dimensions [lat, lon, totnpp, year] per record
#'
#' @examples
#' \dontrun{
#' }
#'
#' @export
read_in_zhang_file <- function(csv_in_file){
  npp_measurements <- read.csv(file = csv_in_file)
  return(npp_measurements[,c("latitude","longitude","totnpp","year")])
}

#' Add modelled NPP to measured zhang data
#'
#' Add the best fitting modelled NPP to the measured zhang data list object 
#' (e.g. obtained from read_in_zhang_file() function)
#'
#' @param zhang_data zhang data list object 
#'        (e.g. obtained from read_in_zhang_file() function)
#' @param modelled_npp array with modelled NPP (dimensions [lat,lon,year])
#'
#' @return modified zhang_data input list, with added colums: measured, modelled
#'
#' @examples
#' \dontrun{
#' }
#'
#' @export
add_modelled_npp_to_zhang_npp_measurements <- function(zhang_data, modelled_npp) {
  npp_years <- as.integer(dimnames(modelled_npp)$year)
  npp_years_length <- length(npp_years)
  min_model_year <- min(npp_years)

  years <- unique(zhang_data$year)
  #      Average,1991,1959-98,1988-90,Potential,1993,1986-93,1982-89,1983-90,1992,1982-96,1988,latest interval,1990-93,1987-97
  from <- c(  51,  91,     59,     88,       NA,  93,     86,     82,     83,  92,     82, 88,              NA,     90,     87) - (min_model_year-1900)
  to <-   c( 122,  91,     98,     90,       NA,  93,     93,     89,     90,  92,     96, 88,              NA,     93,     97) - (min_model_year-1900)
  from[from<1] <- 1
  from[from>npp_years_length] <- npp_years_length
  to[to<1] <- 1
  to[to>npp_years_length] <- npp_years_length
  npp_measured <- zhang_data$totnpp
  npp_modelled <- array(NA,dim=length(npp_measured))
  for (i in 1:length(npp_measured)){
    if (zhang_data$year[i] %in% c("Potential","latest interval") 
        | zhang_data$totnpp[i]<0 
        | is.na(zhang_data$totnpp[i])) {
      npp_measured[i] <- NA
      npp_modelled[i] <- NA
      next
    }
    latPick <- round((zhang_data$latitude[i] - 0.25) * 2) / 2 + 0.25
    lonPick <- round((zhang_data$longitude[i] - 0.25) * 2) / 2 + 0.25
    ind <- which(years == zhang_data$year[i])
    npp_modelled[i] <- mean(NPP[as.character(lonPick),as.character(latPick),from[ind]:to[ind]])
  }
  zhang_data$measured <- npp_measured
  zhang_data$modelled <- npp_modelled
  return(zhang_data)
}

#' Create a scatterplot
#'
#' Create a scatterplot of two data series e.g. measured vs. modelled
#'
#' @param x data to plot on x axis
#' @param y data to plot on y axis
#' @param file character path. file location to save the plot to. 
#'        If not supplied prints to screen. Default: NULL
#' @param title character string printed as title.
#' @param xlab label for x axis
#' @param ylab label for y axis
#'
#' @return None
#'
#' @examples
#' \dontrun{
#' }
#'
#' @export
scatter_plot <- function(x, y, file = NULL, title = "", xlab = "", ylab = ""){
  if (!is.null(file)) png(file, width=7.25, height=3.5, 
                          units="in", res=300, pointsize=6,type="cairo")
  par(mar=c(4.5,4.5,1,1),bty = "o")
  max = max(y,x, na.rm = T)
  plot(x = x, y = y, ylim = c(0,max), xlim = c(0,max), xlab = xlab, 
       ylab = ylab, asp=1, main = title)
  abline(a = 0,b = 1)
  fit <- summary(lm(x~y))
  text(x = 0, y = max*0.9,paste("R^2 =",round(fit$r.squared,4)),adj=c(0,1))
  if (!is.null(file)) dev.off()
}
