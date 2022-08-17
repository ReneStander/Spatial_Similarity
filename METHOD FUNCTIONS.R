##################
# METHOD FUNCTIONS
##################

###############################
### STEP 1: PIXEL IMAGE 
###############################

#####
# GEOSTATISTICAL
#####

geostatistical_pixel <- function(data, window, resolution = 20){
  # make sure data is a data frame and window is of type owin
  data <- as.data.frame(data)
  window <- as.owin(window)
  
  # get the grid and the centroids
  grid <- quadrats(window, nx = resolution)
  centres <- gridcentres(grid, nx = resolution, ny = resolution)
  centres.df <- as.data.frame(centres)
  # determine which centres are inside the window
  centres.in <- as.data.frame(as.ppp(centres, W = window))
  coordinates(centres.in) <- ~ x + y
  
  # krige
  coordinates(data) <- ~ x + y
  vgm <- variogram(log(data$val)~1, locations = data)
  fit <- fit.variogram(vgm, model=vgm(NA, "Sph", NA, NA))
  kriged <- krige(log(data$val) ~ 1, data, centres.in, model=fit)
  kriged.df <- as.data.frame(kriged)
  kriged.df <- kriged.df[,c("x", "y", "var1.pred")]
  colnames(kriged.df) <- c("x", "y", "density")
  
  kriged.im <- as.im(matrix(kriged.df$density, byrow = TRUE, nrow = resolution), W = window)
  return(kriged.im)
}

#####
# LATTICE
#####

lattice_pixel <- function(data, grid_df, wind = enclosed1){
  merged <- merge(grid_df, data, all.x = TRUE)
  merged_subset <- subset(merged, select = -c(lattice))
  sorted <- merged_subset[order(merged_subset$y,merged_subset$x),]
  density <- as.im(matrix(sorted$value, byrow = TRUE, nrow = sqrt(nrow(merged))), W = wind)
  return(density)
}

#####
# POINT PATTERNS
#####

# UNMARKED
unmarked_density <- function(data, gr = 10, wind = win_rec, wind_type = "rectangular"){
  
  # divide the window of the pattern into a grid
  grid <- quadrats(wind, nx = gr)
  
  # take the centre of each grid cell as the spatial unit where the density will be calculated
  grid.centres <- gridcentres(wind, nx = gr, ny = gr)
  # - convert it to a dataframe
  grid.centres <- data.frame(gc_x = grid.centres$x, gc_y = grid.centres$y)
  
  mpp <- ppp(data$x, data$y, window = wind)
  freqmat <- quadratcount(mpp, nx = gr)
  if(wind_type == 'convex'){
    t <- freqmat
    tt <- as.data.frame(t)
    tt$tile <- gsub("Tile ", "", tt$tile)
    tt$tile <- gsub("row ", "", tt$tile)
    tt$tile <- gsub(", col ", " ", tt$tile)
    elems <- unlist( strsplit( tt$tile , " " ) )
    m <- matrix( elems , ncol = 2 , byrow = TRUE )
    tt$row <- m[,1]
    tt$col <- m[,2]
    freqmat <- matrix(NA, nrow = gr, ncol = gr)
    for(i in 1:nrow(tt)){
      freqmat[gr+1-as.numeric(tt$row[[i]]), as.numeric(tt$col[[i]])] <- tt$Freq[i]
    }
  }
  
  grid.centres_x <- as.character(levels(as.factor(grid.centres$gc_x)))
  grid.centres_y <- as.character(levels(as.factor(grid.centres$gc_y)))
  grid.centres$indx <- NA
  grid.centres$indy <- NA
  # give an index to each grid centre coordinate
  for(i in 1:gr){
    d1 <- grid.centres_x[i]
    d2 <- grid.centres_y[i]
    grid.centres$indx[(as.character(grid.centres$gc_x) == d1)] <- i
    grid.centres$indy[(as.character(grid.centres$gc_y) == d2)] <- i
  }
  grid.centres$gc_x <- as.numeric(grid.centres$gc_x)
  grid.centres$gc_y <- as.numeric(grid.centres$gc_y)
  
  # bandwidth of the pattern
  smoothing_bw <- bw.diggle(mpp, hmax = 5)
  bw_matrix <- matrix(c(smoothing_bw, 0, 0, smoothing_bw), nrow = 2)
  
  # create a data frame from the given pattern
  df_use <- data
  # number of points in the pattern
  n <- nrow(df_use)
  # create an empty column in the grid centres data frame that will contain the densities
  grid.centres$density <- NA
  
  # calculate the Nadaraya-Watson smoother for each spatial location
  for(j in 1:nrow(grid.centres)){
    d <- grid.centres[j,1:2]
    d_indy <- grid.centres[j,3]
    d_indx <- grid.centres[j,4]
    if(!is.na(freqmat[d_indx, d_indy])){
      dd <- cbind(d[,1] - df_use[,1], d[,2] - df_use[,2])
      kernel <- apply(dd,1,dmvnorm, mean = c(0,0), sigma = bw_matrix)
      grid.centres$density[j] <- sum(kernel)
    }
  }
  
  ## standardise the density
  # estimate the area under the density and divide each density by that
  # this is so that the area under the used density is equal to 1
  gridcell_size <- area(wind)/((gr^2))
  grid.centres$density[is.na(grid.centres$density)] <- 0
  totalarea_estimate <- sum(grid.centres$density*gridcell_size)
  grid.centres$density <- grid.centres$density/totalarea_estimate
  grid.centres$density[grid.centres$density == 0] <- NA 
  
  density <- as.im(matrix(grid.centres$density, byrow = TRUE, nrow = gr), W = wind)
  return(density)
}

# DISCRETE MARKED
calculate_mode <- function(x) {
  uniqx <- unique(na.omit(x))
  uniqx[which.max(tabulate(match(x, uniqx)))]
}

discrete_marked_density <- function(data, resolution, wind = win_rec){
  data$value <- as.factor(data$value)
  grid.centres <- gridcentres(wind, nx = resolution, ny = resolution)
  grid.centres <- data.frame(gc_x = grid.centres$x, gc_y = grid.centres$y)
  k <- floor(nrow(data)*0.1)
  
  
  grid.centres$dens_val <- as.factor(NA)
  levels(grid.centres$dens_val) <- levels(data$value)
  for(i in 1:nrow(grid.centres)){
    d <- grid.centres[i,1:2]
    
    dist <- as.data.frame(t(cdist(d,data[,1:2])))
    colnames(dist) <- "dist"
    dist$type <- data$value
    
    dist <- dist[order(dist$dist),]
    grid.centres$dens_val[i] <- calculate_mode(dist[1:k,2])[1]
  }
  
  density <- as.im(matrix(grid.centres$dens_val, byrow = TRUE, nrow = resolution), W = wind)
  return(density)
}

# CONTINUOUS MARKED
continuous_marked_density <- function(data, gr = 10, wind = win_rec, wind_type = 'rectangular'){
  
  # divide the window of the pattern into a grid
  grid <- quadrats(wind, nx = gr)
  
  # take the centre of each grid cell as the spatial unit where the density will be calculated
  grid.centres <- gridcentres(wind, nx = gr, ny = gr)
  # - convert it to a dataframe
  grid.centres <- data.frame(gc_x = grid.centres$x, gc_y = grid.centres$y)
  
  mpp <- ppp(data$x, data$y, window = wind, marks = data$value)
  freqmat <- quadratcount(mpp, nx = gr)
  if(wind_type == 'convex'){
    t <- freqmat
    tt <- as.data.frame(t)
    tt$tile <- gsub("Tile ", "", tt$tile)
    tt$tile <- gsub("row ", "", tt$tile)
    tt$tile <- gsub(", col ", " ", tt$tile)
    elems <- unlist( strsplit( tt$tile , " " ) )
    m <- matrix( elems , ncol = 2 , byrow = TRUE )
    tt$row <- m[,1]
    tt$col <- m[,2]
    freqmat <- matrix(NA, nrow = gr, ncol = gr)
    for(i in 1:nrow(tt)){
      freqmat[gr+1-as.numeric(tt$row[[i]]), as.numeric(tt$col[[i]])] <- tt$Freq[i]
    }
  }
  
  grid.centres_x <- as.character(levels(as.factor(grid.centres$gc_x)))
  grid.centres_y <- as.character(levels(as.factor(grid.centres$gc_y)))
  grid.centres$indx <- NA
  grid.centres$indy <- NA
  # give an index to each grid centre coordinate
  for(i in 1:gr){
    d1 <- grid.centres_x[i]
    d2 <- grid.centres_y[i]
    grid.centres$indx[(as.character(grid.centres$gc_x) == d1)] <- i
    grid.centres$indy[(as.character(grid.centres$gc_y) == d2)] <- i
  }
  grid.centres$gc_x <- as.numeric(grid.centres$gc_x)
  grid.centres$gc_y <- as.numeric(grid.centres$gc_y)
  
  # bandwidth of the pattern
  smoothing_bw <- bw.diggle(mpp)
  bw_matrix <- matrix(c(smoothing_bw, 0, 0, smoothing_bw), nrow = 2)
  
  # create a data frame from the given pattern
  df_use <- data
  # number of points in the pattern
  n <- nrow(df_use)
  # create an empty column in the grid centres data frame that will contain the densities
  grid.centres$denom <- NA
  grid.centres$num <- NA
  
  # calculate the Nadaraya-Watson smoother for each spatial location
  for(j in 1:nrow(grid.centres)){
    d <- grid.centres[j,1:2]
    d_indy <- grid.centres[j,3]
    d_indx <- grid.centres[j,4]
    if(!is.na(freqmat[d_indx, d_indy])){
      dd <- cbind(d[,1] - df_use[,1], d[,2] - df_use[,2])
      kernel <- apply(dd,1,dmvnorm, mean = c(0,0), sigma = bw_matrix)
      grid.centres$num[j] <- sum(df_use$value*kernel)
      grid.centres$denom[j] <- sum(kernel)
    }
  }
  
  grid.centres$density <- grid.centres$num/grid.centres$denom
  ## standardise the density
  # estimate the area under the density and divide each density by that
  # this is so that the area under the used density is equal to 1
  gridcell_size <- area(wind)/((gr^2))
  grid.centres$density[is.na(grid.centres$density)] <- 0
  totalarea_estimate <- sum(grid.centres$density*gridcell_size)
  grid.centres$density <- grid.centres$density/totalarea_estimate
  grid.centres$density[grid.centres$density == 0] <- NA 
  
  density <- as.im(matrix(grid.centres$density, byrow = TRUE, nrow = gr), W = wind)
  return(density)
}

###############################
### STEP 2: SIMILARITY MAP 
###############################

SSIM_cont <- function(imx, imy, sl, alpha = 1, beta = 1, gamma = 1){
  # determine the dimensions of the images and make sure they are the same
  dimx <- dim(imx)
  dimy <- dim(imy)
  if(all(dimx) != all(dimy)){
    stop('The input images should have the same dimensions')
  }
  # convert the images to matrices
  imx_mat <- imx$v
  imy_mat <- imy$v
  # calculate the dynamic range
  min_x <- min(imx_mat,na.rm=TRUE)
  max_x <- max(imx_mat,na.rm=TRUE)
  min_y <- min(imy_mat,na.rm=TRUE)
  max_y <- max(imy_mat,na.rm=TRUE)
  # constants
  L <- max(max_x, max_y) - min(min_x, min_y)
  K1 <- 0.01
  K2 <- 0.03
  c1 <- (K1*L)^2
  c2 <- (K2*L)^2
  c3 <- c2/2
  gr <- dimx[1]
  
  # define an empty matrix for the ssim value to be stored to
  res <- matrix(NA, nrow = gr, ncol = gr)
  for(i in 1:gr){ # row indices
    for(j in 1:gr){ # column indices
      if(!is.na(imx_mat[i, j])){
        # calculate the subset of the matrix
        toeachside <- (sl - 1)/2
        toplim <- (i - toeachside)
        botlim <- i + toeachside
        leftlim <- j - toeachside
        rightlim <- j + toeachside
        if(toplim < 1) toplim = 1
        if(botlim > gr) botlim = gr
        if(leftlim < 1) leftlim = 1
        if(rightlim > gr) rightlim = gr
        subx <- imx_mat[(toplim:botlim),(leftlim:rightlim)]
        suby <- imy_mat[(toplim:botlim),(leftlim:rightlim)]
        # mean
        subx_mean <- mean(subx,na.rm=TRUE)
        suby_mean <- mean(suby,na.rm=TRUE)
        # standard deviation
        subx_sd <- sd(subx,na.rm=TRUE)
        suby_sd <- sd(suby,na.rm=TRUE)
        # covariance
        subx<-subx[!is.na(subx)]
        suby<-suby[!is.na(suby)]
        subx_suby_covariance <- cov(as.vector(subx), as.vector(suby))
        # luminance
        lum <- (2*subx_mean*suby_mean + c1)/(subx_mean^2 + suby_mean^2 + c1)
        # contrast
        con <- (2*subx_sd*suby_sd + c2)/(subx_sd^2 + suby_sd^2 + c2)
        12
        # structure
        str <- (subx_suby_covariance + c3)/(subx_sd * suby_sd + c3)
        # ssim
        res[i,j] <- (lum^alpha)*(con^beta)*(str^gamma)
      }
    }
  }
  return(res)
}

ssim_discr <- function(imx, imy){
  # determine the dimensions of the images and make sure they are the same
  dimx <- dim(imx)
  dimy <- dim(imy)
  if(all(dimx) != all(dimy)){
    stop('The input images should have the same dimensions')
  }
  # convert the images to matrices
  imx_mat <- imx$v
  imy_mat <- imy$v
  # direct comparison
  res <- 1*(imx_mat == imy_mat)
  
  return(res)
}


###############################
### STEP 3: GLOBAL INDEX 
###############################

global_index <- function(sim_map){
  ind <- mean(sim_map, na.rm = TRUE)
  return(ind)
}





