#### Spatial Similarity

####
# packages needed
library(matlib)
library(expm)
library(spatstat)
library(mvtnorm)
library(sampling)
####

# Import the functions
source("METHOD FUNCTIONS.R")

########
# Example 1: On unmarked point patterns
########


# Use 'waterstriders' built-in data set in spatstat library
# plot the data
plot(waterstriders)

# From the plot we can see that there are three components. 
# Let us compare component 1 with component 2

# Extract component 1
waterstriders_1 <- waterstriders[[1]]
plot(waterstriders_1, main = "Component 1")
# Extract component 2
waterstriders_2 <- waterstriders[[2]]
plot(waterstriders_2, main = "Component 2")


####
# STEP 1
####

# Pixel representation for 'Component 1'
waterstr_1_density <- unmarked_density(data = as.data.frame(waterstriders_1), # point pattern to compare
                                       gr = 15, # grid size
                                       wind = waterstriders_1$window, # point pattern window
                                       wind_type = "rectangular") # window type

plot(waterstr_1_density)

# Pixel representation for 'Component 2'
waterstr_2_density <- unmarked_density(data = as.data.frame(waterstriders_2), # point pattern to compare
                                       gr = 15, # grid size
                                       wind = waterstriders_2$window, # point pattern window
                                       wind_type = "rectangular") # window type

plot(waterstr_2_density)


####
# STEP 2
####

# Calculate the SSIM values

waterstr_1_2_ssim <- SSIM_cont(imx = waterstr_1_density, # pixel density 1
                               imy = waterstr_2_density, # pixel density 2
                               sl = 3, # size of sliding window
                               alpha = 1, # parameters (default)
                               beta = 1, # parameters (default)
                               gamma = 1) # parameters (default)

plot(as.im(waterstr_1_2_ssim))

####
# STEP 3
####

(result <- global_index(waterstr_1_2_ssim))

# This similarity index is quite low
# therefore we can say that the spatial data sets are not similar


########
# On your own
########

# Replicate the above examples by:
# 1. comparing component 2 with component 3 in the `waterstriders` data.
# 2. comparing component 1 with component 3 in the `waterstriders` data.
# 3. try find your own data sets to compare! You can either simulate some spatial data or use other built-in data sets. 





