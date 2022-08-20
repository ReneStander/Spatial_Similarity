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
# Example 2: On simulated marked point patterns
########

## For this example we are going to simulate two marked point patterns
## The main goal for this example is to simulate a similar point pattern 

# Define a 10x10 square window
window <- square(r = 10)

# Plot the window
plot(window, main = "A 10x10 square window")

# Simulate a random Poisson point pattern with an intensity of 1
set.seed(123) # setting a seed value
firstpp <- rpoispp(lambda = 1, win = window) # simulating the random point pattern
plot(firstpp, main = "Point pattern 1 - only locations", pch = 19, cex = 0.6)

# Let us simulate a mark for each point from a N(20,16) distribution

# first save the point pattern as a data frame
firstpp.df <- as.data.frame(firstpp)
# view the first 6 rows of the data set
head(firstpp.df)

# add the marks 
firstpp.df$mark <- rnorm(nrow(firstpp.df), mean = 20, sd = 4)
# view the first 6 rows of the data set
head(firstpp.df)

# convert back to a point pattern
firstmpp <- as.ppp(firstpp.df, W = window)
# plot the marked point pattern
plot(firstmpp, main = "Marked point pattern 1")

# Let us change 10% of the point's marks in the point pattern
# the marks that we change will be simulated from N(30,16)

# use the data frame we created with point pattern 1
secondpp.df <- firstpp.df

# points to change the marks
(to_change <- sample(1:nrow(firstpp.df), size = floor(nrow(firstpp.df)*0.1), replace = FALSE))
secondpp.df$mark[to_change] <- rnorm(to_change, mean = 30, sd = 4)

# convert to a point pattern
secondmpp <- as.ppp(secondpp.df, W = window)
# plot the marked point pattern
plot(secondmpp, main = "Marked point pattern 2")

# Let us plot the two marked point patterns next to each other
par(mfrow = c(1,2))
plot(firstmpp, main = "Marked point pattern 1")
plot(secondmpp, main = "Marked point pattern 2")


par(mfrow = c(1,1))
####################################################
# Do the similarity test

####
# STEP 1
####

# Pixel representation for `Marked point pattern 1`
mpp_1_density <- continuous_marked_density(data = as.data.frame(firstmpp), # point pattern to compare
                                       gr = 15, # grid size
                                       wind = window, # point pattern window
                                       wind_type = "rectangular") # window type

plot(mpp_1_density)

# Pixel representation for 'Component 2'
mpp_2_density <- continuous_marked_density(data = as.data.frame(secondmpp), # point pattern to compare
                                       gr = 15, # grid size
                                       wind = window, # point pattern window
                                       wind_type = "rectangular") # window type

plot(mpp_2_density)

# plot the two densities next to each other
par(mfrow = c(1,2))
plot(mpp_1_density, main = "Marked point pattern 1 density")
plot(mpp_2_density, main = "Marked point pattern 2 density")


par(mfrow = c(1,1))
####
# STEP 2
####

# Calculate the SSIM values

mpp_1_2_ssim <- SSIM_cont(imx = mpp_1_density, # pixel density 1
                               imy = mpp_2_density, # pixel density 2
                               sl = 3, # size of sliding window
                               alpha = 1, # parameters (default)
                               beta = 1, # parameters (default)
                               gamma = 1) # parameters (default)

plot(as.im(mpp_1_2_ssim))

####
# STEP 3
####

(result <- global_index(mpp_1_2_ssim))


########
# On your own
########

# Explore the following:
# 1. Play around with the simulation, specifically the 
#   a. simulation of the marks
#   b. simulation of the point pattern (try different number of points and intensities)
# 2. Play around with the grid size when creating the pixel image representation
# 3. Play around with the size of the sliding window when calculating the local similarity map

