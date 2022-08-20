# Spatial_Similarity
This repository is used for the code on the Spatial Similarity content. 

A generic similarity test for spatial data is included in this repository. This test is able to handle any type of spatial data, namely geostatistical data, lattice data and point patterns. It consists of three steps and compare two spatial data sets with each other.

* Step 1: Create a pixel image representation
* Step 2: Compare the pixel image representations with the SSIM index - this creates a similarity map
* Step 3: Calculate a global similarity index

The similarity map created in Step 2 can also be used as a local similarity index.

## METHOD FUNCTIONS.R

This script contains all the functions needed to perform the test.

### Step 1

* geostatistical_pixel: creates a pixel image representation for geostatistical data
* lattice_pixel: creates a pixel image representation for lattice data
* unmarked_density: creates a pixel image representation for unmarked point patterns
* discrete_marked_density: creates a pixel image representation for marked point patterns with discrete marks
* continuous_marked_density: creates a pixel image representation for marked point patterns with continuous marks

### Step 2

* SSIM_cont: creates a local similarity map with continuous values
* SSIM_dicr: creates a local similarity map with discrete values

### Step 3

* global_index: calculated the global similarity index

## Example 1.R

This script contains an example where the method is applied to unmarked point patterns. 

## Example 2.R 

This script contains an example where the method is applied to marked point patterns. For this example, the data is simulated. 
