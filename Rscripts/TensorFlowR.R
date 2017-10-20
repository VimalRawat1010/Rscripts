##### Tensor Flow for R
library(devtools)
devtools::install_github("rstudio/tfestimators")
library(tfestimators)
install_tensorflow(version = "1.3.0")


# return an input_fn for a given subset of data
mtcars_input_fn <- function(data) {
  input_fn(data, 
           features = c("disp", "cyl"), 
           response = "mpg")
}

mtcars_input_fn(mpg)