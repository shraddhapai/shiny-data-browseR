# #########################################
# makeColorLighter.R
# An R utility to make an RGB colour lighter. Useful in situations where we want to 
# plot a lighter confidence interval band across the mean trendline of a dataseries.
#
# Requires: hsv2rgb.R also part of Rcolorutil
# 
# Example usage:
# > pal <- c("#D7191C","#FDAE61","#ABDDA4","#2B83BA")
# > y <- sapply(pal,.makeLighter)
# > par(mfrow=c(2,1))
# > barplot(1:length(pal), col=pal,main="Original colours") ## original colours
# > barplot(1:length(y),col=y, main="Lighter colours")
# Copyright 2014 Shraddha Pai <Shraddha.Pai@camh.ca>
#
# LICENSE:
# makeColorLighter.R is part of "Rcolorutil"
#  "Rcolorutil" is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
# 
#  "Rcolorutil" is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
makeColorLighter <- function(rgb_col) {
  source("hsv2rgb.R")
  
  x <- c(
    strtoi(paste("0x",substr(rgb_col,2,3),sep="")),
    strtoi(paste("0x",substr(rgb_col,4,5),sep="")),
    strtoi(paste("0x",substr(rgb_col,6,7),sep=""))
  );
  #cat(sprintf("RGB = { %i, %i, %i} \n", x[1],x[2],x[3]))
  y <- rgb2hsv(x[1],x[2],x[3],max=255)
  #cat("HSV : Before")
  #print(y)
  y[3] <-  min(1,y[3] + 0.4);  
  y[2] <- y[2]/3 # increase brightness, lower saturation
  #cat("HSV: After")
  #print(y)
  z <- hsv2rgb(y[1]*360, y[2],y[3])
  z
  ### (RGB) returns lighter RGB object
}
