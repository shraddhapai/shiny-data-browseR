# #########################################
# hsv2rgb.R
# A simple R utility to convert a colour from HSV space to RGb space
# Implementation based on formula from Wikipedia:
# http:// en.wikipedia.org/wiki/HSL_and_HSV#From_HSV
#
# Example usage:
# hsv2rgb(162,0.25,1.0)
#
# Copyright 2014 Shraddha Pai <Shraddha.Pai@camh.ca>
#  hsv2rgb.R is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
# 
#  hsv2rgb.R is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.

hsv2rgb <- function(h, ##<<(numeric between 0-360) hue
                    s, ##<<(numeric between 0 and 1) saturation
                    v  ##<<(numeric between 0 and 1) value
                    ) { # Converts HSV to RGB. 

hPrime <- h/60;
C <- v*s; # chroma
X <- C * ( 1 - abs((hPrime %% 2) - 1))

i <- floor(hPrime)
vec <- switch(as.character(i),  
  "0"=c(C,X,0),
  "1"=c(X,C,0),
  "2"=c(0,C,X),
  "3"=c(0,X,C),
  "4"=c(X,0,C),
  "5"=c(C,0,X)
);
m <- v - C;
vec <- vec + m;

return(rgb(vec[1],vec[2],vec[3]))
### RGB object
}