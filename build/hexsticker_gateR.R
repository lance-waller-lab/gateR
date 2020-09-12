# ----------------------------------------------------- #
# HexSticker for the {gateR} package
#
# Created by: Ian Buller, Ph.D., M.A. (GitHub: @idblr)
# Created on: September 10, 2020
#
# Recently modified by:
# Recently modified on:
#
# Notes:
# A) Uses the "hexSticker" package
# B) Modified image: https://publicdomainvectors.org/en/free-clipart/Crocodile-in-egg-shell/78578.html
# ----------------------------------------------------- #

# Packages
library(hexSticker)

# Image file
path_image <- "./man/figures/gator.png"

# Create hexSticker
## the alligator palette https://www.color-hex.com/color-palette/33955
## And colors from image https://html-color-codes.info/colors-from-image/
s <- hexSticker::sticker(subplot = path_image,
                         package = "gateR", p_size = 8, p_color = "#577D45",
                         s_x = 1, s_y = 0.8, s_width = 0.8, s_height = 0.8,
                         h_fill = "#ffedb2",
                         h_color = "#577D45",
                         dpi = 1000,
                         filename = "./man/figures/gateR.png",
                         white_around_sticker = F)
# -------------------- END OF CODE -------------------- #
