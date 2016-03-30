
# set the number of colors you need.
num <- 38


myPalette <- colorRampPalette(rev(brewer.pal(num, "Spectral")))
sf <- scale_fill_manual(values = myPalette(num) )
# methods(sf)
# ?scale_fill_manual
# attributes(sf)


myString <- toString(sf$palette(num))
# myString <- gsub('#', '', myString)
myString <- gsub("([#A-Za-z0-9]+)", "'\\1'", myString)
print(myString, quote=FALSE)