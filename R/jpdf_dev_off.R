jpdf_dev_off = function(pdfname, mypattern, ...) {
  dev.off()
  require(magick)
  gpat = paste0(mypattern, ".*\\.jpeg")
  jpegs = list.files(path = tempdir(), pattern = gpat, full.names = TRUE)
  pdfs <- image_convert(image = image_read(jpegs), format = 'pdf')
  image_write(pdfs, pdfname, format = 'pdf', ...)
}
