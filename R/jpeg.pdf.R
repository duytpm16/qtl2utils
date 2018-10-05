jpeg_pdf = function(pdfname, mypattern = "MYTEMPJPEG", ...) {
  fname = paste0(mypattern, "%05d.jpeg")
  gpat = paste0(mypattern, ".*\\.jpeg")
  takeout = list.files(path = tempdir(), pattern = gpat, full.names = TRUE)
  if (length(takeout) > 0)
    file.remove(takeout)
  jpegname = file.path(tempdir(), fname)
  jpeg(jpegname, ...)
  return(list(pdfname = pdfname, mypattern = mypattern))
}
