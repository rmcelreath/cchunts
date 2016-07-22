.onLoad <- function(libname, pkgname) { }

.onAttach <- function(...) {
  theLib <- dirname(system.file(package = "cchunts"))
  pkgdesc <- packageDescription("cchunts", lib.loc = theLib)
  builddate <- gsub(';.*$', '', pkgdesc$Packaged)
  packageStartupMessage(paste("cchunts (Version ", pkgdesc$Version, ")\nType '?cchunts' for list of data sets.", sep = ""))
} 