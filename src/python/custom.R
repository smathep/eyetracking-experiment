# A nice function to load packages.
# It checks if the required libraries are installed. If not installs them quietly and load to workspace.
# At the end it prints all loaded packages(libraries).
load.libraries <- function(libs,libpath) {
  for (package in libs) {
    if (!require(package, character.only = TRUE, quietly = TRUE)) {
      install.packages(package, type="binary",
                       repos = "http://stat.ethz.ch/CRAN/",
                       lib = libpath,
                       dependencies = TRUE)
      library(package, character.only = TRUE)
    }
  }
  print('Loaded packages: ')
  print(.packages())
}
