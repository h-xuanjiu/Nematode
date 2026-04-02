#' @importFrom utils packageVersion
NULL
.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    paste0(
      "Welcome to ", pkgname, " version ", packageVersion(pkgname), "\n",
      "This package is licensed under the GPL (>= 3.0) license.\n",
      "Please cite: He Y., He J., Feng J., et al. (2026). Soil nematode community structure and composition along a pomegranate plantation chronosequence in the Loess Plateau. Applied Soil Ecology, 219, 106820."
    )
  )
}
