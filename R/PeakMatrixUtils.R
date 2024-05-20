#'
#'
#'@importFrom MsCoreUtils closest
#'@export

findMzIndex <- function(se, mz, tol = Inf, ppm = 0, stop.na = TRUE){
    if (ppm) tol <- 0
    .check_mass(se)
    idx <- closest(mz, rowData(se)$mass, tolerance = tol, ppm = ppm,
            duplicates = "closest")
    if (is.na(idx)) {
        if(stop.na) stop("No mass hits found. Try increasing 'tol' or 'ppm'")
        return(numeric())
    }
    idx
}


.check_mass <- function(se){
    if (!"mass" %in% colnames(rowData(se))) {
        stop("SummarizedExperiment must have 'mass' in rowData")
    }
}
