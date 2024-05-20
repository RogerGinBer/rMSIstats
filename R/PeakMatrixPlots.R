#'
#' @importFrom ggplot2 ggplot aes geom_tile vars facet_grid facet_wrap ggtitle
#'
#'
#' @export
plotFeature <- function(se, ft,
                        px = seq_len(ncol(se)),
                        feature = assayNames(se)[1],
                        arrange = c("rows", "cols", "wrap"),
                        scale = c("absolute", "relative"),
                        FUN = identity){
    if (length(ft) != 1) {
        stop("Only one feature can be plotted at the same time. ",
             "length(ft) must be exactly one")
    }
    if(!is.numeric(ft)) {
        stop("ft subindex must be numeric")
    }
    if (!is.character(feature) || length(feature) != 1){
        stop("feature argument must be a character of length 1")
    }
    if (!nrow(se)) stop("SummarizedExperiment has nrow equal to 0")
    if (!is.function(FUN)) stop("FUN must be a valid function")

    .check_mass(se)

    arrange <- match.arg(arrange)
    scale <- match.arg(scale)
    pixelData <- colData(se[,px])
    pixelData[[feature]] <- as.matrix(assay(se[ft, px]))

    pixelData <- FUN(pixelData)

    p <- ggplot(pixelData) +
        geom_tile(aes(x = x, y = y, fill = .data[[feature]])) +
        ggtitle(paste("m/z:", round(rowData(se)$mass[ft], 5)))

    sc <- switch(scale,
                 absolute = "fixed",
                 relative = "free")

    arr <- switch(arrange,
                  rows = facet_grid(rows = vars(sample), scales = sc),
                  cols = facet_grid(cols = vars(sample), scales = sc),
                  wrap = facet_wrap(vars(sample), scales = sc))

    p + arr
}
