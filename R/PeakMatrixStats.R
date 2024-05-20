#'
#' @importFrom HDF5Array loadHDF5SummarizedExperiment
#' @importFrom SummarizedExperiment colData rowData assayNames
#' @export
loadPeakMatrix <- function(folder){
    se <- loadHDF5SummarizedExperiment(folder)
    .validate_peakmat_summ_exp(se)
}

.validate_peakmat_summ_exp <- function(se){
    msg <- c()
    cols <- colnames(colData(se))
    if (!all(c("x", "y") %in% cols)){
        msg <- c(msg, "Missing coordinates data (x,y) in the colData")
    }
    if (!"sample" %in% cols){
        msg <- c(msg, "Missing sample information in the colData")
    }

    feature_info <- colnames(rowData(se))
    if (!"mass" %in% feature_info){
        msg <- c(msg, "Missing m/z (mass) column in the rowData")
    }
    if(!length(assayNames(se))){
        msg <- c(msg, "No assay data detected in the assays slot")
    }
    if (length(msg)){
        stop(paste(c("The following errors have been found:", msg),
                   collpase = "\n"))
    }
    se
}


#### Dimensionality Reduction ####

#' @importFrom fields rdist
matrixReductionOneFile <- function(mat, coords, reductionX = 2, reductionY = 2, sigma = 0.5, tol = 1e-3){
    rangeX <- range(coords$x)
    rangeY <- range(coords$y)
    breaksX <- seq(rangeX[1], rangeX[2],
                   length.out = length(unique(coords$x)) / reductionX)
    breaksY <- seq(rangeY[1], rangeY[2],
                   length.out = length(unique(coords$y)) / reductionY)
    new_mat <- matrix(0, ncol = length(breaksX) * length(breaksY), nrow = nrow(mat))
    new_coords <- expand.grid(breaksX, breaksY)
    colnames(new_coords) <- c("x", "y")
    coords_num <- coords[,c("x", "y")]
    for (i in seq(nrow(new_coords))){
      d <- as.vector(rdist(new_coords[i,,drop = FALSE], coords_num))
      if(all(is.na(d))) next
      dn <- dnorm(d, sd = sigma)
      dn[dn < tol] <- 0
      if(!any(dn > tol)) next
      idx <- which(dn > tol)
      a <- t(dn[idx] %*% t(mat[,idx])) / sum(dn[idx])
      new_mat[,i] <- a
    }
    sums <- colSums(new_mat)

    ## I don't like this
    # if(!all(sums == 0)){
    #     new_mat <- new_mat[, sums > 0]
    #     new_coords <- new_coords[, sums > 0]
    # }

    ## Ugly, TODO let's think how to handle sample outside this function
    new_coords$sample <- coords$sample[1]
    return(list(mat = new_mat, coords = new_coords))
  }

#'
#' @importFrom SummarizedExperiment assay assay<- assays assayNames
#' @importFrom S4Vectors DataFrame
#' @export
#'
reducePixels <- function(se, ...){
    files <- unique(se$sample)
    sample_idx <- lapply(files, function(x){which(se$sample == x)})

    alist <- list()
    anames <- assayNames(se)
    for (an in anames){
        out <- lapply(sample_idx, function(x){
            gc()
            m <- matrixReductionOneFile(
                    as.matrix(assay(se, an)[,x]),
                    colData(se)[x, c("sample", "x", "y")],
                    ...)
        })
        alist[[an]] <- do.call(cbind, lapply(out, function(x) x[[1]]))
        new_coords <- do.call(rbind, lapply(out, function(x) x[[2]]))
    }

    ## Empty assays first to avoid dim mismatch error
    se@assays <- NULL
    se@colData <- DataFrame(new_coords)
    for (an in anames){
        assay(se, an) <- alist[[an]]
    }
    se
}


#### Correlation analysis ####

#'
#' @export
#'
massDifferenceTable <- function(se, ft, maxft = 5e3){
    .check_mass(se)
    if (length(ft) > maxft) {
        stop("Too many features have been selected. ",
             "If you want to proceed, please raise the maxft level.")
    }
    if (!all(ft %in% seq_len(nrow(se)))) {
        stop("Some feature indices cannot be found on the SummarizedExperiment")
    }
    as.matrix(dist(rowData(se[ft,])$mass, method = "manhattan", ))
}


block_correlation <- function(mat, method = "pearson", chunkWidthSize=1){
    method <- match.arg(method)
    Nfeatures = nrow(mat)
    m <- matrix(0, nrow = Nfeatures, ncol=Nfeatures)
    Nchunk = ceiling(Nfeatures / mat@seed@chunkdim[1] / chunkWidthSize)
    chunks <- cut(seq(Nfeatures), breaks = Nchunk, labels = F)
    for(i in unique(chunks)){
        iind <- chunks == i
        xmat <- t(as.matrix(mat[iind,]))
        for(j in unique(chunks)){
            print(paste(i,j))
            jind <- chunks == j
            if(i < j) next
            m_sub <- cor(xmat, t(as.matrix(mat[jind,])), method = method)
            m[iind, jind] <- m_sub
            m[jind, iind] <- t(m_sub) #It's symmetric
        }
        gc()
    }
    return(m)
}

featureCorrelation <- function(se){


}









