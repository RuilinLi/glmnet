#' @export
wls_dense_cpp <- function(alm0, almc, alpha, m, no, ni, x, r, v, intr, ju, vp, cl, 
    nx, thr, maxit, a, aint, g, ia, iy, iz, mm, nino, rsqc, nlp, jerr) {
    result <- .Call("wls_dense", alm0, almc, alpha, m, no, ni, x, r, v, intr, ju, 
        vp, cl, nx, thr, maxit, a, aint, g, ia, iy, iz, mm, nino, rsqc, nlp, jerr)
    names(result) <- c("almc", "m", "no", "ni", "r", "ju", "vp", "cl", "nx", "a", 
        "aint", "g", "ia", "iy", "iz", "mm", "nino", "rsqc", "nlp", "jerr")
    result
}


#' @export
wls_plink_cpp <- function(alm0, almc, alpha, m, no, ni, x, r, v, intr, ju, vp, cl, 
    nx, thr, maxit, a, aint, g, ia, iy, iz, mm, nino, rsqc, nlp, jerr) {
    result <- .Call("wls_plink", alm0, almc, alpha, m, no, ni, x, r, v, intr, ju, 
        vp, cl, nx, thr, maxit, a, aint, g, ia, iy, iz, mm, nino, rsqc, nlp, jerr)
    names(result) <- c("almc", "m", "no", "ni", "r", "ju", "vp", "cl", "nx", "a", 
        "aint", "g", "ia", "iy", "iz", "mm", "nino", "rsqc", "nlp", "jerr")
    result
}

setClass("PlinkMatrix", representation(ptr = "externalptr", samples = "integer", 
    variants = "integer", fname = "character", xm = "numeric", xs = "numeric", xmax = "numeric", 
    xmin = "numeric"), contains = "Matrix")
#' @export
PlinkMatrix <- function(fname, samples, variants, weight = NULL) {
    if (is.unsorted(samples, strictly = TRUE)) {
        stop("The subset indices must be sorted in strictly increasing order.")
    }
    
    if (length(variants) != length(unique(variants))) {
        stop("Must not have duplicated variant index.")
    }
    
    samples <- as.integer(samples)
    variants <- as.integer(variants)
    if (is.null(weight)) {
        weight <- rep(1/length(samples), length(samples))
    } else {
        if (length(samples) != length(weight)) {
            stop("weight must have the same length as the sample subset")
        }
        weight <- weight/(sum(weight))
    }
    
    # info = .Call('PlinkMatrix_info', fname, samples, variants, weight)
    ptr <- .Call("initialize_plinkmatrix_Xptr", fname, samples, variants)
    info <- .Call("PlinkMatrix_info", ptr, weight)
    
    new("PlinkMatrix", ptr = ptr, samples = samples, variants = variants, fname = fname, 
        Dim = c(length(samples), length(variants)), xm = info[[1]], xs = sqrt(info[[2]]), 
        xmax = info[[3]], xmin = info[[4]])
}



#' @export
xptrtest <- function(fname, samples, variants) {
    .Call("initialize_plinkmatrix_Xptr", fname, samples, variants)
}

setGeneric("colmax", function(object) {
    standardGeneric("colmax")
})

setGeneric("colmin", function(object) {
    standardGeneric("colmin")
})

setGeneric("center", function(object, ...) {
    standardGeneric("center")
})

setGeneric("standardize", function(object, ...) {
    standardGeneric("standardize")
})



#' @export
setMethod("colmax", signature(object = "matrix"), function(object) {
    apply(object, 2, max)
})

#' @export
setMethod("colmin", signature(object = "matrix"), function(object) {
    apply(object, 2, min)
})

#' @export
setMethod("colmax", signature(object = "PlinkMatrix"), function(object) {
    attr(object, "xmax", exact = TRUE)
})

#' @export
setMethod("colmin", signature(object = "PlinkMatrix"), function(object) {
    attr(object, "xmin", exact = TRUE)
})


#' @export
setMethod("center", signature(object = "matrix"), function(object, weights) {
    if (length(weights) != nrow(object)) {
        stop("Must provide a weight vector with same number of rows of the matrix")
    }
    xm = apply(x, 2, function(r) weighted.mean(r, weights))
    sweep(object, 2, xm, FUN = "-")
})

#' @export
setMethod("standardize", signature(object = "matrix"), function(object, weights) {
    if (length(weights) != nrow(object)) {
        stop("Must provide a weight vector with same number of rows of the matrix")
    }
    xs <- apply(x, 2, function(r) sqrt(weighted.mean(r^2, weights) - weighted.mean(r, 
        weights)^2))
    sweep(object, 2, xs, FUN = "/")
})


#' @export
setMethod("center", signature(object = "PlinkMatrix"), function(object) {
    .Call("PlinkSetMean", object@ptr, object@xm)
    return(object)
})

#' @export
setMethod("standardize", signature(object = "PlinkMatrix"), function(object) {
    .Call("PlinkSetSd", object@ptr, object@xs)
    return(object)
})

#' @export
setMethod("%*%", signature(x = "PlinkMatrix", y="numeric"), function(x,y) {
  if(length(y) != ncol(x))
  {
    stop('non-conformable arguments.')
  }
  result = rep(0.0, nrow(x))
  .Call('PlinkMultiplyv', x@ptr, y, result)
  return(result)
})