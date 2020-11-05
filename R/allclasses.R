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
    xmin = "numeric", weights = "numeric"), contains = "Matrix")
#' @export
PlinkMatrix <- function(fname, samples, variants, weights = NULL) {
    if (is.unsorted(samples, strictly = TRUE)) {
        stop("The subset indices must be sorted in strictly increasing order.")
    }
    
    if (length(variants) != length(unique(variants))) {
        stop("Must not have duplicated variant index.")
    }
    
    samples <- as.integer(samples)
    variants <- as.integer(variants)
    if (is.null(weights)) {
        weights <- rep(1/length(samples), length(samples))
    } else {
        if (length(samples) != length(weights)) {
            stop("weights must have the same length as the sample subset")
        }
        weights <- weights/(sum(weights))
    }
    
    # info = .Call('PlinkMatrix_info', fname, samples, variants, weight)
    ptr <- .Call("initialize_plinkmatrix_Xptr", fname, samples, variants)
    info <- .Call("PlinkMatrix_info", ptr, weights)
    
    new("PlinkMatrix", ptr = ptr, samples = samples, variants = variants, fname = fname, 
        Dim = c(length(samples), length(variants)), xm = info[[1]], xs = sqrt(info[[2]]), 
        xmax = info[[3]], xmin = info[[4]], weights = weights)
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

setGeneric("center", function(object, weights) {
    standardGeneric("center")
})

setGeneric("standardize", function(object, weights) {
    standardGeneric("standardize")
})

setGeneric("wlsFlex", function(x, ...) {
    standardGeneric("wlsFlex")
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
setMethod("colmax", signature(object = "Matrix"), function(object) {
    apply(object, 2, max)
})

#' @export
setMethod("colmin", signature(object = "Matrix"), function(object) {
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
    xm = apply(object, 2, function(r) weighted.mean(r, weights))
    list(x = sweep(object, 2, xm, FUN = "-"), xm = xm)
})

#' @export
setMethod("standardize", signature(object = "matrix"), function(object, weights) {
    if (length(weights) != nrow(object)) {
        stop("Must provide a weight vector with same number of rows of the matrix")
    }
    xs <- apply(object, 2, function(r) sqrt(weighted.mean(r^2, weights) - weighted.mean(r, 
        weights)^2))
    return(list(x = sweep(object, 2, xs, FUN = "/"), xs = xs))
})

#' @export
setMethod("center", signature(object = "sparseMatrix"), function(object, weights) {
    if (length(weights) != nrow(object)) {
        stop("Must provide a weight vector with same number of rows of the matrix")
    }
    xm = apply(object, 2, function(r) weighted.mean(r, weights))
    attr(object, "xm") <- xm
    list(x = object, xm = xm)
})

#' @export
setMethod("standardize", signature(object = "sparseMatrix"), function(object, weights) {
    if (length(weights) != nrow(object)) {
        stop("Must provide a weight vector with same number of rows of the matrix")
    }
    xs <- apply(object, 2, function(r) sqrt(weighted.mean(r^2, weights) - weighted.mean(r, 
        weights)^2))
    attr(object, "xs") <- xs
    return(list(x = object, xs = xs))
})


#' @export
setMethod("center", signature(object = "PlinkMatrix"), function(object, weights) {
    .Call("PlinkSetMean", object@ptr, object@xm)
    return(list(x = object, xm = object@xm))
})

#' @export
setMethod("standardize", signature(object = "PlinkMatrix"), function(object, weights) {
    .Call("PlinkSetSd", object@ptr, object@xs)
    return(list(x = object, xs = object@xs))
})

#' @export
setMethod("%*%", signature(x = "PlinkMatrix", y = "numeric"), function(x, y) {
    if (length(y) != ncol(x)) {
        stop("non-conformable arguments.")
    }
    result = rep(0, nrow(x))
    .Call("PlinkMultiplyv", x@ptr, y, result)
    return(result)
})

#' @export
setMethod("%*%", signature(x = "numeric", y = "PlinkMatrix"), function(x, y) {
    if (length(x) != nrow(y)) {
        stop("non-conformable arguments.")
    }
    result = numeric(ncol(y))
    .Call("PlinkPreMultiplyv", y@ptr, x, result)
    return(result)
})

#' @export
setMethod("wlsFlex", signature(x = "sparseMatrix"), function(x, alm0, almc, alpha, 
    m, nobs, nvars, r, v, intr, ju, vp, cl, nx, thr, maxit, a, aint, g, ia, iy, iz, 
    mm, nino, rsqc, nlp, jerr) {
    xm <- as.double(attr(x, "xm"))
    xs <- as.double(attr(x, "xs"))
    ix <- as.integer(x@p + 1)
    jx <- as.integer(x@i + 1)
    x <- as.double(x@x)
    wls_fit <- .Fortran("spwls", alm0 = alm0, almc = almc, alpha = alpha, m = m, 
        no = nobs, ni = nvars, x = x, ix = ix, jx = jx, xm = xm, xs = xs, r = r, 
        v = v, intr = intr, ju = ju, vp = vp, cl = cl, nx = nx, thr = thr, maxit = maxit, 
        a = a, aint = aint, g = g, ia = ia, iy = iy, iz = iz, mm = mm, nino = nino, 
        rsqc = rsqc, nlp = nlp, jerr = jerr)
    return(wls_fit)
})

#' @export
setMethod("wlsFlex", signature(x = "matrix"), function(x, alm0, almc, alpha, m, nobs, 
    nvars, r, v, intr, ju, vp, cl, nx, thr, maxit, a, aint, g, ia, iy, iz, mm, nino, 
    rsqc, nlp, jerr) {
    
    wls_fit <- wls_dense_cpp(alm0 = alm0, almc = almc, alpha = alpha, m = m, no = nobs, 
        ni = nvars, x = x, r = r, v = v, intr = intr, ju = ju, vp = vp, cl = cl, 
        nx = nx, thr = thr, maxit = maxit, a = a, aint = aint, g = g, ia = ia, iy = iy, 
        iz = iz, mm = mm, nino = nino, rsqc = rsqc, nlp = nlp, jerr = jerr)
    return(wls_fit)
})

#' @export
setMethod("wlsFlex", signature(x = "PlinkMatrix"), function(x, alm0, almc, alpha, 
    m, nobs, nvars, r, v, intr, ju, vp, cl, nx, thr, maxit, a, aint, g, ia, iy, iz, 
    mm, nino, rsqc, nlp, jerr) {
    wls_fit <- wls_plink_cpp(alm0 = alm0, almc = almc, alpha = alpha, m = m, no = nobs, 
        ni = nvars, x = x@ptr, r = r, v = v, intr = intr, ju = ju, vp = vp, cl = cl, 
        nx = nx, thr = thr, maxit = maxit, a = a, aint = aint, g = g, ia = ia, iy = iy, 
        iz = iz, mm = mm, nino = nino, rsqc = rsqc, nlp = nlp, jerr = jerr)
    return(wls_fit)
})

