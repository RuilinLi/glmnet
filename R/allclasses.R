#' @export
wls_dense_cpp <- function(alm0, almc, alpha, m, no, ni,
                            x, r, v, intr, ju, vp, cl, nx, thr,
                            maxit, a, aint, g, ia, iy, iz, mm,
                            nino, rsqc, nlp, jerr) {
  result = .Call('wls_dense', alm0, almc, alpha, m, no, ni,
      x, r, v, intr, ju, vp, cl, nx, thr,
      maxit, a, aint, g, ia, iy, iz, mm,
      nino, rsqc, nlp, jerr)
  names(result) = c("almc", "m", "no", "ni", "r", "ju", "vp", "cl", "nx",
                          "a", "aint", "g", "ia", "iy", "iz", "mm", "nino",
                          "rsqc", "nlp", "jerr")
 result
}


#' @export
wls_plink_cpp <- function(alm0, almc, alpha, m, no, ni,
                            x, r, v, intr, ju, vp, cl, nx, thr,
                            maxit, a, aint, g, ia, iy, iz, mm,
                            nino, rsqc, nlp, jerr) {
  result = .Call('wls_plink', alm0, almc, alpha, m, no, ni,
      x, r, v, intr, ju, vp, cl, nx, thr,
      maxit, a, aint, g, ia, iy, iz, mm,
      nino, rsqc, nlp, jerr)
  names(result) = c("almc", "m", "no", "ni", "r", "ju", "vp", "cl", "nx",
                          "a", "aint", "g", "ia", "iy", "iz", "mm", "nino",
                          "rsqc", "nlp", "jerr")
 result
}

setClass("PlinkMatrix", representation(samples = "integer", variants="integer", 
                                        fname="character", xm="numeric", xs="numeric", xmax="numeric", xmin="numeric"),
         contains = "Matrix")
#' @export
PlinkMatrix <- function(fname, samples, variants, weight=NULL)
{
  samples =as.integer(sort(unique(samples)))
  if(is.null(weight))
  {
      weight = rep(1.0/length(samples), length(samples))
  }
  info = .Call('PlinkMatrix_info', fname, samples, variants, weight)

  new("PlinkMatrix", samples=samples, variants = variants, 
    fname=fname, Dim=c(length(samples),length(variants)), xm=info[[1]], xs=sqrt(info[[2]]), xmax=info[[3]], xmin=info[[4]])
}



#' @export
xptrtest = function(fname, samples, variants)
{
  .Call('initialize_plinkmatrix_Xptr', fname, samples, variants)
}

setGeneric("colmax", function(object) {
  standardGeneric("colmax")
})

setGeneric("colmin", function(object) {
  standardGeneric("colmin")
})

setMethod("colmax", signature(object = "matrix"), function(object) {
  apply(object, 2, max)
})

setMethod("colmin", signature(object = "matrix"), function(object) {
  apply(object, 2, mind)
})