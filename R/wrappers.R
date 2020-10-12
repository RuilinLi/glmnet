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