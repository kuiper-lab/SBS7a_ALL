library(cellPhyWrapperPlotting)
check_0_mm <- cellPhyWrapperPlotting:::check_0_mm

fit_to_signatures_strict_tree <- function (mut_matrix, signatures, max_delta = 0.01, remove_min = 20, method){
  mut_matrix = check_0_mm(mut_matrix, remove_min)
  contri = fit_to_signatures_strict(mut_matrix, signatures,
                                    max_delta = max_delta,
                                    method = method)$fit_res$contribution
  contri = contri[rowSums(contri) > 0, ]
}
