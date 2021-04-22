#' Test Multiple Collateral Effects
#'
#' Wrapper for `collateral_t_test()`, which quantifies and tests collateral
#' effects in a set of MICs for more than two antibiotics.
#'
#' @param MIC_data data frame or matrix with MIC for strains (rows) and
#' antibiotics (cols, with corresponding colnames)
#' @param MTC_method multiple testing correction method (see
#' `?stats::p.adjust.methods`)
#' @inheritParams collateral_t_test
#'
#' @return Data frame with results for every antibiotic combination, containing
#' the T-value, p-value and corrected p-value (q).
#' @export
#'
#' @examples
#'
#' MIC_test <- data.frame(d1 = c(16, 16, 128, 128, 32, 16, 128, 1, 64, 1),
#'                        d2 = c(32, 32, 16, 32, 8, 32, 8, 8, 8, 8),
#'                        d3 = c(8, 64, 64, 64, 16, 6, 16, 8, 8, 8))
#' collateral_mult_test(MIC_test)

collateral_mult_test <- function(MIC_data, effect_type = "both",
                                 crit_type = "median", criterium = NULL,
                                 MTC_method = "BY") {
  antibiotics <- colnames(MIC_data)
  m <- length(antibiotics)

  t_test_results <- as.data.frame(matrix(NA, nrow = m*(m - 1), ncol = 9))
  names_all <- c("A", "B", "n", "n_high", "tau")
  names_test_successful <- c("mean_high", "mean_low", "t", "p")
  names(t_test_results) <- c(names_all, names_test_successful)

  counter <- 0
  for (dep in 1:m) {
    for (indep in 1:m) {
      if (dep == indep) {
        next
      }

      counter <- counter + 1
      ab <- antibiotics[c(dep, indep)]

      t_test <- collateral_t_test(A = MIC_data[, dep], B = MIC_data[, indep],
                                  effect_type = "both", crit_type = crit_type,
                                  criterium = criterium, warn = FALSE)
      result_i <- data.frame(A = ab[1], B = ab[2],
                             n = sum(lengths(t_test$data)),
                             n_high = length(t_test$data$`A|B = high`),
                             tau = t_test$tau, stringsAsFactors = F)
      result_i[, names_test_successful] <- NA

      #Add results t test if no problems occurred
      if (!is.null(t_test$statistic)) {
        result_i[, names_test_successful] <- data.frame(
          mean_high = t_test$estimate[1],
          mean_low = t_test$estimate[2], t = t_test$statistic,
          p = t_test$p.value, stringsAsFactors = F)
      }
      t_test_results[counter, ] <- result_i
    }
  }

  M <- c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
  if (MTC_method %in% M) {
    t_test_results$q <- stats::p.adjust(t_test_results$p, method = MTC_method)
  }

  t_test_results$effect_size <- with(t_test_results, mean_high - mean_low)
  t_test_results$effect_type <- with(t_test_results, ifelse(t > 0, "CR", "CS"))

  return(t_test_results)
}
