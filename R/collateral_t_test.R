#' Test Collateral Effect
#'
#' Take MIC values of two antibiotics and quantify and test their collateral
#' effect. The quantified effect is the effect of B on A.
#'
#' @param A MIC values of antibiotic A, numeric vector
#' @param B MIC values of antibiotic B, numeric vector, same length as A
#' @param effect_type type of collateral effect to be evaluated, string:
#' `"both"`, `"CR"` or `"CS"`
#' @param crit_type type of dichotomization criterion, string: `"quant"`,
#' `"median"` or `"log2_MIC"`
#' @param criterium value of dichotomization criterium, numeric value (is
#' overwritten by `crit_type = "median"`)
#' @param warn set true to give warning when number of observations is too low
#' for the test.
#'
#' @return List with class `htest` including statistic, parameter, p.value, etc.
#'
#' @export
#'
#' @examples
#' A <- c(16, 16, 128, 128, 32, 16, 128, 1, 64, 1)
#' B <- c(32, 32, 16, 32, 8, 32, 8, 8, 8, 8)
#' collateral_t_test(A, B)
#'
collateral_t_test <- function(A, B, effect_type = "both", crit_type = "median",
                              criterium = NULL, warn = TRUE) {

  #Error handling
  if (length(A) != length(B)) {
    stop("Arguments A and B have different lengths: ", length(A), " and ",
         length(B), ".")
  }

  #Remove NA's
  ind_na <- is.na(A) | is.na(B)
  A <- log2(A[!ind_na])
  B <- log2(B[!ind_na])

  #Translate effect_type to t-test direction
  effect_types <- c("both", "CS", "CR")
  #Error handling
  if (!effect_type %in% effect_types) {
    warning("Argument effect_type is not a valid input, defaults to \"both\"")
    effect_type <- "both"
  }
  direction <- c("two.sided", "less", "greater")[effect_types == effect_type]

  #Calculate dichotomization criterium tau based on criterium type
  if (crit_type == "quant") {
    tau <- stats::quantile(B, criterium)

  } else if (crit_type == "median") {

    tau <- stats::quantile(B, 0.5)
    B_values <- sort(unique(B))
    #Adjust tau to create most equal split
    if (sum(B > tau) < sum(B < tau)) {
      tau <- mean(c(B_values[which(B_values == tau) - 1], tau))
    } else if (sum(B > tau) > sum(B < tau)) {
      tau <- mean(c(B_values[which(B_values == tau) + 1], tau))
    }

  } else {
    tau <- criterium
    if (!crit_type %in% c("quant", "median", "log2_MIC")) {
      warning("crit_type is not specified (correctly), defaults to ",
              "\"log2_MIC\"")
    }
  }

  A_Blow <- A[B < tau]
  A_Bhigh <- A[B >= tau]

  #Error handling
  if (length(A_Blow) < 2 | length(A_Bhigh) < 2) {
    if (warn) {
      warning("Not enough observations for test! Returning list with limited information. Try different dichotomization criterium for test.")
    }
    return(list(A = A, B = B, tau = tau, data = list(`A|B = r` = A_Bhigh,
                                                     `A|B != r` = A_Blow)))
  }

  #Perform t.test
  A_t_test <- stats::t.test(A_Bhigh, A_Blow, var.equal = TRUE,
                            alternative = direction)

  #Adjust names of groups
  names(A_t_test$estimate) <- c("mean of A|B = high", "mean of A|B = low")
  A_t_test$data.name <- "A|B = high and A|B = low"

  #Add data to output
  A_t_test$data <- list(`A|B = high` = A_Bhigh, `A|B = low` = A_Blow)

  #Add tau to object
  A_t_test$tau <- tau

  #Add log2(fold change) tot output
  A_t_test$log2_FC <- A_t_test$estimate[1] - A_t_test$estimate[2]
  names(A_t_test$log2_FC) <- "log2_FC"

  return(A_t_test)
}
