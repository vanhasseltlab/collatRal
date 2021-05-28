#' Translate antibiotic names to abbreviations
#'
#' @param antibiotic_names character vector containing the names to be
#' translated.
#' @param dictionary matrix or data frame with at least two columns containing
#' the names to be translated and the names to which they will be translated.
#' @param which_dictionary_cols character vector with two column names of the
#' from and to translation values, respectively.
#'
#' @return character vector with same length as `antibiotic_names` containing
#' translated names
#' @export
#'
#' @examples
#' dict <- data.frame(antibiotic = c("amikacin", "tobramycin", "colistin"),
#'                    abbreviation = c("AMK", "TOB", "CST"))
#' ab_names <- c("amikacin", "amikacin", "tobramycin")
#' translate_antibiotic(ab_names, dict)
translate_antibiotic <- function(antibiotic_names, dictionary,
                                 which_dictionary_cols = c("antibiotic",
                                                           "abbreviation")) {
  # Prepare dictionary
  dictionary$flat <- make.names(dictionary[, which_dictionary_cols[1]])
  dictionary <- dictionary[!duplicated(dictionary$flat), ]
  rownames(dictionary) <- dictionary$flat

  # Prepare antibiotic names
  ab_names_new <- as.factor(antibiotic_names)
  ab_names_unique <- unique(make.names(levels(ab_names_new)))
  ab_names_unique <- gsub("_", ".", ab_names_unique, fixed = TRUE)
  ab_abbr <- dictionary[ab_names_unique, which_dictionary_cols[2]]

  antibiotic_names <- gsub("_", ".", make.names(antibiotic_names), fixed = TRUE)

  ab_abbr <- dictionary[antibiotic_names, which_dictionary_cols[2]]
  # Warning for antibiotic_names not in the dictionary
  not_available <- which(!antibiotic_names %in% dictionary$flat)
  if (length(not_available) > 0) {
    warning("Not in dictionary: ", paste(unique(antibiotic_names[not_available]),
                                         collapse = ", "))
    ab_abbr[not_available] <- antibiotic_names[not_available]
  }



  return(ab_abbr)
}
