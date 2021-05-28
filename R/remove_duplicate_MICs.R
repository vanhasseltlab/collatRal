#' Remove duplicate MIC measurements from MIC data
#'
#' If multiple MIC measurements are done on the same strain and antibiotic, but
#' this function will summarize MIC value based on `measurement_sign`
#'
#' @param MIC_df a long format data frame with each row an MIC value (`MIC`) for
#' a strain/antibiotic. Contains a `measurement_sign` column.
#' and
#' @param key character (vector) containing column names with the columns over
#' which the duplicates should be removed.
#'
#' @return data frame with same structure as `MIC_df`, with unique rows based on
#' the provided key
#' @export
#'
#' @examples
#' #not yet created
remove_duplicate_MICs <- function(MIC_df, key = NULL) {

  if (!is.null(key)) {
    MIC_df$key <- apply(MIC_df[, key], 1, paste, collapse = "_")
  } else {
    #error handling, if no key available
    if (is.null(MIC_df$key)) {
      stop("No key provided for which measurements are different")
    }
  }

  # two unique measurements for the same key will be summarized:
  #   max if ">" and mean if "<="
  if (length(unique(MIC_df$key)) == nrow(MIC_df)) {
    return(MIC_df)
  }
  find_unique <- ave(MIC_df$key, MIC_df$key, FUN = length) == 1

  # add unique MICs to final data frame
  new_MIC <- MIC_df[find_unique, ]

  search_keys <- unique(MIC_df$key[!find_unique])
  # add rows for keys with duplicates

  MIC_reduced <- as.data.frame(matrix(NA, nrow = length(search_keys)),
                               row.names = search_keys)
  for (i in search_keys) {
    cat("\r", "At (%): ", round(which(search_keys == i)/length(search_keys)*100, 2), ", key = ", i)
    dat <- MIC_df[MIC_df$key == i, ]
    if (nrow(dat) == 1) {
      MIC_reduced[i, ] <- dat
    }
    if (nrow(dat) > 1) {
      if (dat$measurement_sign[1] == dat$measurement_sign[2]) {
        if(dat$measurement_sign[1] == "<=") {
          MIC_i <- mean(dat$MIC)
        }
        if(dat$measurement_sign[1] == ">") {
          MIC_i <- max(dat$MIC)
        }
      } else {
        MIC_i <- mean(dat$MIC[dat$measurement_sign == "<="])
        dat$measurement_sign <- "<="
      }
      dat <- dat[1, ]
      dat$MIC <- MIC_i
      new_MIC <- rbind(new_MIC, dat)
    }
  }
  return(new_MIC)
}
