#' Plot Heatmap of Collateral Effects
#'
#' Plot heatmap of T-values or effect sizes for results from
#' `collateral_mult_test`.
#'
#' @param t_result T test results from wrapper function (to be developed)
#' @param sign_criterium the criterium for which of results are considered
#' significant, numeric value in (0, 1].
#' @param selected_ab optional argument, character vector containing which antibiotics should be included
#' @param t_or_effect plot the effect size or the T value, character `"t"` or `"effect"`
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' MIC_test <- data.frame(d1 = c(16, 16, 128, 128, 32, 16, 128, 1, 64, 1),
#'                        d2 = c(32, 32, 16, 32, 8, 32, 8, 8, 8, 8),
#'                        d3 = c(8, 64, 64, 64, 16, 6, 16, 8, 8, 8))
#' MIC_test_result <- collateral_mult_test(MIC_test)
#' plot_heatmap_CE(MIC_test_result)
#'
#TO-DO: reduce to base and ggplot, remove some arguments, fit to wrapper function output
plot_heatmap_CE <- function(t_result, sign_criterium = 1, selected_ab = NULL,
                            t_or_effect = "effect") {

  t_result[is.na(t_result$t), c("effect_size", "t")] <- 0

  if (t_or_effect == "t") {
    t_result$effect_size <- t_result$t
    effect_label <- "T-value"
  } else {
    effect_label <- expression(log[2]*"(FC)")
  }

  #Set all non significant finding to zero (white in plot)
  if (is.null(t_result$q)) {
    not_significant <- t_result$p > sign_criterium
  } else {
    not_significant <- t_result$q > sign_criterium
  }
  t_result$effect_size[not_significant] <- 0


  # find reciprocal AB combinations
  t_result$reciprocal <- "one-directional"
  combinations <- t(utils::combn(unique(t_result$A), 2))
  for (i in 1:nrow(combinations)) {
    comb <- combinations[i, ]
    ind_comb <- t_result$A %in% comb & t_result$B %in% comb
    df_comb <- t_result[ind_comb, ]
    if (all(df_comb$effect_size > 0) | all(df_comb$effect_size < 0)){
      t_result[ind_comb, "reciprocal"] <- "reciprocal"
    }
  }
  t_result$reciprocal <- factor(t_result$reciprocal,
                                levels = c("one-directional", "reciprocal"))

  #Use only antibiotics from selected_ab argument
  if (!is.null(selected_ab)) {
    only_select <- with(t_result, A %in% selected_ab & B %in% selected_ab)
    t_result <- t_result[only_select, ]
  }

  #define plot colors, limits and labels
  blues <- c("#182450", "#243676", "#283C82", "#2F4286", "#36498A", "#3E508E",
             "#455693", "#4D5D97", "#54649B", "#5B6BA0", "#6371A4", "#6A78A8",
             "#727FAD", "#8893BA", "#9EA7C6", "#CBCFE0")
  oranges <- c("#F54C00", "#F55208", "#F55811", "#F65E1A", "#F66423", "#F66A2B",
               "#F77134", "#F7773D", "#F77D46", "#F8834F", "#F88957", "#F99C72",
               "#FAAE8C", "#FBC1A7", "#FCD3C1", "#FDE6DB")
  limits <- c(-1, 1) * max(abs(t_result$effect_size))
  draw_grid <- function(x) {
    seq(1.5, length(unique(x)) - 0.5, 1)
  }

  # Define variables in environment to circumvent building errors
  A <- B <- effect_size <- reciprocal <- NULL

  # Make plot
  plot_ <- ggplot2::ggplot(t_result, ggplot2::aes(x = B, y = A)) +
    # Add plot layers
    ggplot2::geom_tile(ggplot2::aes(fill = effect_size)) +
    ggplot2::geom_point(ggplot2::aes(shape = reciprocal), size = 5,
                        colour = "white") +
    # Scale aesthetics
    ggplot2::scale_x_discrete(expand = ggplot2::expansion(mult = 0,
                                                          add = rep(0.5, 2))) +
    ggplot2::scale_y_discrete(expand = ggplot2::expansion(mult = 0,
                                                          add = rep(0.5, 2))) +
    ggplot2::scale_fill_gradientn(colours = c(blues, "white", rev(oranges)),
                                  limits = limits) +
    ggplot2::scale_shape_manual(values = c("", "\u2194"), drop = FALSE) +
    ggplot2::labs(x = "Splitting antibiotic (B)",
                  y = paste0("Testing antibiotic (A)"),
                  fill = effect_label, shape = "") +
    # Create custom grid
    ggplot2::geom_vline(xintercept = draw_grid(t_result$A), colour = "grey60") +
    ggplot2::geom_hline(yintercept = draw_grid(t_result$B), colour = "grey60") +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "grey60") +
    # Change theme
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   legend.title = ggplot2::element_text(size = 12),
                   legend.key = ggplot2::element_rect(fill = "grey60"),
                   axis.title = ggplot2::element_text(size = 14),
                   axis.text.x = ggplot2::element_text(angle = 90, hjust = 1,
                                                       vjust = 0.5),
                   aspect.ratio = 1) #+
    # ggplot2::facet_grid(class_A ~ class_B, drop = TRUE, scales = "free")

  return(plot_)
}

#' Plot MIC Distribution
#'
#' WORK IN PROGRESS
#'
#' @param MIC_range range of x-axis (default NULL)
#' @param colors set of two colors for the high and low MIC groups
#' @inheritParams collateral_t_test
#'
#' @return (list with) ggplot object(s) or grid.arrange output
#' @export
#'
#' @examples
#' MIC_test <- data.frame(d1 = c(16, 16, 128, 128, 32, 16, 128, 1, 64, 1),
#'                        d2 = c(32, 32, 16, 32, 8, 32, 8, 8, 8, 8),
#'                        d3 = c(8, 64, 64, 64, 16, 6, 16, 8, 8, 8))
#' MIC_test_result <- collateral_mult_test(MIC_test)
#' plot_histogram_CE(MIC_test$d1, MIC_test$d2)
#'
#'
plot_histogram_CE <- function(A, B, effect_type = "both", crit_type = "median",
                              criterium = NULL, MIC_range = NULL,
                              colors = c("#BFC6B8", "#4A5242")) {
  # Perform t-test to determine tau and the means
  t_result <- collateral_t_test(A, B, effect_type = effect_type,
                                crit_type = crit_type, criterium = criterium,
                                warn = FALSE)
  d <- t_result$tau

  antibiotics <- c(deparse(substitute(A)), deparse(substitute(B)))
  # Remove data.frame name if input includes '$' sign.
  antibiotics <- ifelse(grepl("$", antibiotics, fixed = T),
                        gsub(".*\\$","",antibiotics), antibiotics)

  ind_na <- is.na(A) | is.na(B)
  A <- A[!ind_na]
  B <- B[!ind_na]

  plot_data <- data.frame(A = log2(A), B = log2(B))


  plot_data$Condition <- as.factor(ifelse(plot_data$B >= d,
    paste0(antibiotics[1],"|", antibiotics[2],"	 \u2265 ", round(d, 2)),
    paste0(antibiotics[1],"|", antibiotics[2]," < ", round(d, 2))))

  means <- data.frame(Means = sort(unique(plot_data$Condition)),
                      mean = t_result$estimate[2:1])
  # Define variables in environment to circumvent building errors
  Means <- Condition <- ..count.. <-  NULL

  if (is.null(MIC_range)) {
    MIC_range <- range(plot_data$A) + (0.5 * c(-1, + 1))
  }

  ticks <- floor(seq(MIC_range[1], MIC_range[2] + 2))
  labeling <- c(bquote(hat(mu)[.(as.character(means$Means[1]))]),
                bquote(hat(mu)[.(antibiotics[1])*"|"*
                                 .(paste0(antibiotics[2], " ")) >=
                                 .(paste0(" ", d))]))

  histogram <- ggplot2::ggplot(plot_data,
                               ggplot2::aes(x = A, y = ..count..,
                                            group = Condition,
                                            fill = Condition)) +
    # Add plot layers
    ggplot2::geom_bar(stat = "count", width = 0.6, position = "stack") +
    ggplot2::geom_vline(data = means, ggplot2::aes(xintercept = mean),
                        color = "white") +
    ggplot2::geom_vline(data = means, linetype = 2,
                        ggplot2::aes(xintercept = mean, color = Means)) +
    ggplot2::labs(x = bquote(log[2]*"(MIC) " ~ .(antibiotics[1])),
                  y = "Counts") +
    # Scale aesthetics
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, .05))) +
    ggplot2::scale_x_continuous(breaks = ticks, limits = MIC_range) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::scale_colour_manual(labels = labeling, values = colors) +

    # Change theme
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.minor.x = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.box = "horizontal", legend.direction = "vertical",
      legend.margin = ggplot2::margin(0,0,0,0),
      legend.box.margin = ggplot2::margin(-15,0,0,0),
      legend.background = ggplot2::element_blank()
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(override.aes = list(linetype = 0), order = 1,
                                   title = NULL),
      color = ggplot2::guide_legend(order = 2, title = NULL,
                                     label.theme = ggplot2::element_text(size = 12))
      )

  return(histogram)
}
