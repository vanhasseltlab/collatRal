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
  not_significant <- with(t_result, ifelse(is.null(q), p > sign_criterium,
                                                       q > sign_criterium))
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
  levels(t_result$reciprocal) <- c("one-directional", "reciprocal")

  #Use only antibiotics from selected_ab argument
  if (!is.null(selected_ab)) {
    only_select <- with(t_result, A %in% selected_ab & B %in% selected_ab)
    t_result <- t_result[only_select, ]
  }

  #define plot colors, limits and labels
  blues <- c("#283C82", "#2F4286", "#36498A", "#3E508E", "#455693", "#4D5D97",
             "#54649B", "#5B6BA0", "#6371A4", "#6A78A8", "#727FAD", "#8893BA",
             "#9EA7C6", "#B4BBD3", "#CBCFE0", "#E1E4ED", "#F7F8FA")
  oranges <- c("#F54C00", "#F55208", "#F55811", "#F65E1A", "#F66423", "#F66A2B",
               "#F77134", "#F7773D", "#F77D46", "#F8834F", "#F88957", "#F99C72",
               "#FAAE8C", "#FBC1A7", "#FCD3C1", "#FDE6DB", "#FEF8F6")
  limits <- c(-1, 1) * max(t_result$effect_size)
  draw_grid <- function(x) {
    seq(1.5, length(unique(x)) - 0.5, 1)
  }

  #define variables in environment to circumvent building errors
  A <- B <- effect_size <- reciprocal <- NULL

  #make ggplot
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
    ggplot2::coord_fixed() +
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
                   axis.text.x = ggplot2::element_text(angle = 90))

  return(plot_)
}

#' Plot MIC Distribution
#'
#' WORK IN PROGRESS
#'
#' @param MIC_clean MIC data frame (rows = strains, cols = AB)
#' @param results T test results from wrapper function (currently still needed)
#' @param t_rank in the `results` outcome, which ranking value should be plotted (overwritten by whichAB)
#' @param one_direction plot only A on B  (TRUE) or also plot B on A (FALSE)
#' @param CResponse direction of response either collateral sensitivity (CS) or resistance (CR)
#' @param sign_criterium the criterium for which of results are considered significant
#' @param whichAB pair of antibiotics (names string) you want to plot
#' @param colours pair or colors for B = high and B = low groups
#' @param separate_plots if one_direction = TRUE, whether you want to return two separate plots (in list) or one plot (grid.arrange)
#' @param ran range of MIC values plotted
#' @param t_or_effect Want to plot the highest ranking effect or T value?
#'
#' @return (list with) ggplot object(s) or grid.arrange output
#' @export
#'
#' @examples
#' MIC_test <- data.frame(d1 = c(16, 16, 128, 128, 32, 16, 128, 1, 64, 1),
#'                        d2 = c(32, 32, 16, 32, 8, 32, 8, 8, 8, 8),
#'                        d3 = c(8, 64, 64, 64, 16, 6, 16, 8, 8, 8))
#' MIC_test_result <- collateral_mult_test(MIC_test)
#' plot_histogram_CE(MIC_test, MIC_test_result)
#'
#'
#TO-DO: reduce to base and ggplot, remove most arguments, only include an AB pair
plot_histogram_CE <- function(MIC_clean, results, t_rank = 1,
                              one_direction = TRUE, CResponse = "CS", sign_criterium = 0.15, whichAB = NULL,
                              colours = c("#283c82", "#F46B2D"), separate_plots = NULL, ran = NULL,
                              t_or_effect = "t") {
  return(NULL)
#
#   if (!is.null(whichAB)) {
#     comb <- whichAB
#     sign_dat <- results %>%
#       filter(A == comb[1] & B == comb[2])
#     row.names(sign_dat) <- paste0(sign_dat$A, sign_dat$B)
#     if (whichAB[1] %in% sign_dat$A & whichAB[2] %in% sign_dat$B) {
#       d <- sign_dat[paste0(comb, collapse = ""), "tau"]
#     } else {
#       message("Specified antibiotics are not in data, using highest effect size")
#       whichAB <- NULL
#     }
#   }
#
#
#   if (is.null(whichAB)) {
#     if (t_or_effect == "effect"){
#       results$t <- results$effect_size
#     }
#
#     sign_dat <- results %>%
#       filter(p_BY < sign_criterium & (sign(t) == c(-1, 1)[CResponse == c("CS", "CR")])) %>%
#       arrange(desc(abs(t)))
#     row.names(sign_dat) <- paste0(sign_dat$A, sign_dat$B)
#
#     if (nrow(sign_dat) < 1) {
#       message("No significant ", CResponse, " effect, plotting non significant finding with largest T")
#       if (CResponse == "CS") {
#         sign_dat <- results %>%
#           arrange(t) %>%
#           slice(1)
#       } else {
#         sign_dat <- results %>%
#           arrange(desc(t)) %>%
#           slice(1)
#       }
#       t_rank <- 1
#     }
#
#     if (nrow(sign_dat) < t_rank) {
#       message(paste0("No", t_rank, "th significant ", CResponse, " effect, plotting significant finding with largest T"))
#       t_rank <- 1
#     }
#     comb <- unlist(sign_dat[t_rank, 1:2])
#     d <- sign_dat[t_rank, "tau"]
#   }
#
#
#
#   #for all combinations?
#   dat <- log2(MIC_clean[, comb])
#   dat <- dat[!is.na(dat[, 1]) & !is.na(dat[, 2]), ]
#   names(dat) <- c("A", "B")
#
#   dat$Condition <- as.factor(ifelse(dat$B >= d,
#                                     paste0(comb[1],"|", comb[2]," > ", round(d, 2)),
#                                     paste0(comb[1],"|", comb[2]," < ", round(d, 2))))
#   means <- data.frame(mean = c(mean(dat$A[dat$B < d]), mean(dat$A[dat$B >= d])),
#                       Means = sort((unique(dat$Condition))))
#   change_range <- FALSE
#   if (is.null(ran)){
#     ran <- range(dat$A) + (0.5 * c(-1, + 1))
#     change_range <- TRUE
#   }
#   ticks <- floor(seq(ran[1], ran[2] + 2))
#
#   plotPanel <- function(dat) {
#
#     plotje <- dat %>%
#       ggplot(aes(x = A, y = ..count.., group = Condition, fill = Condition)) +
#       geom_bar(stat = "count", width = 0.6, position = "stack") +
#       labs(x = bquote(log[2]*"(MIC) " ~ .(comb[1])), y = "Counts") +
#       #labs(x = expression(log[2]*"(MIC) "*comb[1]), y = "Counts") +
#       #scale_y_continuous(expand = expansion(mult = c(0, .05)), labels = function(x) sprintf("%.2f", x)) +
#       scale_y_continuous(expand = expansion(mult = c(0, .05))) +
#       scale_x_continuous(breaks = ticks, limits = ran) +
#       geom_vline(data = means, aes(xintercept = mean), colour = "white") +
#       geom_vline(data = means, aes(xintercept = mean, colour = Means), show.legend  = TRUE, linetype = 2) +
#       scale_fill_manual(values = colours) +
#       scale_colour_manual(labels = c(bquote(hat(mu)[.(as.character(means$Means[1]))]),
#                                      #                               bquote(paste(hat(mu)[paste(comb[1], "|", comb[2],>=, d)]))),
#                                      bquote(hat(mu)[.(as.character(means$Means[2]))])),
#                           values = colours) +
#       theme_bw() +
#       theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(), legend.position = "bottom",
#             #     legend.background = element_rect(fill = "transparent", size = 0.5, linetype = "solid", colour = "black"),
#             legend.box = "horizontal", legend.direction = "vertical",
#             legend.margin = margin(0,0,0,0),
#             legend.box.margin = margin(-15,0,0,0),
#             legend.background = element_blank()) +
#       guides(fill = guide_legend(override.aes = list(linetype = 0), order = 1, title = NULL),
#              colour = guide_legend(order = 2, title = NULL, label.theme = element_text(size = 12)))
#
#
#     return(plotje)
#   }
#
#
#   p_A <- plotPanel(dat)
#
#   if (one_direction) {
#     return(p_A)
#   }
#
#   new_dat <- results
#   rownames(new_dat) <- paste(new_dat$A, new_dat$B, sep = "_")
#   comb <- comb[2:1]
#   d <- new_dat[paste(comb, collapse = "_"), "tau"]
#   dat <- log2(MIC_clean[, comb])
#   dat <- dat[!is.na(dat[, 1]) & !is.na(dat[, 2]), ]
#   names(dat) <- c("A", "B")
#   d_min <- max(dat$B[dat$B < d], na.rm = T)
#
#   dat$Condition <- as.factor(ifelse(dat$B >= d,
#                                     paste0(comb[1],"|", comb[2]," > ", round(d, 2)),
#                                     paste0(comb[1],"|", comb[2]," < ", round(d, 2))))
#   means <- data.frame(mean = c(mean(dat$A[dat$B < d]), mean(dat$A[dat$B >= d])),
#                       Means = sort((unique(dat$Condition))))
#
#   if (change_range == TRUE){
#     ran <- range(dat$A) + (0.5 * c(-1, + 1))
#   }
#   ticks <- floor(seq(ran[1], ran[2] + 2))
#   p_B <- plotPanel(dat)
#
#
#
#   if (is.null(separate_plots)) {
#     p_A <- p_A + theme(axis.title.y = element_blank())
#     p_B <- p_B + theme(axis.title.y = element_blank())
#     return(grid.arrange(p_A, p_B, ncol = 2, left = textGrob("Counts", hjust = -0.45, rot = 90)))
#   } else {
#     return(list(p_A, p_B + ylab("  ")))
#   }
}
