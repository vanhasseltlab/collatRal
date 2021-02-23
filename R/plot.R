#' Plot Heatmap of Collateral Effects
#'
#' WORK IN PROGRESS
#'
#' @param t_result T test results from wrapper function (to be developed)
#' @param FDR_crit the criterium for which of results are considered significant
#' @param species name of species can be added to put in title of plot
#' @param selectedAB optional argument, character vector containing which antibiotics should be included
#' @param t_or_effect plot the effect size or the T value, character `"t"` or `"effect"`
#'
#' @return
#' @export
#'
#' @examples
#'
#'
#TO-DO: reduce to base and ggplot, remove some arguments, fit to wrapper function output
plot_heatmap_CE <- function(t_result, FDR_crit = 0.15, species,
  selectedAB = NULL, t_or_effect = "effect") {

  t_result[is.na(t_result$t), c("effect_size", "t")] <- 0

  if (t_or_effect == "t") {
    t_result$effect_size <- t_result$t
    effect_label <- "T-value"
  } else {
    effect_label <- expression(log[2]*"(FC)")
  }


  t_result$effect_size[t_result$p_BY > FDR_crit] <- 0


  # find reciprocal AB combinations
  t_result$reciprocal <- "one-directional"
  combinations <- t(combn(unique(t_result$A), 2))

  for (i in 1:nrow(combinations)) {
    comb <- combinations[i, ]
    ind_comb <- t_result$A %in% comb & t_result$B %in% comb
    df_comb <- t_result[ind_comb, ]
    if (all(df_comb$effect_size > 0) | all(df_comb$effect_size < 0)){
      t_result[ind_comb, "reciprocal"] <- "reciprocal"
    }
  }
  levels(t_result$reciprocal) <- c("one-directional", "reciprocal")

  if (!is.null(selectedAB)) {
    t_result <- t_result %>%
      filter(A %in% selectedAB & B %in% selectedAB)
  }
  bl <- colorRampPalette(c("#283c82", "white"))(30) [c(1:10, seq(11, 30, by = 3))]
  re <- colorRampPalette(c("#f54c00", "white"))(30)[c(1:10, seq(11, 30, by = 3))]

  limits <- c(-1, 1)*max(t_result$effect_size)

  ll <- c("Collateral Sensitivity", "Collateral Resistance") # labels.

  t_result$Direction <- ifelse(t_result$effect_size < 0, ll[1], ll[2] )


  plot_ <- ggplot(t_result, aes(x = B, y = A)) +
    geom_tile(aes(fill = effect_size)) +
    scale_fill_gradientn(colours = c(bl, "white", rev(re)), limits = limits) +
    scale_x_discrete(expand = expansion(mult = 0, add = rep(0.5, 2))) +
    scale_y_discrete(expand = expansion(mult = 0, add = rep(0.5, 2))) +
    geom_point(aes(shape = reciprocal), size = 5, colour = "white") +
    scale_shape_manual(values = c("", "\u2194"), drop = FALSE) +
    coord_fixed() +
    theme(axis.text.x = element_text(angle = 90)) +
    labs(x = "Splitting antibiotic (B)", y = paste0("Testing antibiotic (A)"),
         fill = effect_label, shape = "",
         title = bquote("Significant collateral responses"~italic(.(species)))) +
    geom_vline(xintercept = seq(1.5, length(unique(t_result$A)) - 0.5, 1), colour = "grey60") +
    geom_hline(yintercept = seq(1.5, length(unique(t_result$B)) - 0.5, 1), colour = "grey60") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.title = element_text(size = 12), legend.key = element_rect(fill = "grey60"),
          axis.title = element_text(size = 14)) +
    geom_abline(slope = 1, intercept = 0, colour = "grey60")

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
#' @param FDR_crit the criterium for which of results are considered significant
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
#'
#'
#'
#TO-DO: reduce to base and ggplot, remove most arguments, only include an AB pair
plot_histogram_CE <- function(MIC_clean, results, t_rank = 1,
  one_direction = TRUE, CResponse = "CS", FDR_crit = 0.15, whichAB = NULL,
  colours = c("#283c82", "#F46B2D"), separate_plots = NULL, ran = NULL,
  t_or_effect = "t") {

  if (!is.null(whichAB)) {
    comb <- whichAB
    sign_dat <- results %>%
      filter(A == comb[1] & B == comb[2])
    row.names(sign_dat) <- paste0(sign_dat$A, sign_dat$B)
    if (whichAB[1] %in% sign_dat$A & whichAB[2] %in% sign_dat$B) {
      d <- sign_dat[paste0(comb, collapse = ""), "tau"]
    } else {
      message("Specified antibiotics are not in data, using highest effect size")
      whichAB <- NULL
    }
  }


  if (is.null(whichAB)) {
    if (t_or_effect == "effect"){
      results$t <- results$effect_size
    }

    sign_dat <- results %>%
      filter(p_BY < FDR_crit & (sign(t) == c(-1, 1)[CResponse == c("CS", "CR")])) %>%
      arrange(desc(abs(t)))
    row.names(sign_dat) <- paste0(sign_dat$A, sign_dat$B)

    if (nrow(sign_dat) < 1) {
      message("No significant ", CResponse, " effect, plotting non significant finding with largest T")
      if (CResponse == "CS") {
        sign_dat <- results %>%
          arrange(t) %>%
          slice(1)
      } else {
        sign_dat <- results %>%
          arrange(desc(t)) %>%
          slice(1)
      }
      t_rank <- 1
    }

    if (nrow(sign_dat) < t_rank) {
      message(paste0("No", t_rank, "th significant ", CResponse, " effect, plotting significant finding with largest T"))
      t_rank <- 1
    }
    comb <- unlist(sign_dat[t_rank, 1:2])
    d <- sign_dat[t_rank, "tau"]
  }



  #for all combinations?
  dat <- log2(MIC_clean[, comb])
  dat <- dat[!is.na(dat[, 1]) & !is.na(dat[, 2]), ]
  names(dat) <- c("A", "B")

  dat$Condition <- as.factor(ifelse(dat$B >= d,
                                    paste0(comb[1],"|", comb[2]," > ", round(d, 2)),
                                    paste0(comb[1],"|", comb[2]," < ", round(d, 2))))
  means <- data.frame(mean = c(mean(dat$A[dat$B < d]), mean(dat$A[dat$B >= d])),
                      Means = sort((unique(dat$Condition))))
  change_range <- FALSE
  if (is.null(ran)){
    ran <- range(dat$A) + (0.5 * c(-1, + 1))
    change_range <- TRUE
  }
  ticks <- floor(seq(ran[1], ran[2] + 2))

  plotPanel <- function(dat) {

    plotje <- dat %>%
      ggplot(aes(x = A, y = ..count.., group = Condition, fill = Condition)) +
      geom_bar(stat = "count", width = 0.6, position = "stack") +
      labs(x = bquote(log[2]*"(MIC) " ~ .(comb[1])), y = "Counts") +
      #labs(x = expression(log[2]*"(MIC) "*comb[1]), y = "Counts") +
      #scale_y_continuous(expand = expansion(mult = c(0, .05)), labels = function(x) sprintf("%.2f", x)) +
      scale_y_continuous(expand = expansion(mult = c(0, .05))) +
      scale_x_continuous(breaks = ticks, limits = ran) +
      geom_vline(data = means, aes(xintercept = mean), colour = "white") +
      geom_vline(data = means, aes(xintercept = mean, colour = Means), show.legend  = TRUE, linetype = 2) +
      scale_fill_manual(values = colours) +
      scale_colour_manual(labels = c(bquote(hat(mu)[.(as.character(means$Means[1]))]),
                                     #                               bquote(paste(hat(mu)[paste(comb[1], "|", comb[2],>=, d)]))),
                                     bquote(hat(mu)[.(as.character(means$Means[2]))])),
                          values = colours) +
      theme_bw() +
      theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(), legend.position = "bottom",
            #     legend.background = element_rect(fill = "transparent", size = 0.5, linetype = "solid", colour = "black"),
            legend.box = "horizontal", legend.direction = "vertical",
            legend.margin = margin(0,0,0,0),
            legend.box.margin = margin(-15,0,0,0),
            legend.background = element_blank()) +
      guides(fill = guide_legend(override.aes = list(linetype = 0), order = 1, title = NULL),
             colour = guide_legend(order = 2, title = NULL, label.theme = element_text(size = 12)))


    return(plotje)
  }


  p_A <- plotPanel(dat)

  if (one_direction) {
    return(p_A)
  }

  new_dat <- results
  rownames(new_dat) <- paste(new_dat$A, new_dat$B, sep = "_")
  comb <- comb[2:1]
  d <- new_dat[paste(comb, collapse = "_"), "tau"]
  dat <- log2(MIC_clean[, comb])
  dat <- dat[!is.na(dat[, 1]) & !is.na(dat[, 2]), ]
  names(dat) <- c("A", "B")
  d_min <- max(dat$B[dat$B < d], na.rm = T)

  dat$Condition <- as.factor(ifelse(dat$B >= d,
                                    paste0(comb[1],"|", comb[2]," > ", round(d, 2)),
                                    paste0(comb[1],"|", comb[2]," < ", round(d, 2))))
  means <- data.frame(mean = c(mean(dat$A[dat$B < d]), mean(dat$A[dat$B >= d])),
                      Means = sort((unique(dat$Condition))))

  if (change_range == TRUE){
    ran <- range(dat$A) + (0.5 * c(-1, + 1))
  }
  ticks <- floor(seq(ran[1], ran[2] + 2))
  p_B <- plotPanel(dat)



  if (is.null(separate_plots)) {
    p_A <- p_A + theme(axis.title.y = element_blank())
    p_B <- p_B + theme(axis.title.y = element_blank())
    return(grid.arrange(p_A, p_B, ncol = 2, left = textGrob("Counts", hjust = -0.45, rot = 90)))
  } else {
    return(list(p_A, p_B + ylab("  ")))
  }
}
