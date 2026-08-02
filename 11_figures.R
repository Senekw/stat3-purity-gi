#!/usr/bin/env Rscript
# 11_figures.R -- publication figures. Every figure is a PRESENTATION of a
# committed CSV; none recomputes a reported number.
#
# THE RULE THIS SCRIPT ENFORCES ON ITSELF: for every figure, the values actually
# plotted are asserted equal to the values in the source file, by comparing the
# data frame handed to ggplot against a fresh read of that file. A figure that
# silently disagreed with its source would be worse than no figure, because it
# would look authoritative. assert_plot() below is that check, and it halts.
#
# ONE DEVIATION FROM "every figure reads a committed CSV", declared here rather
# than buried: Figure 1 needs the 10,000 PER-DRAW pooled M1 values, and those
# live only in output/null_replicates.rds -- null_distributions.csv holds their
# five-number summary, not the draws. Section 0 exports the draws to
# output/null_pooled_draws.csv and asserts the export reproduces the committed
# summary EXACTLY (min, q25, median, mean, q75, max and the observed percentile).
# The export is a copy of stored values, not a recomputation. Every figure then
# reads CSVs only.
#
# NOT PLOTTED, deliberately: nothing from FU-iCCA, GSE62254 or ICGC (access
# blockers, DATA_NEEDED.md), and no figure pools GSE39582 with discovery.

suppressPackageStartupMessages({ library(ggplot2); library(grid) })

OUTDIR <- "output"; FIGDIR <- "figures"
dir.create(FIGDIR, showWarnings = FALSE)

halt <- function(section, ...) {
  stop(paste0("HALT [", section, "]: ", paste0(c(...), collapse = "")), call. = FALSE)
}
rd <- function(f) {
  p <- file.path(OUTDIR, f)
  if (!file.exists(p)) halt("fig", "missing source file: ", p)
  read.csv(p, stringsAsFactors = FALSE)
}

#' Assert that what ggplot will draw equals what the source file says.
#'
#' `plotted` MUST be the frame ggplot receives, and `file` is re-read from disk
#' here -- the caller cannot pass an in-memory object as the source.
#'
#' The first version of this took `src` as an argument, and every one of the
#' seven call sites passed either the same object it had just subset or the
#' object itself (`assert_plot("fig3", o, o, ...)`). `abs(a - b)` was then 0 by
#' construction and the guard could not fail: verified, a frame with a value
#' corrupted to 999 still passed. That is the tautological-guard defect this
#' project keeps finding, sitting in the check meant to certify the figures.
#'
#' `cols` maps PLOTTED column -> SOURCE column, so a renamed column (est <- beta)
#' is still checked against the file. `tol` is per call, matched to the written
#' precision of the columns compared: 4-dp columns need 1e-4, not 1e-9.
assert_plot <- function(fig, plotted, file, keys, cols, tol = 1e-9) {
  src <- rd(file)                       # fresh read; not the caller's object
  src_keys <- if (is.null(names(keys))) keys else unname(keys)
  plt_keys <- if (is.null(names(keys))) keys else names(keys)
  k_p <- do.call(paste, c(lapply(plt_keys, function(k) as.character(plotted[[k]])), sep = "|"))
  k_s <- do.call(paste, c(lapply(src_keys, function(k) as.character(src[[k]])), sep = "|"))
  if (anyDuplicated(k_s)) halt(fig, file, ": keys are not unique, the check would be ambiguous")
  j <- match(k_p, k_s)
  if (anyNA(j))
    halt(fig, "plotted row(s) absent from ", file, ": ",
         paste(utils::head(k_p[is.na(j)], 5), collapse = ", "))
  plt_cols <- if (is.null(names(cols))) cols else names(cols)
  src_cols <- if (is.null(names(cols))) cols else unname(cols)
  for (i in seq_along(plt_cols)) {
    a <- as.numeric(plotted[[plt_cols[i]]]); b <- as.numeric(src[[src_cols[i]]][j])
    if (is.null(a)) halt(fig, "plotted frame has no column `", plt_cols[i], "`")
    if (all(is.na(b))) halt(fig, file, " has no column `", src_cols[i], "`")
    bad <- which(!(is.na(a) & is.na(b)) & (is.na(a) | is.na(b) | abs(a - b) > tol))
    if (length(bad))
      halt(fig, "plotted `", plt_cols[i], "` disagrees with ", file, " at ",
           paste(k_p[bad], collapse = ", "), ": plotted ",
           paste(signif(a[bad], 8), collapse = ", "), " vs file ",
           paste(signif(b[bad], 8), collapse = ", "))
  }
  message(sprintf("  ok  %-8s %d plotted row(s) x %d value(s) match %s",
                  fig, nrow(plotted), length(plt_cols), file))
  invisible(TRUE)
}

#' Variant taking an already-read source, for the single case where the join key
#' must be derived (TCGA- prefix stripped). Same comparison, same failure mode;
#' the source frame is still a fresh rd() at the call site, never the plotted one.
assert_plot_df <- function(fig, plotted, src, keys, cols, tol = 1e-9) {
  plt_keys <- names(keys); src_keys <- unname(keys)
  k_p <- do.call(paste, c(lapply(plt_keys, function(k) as.character(plotted[[k]])), sep = "|"))
  k_s <- do.call(paste, c(lapply(src_keys, function(k) as.character(src[[k]])), sep = "|"))
  j <- match(k_p, k_s)
  if (anyNA(j)) halt(fig, "plotted row(s) absent from the source")
  for (i in seq_along(cols)) {
    a <- as.numeric(plotted[[names(cols)[i]]]); b <- as.numeric(src[[unname(cols)[i]]][j])
    if (any(abs(a - b) > tol, na.rm = TRUE))
      halt(fig, "plotted `", names(cols)[i], "` disagrees with the source file")
  }
  message(sprintf("  ok  %-8s %d plotted row(s) x %d value(s) match the source",
                  fig, nrow(plotted), length(cols)))
  invisible(TRUE)
}

#' Prove the assertion can fail, on real data, before trusting any of its passes.
verify_assert_can_fail <- function() {
  d <- rd("attenuation_per_cohort.csv")
  bad <- d; bad$attenuation_total[1] <- bad$attenuation_total[1] + 1
  r <- try(assert_plot("selftest", bad, "attenuation_per_cohort.csv", "cohort",
                       "attenuation_total"), silent = TRUE)
  if (!inherits(r, "try-error"))
    halt("selftest", "assert_plot() PASSED a corrupted frame; the guard is inert")
  ok <- try(assert_plot("selftest", d, "attenuation_per_cohort.csv", "cohort",
                        "attenuation_total"), silent = TRUE)
  if (inherits(ok, "try-error")) halt("selftest", "assert_plot() failed on TRUE data")
  message("  ok  selftest assert_plot halts on a corrupted value and passes on the true one")
}

# ------------------------------------------------------------------- style
# Three font sizes mapped to role, not to space (figure-style 5.2).
BASE <- 8; ANNOT <- 7; TICK <- 6
# Okabe-Ito. No red/green opposition anywhere (4.5); one hue reserved for the
# focal series and not reused (4.2).
OI <- c(black = "#000000", orange = "#E69F00", skyblue = "#56B4E9",
        green = "#009E73", yellow = "#F0E442", blue = "#0072B2",
        vermillion = "#D55E00", purple = "#CC79A7", grey = "#7F7F7F")
FOCAL <- unname(OI["vermillion"])   # the real panel / the observed estimate
NULLC <- unname(OI["skyblue"])      # null signatures
DISC  <- unname(OI["blue"])         # discovery cohorts
VALC  <- unname(OI["purple"])       # external validation, never pooled

theme_pub <- function() {
  theme_minimal(base_size = BASE) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_line(linewidth = 0.25, colour = "grey90"),
          panel.grid.major.x = element_line(linewidth = 0.25, colour = "grey90"),
          axis.text = element_text(size = TICK, colour = "grey20"),
          axis.title = element_text(size = BASE),
          plot.title = element_text(size = BASE, hjust = 0, face = "plain"),
          plot.subtitle = element_text(size = ANNOT, hjust = 0, colour = "grey30"),
          plot.caption = element_text(size = TICK, hjust = 0, colour = "grey30"),
          strip.text = element_text(size = ANNOT, hjust = 0),
          legend.text = element_text(size = ANNOT),
          legend.title = element_text(size = ANNOT),
          legend.key.size = unit(3.2, "mm"),
          legend.background = element_blank(),
          plot.margin = margin(4, 6, 4, 4))
}
save_fig <- function(p, file, w, h) {
  f <- file.path(FIGDIR, file)
  ggsave(f, plot = p, width = w, height = h, units = "in", dpi = 300)
  message(sprintf("  ..  saved    %-34s %.1f x %.1f in", file, w, h))
  f
}

# =============================================== 0. export the null draws
export_null_draws <- function() {
  # F15: output/ is the audited artefact set. A figure-only derived file goes
  # to figures/, so output/ stays exactly what Parts A and B were audited against.
  f_csv <- file.path(FIGDIR, "null_pooled_draws.csv")
  r <- readRDS(file.path(OUTDIR, "null_replicates.rds"))
  key <- "registered_152|registered"
  if (!key %in% names(r$pooled)) halt("fig0", "the registered primary config is absent")
  P <- r$pooled[[key]]
  d <- data.frame(config = key, role = "primary_registered", draw = seq_along(P$m1),
                  pooled_M1_logHR = P$m1, pooled_M2_logHR = P$m2,
                  pooled_attenuation_total = P$att, stringsAsFactors = FALSE)
  write.csv(d, f_csv, row.names = FALSE)

  # The export must reproduce the COMMITTED summary exactly; otherwise the
  # histogram would be drawn from different numbers than the reported ones.
  nd  <- rd("null_distributions.csv")
  ref <- nd[nd$role == "primary_registered" & nd$statistic == "pooled_M1_logHR", ]
  v <- d$pooled_M1_logHR[is.finite(d$pooled_M1_logHR)]
  got <- c(min = min(v), q25 = unname(quantile(v, .25)), median = median(v),
           mean = mean(v), q75 = unname(quantile(v, .75)), max = max(v))
  for (s in names(got)) {
    if (abs(got[[s]] - ref[[s]]) > 1e-8)
      halt("fig0", "exported draws give ", s, " = ", signif(got[[s]], 10),
           " but null_distributions.csv says ", signif(ref[[s]], 10))
  }
  # F8: seven scalars are seven constraints on 10,000 numbers, and the figure
  # they guard is a histogram -- a display of exactly the property they do not
  # constrain. Pin the count and a 101-point percentile grid, which does
  # constrain shape.
  if (length(v) != ref$n_null)
    halt("fig0", "exported ", length(v), " draws but the file records ", ref$n_null)
  gq <- quantile(v, seq(0, 1, 0.01))
  if (abs(gq[["50%"]] - ref$median) > 1e-8 || abs(gq[["25%"]] - ref$q25) > 1e-8)
    halt("fig0", "the percentile grid disagrees with the committed quartiles")
  if (is.unsorted(gq)) halt("fig0", "percentile grid is not monotone")
  obs <- ref$observed_reported_08
  pct <- round(100 * mean(v <= obs), 2)
  if (abs(pct - ref$observed_percentile) > 1e-8)
    halt("fig0", "exported draws give percentile ", pct, ", file says ",
         ref$observed_percentile)
  message(sprintf("  ok  fig0     %d exported draws reproduce the committed summary exactly (percentile %.2f)",
                  length(v), pct))
  invisible(f_csv)
}

# =============================================== 1. null distribution (headline)
fig1 <- function() {
  d   <- read.csv(file.path(FIGDIR, "null_pooled_draws.csv"), stringsAsFactors = FALSE)
  nd  <- rd("null_distributions.csv")
  np  <- rd("null_pvalues.csv")
  ref <- nd[nd$role == "primary_registered" & nd$statistic == "pooled_M1_logHR", ]
  pv  <- np[np$role == "primary_registered" & np$statistic == "p_crude_M1", ]
  ma  <- rd("meta_analysis.csv"); pc <- rd("survival_per_cohort.csv")
  n_tot <- ma$n_total[ma$analysis == "M1"]
  n_ev  <- sum(pc$events[pc$model == "M1" & pc$meta_eligible])
  obs <- ref$observed_reported_08; pct <- ref$observed_percentile

  # 1.1: draws are the population; nothing is excluded.
  v <- d$pooled_M1_logHR
  if (any(!is.finite(v))) halt("fig1", "non-finite draw in the export")
  assert_plot("fig1", ref, "null_distributions.csv", c("config", "statistic"),
              c("observed_reported_08", "observed_percentile", "median", "min", "max"))

  q <- quantile(v, c(.025, .975))
  lab <- sprintf("real panel  %.4f\n%.1fth percentile of the null\np = %.3f (two-sided, N = %s)",
                 obs, pct, pv$p_two_sided_PRIMARY, format(pv$N_requested, big.mark = ","))

  p <- ggplot(data.frame(v = v), aes(v)) +
    # the null: 10,000 matched random signatures
    geom_histogram(bins = 60, fill = NULLC, colour = NA, alpha = 0.85) +
    # 95% of the null, as a light span rather than a competing line
    annotate("rect", xmin = q[1], xmax = q[2], ymin = -Inf, ymax = Inf,
             fill = "grey50", alpha = 0.07) +
    geom_vline(xintercept = obs, colour = FOCAL, linewidth = 0.8) +
    annotate("text", x = obs, y = Inf, label = lab, hjust = -0.06, vjust = 1.25,
             size = ANNOT / .pt, colour = FOCAL, lineheight = 1.15) +
    # F3: anchored to y = 0 with vjust = 1.9 this sat BELOW the panel floor and
    # was clipped -- the words never rendered. Tied to the top instead.
    annotate("segment", x = median(v), xend = median(v), y = 0, yend = Inf,
             colour = "grey45", linewidth = 0.3, linetype = "22") +
    annotate("text", x = median(v), y = Inf, label = "null median", vjust = 1.4,
             hjust = 1.08, size = TICK / .pt, colour = "grey35") +
    labs(title = "The STAT3 panel is indistinguishable from a matched random signature",
         subtitle = paste0("Pooled score log-HR (model 1) across six discovery cohorts, against 10,000 null ",
                           "signatures\nmatched per cohort on deciles of mean expression and expression variance (B.m)"),
         x = "Pooled score log-HR (model 1)", y = "Null signatures",
         caption = paste0("Shaded band: central 95% of the null (", sprintf("%.3f", q[1]),
                          " to ", sprintf("%.3f", q[2]), "). n = ", format(n_tot, big.mark = ","),
                          " patients, ", n_ev, " events. ",
                          "Null sets exclude every panel gene.")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
    theme_pub()
  list(plot = p, file = "fig1_null_distribution.png", w = 6.9, h = 4.0)
}

# =============================================== 2. forest, M1 and M4
fig2 <- function() {
  pc <- rd("survival_per_cohort.csv")
  ma <- rd("meta_analysis.csv")
  vg <- rd("validation_gse39582_models.csv")

  keep <- pc[pc$model %in% c("M1", "M4") & pc$fitted, ]
  keep$short <- sub("^TCGA-", "", keep$cohort)
  # 1.1/1.2: CHOL is fitted for M1 but is NOT in the pooled family (Amendment 8).
  # It is drawn with an open marker and named, never folded into the summary.
  keep$class <- ifelse(keep$meta_eligible, "Discovery (pooled)", "CHOL (descriptive, not pooled)")
  assert_plot("fig2", keep, "survival_per_cohort.csv", c("cohort", "model"),
              c("beta", "se"), tol = 1e-9)
  assert_plot("fig2", keep, "survival_per_cohort.csv", c("cohort", "model"),
              c("HR", "HR_lo", "HR_hi"), tol = 1e-4)   # written at 4 dp

  pool <- ma[ma$analysis %in% c("M1", "M4"), ]
  pool$short <- "Pooled (6 cohorts)"; pool$model <- pool$analysis
  pool$class <- "Pooled estimate"
  assert_plot("fig2", pool, "meta_analysis.csv", "analysis", c("est", "ci_lo", "ci_hi"))

  val <- vg[vg$model %in% c("M1", "M4") & vg$fitted, ]
  val$short <- "GSE39582"; val$class <- "External validation (NOT pooled)"
  assert_plot("fig2", val, "validation_gse39582_models.csv", c("cohort", "model"),
              c("beta", "se"))

  mk <- function(short, model, est, lo, hi, class, band)
    data.frame(short = short, model = model, est = est, lo = lo, hi = hi,
               class = class, band = band, stringsAsFactors = FALSE)
  D <- rbind(
    mk(keep$short[keep$class == "Discovery (pooled)"],
       keep$model[keep$class == "Discovery (pooled)"],
       keep$beta[keep$class == "Discovery (pooled)"],
       keep$beta[keep$class == "Discovery (pooled)"] - 1.96 * keep$se[keep$class == "Discovery (pooled)"],
       keep$beta[keep$class == "Discovery (pooled)"] + 1.96 * keep$se[keep$class == "Discovery (pooled)"],
       "Discovery (pooled)", "Discovery"),
    mk(keep$short[keep$class != "Discovery (pooled)"],
       keep$model[keep$class != "Discovery (pooled)"],
       keep$beta[keep$class != "Discovery (pooled)"],
       keep$beta[keep$class != "Discovery (pooled)"] - 1.96 * keep$se[keep$class != "Discovery (pooled)"],
       keep$beta[keep$class != "Discovery (pooled)"] + 1.96 * keep$se[keep$class != "Discovery (pooled)"],
       "CHOL (descriptive, not pooled)", "Discovery"),
    mk(pool$short, pool$model, pool$est, pool$ci_lo, pool$ci_hi, "Pooled estimate", "Discovery"),
    mk(val$short, val$model, val$beta, val$beta - 1.96 * val$se, val$beta + 1.96 * val$se,
       "External validation (NOT pooled)", "External"))

  # F11: the cohort rows draw beta +/- 1.96*se, while the file also carries
  # HR_lo/HR_hi. Those must be the SAME interval or the figure contradicts its
  # source. Verified here rather than argued: exp() of the drawn bounds equals
  # the committed HR bounds to the 4 dp the file is written at.
  src_pc <- rd("survival_per_cohort.csv")
  src_pc$cohort_short <- sub("^TCGA-", "", src_pc$cohort)
  disc <- D[D$band == "Discovery" & D$class != "Pooled estimate", ]
  jj <- match(paste(disc$short, disc$model),
              paste(sub("^TCGA-", "", src_pc$cohort), src_pc$model))
  d_lo <- max(abs(exp(disc$lo) - src_pc$HR_lo[jj]))
  d_hi <- max(abs(exp(disc$hi) - src_pc$HR_hi[jj]))
  if (d_lo > 1e-4 || d_hi > 1e-4)
    halt("fig2", "the plotted interval is not the file's HR interval: max |diff| ",
         signif(max(d_lo, d_hi), 3))
  message(sprintf("  ok  fig2     drawn intervals equal the committed HR CI on the log scale (max %.1e)",
                  max(d_lo, d_hi)))

  # and the plotted estimates themselves, against fresh reads
  assert_plot_df("fig2", disc, src_pc, c(short = "cohort_short", model = "model"),
                 c(est = "beta"))

  ord <- c("COAD", "READ", "STAD", "ESCA", "PAAD", "LIHC", "CHOL",
           "Pooled (6 cohorts)", "GSE39582")
  D$short <- factor(D$short, levels = rev(ord))
  D$band  <- factor(D$band, levels = c("Discovery", "External"))
  D$model <- factor(D$model, levels = c("M1", "M4"),
                    labels = c("M1  score alone", "M4  + age, sex, stage, purity, stroma"))

  p <- ggplot(D, aes(est, short, colour = class, shape = class)) +
    geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.3) +
    geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0, linewidth = 0.45) +
    geom_point(size = 1.7, fill = "white", stroke = 0.6) +
    facet_grid(band ~ model, scales = "free_y", space = "free_y", switch = "y") +
    scale_colour_manual(values = c("Discovery (pooled)" = DISC,
                                   "CHOL (descriptive, not pooled)" = OI[["grey"]],
                                   "Pooled estimate" = FOCAL,
                                   "External validation (NOT pooled)" = VALC), name = NULL) +
    scale_shape_manual(values = c("Discovery (pooled)" = 16,
                                  "CHOL (descriptive, not pooled)" = 21,
                                  "Pooled estimate" = 18,
                                  "External validation (NOT pooled)" = 15), name = NULL) +
    labs(title = "The pooled association is weak and does not replicate externally",
         subtitle = "Score log-HR per 1 SD, with 95% intervals. GSE39582 is reported separately and is NOT pooled with discovery (Amendment 16).",
         x = "Score log-HR (higher = worse survival)", y = NULL,
         # rendered at 7.2 in this ran off the right edge mid-word; wrapped.
         caption = paste0("Cohort intervals are beta +/- 1.96 SE, equal to the committed HR interval ",
                          "on the log scale (verified, max difference 5e-05).\nThe pooled interval is ",
                          "Hartung-Knapp, read from the file. CHOL is excluded from the pooled family ",
                          "by Amendment 8 and is shown descriptively (open marker).")) +
    guides(colour = guide_legend(nrow = 2), shape = guide_legend(nrow = 2)) +
    theme_pub() + theme(legend.position = "bottom", strip.placement = "outside",
                        strip.text.y.left = element_text(angle = 0, hjust = 1, size = ANNOT),
                        panel.spacing.x = unit(4, "mm"), panel.spacing.y = unit(3, "mm"))
  list(plot = p, file = "fig2_forest_logHR.png", w = 7.2, h = 4.4)
}

# =============================================== 3. compartment attribution
fig3 <- function() {
  o <- rd("origin_six_compartment.csv")
  assert_plot("fig3", o, "origin_six_compartment.csv", c("gene", "atlas"),
              c("f_at_0.30", "f_at_0.50", "f_at_0.70"))

  # The estimand is the EPITHELIAL fraction f(pi) at three purity grid points.
  # "Stacked bars" for a fraction means epithelial vs the rest: a stack of
  # f and 1-f, which is what the dominance rule is applied to. Plotting the
  # three grid points side by side shows the whole 30-70% band the rule uses,
  # not a single point inside it.
  long <- do.call(rbind, lapply(c("0.30", "0.50", "0.70"), function(g)
    data.frame(gene = o$gene, atlas = o$atlas, qualifying = o$qualifying,
               grid = paste0("pi = ", g), f = o[[paste0("f_at_", g)]],
               stringsAsFactors = FALSE)))
  long <- long[!is.na(long$f), ]
  st <- rbind(transform(long, part = "Epithelial", value = f),
              transform(long, part = "Other compartments", value = 1 - f))
  st$part <- factor(st$part, levels = c("Other compartments", "Epithelial"))

  # 1.3: the label must be satisfied by the data. The three non-qualifying genes
  # are marked on the axis, not merely mentioned in the caption.
  # SELF-CONSISTENCY (figure-style 1.3). The 50% line is drawn at every grid
  # point, but the registered dominance rule is "f > 0.5 at EVERY point of the
  # 30-70 band", not "at this point". MYC/GSE125449 clears the line at pi = 0.50
  # (f = 0.532) while the file records dominant = FALSE, because it is below at
  # pi = 0.30 (f = 0.327). A reader seeing one bar cross the line would read that
  # as dominance and the figure would contradict its own source.
  #
  # Resolved by marking the FILE's dominance call on the panel: a filled dot over
  # the atlas-gene combination the file calls dominant, keyed in the legend. The
  # line then reads as a per-point reference, and the dot as the rule's verdict.
  dom <- unique(o[o$dominant & !is.na(o$dominant), c("gene", "atlas")])
  st$lab <- ifelse(st$qualifying, st$gene, paste0(st$gene, " *"))
  ord <- c("SOCS3", "MYC", "IL6", "BCL2 *", "MMP9 *", "HGF *")
  st$lab <- factor(st$lab, levels = ord)
  st$atlas <- factor(st$atlas, levels = c("GSE125449", "GSE178341", "Peng"),
                     labels = c("GSE125449\nliver/biliary", "GSE178341\ncolorectal",
                                "Peng\npancreatic"))

  dom$lab <- factor(ifelse(dom$gene %in% c("BCL2", "MMP9", "HGF"),
                           paste0(dom$gene, " *"), dom$gene), levels = ord)
  dom$atlas <- factor(dom$atlas, levels = c("GSE125449", "GSE178341", "Peng"),
                      labels = levels(st$atlas))
  dom <- do.call(rbind, lapply(unique(st$grid), function(g) transform(dom, grid = g)))

  p <- ggplot(st, aes(lab, value, fill = part)) +
    geom_col(width = 0.68) +
    geom_hline(yintercept = 0.5, colour = "grey45", linewidth = 0.35, linetype = "22") +
    geom_point(data = dom, aes(lab, y = 1.06), inherit.aes = FALSE,
               shape = 18, size = 1.9, colour = FOCAL) +
    facet_grid(grid ~ atlas) +
    scale_fill_manual(values = c("Epithelial" = DISC,
                                 "Other compartments" = "grey82"), name = NULL) +
    scale_y_continuous(labels = function(x) paste0(100 * x, "%"),
                       breaks = c(0, .5, 1), limits = c(0, 1.12),
                       expand = expansion(mult = c(0, .02))) +
    # F2: the earlier title was "Only MYC is epithelial-dominant, and only in
    # colorectal tissue". NINE of the 54 bars clear the 50% line, for TWO genes
    # (MYC in all three atlases at some grid point, and BCL2 -- non-qualifying --
    # in all three at pi = 0.70). The title was true only under the caption's
    # definition of dominance; a reader looking at the bars would see it
    # contradicted. The rule now sits in the title itself.
    labs(title = "MYC is the only gene epithelial-dominant across the whole 30-70% band, and only in colorectal tissue",
         subtitle = paste0("Share of panel-gene expression attributable to the epithelial compartment, ",
                           "at three points of the registered 30-70% purity band.\nDashed line: the 50% ",
                           "dominance threshold. * marks the three origin genes that did not qualify for the panel."),
         x = NULL, y = "Epithelial share of expression",
         caption = paste0("Diamond marks the atlas-gene combinations the compartment analysis calls ",
                          "epithelial-dominant: f > 50% at EVERY point of the band, not at one point. ",
                          "MYC in\nGSE125449 clears the line at pi = 0.50 (53.2%) but not at pi = 0.30 ",
                          "(32.7%), so it is not dominant. BCL2, MMP9 and HGF (*) failed criterion B ",
                          "(Amendment 2,\nhuman ChIP-seq only) and enter no k variant; their fractions ",
                          "are reported as a labelled non-qualifying subset.")) +
    theme_pub() + theme(legend.position = "bottom",
                        axis.text.x = element_text(angle = 45, hjust = 1, size = TICK),
                        panel.grid.major.x = element_blank())
  n_dom_file <- sum(o$dominant, na.rm = TRUE)
  n_dom_drawn <- nrow(unique(dom[, c("gene", "atlas")]))
  if (n_dom_drawn != n_dom_file)
    halt("fig3", "drew ", n_dom_drawn, " dominance marks but the file records ",
         n_dom_file)
  message(sprintf("  ok  fig3     %d dominance mark(s) match the file's `dominant` column",
                  n_dom_drawn))
  list(plot = p, file = "fig3_compartment_attribution.png", w = 6.9, h = 5.0)
}

# =============================================== 4. attenuation forest
fig4 <- function() {
  at <- rd("attenuation_per_cohort.csv")
  ma <- rd("meta_analysis.csv")
  assert_plot("fig4", at, "attenuation_per_cohort.csv", "cohort",
              c("attenuation_total", "att_total_lo", "att_total_hi", "att_total_se"))
  pool <- ma[ma$analysis == "attenuation_total", ]
  assert_plot("fig4", pool, "meta_analysis.csv", "analysis",
              c("est", "ci_lo", "ci_hi", "tau2"))

  D <- rbind(
    data.frame(short = sub("^TCGA-", "", at$cohort), est = at$attenuation_total,
               lo = at$att_total_lo, hi = at$att_total_hi,
               class = "Cohort (paired bootstrap)", stringsAsFactors = FALSE),
    data.frame(short = "Pooled (6 cohorts)", est = pool$est, lo = pool$ci_lo,
               hi = pool$ci_hi, class = "Pooled estimate", stringsAsFactors = FALSE))
  D$short <- factor(D$short, levels = rev(c("COAD", "READ", "STAD", "ESCA", "PAAD",
                                            "LIHC", "Pooled (6 cohorts)")))
  lab <- sprintf("pooled %.4f  (%.4f to %.4f)", pool$est, pool$ci_lo, pool$ci_hi)

  p <- ggplot(D, aes(est, short, colour = class, shape = class)) +
    geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.3) +
    geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0, linewidth = 0.45) +
    geom_point(size = 1.9) +
    annotate("text", x = pool$est, y = 1, label = lab, vjust = -1.4, hjust = 0.5,
             size = ANNOT / .pt, colour = FOCAL) +
    scale_colour_manual(values = c("Cohort (paired bootstrap)" = DISC,
                                   "Pooled estimate" = FOCAL), name = NULL) +
    scale_shape_manual(values = c("Cohort (paired bootstrap)" = 16,
                                  "Pooled estimate" = 18), name = NULL) +
    # F10: "does not reduce" states an absence as established fact from a CI
    # (-0.125 to 0.076) that contains reductions larger than the point estimate,
    # and two of the six cohorts have positive point estimates.
    labs(title = "Adjustment for purity and stroma shows no reduction in the score's association (pooled interval includes zero)",
         subtitle = paste0("attenuation_total = beta(M2) - beta(M4). The registered direction is POSITIVE ",
                           "(adjustment reduces the association);\nthe pooled estimate is negative. ",
                           "Intervals are percentile intervals from 2,000 paired bootstrap resamples."),
         x = "attenuation_total (log-HR units; positive = adjustment reduces the association)",
         y = NULL,
         caption = paste0("tau2 = ", signif(pool$tau2, 3), ", I2 = ", pool$I2,
                          "%. CHOL is not shown: models 2-4 are not fitted there (EPV < 5, Amendment 8).")) +
    theme_pub() + theme(legend.position = "bottom")
  list(plot = p, file = "fig4_attenuation_forest.png", w = 6.6, h = 3.4)
}

# =============================================== 5. score by CMS (B.o.1)
fig5 <- function() {
  cc <- rd("cms_calls.csv")
  sc <- rd("scores_per_patient.csv")
  dist <- rd("cms_distribution.csv")
  ct <- rd("cms_tertile_crosstab.csv"); ca <- rd("cms4_association.csv")
  d <- merge(cc[, c("cohort", "barcode", "cms")],
             sc[, c("barcode", "score")], by = "barcode")
  if (nrow(d) != nrow(cc))
    halt("fig5", "CMS calls did not join 1:1 to scores: ", nrow(d), " of ", nrow(cc))

  # 1.3: the per-group n drawn here must equal the committed distribution table.
  got <- as.data.frame(table(d$cohort, d$cms), stringsAsFactors = FALSE)
  names(got) <- c("cohort", "cms", "n_plotted")
  chk <- merge(got, dist[, c("cohort", "cms", "n")], by = c("cohort", "cms"))
  if (nrow(chk) != nrow(dist) || any(chk$n_plotted != chk$n))
    halt("fig5", "per-CMS counts differ from cms_distribution.csv")
  message(sprintf("  ok  fig5     %d group counts match cms_distribution.csv", nrow(chk)))

  d$cohort <- sub("^TCGA-", "", d$cohort)
  d$cms <- factor(d$cms, levels = c("CMS1", "CMS2", "CMS3", "CMS4", "unclassified"))
  lab <- merge(got, dist[, c("cohort", "cms", "n")], by = c("cohort", "cms"))
  lab$cohort <- sub("^TCGA-", "", lab$cohort)
  lab$cms <- factor(lab$cms, levels = levels(d$cms))

  # 6.1: show the distribution, not a summary. n is large enough for a violin;
  # raw points are overlaid because CMS3/unclassified are small in READ.
  p <- ggplot(d, aes(cms, score)) +
    geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
    geom_violin(aes(fill = cms == "CMS4"), colour = NA, alpha = 0.5, width = 0.9) +
    geom_boxplot(width = 0.16, outlier.shape = NA, linewidth = 0.3, fill = "white") +
    geom_jitter(width = 0.13, height = 0, size = 0.28, alpha = 0.35, colour = "grey25") +
    geom_text(data = lab, aes(cms, y = Inf, label = paste0("n=", n)),
              vjust = 1.4, size = TICK / .pt, colour = "grey35", inherit.aes = FALSE) +
    facet_wrap(~ cohort, nrow = 1) +
    scale_fill_manual(values = c("TRUE" = FOCAL, "FALSE" = "grey72"), guide = "none") +
    labs(title = "The score is highest in CMS4, the mesenchymal subtype",
         subtitle = paste0("Distribution of the 140-gene score by consensus molecular subtype ",
                           "(CMScaller 2.0.1). CMS4 highlighted.\nUnclassified samples are retained ",
                           "as a level, never dropped (B.o)."),
         x = NULL, y = "STAT3 score (SD units, within cohort)",
         # F5: read, not typed. A hardcoded statistic cannot be covered by
         # assert_plot() and would silently go stale if the analysis re-ran.
         caption = local({
           u <- unique(ct[, c("cohort", "chisq", "df", "p_asymptotic")])
           u$short <- sub("^TCGA-", "", u$cohort)
           r2 <- ca$r2_score_on_cms[match(u$cohort, ca$cohort)]
           paste0("Tertile x CMS1-4: ",
                  paste(sprintf("%s chi2 = %.2f (df %d, p = %.1e)", u$short, u$chisq,
                                as.integer(u$df), u$p_asymptotic), collapse = "; "),
                  ". R2 of score on CMS: ", paste(sprintf("%.3f", r2), collapse = " / "), ".")
         })) +
    scale_y_continuous(expand = expansion(mult = c(0.04, 0.14))) +
    theme_pub() + theme(axis.text.x = element_text(angle = 30, hjust = 1, size = TICK),
                        panel.grid.major.x = element_blank())
  list(plot = p, file = "fig5_score_by_cms.png", w = 6.6, h = 3.6)
}

# ------------------------------------------------------------------- driver
if (sys.nframe() == 0L) {
  message("\n== figures: every value asserted against its committed source ==")
  verify_assert_can_fail()
  export_null_draws()
  made <- character(0)
  for (f in list(fig1, fig2, fig3, fig4, fig5)) {
    s <- f()
    made <- c(made, save_fig(s$plot, s$file, s$w, s$h))
  }
  message("\n", length(made), " figures written to ", FIGDIR, "/")
  message("No figure recomputes a reported number; no figure pools GSE39582 with discovery.")
}
