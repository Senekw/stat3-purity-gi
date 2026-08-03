#!/usr/bin/env Rscript
# 11b_figures_refined.R -- refined figure set. PRESENTATION ONLY.
#
# 11_figures.R is NOT modified and its outputs in figures/ are NOT touched.
# Everything here is written to figures/refined/.
#
# No reported number changes. Every plotted value is read from a committed CSV
# and asserted against a FRESH read of that file by assert_plot(), which is
# proven able to fail before any figure is drawn.
#
# TWO DISCREPANCIES between the run order and the committed record, resolved in
# favour of the files and annotated where they appear:
#   - the run order says the null is "centred at 0.118". The committed null
#     median is 0.117486 and the mean 0.117481; both round to 0.117. Panel A
#     reads the value from null_distributions.csv and prints 0.117.
#   - the run order names exploratory_crc_pooled.csv; the committed file is
#     output/exploratory_colorectal_pooled.csv. Its HR 0.9766 (0.8756-1.0892)
#     matches the quoted "0.977 (0.876-1.089)" exactly -- the WALD interval.

suppressPackageStartupMessages({ library(ggplot2); library(grid); library(cowplot) })

OUTDIR <- "output"
FIGDIR <- "figures/refined"
TAG    <- "EXPLORATORY_POSTHOC"

halt <- function(fig, ...) stop(paste0("HALT [", fig, "]: ", paste0(c(...), collapse="")), call.=FALSE)
rd   <- function(f) {
  p <- file.path(OUTDIR, f)
  if (!file.exists(p)) halt("read", "missing committed source: ", p)
  read.csv(p, stringsAsFactors = FALSE)
}

# ----------------------------------------------------------------- palette
# Okabe-Ito. ACCENT is reserved for the real panel / observed value and is used
# for nothing else in any figure (figure-style 4.1/4.2).
OI <- c(black="#000000", orange="#E69F00", skyblue="#56B4E9", green="#009E73",
        yellow="#F0E442", blue="#0072B2", vermillion="#D55E00", purple="#CC79A7")
ACCENT   <- unname(OI["vermillion"])   # the panel / observed value. Nothing else.
DISC     <- unname(OI["blue"])         # discovery cohorts
EXTERNAL <- unname(OI["purple"])       # external validation, never pooled
NEUTRAL  <- "grey72"                   # nulls and reference material
INK      <- "grey15"
MUTED    <- "grey40"

# --------------------------------------------------------------- typography
# Three roles, one family (figure-style 5.2): base 10, annotation 9, tick 9,
# title 13 bold, caption 8. Titles and captions left-aligned to the panel.
PT_TITLE <- 13; PT_SUB <- 10; PT_AXIS <- 10; PT_TICK <- 9; PT_ANNOT <- 9; PT_CAP <- 8

theme_pub <- function(grid_y = TRUE, grid_x = FALSE) {
  theme_minimal(base_size = PT_AXIS, base_family = "sans") +
    theme(
      plot.title       = element_text(size = PT_TITLE, face = "bold", colour = INK,
                                      hjust = 0, margin = margin(b = 3)),
      plot.subtitle    = element_text(size = PT_SUB, colour = MUTED, hjust = 0,
                                      margin = margin(b = 8), lineheight = 1.15),
      plot.caption     = element_text(size = PT_CAP, colour = MUTED, hjust = 0,
                                      margin = margin(t = 8), lineheight = 1.15),
      plot.caption.position = "plot", plot.title.position = "plot",
      axis.title       = element_text(size = PT_AXIS, colour = INK),
      axis.text        = element_text(size = PT_TICK, colour = MUTED),
      legend.text      = element_text(size = PT_ANNOT, colour = INK),
      legend.title     = element_blank(),
      legend.key.size  = unit(3.6, "mm"),
      legend.margin    = margin(0, 0, 0, 0),
      strip.text       = element_text(size = PT_ANNOT, colour = INK, hjust = 0),
      panel.grid.major.y = if (grid_y) element_line(colour = "grey92", linewidth = 0.3) else element_blank(),
      panel.grid.major.x = if (grid_x) element_line(colour = "grey92", linewidth = 0.3) else element_blank(),
      panel.grid.minor = element_blank(),
      panel.border     = element_blank(),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background  = element_rect(fill = "white", colour = NA))
}

#' Assert that what ggplot will draw equals what the source file says.
#' `file` is re-read HERE; the caller cannot pass an in-memory object as source.
assert_plot <- function(fig, plotted, file, keys, cols, tol = 1e-9) {
  src <- rd(file)
  pk <- if (is.null(names(keys))) keys else names(keys)
  sk <- if (is.null(names(keys))) keys else unname(keys)
  k_p <- do.call(paste, c(lapply(pk, function(k) as.character(plotted[[k]])), sep="|"))
  k_s <- do.call(paste, c(lapply(sk, function(k) as.character(src[[k]])), sep="|"))
  if (anyDuplicated(k_s)) halt(fig, file, ": source keys are not unique")
  j <- match(k_p, k_s)
  if (anyNA(j)) halt(fig, "plotted row(s) absent from ", file, ": ",
                     paste(utils::head(k_p[is.na(j)], 5), collapse=", "))
  pc <- if (is.null(names(cols))) cols else names(cols)
  sc <- if (is.null(names(cols))) cols else unname(cols)
  for (i in seq_along(pc)) {
    a <- as.numeric(plotted[[pc[i]]]); b <- as.numeric(src[[sc[i]]][j])
    if (all(is.na(b))) halt(fig, file, " has no column `", sc[i], "`")
    bad <- which(!(is.na(a) & is.na(b)) & (is.na(a) | is.na(b) | abs(a-b) > tol))
    if (length(bad))
      halt(fig, "plotted `", pc[i], "` disagrees with ", file, " at ",
           paste(k_p[bad], collapse=", "), ": ", paste(signif(a[bad],8), collapse=","),
           " vs ", paste(signif(b[bad],8), collapse=","))
  }
  message(sprintf("  ok  %-14s %3d row(s) x %d value(s) match %s", fig, nrow(plotted), length(pc), file))
  invisible(TRUE)
}

#' Prove the assertion can fail on real data before trusting any of its passes.
verify_assert_can_fail <- function() {
  d <- rd("attenuation_per_cohort.csv"); bad <- d
  bad$attenuation_total[1] <- bad$attenuation_total[1] + 1
  if (!inherits(try(assert_plot("selftest", bad, "attenuation_per_cohort.csv",
                                "cohort", "attenuation_total"), silent=TRUE), "try-error"))
    halt("selftest", "assert_plot() PASSED a corrupted frame; the guard is inert")
  if (inherits(try(assert_plot("selftest", d, "attenuation_per_cohort.csv",
                               "cohort", "attenuation_total"), silent=TRUE), "try-error"))
    halt("selftest", "assert_plot() failed on TRUE data")
  message("  ok  selftest      halts on a corrupted value, passes on the true one")
}

MM <- function(mm) mm / 25.4
save_fig <- function(p, file, w_mm, h_mm) {
  dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)
  f <- file.path(FIGDIR, file)
  ggsave(f, plot = p, width = MM(w_mm), height = MM(h_mm), dpi = 300, bg = "white")
  message(sprintf("  ..  saved         %-42s %.0f x %.0f mm", file, w_mm, h_mm))
  f
}

# ============================================================ panel builders
# Each returns a ggplot. Titles/subtitles are set by the caller for the main
# figure (where panel tags carry the role) and inside for the standalone.

#' Export the 10,000 per-draw pooled M1 values. They exist only in the .rds;
#' null_distributions.csv holds their summary. The export is asserted to
#' reproduce that summary before any figure reads it.
export_null_draws <- function() {
  f <- file.path(FIGDIR, "null_pooled_draws.csv")
  r <- readRDS(file.path(OUTDIR, "null_replicates.rds"))
  P <- r$pooled[["registered_152|registered"]]
  if (is.null(P)) halt("fig0", "the registered primary configuration is absent")
  d <- data.frame(draw = seq_along(P$m1), pooled_M1_logHR = P$m1)
  dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)
  write.csv(d, f, row.names = FALSE)
  nd  <- rd("null_distributions.csv")
  ref <- nd[nd$role=="primary_registered" & nd$statistic=="pooled_M1_logHR", ]
  v <- d$pooled_M1_logHR
  if (length(v) != ref$n_null) halt("fig0", "exported ", length(v), " draws, file says ", ref$n_null)
  got <- c(min=min(v), q25=unname(quantile(v,.25)), median=median(v), mean=mean(v),
           q75=unname(quantile(v,.75)), max=max(v))
  for (s in names(got)) if (abs(got[[s]] - ref[[s]]) > 1e-8)
    halt("fig0", "exported draws give ", s, " = ", signif(got[[s]],10),
         " but the committed summary says ", signif(ref[[s]],10))
  if (abs(round(100*mean(v <= ref$observed_reported_08), 2) - ref$observed_percentile) > 1e-8)
    halt("fig0", "exported draws do not reproduce the committed percentile")
  gq <- quantile(v, seq(0,1,0.01))
  if (is.unsorted(gq)) halt("fig0", "percentile grid is not monotone")
  message(sprintf("  ok  fig0           %d draws reproduce the committed summary (percentile %.2f)",
                  length(v), ref$observed_percentile))
  invisible(f)
}

panel_null <- function(standalone = FALSE) {
  d  <- read.csv(file.path(FIGDIR, "null_pooled_draws.csv"), stringsAsFactors = FALSE)
  nd <- rd("null_distributions.csv"); np <- rd("null_pvalues.csv")
  ma <- rd("meta_analysis.csv");      pc <- rd("survival_per_cohort.csv")
  ref <- nd[nd$role=="primary_registered" & nd$statistic=="pooled_M1_logHR", ]
  pv  <- np[np$role=="primary_registered" & np$statistic=="p_crude_M1", ]
  assert_plot("figA-null", ref, "null_distributions.csv", c("config","statistic"),
              c("observed_reported_08","observed_percentile","median","min","max"))
  v <- d$pooled_M1_logHR
  if (any(!is.finite(v))) halt("figA", "non-finite draw in the export")
  obs <- ref$observed_reported_08; med <- ref$median
  q <- quantile(v, c(.025,.975))
  # These two reach the SUBTITLE, so they are reported numbers. n_total is read
  # from a committed column; the event total has NO committed column, so its
  # per-cohort components are asserted instead of the sum being trusted.
  n_tot <- ma$n_total[ma$analysis=="M1"]
  ev <- pc[pc$model=="M1" & pc$meta_eligible, ]
  assert_plot("figA-events", ev, "survival_per_cohort.csv", c("cohort","model"), c("events","n"))
  n_ev <- sum(ev$events)
  if (sum(ev$n) != n_tot)
    halt("figA", "per-cohort n sums to ", sum(ev$n), " but meta_analysis.csv says ", n_tot)

  p <- ggplot(data.frame(v=v), aes(v)) +
    annotate("rect", xmin=q[1], xmax=q[2], ymin=-Inf, ymax=Inf, fill="grey50", alpha=0.07) +
    geom_histogram(bins=60, fill=NEUTRAL, colour=NA) +
    # the null's centre, read from the file -- the run order said 0.118; it is 0.117
    annotate("segment", x=med, xend=med, y=0, yend=Inf, colour=MUTED,
             linewidth=0.3, linetype="22") +
    annotate("text", x=med, y=Inf, hjust=1.06, vjust=1.5, size=PT_TICK/.pt, colour=MUTED,
             label=sprintf("null centred at %.3f,\nnot at zero", med)) +
    geom_vline(xintercept=obs, colour=ACCENT, linewidth=0.9) +
    annotate("text", x=obs, y=Inf, hjust=-0.07, vjust=1.5, size=PT_ANNOT/.pt,
             colour=ACCENT, lineheight=1.2,
             label=sprintf("real panel  %.4f\n%.1fth percentile\np = %.3f", obs,
                           ref$observed_percentile, pv$p_two_sided_PRIMARY)) +
    scale_y_continuous(expand=expansion(mult=c(0,0.20))) +
    labs(x="Pooled score log-HR (model 1)", y="Null signatures") +
    theme_pub()
  if (standalone)
    p <- p + labs(title="The panel is indistinguishable from a matched random signature",
                  subtitle=sprintf(paste0("Pooled score log-HR (model 1), %s patients and %d events ",
                                          "across six discovery cohorts,\nagainst %s null signatures ",
                                          "matched per cohort on deciles of mean expression and variance (B.m)"),
                                   format(n_tot, big.mark=","), n_ev,
                                   format(pv$N_requested, big.mark=",")),
                  caption=sprintf(paste0("Shaded band: central 95%% of the null (%.3f to %.3f). ",
                                         "Null sets exclude every panel gene."), q[1], q[2]))
  p
}

panel_forest <- function(models = c("M1","M4"), standalone = FALSE, compact = FALSE) {
  pc <- rd("survival_per_cohort.csv"); ma <- rd("meta_analysis.csv")
  vg <- rd("validation_gse39582_models.csv")
  keep <- pc[pc$model %in% models & pc$fitted, ]
  keep$short <- sub("^TCGA-","",keep$cohort)
  keep$class <- ifelse(keep$meta_eligible, "Discovery (pooled)", "CHOL (descriptive, not pooled)")
  assert_plot("figB-cohorts", keep, "survival_per_cohort.csv", c("cohort","model"), c("beta","se"))
  pool <- ma[ma$analysis %in% models, ]
  assert_plot("figB-pooled", pool, "meta_analysis.csv", "analysis", c("est","ci_lo","ci_hi"))
  val <- vg[vg$model %in% models & vg$fitted, ]
  assert_plot("figB-external", val, "validation_gse39582_models.csv", c("cohort","model"), c("beta","se"))

  mk <- function(sh, mo, e, lo, hi, cl, bd)
    data.frame(short=sh, model=mo, est=e, lo=lo, hi=hi, class=cl, band=bd, stringsAsFactors=FALSE)
  D <- rbind(
    mk(keep$short, keep$model, keep$beta, keep$beta-1.96*keep$se, keep$beta+1.96*keep$se,
       keep$class, "Discovery"),
    mk(rep("Pooled (6 cohorts)", nrow(pool)), pool$analysis, pool$est, pool$ci_lo, pool$ci_hi,
       "Pooled estimate", "Discovery"),
    mk(rep("GSE39582", nrow(val)), val$model, val$beta, val$beta-1.96*val$se,
       val$beta+1.96*val$se, "External validation (NOT pooled)", "External"))
  # the drawn interval IS the file's interval, on the log scale
  src <- pc; src$short <- sub("^TCGA-","",src$cohort)
  dd <- D[D$band=="Discovery" & D$class!="Pooled estimate", ]
  j <- match(paste(dd$short, dd$model), paste(src$short, src$model))
  if (max(abs(exp(dd$lo)-src$HR_lo[j])) > 1e-4 || max(abs(exp(dd$hi)-src$HR_hi[j])) > 1e-4)
    halt("figB", "the plotted interval is not the committed HR interval")

  ord <- c("COAD","READ","STAD","ESCA","PAAD","LIHC","CHOL","Pooled (6 cohorts)","GSE39582")
  D$short <- factor(D$short, levels=rev(ord))
  D$band  <- factor(D$band, levels=c("Discovery","External"))
  D$model <- factor(D$model, levels=c("M1","M4"),
                    labels=c("M1  score alone","M4  + age, sex, stage, purity, stroma"))
  cols <- c("Discovery (pooled)"=DISC, "CHOL (descriptive, not pooled)"=NEUTRAL,
            "Pooled estimate"=ACCENT, "External validation (NOT pooled)"=EXTERNAL)
  shp  <- c("Discovery (pooled)"=16, "CHOL (descriptive, not pooled)"=21,
            "Pooled estimate"=18, "External validation (NOT pooled)"=15)
  p <- ggplot(D, aes(est, short, colour=class, shape=class)) +
    geom_vline(xintercept=0, colour="grey60", linewidth=0.3) +
    geom_errorbar(aes(xmin=lo, xmax=hi), orientation="y", width=0, linewidth=0.45) +
    geom_point(size=1.9, fill="white", stroke=0.6) +
    scale_colour_manual(values=cols) + scale_shape_manual(values=shp) +
    guides(colour=guide_legend(nrow=2, byrow=TRUE), shape=guide_legend(nrow=2, byrow=TRUE)) +
    labs(x=if (compact) "Score log-HR per 1 SD" else
             "Score log-HR per 1 SD  (higher = worse survival)", y=NULL) +
    theme_pub(grid_y=FALSE, grid_x=TRUE) +
    theme(legend.position="bottom", panel.spacing.y=unit(3,"mm"),
          # without strip.placement="outside" the band label renders INSIDE the
          # panel, on top of the cohort tick labels (caught in the render).
          strip.placement="outside",
          strip.text.y.left=element_text(angle=0, hjust=1, size=PT_ANNOT))
  p <- if (length(models) > 1)
    p + facet_grid(band ~ model, scales="free_y", space="free_y", switch="y")
  else
    p + facet_grid(band ~ ., scales="free_y", space="free_y", switch="y")
  # In the composed figure the band strip label overran the panel width; the
  # gap plus the legend already distinguish external from discovery there.
  if (compact) p <- p + theme(strip.text.y.left=element_blank(),
                              legend.direction="vertical",
                              legend.box.margin=margin(t=-2)) +
                        guides(colour=guide_legend(ncol=1), shape=guide_legend(ncol=1))
  has_chol <- any(D$class == "CHOL (descriptive, not pooled)")
  if (standalone)
    p <- p + labs(title=if (identical(models, "M4"))
                    "Adjusted for purity and stroma, the association is no weaker"
                  else "The pooled association is weak and does not replicate externally",
                  subtitle=paste0("Score log-HR per 1 SD with 95% intervals. GSE39582 is reported separately and is\n",
                                  "NOT pooled with discovery (Amendment 16)",
                                  if (identical(models, "M4"))
                                    ". Model 4 adds age, sex, stage, purity and stroma." else "."),
                  # the CHOL clause is stated ONLY when a CHOL marker is drawn:
                  # models 2-4 are not fitted there, so the M4-only supplement
                  # has no such marker to describe.
                  caption=paste0("Cohort intervals are beta +/- 1.96 SE, equal to the committed HR interval on the log\n",
                                 "scale. The pooled interval is Hartung-Knapp, read from the file. ",
                                 if (has_chol)
                                   paste0("CHOL is excluded from\nthe pooled family by Amendment 8 and is shown ",
                                          "descriptively (open marker).")
                                 else
                                   paste0("CHOL is absent:\nmodels 2-4 are not fitted there (EPV < 5, Amendment 8).")))
  p
}

panel_fuicca <- function(standalone = FALSE) {
  cr <- rd("fuicca_correlations.csv"); pp <- rd("fuicca_per_patient.csv")
  pr <- cr[cr$comparison=="score_140 vs STAT3:Y705" & cr$method=="pearson", ]
  sp <- cr[cr$comparison=="score_140 vs STAT3:Y705" & cr$method=="spearman", ]
  assert_plot("figC-fuicca", pr, "fuicca_correlations.csv", c("comparison","method"),
              c("n","r","ci_lo","ci_hi"))
  d <- pp[is.finite(pp$STAT3_Y705) & is.finite(pp$score_140), ]
  if (nrow(d) != pr$n) halt("figC", "scatter would plot ", nrow(d), " points, the file says ", pr$n)
  p <- ggplot(d, aes(score_140, STAT3_Y705)) +
    geom_point(shape=21, fill=paste0(DISC,"55"), colour=DISC, size=1.5, stroke=0.35) +
    geom_smooth(method="lm", formula=y~x, se=FALSE, colour=ACCENT, linewidth=0.8) +
    annotate("text", x=-Inf, y=Inf, hjust=-0.06, vjust=1.4, size=PT_ANNOT/.pt, colour=INK,
             lineheight=1.25,
             label=sprintf("Pearson  r = %.3f  (%.3f to %.3f)\nSpearman r = %.3f  (%.3f to %.3f)\nn = %d",
                           pr$r, pr$ci_lo, pr$ci_hi, sp$r, sp$ci_lo, sp$ci_hi, pr$n)) +
    scale_y_continuous(expand=expansion(mult=c(0.05,0.16))) +
    labs(x="STAT3 score, 140 genes (SD units, within cohort)",
         y="STAT3 pY705 (median-normalised log2)") +
    theme_pub(grid_y=TRUE, grid_x=TRUE)
  if (standalone)
    p <- p + labs(title="The score tracks measured STAT3 pY705 in FU-iCCA",
                  subtitle=sprintf(paste0("Amendment 16 PRIMARY validation, n = %d of the 208 patients with both assays;\n",
                                          "94 excluded for a missing Y705 value, none imputed."), pr$n),
                  caption="Line: ordinary least squares. Intrahepatic cholangiocarcinoma, not one of the six discovery cohorts.")
  p
}

# ===================================================== new: colorectal pooling
fig_crc <- function() {
  # the run order names exploratory_crc_pooled.csv; the committed file is this one
  po <- rd("exploratory_colorectal_pooled.csv"); ip <- rd("exploratory_colorectal_inputs.csv")
  assert_plot("figCRC-pool", po, "exploratory_colorectal_pooled.csv", "analysis",
              c("est","HR","HR_lo_wald","HR_hi_wald","I2","tau2"))
  assert_plot("figCRC-inputs", ip, "exploratory_colorectal_inputs.csv", "cohort",
              c("beta","se","n","events"))
  # the subtitle quotes n_total and events_total, so check they are the sums
  if (po$n_total != sum(ip$n) || po$events_total != sum(ip$events))
    halt("figCRC", "pooled n/events do not equal the sum of the inputs")
  ip$short <- sub("^TCGA-","",ip$cohort)
  D <- rbind(
    data.frame(short=ip$short, est=exp(ip$beta), lo=exp(ip$beta-1.96*ip$se),
               hi=exp(ip$beta+1.96*ip$se),
               class=ifelse(ip$cohort=="GSE39582","External (NOT pooled with discovery)","Discovery cohort"),
               lab=sprintf("%d / %d", ip$events, ip$n), stringsAsFactors=FALSE),
    data.frame(short="Pooled (3 cohorts)", est=po$HR, lo=po$HR_lo_wald, hi=po$HR_hi_wald,
               class="Pooled (exploratory)",
               lab=sprintf("%d / %d", po$events_total, po$n_total), stringsAsFactors=FALSE))
  D$short <- factor(D$short, levels=rev(c("COAD","READ","GSE39582","Pooled (3 cohorts)")))
  p <- ggplot(D, aes(est, short, colour=class, shape=class)) +
    geom_vline(xintercept=1, colour="grey60", linewidth=0.3) +
    # the compatible ceiling. On a DISCRETE y axis, numeric y in annotate() is a
    # continuous-scale error -- geom_vline spans the panel without needing one.
    geom_vline(xintercept=po$HR_hi_wald, colour=ACCENT, linewidth=0.4, linetype="22") +
    annotate("text", x=po$HR_hi_wald, y=Inf, hjust=-0.07, vjust=1.25, colour=ACCENT,
             size=PT_TICK/.pt, lineheight=1.15,
             label=sprintf("compatible ceiling:\nHR %.3f per SD", po$HR_hi_wald)) +
    geom_errorbar(aes(xmin=lo, xmax=hi), orientation="y", width=0, linewidth=0.45) +
    geom_point(size=2.0) +
    # events/n in the right margin. The first draft passed BOTH a positional
    # `aes()` and a named `mapping=`, which ggplot read as a `stat` -- fixed to a
    # single mapping.
    geom_text(data=D, mapping=aes(x=Inf, y=short, label=lab), hjust=1.06,
              size=PT_TICK/.pt, colour=MUTED, inherit.aes=FALSE) +
    scale_x_continuous(limits=c(0.58, 1.78), breaks=c(0.6,0.8,1.0,1.2,1.4)) +
    scale_y_discrete(expand=expansion(add=c(0.55, 1.15))) +
    scale_colour_manual(values=c("Discovery cohort"=DISC,
                                 "External (NOT pooled with discovery)"=EXTERNAL,
                                 "Pooled (exploratory)"=ACCENT)) +
    scale_shape_manual(values=c("Discovery cohort"=16,
                                "External (NOT pooled with discovery)"=15,
                                "Pooled (exploratory)"=18)) +
    # The first title said "excludes ... above 9%" while the caption said a bound
    # is "not a proof of absence" -- the figure contradicted itself. "compatible
    # with at most" is what a 95% upper bound licenses.
    labs(title="Pooled across three colorectal cohorts, the data are compatible with at most a 9% per-SD hazard increase",
         subtitle=sprintf(paste0("EXPLORATORY, POST-HOC. Score log-HR per 1 SD, model 1, ",
                                 "%s patients and %d events.\nPooled HR %.3f (95%% CI %.3f to %.3f), ",
                                 "I2 = %.1f%%, tau2 = %.3g, Q p = %.3f."),
                          format(po$n_total, big.mark=","), po$events_total,
                          po$HR, po$HR_lo_wald, po$HR_hi_wald, po$I2, po$tau2, po$Q_p),
         x="Hazard ratio per 1 SD of score  (right = worse survival)", y=NULL,
         caption=paste0("EXPLORATORY, POST-HOC -- NOT the registered discovery meta-analysis, which pools six cohorts to HR 1.129.\n",
                        "GSE39582 is an EXTERNAL validation cohort; Amendment 16 states validation cohorts are not meta-analysed\n",
                        "with discovery. The three inputs are not exchangeable: COAD and READ are RNA-seq, GSE39582 is GPL570\n",
                        "microarray; READ's registered endpoint is PFI (Amendment 7), the others OS. Right margin: events / n.\n",
                        "Wald interval quoted (at k = 3 with tau2 = 0 the Hartung-Knapp interval is narrower). The ceiling is the\n",
                        "largest per-SD HR compatible with these data at alpha = 0.05, not a proof of absence.")) +
    theme_pub(grid_y=FALSE, grid_x=TRUE) +
    theme(legend.position="bottom")
  list(plot=p, file="figure_crc_pool.png", w=180, h=110)
}

# ================================================== new: RPPA concordance grid
fig_rppa <- function() {
  cr <- rd("exploratory_rppa_control_correlations.csv")
  po <- rd("exploratory_rppa_control_pooled.csv")
  pe <- cr[cr$method=="pearson", ]
  assert_plot("figRPPA", pe, "exploratory_rppa_control_correlations.csv",
              c("cohort","comparison","method"), c("n","r","ci_lo","ci_hi"))
  assert_plot("figRPPA-pool", po, "exploratory_rppa_control_pooled.csv", "comparison",
              c("r_pooled","ci_lo_wald","ci_hi_wald","I2"))
  pe$short <- sub("^TCGA-","",pe$cohort)
  keep <- c("1. total STAT3 RPPA vs STAT3 mRNA","2. score_140 vs STAT3 mRNA",
            "3. total STAT3 RPPA vs STAT3_pY705 RPPA","4. score_140 vs STAT3_pY705 RPPA (from 13)")
  pe <- pe[pe$comparison %in% keep, ]
  lab <- c("1. total STAT3 protein\nvs STAT3 mRNA\n(cross-assay control)",
           "2. score vs STAT3 mRNA\n(panel coherence)",
           "3. total STAT3 protein vs pY705\n(SAME lysate, same array)",
           "4. score vs pY705\n(the quantity under test)")
  names(lab) <- keep
  # po carries a 2b (leave-STAT3-out) row that is NOT one of the four panels.
  # The first draft built the strip labels before filtering, so 2b became an
  # NA-labelled fifth facet with a spurious pooled line drawn in it.
  po <- po[po$comparison %in% keep, ]
  if (nrow(po) != 4L) halt("figRPPA", "expected 4 pooled rows, got ", nrow(po))
  stat_lab <- setNames(sprintf("%s\npooled r = %.3f (%.3f to %.3f), I2 %.0f%%",
                               lab[match(po$comparison, keep)], po$r_pooled,
                               po$ci_lo_wald, po$ci_hi_wald, po$I2), po$comparison)
  if (anyNA(stat_lab)) halt("figRPPA", "a strip label is NA")
  pe$facet <- factor(stat_lab[pe$comparison], levels=stat_lab[keep])
  po$facet <- factor(stat_lab[po$comparison], levels=stat_lab[keep])
  if (anyNA(pe$facet) || anyNA(po$facet)) halt("figRPPA", "a row has no facet")
  if (nlevels(pe$facet) != 4L) halt("figRPPA", "expected 4 facets, got ", nlevels(pe$facet))
  pe$short <- factor(pe$short, levels=rev(c("COAD","READ","STAD","ESCA","PAAD","LIHC")))
  # comparison 3 is the reference-standard limitation: it alone gets the accent
  pe$focal <- pe$comparison == keep[3]; po$focal <- po$comparison == keep[3]

  p <- ggplot(pe, aes(r, short)) +
    geom_vline(xintercept=0, colour="grey60", linewidth=0.3) +
    geom_rect(data=po, inherit.aes=FALSE,
              aes(xmin=ci_lo_wald, xmax=ci_hi_wald, ymin=-Inf, ymax=Inf, fill=focal),
              alpha=0.13) +
    geom_vline(data=po, aes(xintercept=r_pooled, colour=focal), linewidth=0.6) +
    geom_errorbar(aes(xmin=ci_lo, xmax=ci_hi, colour=focal), orientation="y",
                  width=0, linewidth=0.45) +
    geom_point(aes(colour=focal), size=1.7) +
    # the pooled statistic used to sit inside the panel and clipped at the right
    # on every facet; it now rides in the strip label, which has full width.
    geom_point(data=po, inherit.aes=FALSE,
               aes(x=r_pooled, y=Inf, colour=focal), shape=18, size=2.2,
               show.legend=FALSE) +
    facet_wrap(~ facet, nrow=2) +
    scale_colour_manual(values=c("FALSE"=DISC, "TRUE"=ACCENT), guide="none") +
    scale_fill_manual(values=c("FALSE"="grey55", "TRUE"=ACCENT), guide="none") +
    # three widely-spaced ticks: five 5-character labels collided at 40 mm/facet
    scale_x_continuous(limits=c(-0.45,1.0), breaks=c(-0.25,0,0.25,0.5,0.75)) +
    scale_y_discrete(expand=expansion(add=c(0.55,0.9))) +
    # NOT "disagrees with": the pooled r is 0.139 with a Wald interval of -0.064
    # to 0.331, which SPANS ZERO. An interval spanning zero cannot establish
    # disagreement; it fails to establish agreement. The title now says what the
    # interval supports.
    labs(title="The RPPA phosphosite shows no established agreement with total STAT3 on the same lysate",
         subtitle=paste0("EXPLORATORY, POST-HOC. Pearson r with 95% intervals per discovery cohort; band, line and\n",
                         "diamond are the pooled random-effects estimate. Every cohort is shown because I2 is 60-92%."),
         x="Pearson correlation", y=NULL,
         caption=paste0("EXPLORATORY, POST-HOC. Comparison 3 (accent) is the reference-standard limitation: both antibodies\n",
                        "come from the SAME lysate on the same array, normalised together, so shared loading INFLATES it -- and\n",
                        "it is still only 0.139, with I2 92% and an interval spanning zero. Comparison 1 is the only genuinely\n",
                        "cross-assay test. Comparison 2 is partly a self-correlation: STAT3 is 1 of the 140 scoring genes (0.71%\n",
                        "of the mean), and dropping it moves the pooled r by 0.008. Comparisons 1 and 3 rest on 19 fewer\n",
                        "patients -- total STAT3 is absent from 19 of the 1,278 primary-tumour files.")) +
    theme_pub(grid_y=FALSE, grid_x=TRUE) +
    theme(panel.spacing.x=unit(3,"mm"),
          panel.spacing.y=unit(5,"mm"),
          strip.text=element_text(size=PT_TICK, lineheight=1.15, hjust=0))
  list(plot=p, file="figure_rppa_concordance.png", w=180, h=150)
}

# ============================== compartment attribution: two forms to compare
# Shared data prep, so the two renderings cannot diverge in what they plot.
compartment_data <- function() {
  o <- rd("origin_six_compartment.csv")
  assert_plot("fig3-data", o, "origin_six_compartment.csv", c("gene","atlas"),
              c("f_at_0.30","f_at_0.50","f_at_0.70"))
  long <- do.call(rbind, lapply(c("0.30","0.50","0.70"), function(g)
    data.frame(gene=o$gene, atlas=o$atlas, qualifying=o$qualifying,
               grid=paste0("pi = ", g), f=o[[paste0("f_at_", g)]],
               dominant=o$dominant, stringsAsFactors=FALSE)))
  long <- long[!is.na(long$f), ]
  long$lab <- ifelse(long$qualifying, long$gene, paste0(long$gene, " *"))
  ord <- c("SOCS3","MYC","IL6","BCL2 *","MMP9 *","HGF *")
  long$lab <- factor(long$lab, levels=ord)
  long$atlas <- factor(long$atlas, levels=c("GSE125449","GSE178341","Peng"),
                       labels=c("GSE125449 liver/biliary","GSE178341 colorectal","Peng pancreatic"))
  list(long=long, n_dom=sum(o$dominant, na.rm=TRUE))
}

CAP3 <- paste0("The mark shows which atlas-gene combinations the compartment analysis calls epithelial-dominant:\n",
               "f > 50% at EVERY point of the band, not at one point. Nine of the 54 values clear 50% somewhere -- MYC\n",
               "in all three atlases and BCL2 in all three at pi = 0.70 -- but only MYC in GSE178341 clears it throughout.\n",
               "BCL2, MMP9 and HGF (*) failed criterion B (Amendment 2, human ChIP-seq only) and enter no k variant;\n",
               "their fractions are reported as a labelled non-qualifying subset (prespecification section 4).")
TITLE3 <- paste0("MYC is the only gene epithelial-dominant across the whole\n30-70% band, and only in colorectal tissue")

fig3_bars <- function() {
  cd <- compartment_data(); long <- cd$long
  # rename BEFORE deriving st and dom: renaming only one of them put the
  # dominance marks in a spurious fourth facet (caught in the render).
  levels(long$atlas) <- sub(" ", "\n", levels(long$atlas))
  st <- rbind(transform(long, part="Epithelial", value=f),
              transform(long, part="Other compartments", value=1-f))
  st$part <- factor(st$part, levels=c("Other compartments","Epithelial"))
  dom <- unique(long[long$dominant & !is.na(long$dominant), c("lab","atlas")])
  dom <- do.call(rbind, lapply(unique(long$grid), function(g) transform(dom, grid=g)))
  if (nrow(unique(dom[,c("lab","atlas")])) != cd$n_dom)
    halt("fig3-bars", "dominance marks do not match the file's `dominant` column")
  if (!all(levels(droplevels(dom$atlas)) %in% levels(st$atlas)) ||
      nlevels(st$atlas) != 3L)
    halt("fig3-bars", "the dominance frame does not share the bar frame's atlas levels")
  p <- ggplot(st, aes(lab, value, fill=part)) +
    geom_col(width=0.68) +
    geom_hline(yintercept=0.5, colour=MUTED, linewidth=0.35, linetype="22") +
    geom_point(data=dom, aes(lab, y=1.06), inherit.aes=FALSE, shape=18, size=1.9, colour=ACCENT) +
    facet_grid(grid ~ atlas) +
    scale_fill_manual(values=c("Epithelial"=DISC, "Other compartments"="grey86")) +
    scale_y_continuous(labels=function(x) paste0(100*x,"%"), breaks=c(0,.5,1),
                       limits=c(0,1.12), expand=expansion(mult=c(0,0.02))) +
    labs(title=TITLE3,
         subtitle=paste0("Share of panel-gene expression attributable to the epithelial compartment, at three points\n",
                         "of the registered 30-70% purity band. Diamond above a bar marks a dominance call."),
         x=NULL, y="Epithelial share of expression", caption=CAP3) +
    theme_pub(grid_y=TRUE) +
    theme(legend.position="bottom", panel.grid.major.x=element_blank(),
          axis.text.x=element_text(angle=45, hjust=1, size=PT_TICK))
  list(plot=p, file="fig3_compartment_bars.png", w=180, h=160)
}

fig3_heat <- function() {
  cd <- compartment_data(); long <- cd$long
  dom <- long[long$dominant & !is.na(long$dominant), ]
  if (nrow(unique(dom[,c("lab","atlas")])) != cd$n_dom)
    halt("fig3-heat", "dominance marks do not match the file's `dominant` column")
  # 6.5: 54 cells, so print the value in every cell.
  # 4.4: diverging map centred on the SEMANTIC zero -- the 50% dominance
  # threshold -- not the data midpoint.
  p <- ggplot(long, aes(atlas, lab, fill=f)) +
    geom_tile(colour="white", linewidth=0.6) +
    geom_text(aes(label=sprintf("%.0f", 100*f),
                  colour=abs(f-0.5) > 0.32), size=PT_TICK/.pt, show.legend=FALSE) +
    # A diamond centred on the cell printed ON TOP of the value it marks (caught
    # in the render). A cell outline conveys the same call and occludes nothing.
    geom_tile(data=dom, aes(atlas, lab), inherit.aes=FALSE, fill=NA,
              colour=ACCENT, linewidth=1.1) +
    facet_wrap(~ grid, nrow=1) +
    scale_fill_gradient2(low="#F0F0F0", mid="#C9DDEE", high=DISC, midpoint=0.5,
                         limits=c(0,1), labels=function(x) paste0(100*x,"%"),
                         name="Epithelial share") +
    scale_colour_manual(values=c("FALSE"=INK, "TRUE"="white")) +
    labs(title=TITLE3,
         subtitle=paste0("Epithelial share of expression (%) per gene per atlas, at three points of the registered 30-70%\n",
                         "purity band. Colour is centred on the 50% dominance threshold; a boxed cell is a dominance call."),
         x=NULL, y=NULL, caption=CAP3) +
    theme_pub(grid_y=FALSE) +
    theme(legend.position="bottom", legend.key.width=unit(10,"mm"),
          # full atlas names collided at three-across; angled and single-line
          axis.text.x=element_text(size=PT_TICK, angle=35, hjust=1),
          # the leftmost angled label overran the canvas at margin 0
          plot.margin=margin(5, 6, 5, 16),
          panel.spacing.x=unit(3,"mm"))
  list(plot=p, file="fig3_compartment_heatmap.png", w=180, h=150)
}

fig4_atten <- function() {
  at <- rd("attenuation_per_cohort.csv"); ma <- rd("meta_analysis.csv")
  assert_plot("fig4", at, "attenuation_per_cohort.csv", "cohort",
              c("attenuation_total","att_total_lo","att_total_hi"))
  pool <- ma[ma$analysis=="attenuation_total", ]
  assert_plot("fig4-pool", pool, "meta_analysis.csv", "analysis",
              c("est","ci_lo","ci_hi","tau2","I2"))
  D <- rbind(
    data.frame(short=sub("^TCGA-","",at$cohort), est=at$attenuation_total,
               lo=at$att_total_lo, hi=at$att_total_hi, class="Cohort (paired bootstrap)",
               stringsAsFactors=FALSE),
    data.frame(short="Pooled (6 cohorts)", est=pool$est, lo=pool$ci_lo, hi=pool$ci_hi,
               class="Pooled estimate", stringsAsFactors=FALSE))
  D$short <- factor(D$short, levels=rev(c("COAD","READ","STAD","ESCA","PAAD","LIHC","Pooled (6 cohorts)")))
  ggplot(D, aes(est, short, colour=class, shape=class)) +
    geom_vline(xintercept=0, colour="grey60", linewidth=0.3) +
    geom_errorbar(aes(xmin=lo, xmax=hi), orientation="y", width=0, linewidth=0.45) +
    geom_point(size=1.9) +
    # discrete y: anchor to the factor level by name, not by a numeric index
    geom_text(data=data.frame(x=pool$est, y="Pooled (6 cohorts)",
                              l=sprintf("pooled %.4f  (%.4f to %.4f)",
                                        pool$est, pool$ci_lo, pool$ci_hi)),
              mapping=aes(x=x, y=y, label=l), vjust=2.1, hjust=0.5,
              size=PT_ANNOT/.pt, colour=ACCENT, inherit.aes=FALSE) +
    scale_y_discrete(expand=expansion(add=c(0.95, 0.6))) +
    scale_colour_manual(values=c("Cohort (paired bootstrap)"=DISC, "Pooled estimate"=ACCENT)) +
    scale_shape_manual(values=c("Cohort (paired bootstrap)"=16, "Pooled estimate"=18)) +
    labs(title="Adjustment for purity and stroma shows no reduction in the score's association",
         subtitle=paste0("attenuation_total = beta(M2) - beta(M4). The registered direction is POSITIVE -- adjustment\n",
                         "reduces the association. The pooled estimate is negative."),
         x="attenuation_total  (log-HR units; positive = adjustment reduces the association)", y=NULL,
         caption=sprintf(paste0("Percentile intervals from 2,000 paired bootstrap resamples. tau2 = %.4f, I2 = %.1f%%.\n",
                                "CHOL is not shown: models 2-4 are not fitted there (EPV < 5, Amendment 8). The interval\n",
                                "contains reductions larger than the point estimate."), pool$tau2, pool$I2)) +
    theme_pub(grid_y=FALSE, grid_x=TRUE) + theme(legend.position="bottom")
}

fig5_cms <- function() {
  cc <- rd("cms_calls.csv"); sc <- rd("scores_per_patient.csv"); dist <- rd("cms_distribution.csv")
  ct <- rd("cms_tertile_crosstab.csv"); ca <- rd("cms4_association.csv")
  d <- merge(cc[,c("cohort","barcode","cms")], sc[,c("barcode","score")], by="barcode")
  if (nrow(d) != nrow(cc)) halt("fig5", "CMS calls did not join 1:1 to scores")
  got <- as.data.frame(table(d$cohort, d$cms), stringsAsFactors=FALSE)
  names(got) <- c("cohort","cms","n_plotted")
  chk <- merge(got, dist[,c("cohort","cms","n")], by=c("cohort","cms"))
  if (nrow(chk) != nrow(dist) || any(chk$n_plotted != chk$n))
    halt("fig5", "per-CMS counts differ from cms_distribution.csv")
  message(sprintf("  ok  fig5           %d group counts match cms_distribution.csv", nrow(chk)))
  d$cohort <- sub("^TCGA-","",d$cohort)
  d$cms <- factor(d$cms, levels=c("CMS1","CMS2","CMS3","CMS4","unclassified"))
  lab <- got; lab$cohort <- sub("^TCGA-","",lab$cohort)
  lab$cms <- factor(lab$cms, levels=levels(d$cms))
  # chi2/df/p and R2 are interpolated into the CAPTION, so they are asserted too
  # the test statistic is repeated on all 12 tertile x cms rows of each cohort,
  # so `cohort` is not a unique key -- assert on the full key, then collapse only
  # after confirming the statistic really is constant within cohort.
  assert_plot("fig5-chisq", ct, "cms_tertile_crosstab.csv", c("cohort","tertile","cms"),
              c("chisq","df","p_asymptotic","n"))
  nu <- vapply(split(ct$chisq, ct$cohort), function(x) length(unique(x)), integer(1))
  if (any(nu != 1L)) halt("fig5", "chi-square is not constant within cohort")
  u  <- unique(ct[,c("cohort","chisq","df","p_asymptotic")]); u$short <- sub("^TCGA-","",u$cohort)
  r2df <- ca[match(u$cohort, ca$cohort), ]
  assert_plot("fig5-r2", r2df, "cms4_association.csv", "cohort", "r2_score_on_cms")
  r2 <- r2df$r2_score_on_cms
  ggplot(d, aes(cms, score)) +
    geom_hline(yintercept=0, colour="grey80", linewidth=0.3) +
    geom_violin(aes(fill=cms=="CMS4"), colour=NA, alpha=0.55, width=0.9) +
    geom_boxplot(width=0.15, outlier.shape=NA, linewidth=0.3, fill="white", colour=INK) +
    geom_jitter(width=0.13, height=0, size=0.3, alpha=0.32, colour=MUTED) +
    geom_text(data=lab, aes(cms, y=Inf, label=paste0("n=", n_plotted)), vjust=1.4,
              size=PT_TICK/.pt, colour=MUTED, inherit.aes=FALSE) +
    facet_wrap(~ cohort, nrow=1) +
    scale_fill_manual(values=c("TRUE"=ACCENT, "FALSE"="grey80"), guide="none") +
    scale_y_continuous(expand=expansion(mult=c(0.04,0.14))) +
    labs(title="The score is highest in CMS4, the mesenchymal subtype",
         subtitle=paste0("Distribution of the 140-gene score by consensus molecular subtype (CMScaller 2.0.1).\n",
                         "Unclassified samples are retained as a level, never dropped (B.o)."),
         x=NULL, y="STAT3 score (SD units, within cohort)",
         caption=paste0("Tertile x CMS1-4: ",
                        paste(sprintf("%s chi2 = %.2f (df %d, p = %.1e)", u$short, u$chisq,
                                      as.integer(u$df), u$p_asymptotic), collapse=";\n"),
                        ".\nR2 of score on CMS: ", paste(sprintf("%.3f", r2), collapse=" / "),
                        ". Only COAD and READ are CMS-classifiable (B.o.1).")) +
    theme_pub(grid_y=TRUE) +
    theme(panel.grid.major.x=element_blank(),
          axis.text.x=element_text(angle=30, hjust=1, size=PT_TICK))
}

# -------------------------------------------------------------------- driver
if (sys.nframe() == 0L) {
  message("\n== refined figures -> figures/refined/ (presentation only) ==")
  if (!file.exists("11b_figures_refined.R")) halt("driver", "run from the repository root")
  dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)
  verify_assert_can_fail()
  export_null_draws()

  # ---- the new main figure: A null, B forest (M1 only), C FU-iCCA scatter ----
  A <- panel_null(); B <- panel_forest("M1", compact=TRUE); C <- panel_fuicca()
  # Panel C spans the full width: as a half-width panel it left the bottom row
  # half empty (figure-style 3.5, fill the box). Panel tags are added with a
  # left inset so they cannot collide with the y-axis title, which happened to
  # "C" in the first render.
  top <- plot_grid(A, NULL, B, nrow = 1, rel_widths = c(1, 0.05, 1.05),
                   labels = c("A", "", "B"), label_size = 12, label_fontface = "bold",
                   label_x = 0.005, hjust = 0, vjust = 1.02)
  # "C" is drawn in its own thin strip ABOVE the panel: as a plot_grid label it
  # landed on top of the y-axis title (caught in the render).
  tagC <- ggdraw() + draw_label("C", fontface = "bold", size = 12, colour = INK,
                                hjust = 0, x = 0.008, y = 0.30)
  bot <- plot_grid(tagC, C + theme(plot.margin = margin(2, 4, 2, 6)),
                   ncol = 1, rel_heights = c(0.07, 0.93))
  # title/subtitle wrapped to the 180 mm width -- both clipped at the right edge
  # in the first render.
  # The title block was a separate ggdraw strip and kept clipping at the right
  # and colliding with the panels below. Drawn as a titled wrapper plot instead,
  # so ggplot's own layout reserves the space and wraps nothing off-canvas.
  ttl <- ggdraw() + theme(plot.margin = margin(0,0,0,0)) +
    draw_label(paste0("A preregistered STAT3 score tracks its target yet performs at the\n",
                      "median of matched random signatures"),
               fontface = "bold", size = PT_TITLE, colour = INK, hjust = 0, vjust = 1,
               x = 0.008, y = 0.95, lineheight = 1.2) +
    draw_label(paste0("A: pooled score log-HR (model 1) against 10,000 matched null signatures.\n",
                      "B: per-cohort and pooled log-HR per 1 SD; GSE39582 is external and never pooled.\n",
                      "C: the score vs measured STAT3 pY705 in FU-iCCA. Model 4 is in the supplement."),
               size = PT_SUB, colour = MUTED, hjust = 0, vjust = 1,
               x = 0.008, y = 0.52, lineheight = 1.35)
  main <- plot_grid(ttl, top, bot, ncol = 1, rel_heights = c(0.13, 0.45, 0.42))
  save_fig(main, "figure_main.png", 180, 225)

  # ---- the two new figures ----
  for (f in list(fig_crc(), fig_rppa())) save_fig(f$plot, f$file, f$w, f$h)

  # ---- refined re-renders of the existing set ----
  save_fig(panel_null(standalone=TRUE),   "fig1_null_distribution.png", 180, 105)
  save_fig(panel_forest(c("M1","M4"), standalone=TRUE), "fig2_forest_logHR.png", 180, 125)
  for (f in list(fig3_bars(), fig3_heat())) save_fig(f$plot, f$file, f$w, f$h)
  save_fig(fig4_atten(), "fig4_attenuation_forest.png", 180, 110)
  save_fig(fig5_cms(),   "fig5_score_by_cms.png", 180, 115)
  save_fig(panel_fuicca(standalone=TRUE), "fig6_fuicca_y705_scatter.png", 180, 120)
  # supplement: the M4 panel dropped from the main figure
  save_fig(panel_forest("M4", standalone=TRUE), "figS_forest_M4.png", 180, 118)

  message("\n", length(list.files(FIGDIR, pattern="\\.png$")), " figures in ", FIGDIR)
  message("11_figures.R untouched; figures/ untouched. No reported number changed.")
}
