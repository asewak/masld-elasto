################################################################################
# SWE-based score: prespecified covariates, logistic regression (LR) as
# placeholder (mmlt to be added); bootstrap optimism
# Needs d1 and res from 02_data_preparation.R; adds the sc_* objects to res
################################################################################

# AST and ALT are right-skewed: log scale
d1$lalt <- log(d1$alt)
d1$last <- log(d1$ast)

covs <- c("swe", "lalt", "last", "plt", "bmi", "age")
if (!exists("n_boot")) n_boot <- 1000  # set smaller before sourcing to test

# Complete cases on the covariates
d2 <- d1[complete.cases(d1[covs]), ]
d2$D <- factor(d2$fmash, levels = 0:1)
y_sc <- d2$fmash

# AUC via ranks (ties count 1/2), score higher in disease
auc_rank <- function(lp, y) {
  r <- rank(lp); n1 <- sum(y == 1); n0 <- sum(y == 0)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# LR: linear predictor as score
fit_lr <- function(d) glm(reformulate(covs, "D"), data = d, family = binomial)
score_lr <- function(f, d) as.numeric(predict(f, newdata = d, type = "link"))

# mmlt: Box-Cox type marginal model per covariate with shift and scale
# depending on D, joint Gaussian dependence between covariates (covariance
# does not depend on D). Score = log-likelihood ratio log f(x | D=1) / f(x | D=0).
fit_mm <- function(d) {
  ms <- lapply(paste(covs, "~ D | D"), BoxCox, data = d, scale_shift = TRUE,
               prob = c(0, 1))
  do.call(mmlt, c(ms, data = list(d)))
}
score_mm <- function(f, d) {
  d0 <- d1_ <- d
  d0$D <- factor(0, levels = 0:1)
  d1_$D <- factor(1, levels = 0:1)
  predict(f, newdata = d1_, type = "density", log = TRUE) -
    predict(f, newdata = d0, type = "density", log = TRUE)
}

# Bootstrap currently for LR only (mmlt bootstrap to be run on the cluster)
models <- list(
  lr = list(fit = fit_lr, score = score_lr)
  # , mmlt = list(fit = fit_mm, score = score_mm)
)

# Apparent fits
fits <- lapply(models, function(m) m$fit(d2))
lp_app <- mapply(function(m, f) m$score(f, d2), models, fits, SIMPLIFY = FALSE)
auc_app <- sapply(lp_app, auc_rank, y = y_sc)
auc_ci <- sapply(lp_app, function(lp)
  as.numeric(ci.auc(roc(y_sc, lp, direction = "<", quiet = TRUE), method = "delong")))

# Bootstrap: refit each model in each resample (same resamples for all models);
# optimism = AUC in the resample minus AUC of that fit in the original data
RNGkind("L'Ecuyer-CMRG")
set.seed(2026)
boot_one <- function(b) {
  i <- sample(nrow(d2), replace = TRUE)
  out <- sapply(models, function(m) {
    r <- try({
      f <- m$fit(d2[i, ])
      c(app = auc_rank(m$score(f, d2[i, ]), y_sc[i]),
        test = auc_rank(m$score(f, d2), y_sc))
    }, silent = TRUE)
    if (inherits(r, "try-error")) c(app = NA_real_, test = NA_real_) else r
  })
  out  # rows app/test, columns models
}
boot <- mclapply(seq_len(n_boot), boot_one, mc.cores = 8, mc.set.seed = TRUE)
stopifnot(all(sapply(boot, is.matrix)))
boot <- simplify2array(boot)  # app/test x model x resample
# Failed fits or non-finite scores (NA AUC) are dropped, counted per model
boot_ok <- apply(boot, c(2, 3), function(z) all(is.finite(z)))
n_boot_ok <- rowSums(boot_ok)
optimism <- sapply(colnames(boot), function(m) {
  ok <- boot_ok[m, ]
  mean(boot["app", m, ok] - boot["test", m, ok])
})

# Coefficients of the LR score (original scale)
score_tab <- data.frame(variable = names(coef(fits$lr))[-1], coef = coef(fits$lr)[-1],
                        or = exp(coef(fits$lr)[-1]))

# Apparent AUC (DeLong 95% CI) of single biomarkers, orientation chosen
# automatically (AUC >= 0.5); each on its available cases. FIB-4 and FAST are
# established scores shown for comparison.
single <- c(covs, "fib4", "fast")
single_tab <- do.call(rbind, lapply(single, function(v) {
  ok <- !is.na(d1[[v]])
  r <- roc(d1$fmash[ok], d1[[v]][ok], direction = "auto", quiet = TRUE)
  ci <- as.numeric(ci.auc(r, method = "delong"))
  data.frame(variable = v, n = sum(ok), auc = ci[2], lo = ci[1], hi = ci[3],
             higher_in_disease = r$direction == "<")
}))

################################################################################
# Apparent cut-offs: 90% sensitivity, 90% specificity, Youden optimum
# Score >= cut-off is test positive. Candidate cut-offs are the observed values.
################################################################################

# mmlt: apparent fit (also used for the coefficient table below)
fit_mmlt <- fit_mm(d2)

# FAST: Excel column (= logistic of the code.R formula, checked identical)
cut_scores <- list(
  "LR score (probability)" = plogis(lp_app$lr),
  "mmlt (log LR)"          = score_mm(fit_mmlt, d2),
  "FAST"                   = d2$fast,
  "SWE"                    = d2$swe,
  "Fibroscan (TE)"         = d2$fst
)

# Accuracy measures at cut-off ct
perf_at <- function(s, y, ct) {
  tp <- sum(s >= ct & y == 1); fp <- sum(s >= ct & y == 0)
  fn <- sum(s <  ct & y == 1); tn <- sum(s <  ct & y == 0)
  se <- tp / (tp + fn); sp <- tn / (tn + fp)
  c(cutoff = ct, sens = se, spec = sp, ppv = tp / (tp + fp), npv = tn / (tn + fn),
    acc = (tp + tn) / length(y), youden = se + sp - 1,
    lr_pos = se / (1 - sp), lr_neg = (1 - se) / sp)
}

cut_one <- function(s, y) {
  cand <- sort(unique(s))
  se <- sapply(cand, function(ct) mean(s[y == 1] >= ct))
  sp <- sapply(cand, function(ct) mean(s[y == 0] <  ct))
  # highest cut-off with sens >= 90%; lowest cut-off with spec >= 90%;
  # Youden: first maximum (lowest cut-off) if tied
  ct <- c("Sens >= 90%" = max(cand[se >= 0.9]),
          "Spec >= 90%" = min(cand[sp >= 0.9]),
          "Youden"      = cand[which.max(se + sp - 1)])
  t(sapply(ct, function(z) perf_at(s, y, z)))
}

cut_tab <- do.call(rbind, lapply(names(cut_scores), function(m)
  data.frame(method = m, criterion = c("Sens >= 90%", "Spec >= 90%", "Youden"),
             cut_one(cut_scores[[m]], y_sc), row.names = NULL)))

# FAST: published cut-offs replace the data-driven ones (rule-out < 0.35,
# rule-in >= 0.67); the Youden row stays data-driven
fast_cut <- c("Sens >= 90%" = 0.35, "Spec >= 90%" = 0.67)
for (k in names(fast_cut))
  cut_tab[cut_tab$method == "FAST" & cut_tab$criterion == k, -(1:2)] <-
    perf_at(cut_scores[["FAST"]], y_sc, fast_cut[[k]])

# Apparent AUC with DeLong 95% CI (all scores are higher in disease)
auc_tab <- do.call(rbind, lapply(names(cut_scores), function(m) {
  ci <- as.numeric(ci.auc(roc(y_sc, cut_scores[[m]], direction = "<", quiet = TRUE),
                          method = "delong"))
  data.frame(method = m, auc = ci[2], lo = ci[1], hi = ci[3])
}))

# Rule-out: score < cut-off for 90% sens; rule-in: score >= cut-off for 90% spec
zone_tab <- do.call(rbind, lapply(names(cut_scores), function(m) {
  s <- cut_scores[[m]]
  lo <- cut_tab$cutoff[cut_tab$method == m & cut_tab$criterion == "Sens >= 90%"]
  hi <- cut_tab$cutoff[cut_tab$method == m & cut_tab$criterion == "Spec >= 90%"]
  data.frame(method = m, cut_out = lo, cut_in = hi,
             rule_out = mean(s < lo), grey = mean(s >= lo & s < hi),
             rule_in = mean(s >= hi))
}))

################################################################################
# mmlt coefficients: shift and scale of D=1 per covariate, correlation matrix
# (Bernstein basis coefficients of the marginal transformations are not shown)
################################################################################

mm_cf <- coef(fit_mmlt)
mm_se <- sqrt(diag(vcov(fit_mmlt)))
stopifnot(all(names(mm_cf) %in% names(mm_se)))
mm_se <- mm_se[names(mm_cf)]

mm_marg <- do.call(rbind, lapply(covs, function(v) {
  k <- paste0(v, c(".D1", ".scl_D1"))
  data.frame(variable = v,
             shift = mm_cf[k[1]], shift_se = mm_se[k[1]],
             scale = mm_cf[k[2]], scale_se = mm_se[k[2]], row.names = NULL)
}))
mm_marg$shift_p <- 2 * pnorm(-abs(mm_marg$shift / mm_marg$shift_se))
mm_marg$scale_p <- 2 * pnorm(-abs(mm_marg$scale / mm_marg$scale_se))

# Correlation matrix of the joint Gaussian model (on the transformed scale)
# coef() gives the packed lower triangle, named "<row>.<column>"
cr <- unclass(drop(coef(fit_mmlt, type = "Corr")))
stopifnot(length(cr) == length(covs) * (length(covs) - 1) / 2)
mm_corr <- diag(length(covs)); dimnames(mm_corr) <- list(covs, covs)
for (k in names(cr)) {
  rc <- strsplit(k, ".", fixed = TRUE)[[1]]
  mm_corr[rc[1], rc[2]] <- mm_corr[rc[2], rc[1]] <- cr[[k]]
}

################################################################################
# Empirical ROC curves (base R)
################################################################################

# Points at all observed cut-offs, score >= cut-off is test positive
roc_xy <- function(s, y) {
  ct <- c(Inf, sort(unique(s), decreasing = TRUE))
  data.frame(fpr = sapply(ct, function(z) mean(s[y == 0] >= z)),
             tpr = sapply(ct, function(z) mean(s[y == 1] >= z)))
}

# scores: named list of scores (higher in disease); AUC in the legend
roc_plot <- function(scores, y, cols, main = "") {
  plot(0:1, 0:1, type = "n", asp = 1, xlab = "1 - Specificity", ylab = "Sensitivity",
       main = main)
  abline(0, 1, col = "grey70", lty = 2)
  for (i in seq_along(scores)) {
    r <- roc_xy(scores[[i]], y)
    lines(r$fpr, r$tpr, type = "s", col = cols[i], lwd = 2)  # step function (as geom_step)
  }
  legend("bottomright", col = cols, lwd = 2, bty = "n", cex = 0.8,
         legend = sprintf("%s (AUC %.2f)", names(scores),
                          sapply(scores, auc_rank, y = y)))
}

pdf(paths$f("roc_fmash.pdf"), width = 5, height = 5)
roc_plot(list("SWE-based score (LR)"   = cut_scores[[1]],
              "SWE-based score (mmlt)" = cut_scores[[2]],
              "FAST"                   = cut_scores[["FAST"]]),
         y_sc, cols = c("#1b9e77", "#d95f02", "#7570b3"))
dev.off()

################################################################################
# Objects for the report
################################################################################

res <- c(res, list(
  sc_cut = cut_tab,
  sc_auc_cmp = auc_tab,
  sc_zone = zone_tab,
  sc_n = nrow(d2),
  sc_events = sum(y_sc),
  sc_ids_miss = d1$id[!complete.cases(d1[covs])],
  # per model (lr; mmlt once enabled)
  sc_auc_app = auc_app,
  sc_auc_ci = auc_ci[c(1, 3), , drop = FALSE],
  sc_optimism = optimism,
  sc_auc_corr = auc_app - optimism,
  sc_boot_ok = n_boot_ok,
  sc_n_boot = n_boot,
  sc_tab = score_tab,
  sc_lr_fit = fits$lr,
  sc_mm_marg = mm_marg,
  sc_mm_corr = mm_corr,
  sc_single = single_tab
))
