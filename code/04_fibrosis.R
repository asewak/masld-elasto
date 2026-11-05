################################################################################
# Liver fibrosis (independent of NAS): SWE and Fibroscan for >= F2, >= F3, F4
# Needs dat, res (02), perf_at and cut_one (03); adds fb_* objects to res
################################################################################

# Patients with a fibrosis stage and both elastography values (same patients
# for both methods, so the AUCs are comparable)
d3 <- dat[!is.na(dat$fib), ]
fb_miss <- d3$id[is.na(d3$swe) | is.na(d3$fst)]
d3 <- d3[!is.na(d3$swe) & !is.na(d3$fst), ]

fib_out <- list(
  "F0/F1 vs F2-F4" = "F2",
  "F0-F2 vs F3-F4" = "F3",
  "F0-F3 vs F4"    = "F4"
)
fb_scores <- list("SWE" = d3$swe, "Fibroscan (TE)" = d3$fst)

pdf(paths$f("roc_fibrosis.pdf"), width = 10, height = 3.8)
par(mfrow = c(1, 3), mar = c(4, 4, 2, 1))
fb_res <- lapply(fib_out, function(st) {
  y <- as.integer(as.integer(d3$fib) >= match(st, levels(d3$fib)))
  roc_plot(fb_scores, y, cols = c("#0072B2", "#D55E00"), main = names(fib_out)[fib_out == st])
  cut <- do.call(rbind, lapply(names(fb_scores), function(m)
    data.frame(method = m, criterion = c("Sens >= 90%", "Spec >= 90%", "Youden"),
               cut_one(fb_scores[[m]], y), row.names = NULL)))
  auc <- do.call(rbind, lapply(names(fb_scores), function(m) {
    ci <- as.numeric(ci.auc(roc(y, fb_scores[[m]], direction = "<", quiet = TRUE),
                            method = "delong"))
    data.frame(method = m, auc = ci[2], lo = ci[1], hi = ci[3])
  }))
  list(cut = cut, auc = auc, n = length(y), events = sum(y))
})
dev.off()

res <- c(res, list(
  fb_res = fb_res,
  fb_n = nrow(d3),
  fb_ids_miss = fb_miss,
  fb_stage_n = table(d3$fib)
))
