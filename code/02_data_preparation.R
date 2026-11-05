################################################################################
# Data preparation
# Input:  raw (from 01_import.R)
# Output: dat (one row per patient), res (objects used by the report)
################################################################################

# Excel text such as "k.A." becomes NA
to_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

fib_map <- c(
  "0"  = "F0",
  "1a" = "F1", "1" = "F1", "1b" = "F1", "1B" = "F1", "1c" = "F1",
  "2"  = "F2",
  "3"  = "F3",
  "4"  = "F4"
)

################################################################################
# Patient-level table
################################################################################

# Patient rows have a PID; legend rows below the data do not
raw <- raw[!is.na(raw[["PID"]]), ]

dat <- with(raw, data.frame(
  id    = StudyID,
  pid   = PID,
  # clinical
  sex   = factor(c(m = "Male", w = "Female")[trimws(Sex)], levels = c("Male", "Female")),
  age   = to_num(Age),
  bmi   = to_num(BMI),
  dm    = factor(c("1" = "Yes", "2" = "No")[as.character(`Diabetes 1=ja/2=nein`)],
                 levels = c("No", "Yes")),
  hba1c = to_num(`HbA1c (%)`),
  # 1 (with any text) = yes, 2 = no, 0 / k.a. / empty = NA
  smoker = factor(ifelse(grepl("^1", `Smoker (1=ja / 2=nein)`), "Yes",
                  ifelse(trimws(`Smoker (1=ja / 2=nein)`) %in% "2", "No", NA)),
                  levels = c("No", "Yes")),
  # serological
  ast   = to_num(`AST/GOT (U/l)`),
  alt   = to_num(`ALT/GPT (U/l)`),
  ggt   = to_num(`GGT (U/l)`),
  alp   = to_num(`ALP (alk. Phosphatase) (U/l)`),
  bili  = to_num(`Bilirubin (ymol/l)`),
  alb   = to_num(`Albumin (g/L)`),
  crea  = to_num(`Kreatinin (ymol/L)`),
  na    = to_num(`Natrium (mmol/L)`),
  plt   = to_num(`Platelets count (G/l)`),
  inr   = to_num(INR),
  # existing non-invasive scores
  fib4  = to_num(`FIB 4`),
  aar   = to_num(`AST/ALT ratio`),
  fast  = to_num(`FAST Score`),
  # elastography
  swe   = to_num(`SWE Stiffness kPa median`),
  fst   = to_num(`Fibroscan Stiffness (kPa)`),
  cap   = to_num(`CAP median (dB/m)`),
  # histology
  steat = to_num(`Steatosis (NAS)`),
  lob   = to_num(`Lobular inflammation`),
  ball  = to_num(Ballooning),
  nas   = to_num(`NAS score (NASH grading)`),
  fib_raw = trimws(`Fibrosis Staging`),
  cirr_clin = !is.na(`Klinisch Zirrhose keine Histo zu Gruppe F4 zählen`),
  stringsAsFactors = FALSE
))

# NAS recomputed from its components, to flag inconsistent entries
dat$nas_sum <- dat$steat + dat$lob + dat$ball

# "keine Biopsie" is not in the map and stays NA
dat$fib <- factor(unname(fib_map[dat$fib_raw]), levels = c("F0", "F1", "F2", "F3", "F4"))

################################################################################
# Endpoint: fibrotic MASH (NAS >= 4 & fibrosis >= F2) vs the rest
################################################################################

# NA if NAS or fibrosis stage is missing (%in% would turn a missing stage into FALSE)
dat$fmash <- as.integer(dat$nas >= 4 & dat$fib %in% c("F2", "F3", "F4"))
dat$fmash[is.na(dat$nas) | is.na(dat$fib)] <- NA
dat$disease <- factor(dat$fmash, levels = c(0, 1),
                      labels = c("No fibrotic MASH", "Fibrotic MASH"))

################################################################################
# Table 1: clinical, serological and histological data
################################################################################

# Patients without a defined endpoint are excluded
d1 <- dat[!is.na(dat$fmash), ]

tbl1 <- d1 %>%
  select(disease,
         age, sex, bmi, dm, hba1c, smoker,
         ast, alt, ggt, alp, bili, alb, crea, na, plt, inr,
         fib4, aar, fast,
         steat, lob, ball, nas, fib) %>%
  tbl_summary(
    by = disease,
    type = list(c(steat, lob, ball, nas) ~ "categorical", c(inr, hba1c) ~ "continuous"),
    statistic = list(all_continuous() ~ "{median} ({p25}, {p75})",
                     all_categorical() ~ "{n} ({p}%)"),
    digits = list(all_continuous() ~ 1, all_categorical() ~ c(0, 1)),
    label = list(
      age ~ "Age, years", sex ~ "Sex", bmi ~ "BMI, kg/m2",
      dm ~ "Diabetes", hba1c ~ "HbA1c, %", smoker ~ "Smoker (current or former)",
      ast ~ "AST, U/l", alt ~ "ALT, U/l", ggt ~ "GGT, U/l", alp ~ "ALP, U/l",
      bili ~ "Bilirubin, µmol/l", alb ~ "Albumin, g/l",
      crea ~ "Creatinine, µmol/l", na ~ "Sodium, mmol/l",
      plt ~ "Platelets, G/l", inr ~ "INR",
      fib4 ~ "FIB-4", aar ~ "AST/ALT ratio", fast ~ "FAST score",
      steat ~ "Steatosis grade", lob ~ "Lobular inflammation grade",
      ball ~ "Ballooning grade", nas ~ "NAS", fib ~ "Fibrosis stage"),
    missing = "ifany",
    missing_text = "Missing"
  ) %>%
  add_overall(last = FALSE) %>%
  modify_header(label ~ "**Characteristic**")

################################################################################
# Box plots of SWE and Fibroscan by outcome definition
# Each outcome: patients with the outcome defined and both SWE and Fibroscan
################################################################################

box_out <- list(
  "Fibrotic MASH"  = dat$fmash,
  "F0/F1 vs F2-F4" = as.integer(as.integer(dat$fib) >= 3),
  "F0-F2 vs F3-F4" = as.integer(as.integer(dat$fib) >= 4),
  "F0-F3 vs F4"    = as.integer(as.integer(dat$fib) >= 5)
)
box_dat <- do.call(rbind, lapply(names(box_out), function(o) {
  ok <- !is.na(box_out[[o]]) & !is.na(dat$swe) & !is.na(dat$fst)
  y <- factor(box_out[[o]][ok], levels = 0:1, labels = c("No", "Yes"))
  data.frame(outcome = o, y = y,
             ylab = sprintf("%s\n(n = %d)", as.character(y), ave(seq_along(y), y, FUN = length)),
             swe = dat$swe[ok], fst = dat$fst[ok])
}))
box_dat <- data.frame(box_dat[c("outcome", "y", "ylab")],
  method = rep(c("SWE", "Fibroscan (TE)"), each = nrow(box_dat)),
  kpa = c(box_dat$swe, box_dat$fst))
box_dat$outcome <- factor(box_dat$outcome, levels = names(box_out))
box_dat$method <- factor(box_dat$method, levels = c("SWE", "Fibroscan (TE)"))

p_box <- ggplot(box_dat, aes(ylab, kpa, fill = method)) +
  geom_boxplot(show.legend = FALSE) +
  facet_grid(method ~ outcome, scales = "free_x") +
  scale_fill_manual(values = c("SWE" = "#0072B2", "Fibroscan (TE)" = "#D55E00")) +
  labs(x = "Outcome", y = "Stiffness, kPa") +
  theme_bw()
ggsave(paths$f("box_outcomes.pdf"), p_box, width = 10, height = 6)

################################################################################
# Objects for the report
################################################################################

mm <-dat[which(dat$nas != dat$nas_sum), ]
excl <- dat[is.na(dat$fmash), ]

res <- list(
  tbl1 = tbl1,
  # no SWE or Fibroscan value: left out of the box plots
  box_ids_miss = dat$id[!is.na(dat$fib) & (is.na(dat$swe) | is.na(dat$fst))],
  n_total = nrow(dat),
  n_class = nrow(d1),
  n_unclass = nrow(excl),
  ids_nobx = excl$id[excl$fib_raw %in% "keine Biopsie"],
  n_nobx_cirr = sum(excl$fib_raw %in% "keine Biopsie" & excl$cirr_clin),
  n_fmash = sum(dat$fmash == 1, na.rm = TRUE),
  # NAS entered differs from steatosis + inflammation + ballooning
  nas_mismatch = sprintf("%s (NAS %g, components sum to %g)", mm$id, mm$nas, mm$nas_sum),
  # k.a. / 0 smoker entries
  n_smoker_na = sum(is.na(d1$smoker))
)
