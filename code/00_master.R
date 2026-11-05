################################################################################
# Master script: packages, paths, runs the pipeline (01 to 04)
# Run from anywhere: Rscript code/00_master.R (or `make res` in report/)
################################################################################

library("readxl")
library("dplyr")
library("gtsummary")
library("ggplot2")
library("tram")
library("pROC")
library("parallel")

# Project root = parent of the folder holding this script
if (!exists("root_dir")) {
  args <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", args[grepl("^--file=", args)])
  root_dir <- if (length(f)) normalizePath(file.path(dirname(f), "..")) else {
    if (basename(getwd()) %in% c("code", "report")) normalizePath("..") else getwd()
  }
}

paths <- list(
  data   = file.path(root_dir, "data"),
  code   = file.path(root_dir, "code"),
  report = file.path(root_dir, "report"),
  fig    = file.path(root_dir, "report", "figures")
)
paths$d <- function(...) file.path(paths$data, ...)
paths$f <- function(...) file.path(paths$fig, ...)
dir.create(paths$fig, showWarnings = FALSE, recursive = TRUE)

source(file.path(paths$code, "01_import.R"))
source(file.path(paths$code, "02_data_preparation.R"))
source(file.path(paths$code, "03_models.R"))
source(file.path(paths$code, "04_fibrosis.R"))

saveRDS(res, paths$d("res.rds"))
