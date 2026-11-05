################################################################################
# Import the Excel masterfile (unmodified)
# Output: raw (one row per Excel row, original column names)
################################################################################

xlsx_file <- paths$d("Masterfile MASLD ELASTO PZ 20261010.xlsx")

raw <- suppressMessages(read_excel(xlsx_file, sheet = "Tabelle1"))
raw <- as.data.frame(raw, check.names = FALSE)
