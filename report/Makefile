# ============================================================
# Study: 
# Purpose: knit analysis.rnw and compile report.tex to PDF
# ============================================================

REPORT        = report
SWEAVE_SRC    = analysis.rnw
DATE          := $(shell date +%Y-%m-%d)
OUTDIR        = .
RES_RDS       = ../data/res.rds
FINAL_PDF     = $(OUTDIR)/MASLD_AnalysisReport_$(DATE).pdf

.PHONY: all initial clean open close res

all: $(FINAL_PDF)

# --- Full report -------------------------------------------------------------

# Step 0: pipeline results (data/res.rds); rebuilt when the R scripts change
# After changing the Excel file run `make res` (file name has spaces, so not tracked here)
$(RES_RDS): $(wildcard ../code/*.r ../code/*.R)
	cd ../code && Rscript 00_master.R

res:
	cd ../code && Rscript 00_master.R

# Step 1: knitr — generates analysis.tex (figures are written to figures/)
analysis.tex: $(SWEAVE_SRC) $(RES_RDS)
	echo "library(\"knitr\"); knit(\"$(SWEAVE_SRC)\");" | R --vanilla

# Step 2: Compile LaTeX (3 passes) + BibTeX
$(REPORT).pdf: $(REPORT).tex analysis.tex
	pdflatex -interaction batchmode $(REPORT).tex
	bibtex $(REPORT) || true
	pdflatex -interaction batchmode $(REPORT).tex
	pdflatex -interaction batchmode $(REPORT).tex

# Step 3: Copy to output directory with datestamped name
$(FINAL_PDF): $(REPORT).pdf
	cp $(REPORT).pdf $@


# --- Open / close PDF in Windows (from WSL) ---------------------------------
open: $(FINAL_PDF)
	cmd.exe /c start "" "$$(wslpath -w $(FINAL_PDF))"

close:
	-cmd.exe /c "taskkill /F /IM Acrobat.exe /T" 2>/dev/null
	-cmd.exe /c "taskkill /F /IM AcroRd32.exe /T" 2>/dev/null

# --- Clean -------------------------------------------------------------------
clean: close
	rm -f $(REPORT).aux $(REPORT).bbl $(REPORT).blg \
	      $(REPORT).out $(REPORT).toc $(REPORT).log \
	      analysis.tex
	rm -rf figures
	find . -maxdepth 1 -name "*.pdf" \
	  ! -name "ub_16pt-cmyk.pdf" \
	  ! -name "ub_16pt-bl.pdf" \
	  -delete
	find . -maxdepth 1 -name "*.tmp.*" -delete
	find .. -maxdepth 1 -name ".aider*" -delete
