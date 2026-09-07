LATEXMK ?= latexmk
BUILD_DIR := build
TARGET := main.pdf

.PHONY: all clean

all: $(TARGET)

$(TARGET): main.tex references.bib figures/film-thickness-vs-ca.pdf data/baseline-film.csv
	mkdir -p $(BUILD_DIR)
	$(LATEXMK) -pdf -interaction=nonstopmode -halt-on-error -outdir=$(BUILD_DIR) main.tex
	cp $(BUILD_DIR)/$(TARGET) $(TARGET)

clean:
	$(LATEXMK) -C -outdir=$(BUILD_DIR) main.tex
	rm -f $(TARGET)
	rmdir $(BUILD_DIR) 2>/dev/null || true
