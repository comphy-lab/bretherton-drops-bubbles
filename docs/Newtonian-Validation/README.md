# Newtonian bubble report

Scientific report on the Newtonian bubble calculations and their comparison
with reference solutions.

| File or directory | Contents |
|---|---|
| `main.tex` | Report source |
| `references.bib` | Bibliography |
| `main.pdf` | Compiled report |
| `figures/` | Vector/PDF and PNG figures, with their plotting script |
| `data/baseline-film.csv` | Compact data used in the film comparison |
| `Makefile` | LaTeX build |

## Compile

Requires LaTeX with `latexmk`, BibTeX, `natbib`, `graphicx`, `booktabs`,
`amsmath`, `amssymb`, `geometry` and `hyperref`.

From the repository root:

```bash
make -C docs/Newtonian-Validation
```

This produces `docs/Newtonian-Validation/main.pdf`; auxiliary files remain
under `build/`. Rebuild after editing the LaTeX, bibliography or figures.

## Regenerate the figure

Requires Python with NumPy and Matplotlib, plus LaTeX and `dvipng`:

```bash
python3 docs/Newtonian-Validation/figures/plot_film.py
make -C docs/Newtonian-Validation
```

The script reads the compact CSV and writes the PDF and PNG in `figures/`.
It does not require simulation dumps.
