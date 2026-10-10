# covid_disinfo

Analysis code for:

**Typologies of COVID-19 misinformation beliefs in Canada: a latent class and latent transition analysis of the EKOS COVID-19 Monitor**

Fisman DN, Greer AL, Tuite AR. Submitted to *CMAJ*.

## Overview

This repository contains Python code for:

1. **Latent class analysis (LCA)** of four COVID-19 misinformation items using a custom EM algorithm with 100 random starts and 3-category item coding (Correct / Don't Know / Incorrect).
2. **Posterior-weighted baseline characterization** of latent classes (Table 3), including demographics, vaccination status, institutional trust, media reliance, and concern about disinformation.
3. **Posterior-weighted Wald tests** (chi-square for categorical variables, F-test/ANOVA for continuous variables) comparing characteristics across classes.
4. **Formatted Word table generation** for the manuscript.

## Scripts

| Script | Description |
|--------|-------------|
| `lca_table3.py` | Fits 4-class LCA, computes posterior-weighted estimates and Wald tests, outputs CSV |
| `make_table3_docx.py` | Generates formatted Word table from CSV output |

## Usage

```bash
# Run the LCA and generate the CSV
python lca_table3.py path/to/05021W27_client_data_file.sav --output table3.csv

# Generate the formatted Word table
python make_table3_docx.py --input table3.csv --output Table_3.docx
```

## Dependencies

- Python 3.8+
- numpy
- pandas
- scipy
- pyreadstat
- python-docx

Install with:

```bash
pip install numpy pandas scipy pyreadstat python-docx
```

## Data

The EKOS COVID-19 Monitor data are not publicly available. The data were collected by EKOS Research Associates and are available from the authors on reasonable request, subject to the data sharing agreement with EKOS Research.

## Items

The four common misinformation items (Q34A-D), all false statements:

| Item | Statement |
|------|-----------|
| A | The number of deaths due to COVID-19 has been greatly exaggerated |
| B | The government is hiding deaths caused by vaccines |
| C | COVID-19 vaccines cause infertility |
| D | COVID-19 vaccines alter human DNA |

Response coding: Strongly disagree/Disagree = Correct (0); Don't know = DK (1); Agree/Strongly agree = Incorrect (2).

## Latent classes

The 4-class solution identifies:

- **Well-Informed** (~71%): High P(Correct) across all items
- **Selectively Skeptical** (~9%): Mixed responses, moderate P(Correct)
- **Uncertain/Disengaged** (~15%): Elevated P(Don't Know)
- **Misinformed** (~5%): High P(Incorrect), especially for vaccine items

## License

This code is provided for reproducibility of the analyses described in the manuscript. Please cite the paper if you use or adapt this code.
