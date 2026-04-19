# GC-MS Cuticular Compound Pipeline

This folder is a starter processing layer for the GC-MS Cuticular Compound Lab.

The first script normalizes exported peak tables into the CSV shape that
`gcms_cuticle_compound_lab.html` can open. It also adds conservative class hints
from EI mass-spectrum fragments.

## Input table

Use CSV or TSV with any of these columns:

- `rt`, `retention_time`, `retention time`, or `time`
- `area`, `intensity`, `height`, or `abundance`
- `name`, `compound`, or `candidate`
- `formula`
- `class`, `compound_class`, or `type`
- `match`, `score`, or `similarity`
- `spectrum`, `ions`, or `mass_spectrum`
- `notes` or `comment`

The `spectrum` column should contain ion pairs like:

```text
57:100 71:66 85:42 99:24 113:12 408:7
```

## Run

```bash
python3 gcms_pipeline/process_peak_table.py input.csv reviewed_input.csv
```

Then open `gcms_cuticle_compound_lab.html` and upload `reviewed_input.csv`.

## Agilent export options

You can use this project at three levels.

### 1. Quick MassHunter or ChemStation report export

Export a peak, compound, or library-search report as CSV, TXT, or XLSX. This is
the best first option because it uses the workflow the lab already trusts.

Try to include:

- retention time
- area or height
- candidate compound name
- library match score or quality
- formula, if available
- class or compound type, if available

### 2. Spectrum-rich Agilent export

If the report template can include apex spectra, major ions, or target ions,
include them. The review interface can use ion pairs such as `57:100 71:66` to
show spectra and generate class hints.

This is the most useful option for explaining why a peak looks like an alkane,
methyl-branched alkane, alkene, ester, or unresolved oxygenated compound.

### 3. Raw `.D` conversion

Agilent raw GC-MS runs are often stored as `.D` folders. A later version of this
pipeline can convert those folders to open formats such as mzML or NetCDF, then
do chromatogram reading, peak detection, deconvolution, library search, and
retention-index calculation before the review step.

This is the most powerful path, but it is not required to start using the
interface.

## Identification caution

The class hints are deliberately conservative. They can help flag likely
alkanes, methyl-branched alkanes, alkenes, esters, or oxygenated/unknown
compounds, but exact structural calls still need retention indices, standards,
derivatization when appropriate, and biological context.
