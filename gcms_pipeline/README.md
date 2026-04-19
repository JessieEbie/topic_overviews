# GC-MS Cuticular Compound Pipeline

This folder is a starter processing layer for the GC-MS Cuticular Compound Lab.

The first script normalizes exported peak tables into the CSV shape that
`gcms_cuticle_compound_lab.html` can open. It also adds conservative class hints
from EI mass-spectrum fragments.

There is now also a raw-run script that turns an Agilent `.D` folder or zipped
`.D` folder into a workstation bundle JSON that the browser app can reopen.

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

## Raw Agilent bundle

For a stronger local-processing path, build a workstation bundle from a raw
Agilent run:

```bash
python3 gcms_pipeline/process_agilent_run.py /path/to/run.D processed_run.json
```

or

```bash
python3 gcms_pipeline/process_agilent_run.py /path/to/run.D.zip processed_run.json
```

Then open `gcms_cuticle_compound_lab.html` and use **Load project or bundle** to
open `processed_run.json`.

The bundle includes:

- TIC trace
- first-pass detected peaks
- apex-window averaged spectra
- run metadata from `runstart.txt`, `acqmeth.txt`, and `PRE_POST.INI` when present

This is the beginning of the local analysis core layer. It is still more modest
than OpenChrom, but it is a better foundation than browser-only raw parsing.

## Agilent export options

You can use this project at three levels.

### Which file should I try first?

From the Agilent-style file list, start with:

- `RESULTS.CSV`

That is the most likely file to contain the peak or library-search table this
prototype can read directly.

These files are useful later, but are not the right first upload for the browser
interface:

- `DATA.MS`: raw or semi-raw MS data for a future conversion/parsing path
- `FID1A.CH`: chromatogram channel data
- `EDFIL.DAT`, `edfil.val`: instrument/software data files
- `acqmeth.txt`, `rteres.txt`, `runstart.txt`: method or run metadata
- `cnorm.ini`, `pre_post.ini`: configuration files

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

The new `process_agilent_run.py` script is a first step toward that stronger
local path, even though it is not yet doing full deconvolution or library
search.

## Identification caution

The class hints are deliberately conservative. They can help flag likely
alkanes, methyl-branched alkanes, alkenes, esters, or oxygenated/unknown
compounds, but exact structural calls still need retention indices, standards,
derivatization when appropriate, and biological context.
