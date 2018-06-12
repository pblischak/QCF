# Change log

## v0.2.0

**New output option for `qcf`**

  - Can now print raw, gene-level quartet values for downstream
    processing and per-gene bootstrapping (add `--printRaw` flag).

**Per-gene boostrapping for QCF confidence intervals**

  - Using the raw output from a `qcf` run, we have added a new Python script
    (`qcf_boot.py`) that can perform bootstrap resampling of gene-level quartets
    to calculate confidence intervals on the species-level QCF values.

**Independent runs across gene sets**

  - For large multilocus data sets, analyzing smaller sets of genes in parallel
    can now be done to make analyses faster. These independent runs are combined
    and analyzed with the `qcf_boot.py` script to calculate QCF values and
    confidence intervals.
