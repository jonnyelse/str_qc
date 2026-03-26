# STR QC Pipeline (DNAnexus / UK Biobank RAP)

This repository contains two DNAnexus applets for quality control of short tandem repeat (STR) data derived from UK Biobank RAP.

## Overview

The pipeline performs per-chromosome STR quality control, followed by merging and filtering to produce a final QCed VCF.

### Applets

* **str_qc_master**
  Orchestrates chromosome-level processing, merges outputs, applies sample and variant filtering, and generates the final QCed dataset.

* **str_qc_eh**
  Performs per-chromosome QC, including:

  * Genotype-level filtering (AAP, NTSR, NLC, MADSP)
  * Variant missingness filtering
  * Homozygosity-based filtering (HWE approximation)
  * ADSP-based QC comparisons

## Environment

These applets are designed to run within the **UK Biobank Research Analysis Platform (RAP)** on DNAnexus.

## Requirements

* DNAnexus platform (`dx-toolkit`)
* `bcftools`
* R (with required packages listed in `extraQC.R`)

## Usage

Build the applets:

```
dx build str_qc_master/
dx build str_qc_eh/
```

Run via DNAnexus using appropriate inputs.

## Data Privacy

No patient-level or identifiable data is included in this repository.

All processing is intended to occur within the secure RAP environment.
Do not export or commit any generated data files (e.g. VCFs, TSVs, sample lists).

## Notes

* Outputs and intermediate files are excluded via `.gitignore`
* Scripts are provided for reproducibility and transparency of QC steps

## Author

Jonny Else
