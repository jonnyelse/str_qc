# str_qc_eh

`str_qc_eh` is a DNAnexus applet for the UK Biobank Research Analysis Platform (RAP) that performs chromosome-level quality control of ExpansionHunter Short Tandem Repeat (STR) genotype calls.

The applet applies genotype-level quality filters, identifies poorly performing STR loci using multiple complementary quality-control metrics, and produces filtered chromosome-level VCFs together with QC reports.

This repository contains the implementation used for the STR quality-control analyses described in my PhD thesis. The implementation has been preserved for computational reproducibility

---

# Overview

The applet performs chromosome-level STR quality control in two stages.

First, low-confidence individual genotype calls are removed using ExpansionHunter quality metrics.

Second, variant-level quality-control statistics are calculated to identify loci showing evidence of systematic genotyping artefacts.

The resulting chromosome-level VCFs and QC summaries are subsequently combined by the `str_qc_master` applet to generate the final analysis-ready STR dataset.

---

# Quality-control workflow

For each chromosome the applet:

1. extracts chromosome-specific STR genotypes;
2. removes individual genotype calls failing quality thresholds;
3. calculates variant missingness;
4. evaluates deviations from expected homozygosity;
5. identifies loci showing abnormal Allelic Difference in Spanning Reads (ADSP);
6. generates chromosome-level QC summaries;
7. outputs a filtered chromosome VCF.

---


# Genotype-level quality control

Individual genotype calls are set to missing if they fail one or more user-defined thresholds:

* **AAP (Average Absolute Purity)** below threshold
* **NTSR (Number of Third Species Reads)** above threshold
* **NLC (Normalised Locus Coverage)** below threshold
* **MADSP (Minimum Allelic Depth of Spanning Reads)** below threshold

These filters remove low-confidence genotype calls.

---

# Variant-level quality control

Following genotype filtering, the applet evaluates each STR locus using several complementary metrics.

## Variant missingness

Variants exceeding the specified proportion of missing genotypes are flagged for exclusion.

## Homozygosity test

Observed homozygote frequencies are compared with Hardy-Weinberg expectations for each allele using exact binomial tests.

For each STR locus the smallest homozygosity p-value is retained, and loci within the user-defined lowest proportion are flagged as potential genotyping artefacts.

## Allelic Depth of Spanning Reads (ADSP)

For each allele, spanning-read support is compared between homozygous and heterozygous genotypes using Welch's t-test.

For each locus the smallest comparison p-value is retained, and loci within the user-defined lowest proportion are flagged for potential genotyping bias.

This metric was developed to identify loci where spanning-read support differs systematically between genotype classes, suggesting possible genotyping error.

---

# Inputs

The applet requires:

* genome-wide ExpansionHunter VCF
* corresponding VCF index
* chromosome identifier
* output prefix
* quality-control thresholds for:

  * AAP
  * NTSR
  * NLC
  * MADSP
  * variant missingness
  * homozygosity
  * ADSP comparison

---

# Outputs

The applet produces:

## Filtered chromosome VCF

A bgzipped chromosome VCF in which low-quality genotype calls have been set to missing.

## Variant QC summary

A table listing loci identified by each quality-control procedure.

## High missingness variants

Variants exceeding the missingness threshold.

## Homozygosity failures

Variants identified by the homozygosity test.

## ADSP failures

Variants identified by the ADSP comparison procedure.

These outputs are intended for downstream aggregation by the `str_qc_master` applet.

---


# Dependencies

The workflow uses:

* bcftools
* tabix
* R
* data.table
* dplyr
* tidyr
* ggplot2
* ggrepel
* psych
* additional supporting R packages

---

# Reproducibility

This repository contains the implementation used to generate the chromosome-level STR quality-control results presented in my PhD thesis. The code has been retained to document the computational workflow used during the study.

---

# Author

Jonny Else




