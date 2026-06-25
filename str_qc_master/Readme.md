# str_qc_master

`str_qc_master` is a DNAnexus applet for the UK Biobank Research Analysis Platform (RAP) that performs quality control for genome-wide Short Tandem Repeat (STR) genotype datasets.

The applet coordinates chromosome-level QC jobs, combines their outputs, removes low-quality samples and variants, and generates a final filtered STR VCF suitable for downstream genome-wide association analyses.

This repository contains the implementation used in my PhD thesis. The code has been preserved to ensure computational reproducibility

---

# Overview

The applet acts as the master controller for the STR quality-control pipeline.

Rather than performing all quality-control steps itself, it launches chromosome-specific QC applets in parallel, waits for their completion, combines the chromosome-level outputs, applies sample- and variant-level filters, and produces a single genome-wide STR VCF.

The workflow was developed for large whole-genome sequencing datasets within the UK Biobank Research Analysis Platform.

---

# Quality-control workflow

The pipeline performs the following stages:

1. Launch chromosome-specific STR QC jobs in parallel.
2. Calculate genome-wide Genotype Average Allele Probability (GAAP) for every sample.
3. Remove samples failing SNP QC or GAAP thresholds.
4. Collect chromosome-level QC summaries.
5. Create a genome-wide list of STR loci failing QC.
6. Concatenate chromosome VCFs.
7. Remove excluded samples.
8. Remove excluded STR loci.
9. Calculate per-sample STR missingness.
10. Remove samples exceeding the missingness threshold.
11. Produce a final indexed genome-wide VCF.

---


---

# Quality-control criteria

The applet coordinates several complementary quality-control metrics, including:

* Average absolute purity (AAP)
* Genome-wide AAP (GAAP)
* Normalised Locus Coverage (NLC)
* Number of Third Species Reads (NTSR)
* Minimum Allelic Depth of Spanning reads (MADSP)
* Variant missingness
* Sample missingness
* Homozygosity-based filtering
* Allelic Depth of Spanning Reads (ADSP) comparison
* External SNP quality-control exclusions
* Removal of multi-STR catalogue entries

Thresholds for these metrics are supplied by the user when launching the applet.

---

# Inputs

The applet requires:

* genome-wide STR VCF
* corresponding VCF index
* chromosome list
* SNP QC exclusion list
* list of multi-STR catalogue entries
* quality-control thresholds including:

  * AAP
  * GAAP
  * NLC
  * NTSR
  * MADSP
  * variant missingness
  * sample missingness
  * homozygosity threshold
  * ADSP comparison threshold

---

# Outputs

The applet produces four principal outputs.

## Final STR dataset

A bgzipped, indexed genome-wide STR VCF containing all retained samples and loci after quality control.

## Variant exclusion list

A text file listing every STR locus removed during quality control together with the reason for exclusion.

## Sample exclusion list

A text file listing every removed sample together with the corresponding exclusion criterion (for example SNP QC, low GAAP or excessive missingness).

## VCF index

A CSI index for the final VCF.

---

# Parallel execution

Chromosome-specific quality-control jobs are submitted independently through DNAnexus.

Once all chromosome jobs have completed successfully, the master applet:

* downloads chromosome outputs;
* concatenates chromosome VCFs;
* performs genome-wide filtering;
* writes the final analysis-ready STR dataset.

This parallel design substantially reduces total runtime on large whole-genome sequencing datasets.



---

# Dependencies

The workflow relies on:

* DNAnexus
* bcftools
* tabix
* jq
* Bash
* chromosome-level `str_qc_eh` applets

---

# Reproducibility

This repository contains the implementation used to generate the final quality-controlled STR dataset analysed in my PhD thesis. The code has been preserved to document the computational workflow used during the study.

---

# Author

Jonny Else


