# Bismark Bisulfite-Seq Alignment & Methylation Calling Pipeline

An optimized bioinformatics workflow for the time-efficient analysis of Bisulfite-Seq (BS-Seq) data, featuring reference genome preparation, read mapping, context-specific methylation extraction, and HTML reporting.

## Project Objective
To map bisulfite-treated single-end sequencing reads and perform quantitative cytosine methylation calling against the human reference genome (`GRCh37`) using sample data derived from human embryonic stem (hES) cells.

## Resources & System Dependencies
*   **Core Toolset:** [Bismark](https://github.com) 
*   **Documentation:** [Bismark User Guide Reference](https://rawgit.com/FelixKrueger/Bismark/master/Docs/Bismark_User_Guide.html)
*   **Underlying Aligner:** Bowtie2 (or Hisat2)

### Input Requirements
*   **Reference Genome Folder:** `human_GRCh37/` (Contains unmodified genomic FASTA files)
*   **Sequencing Data:** `test_data.fastq` (10,000 single-end shotgun BS reads from human ES cells)

---

## Step-by-Step Workflow Execution

### Step 1: Bisulfite Genome Preparation (`bismark_genome_preparation`)
This step prepares the reference genome by performing in silico bisulfite conversions ($C \rightarrow T$ and $G \rightarrow A$) and creating parallel genomic indexes using `bowtie2-build`.

```bash
# Usage: bismark_genome_preparation [options] <path_to_genome_folder>
bismark_genome_preparation human_GRCh37/
```
*   **Output:** Two dedicated bisulfite genome subdirectories generated directly within `human_GRCh37/`.

### Step 2: Alignment & Core Processing (`bismark`)
Performs the actual mapping of bisulfite reads to the newly prepared indexes and tracks directional methylation call strings.

```bash
# Usage: bismark [options] --genome <genome_folder> {-1 <mates1> -2 <mates2> | <singles>}
bismark human_GRCh37/ test_data.fastq
```
*   **Output Files:**
    1.  `test_data_bismark_bt2.bam`: Compressed binary alignment file containing all mappings and raw methylation strings.
    2.  `test_data_bismark_SE_report.txt`: Text file logging alignment rates and overall methylation metrics.

### Step 3: Methylation Extraction (`bismark_methylation_extractor`)
Operates on the BAM file to dissect genomic coordinates and separate every individual analyzed cytosine into distinct sequence context tracks.

```bash
# Usage: bismark_methylation_extractor [options] <filenames>
bismark_methylation_extractor test_data_bismark_bt2.bam
```
*   **Output Strands:** Methylated cytosines are flagged as forward reads (`+`), and unmethylated cytosines are logged as reverse reads (`-`) across three separate files:
    *   `CpG_context_test_data_bismark_bt2.txt`
    *   `CHG_context_test_data_bismark_bt2.txt`
    *   `CHH_context_test_data_bismark_bt2.txt`

### Step 4: Graphical Report Generation (`bismark2report`)
Automates data compilation by auto-detecting all alignment, splitting, extraction, and M-bias reports within the current working folder to build an interactive summary dashboard.

```bash
# Run within the directory containing your step 2 and step 3 output files
bismark2report
```
*   **Output:** An interactive, single-file graphical HTML quality report page (`.html`) suitable for web browser evaluation.

---

## Summary of Biological Context
*   **Sample Source:** Human embryonic stem cells (hES cells) — Pluripotent cells captured from the inner cell mass of 3–5 day old blastocyst-stage embryos.
*   **Downstream Applications:** The resulting sequence-context files (`CpG`, `CHG`, `CHH`) provide the necessary single-base resolution required to map global epigenetic patterns, pluripotency markers, and cellular differentiation footprints.
A
