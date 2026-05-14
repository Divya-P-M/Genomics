# Hi-C Chromatin Contact Mapping & Compartment Analysis Pipeline

An end-to-end structural genomics workflow for processing chromosome-conformation-capture (Hi-C) sequencing data to construct spatial chromatin contact maps, discover Topologically Associating Domains (TADs), and determine A/B genomic compartments.

## Project Objective
To analyze paired-end human adrenal tissue sequence reads (Sample `GSM2322539`, Illumina HiSeq 2000 platform) restricted to chromosomes 18 and 19. The pipeline maps spatial proximity interactions to resolve local 3D architecture, insulation boundaries, and active/inactive structural compartments.

## Resources & Environment Setup

### Server Environment
*   **Host System:** Remote Laboratory `mgu` Server
*   **Network IP Address:** `10.20.9.99`
*   **Shared User Credentials:** Username: `mgu` | Password: `genomeinfo`
*   **Workspace Allocation:** Individual tracking directories located inside designated `student*` workspace folders.

### Software Stack
*   **Primary Framework:** [FAN-C (Framework for Analysis of Hi-C)](https://readthedocs.io)
*   **Core Utility:** [Samtools](https://htslib.org)
*   **Plotting Engine:** [fancplot](https://pypi.org) (Requires local desktop setup via: `pip install fanc`)

---

## Step-by-Step Workflow Execution

The analysis maps raw data through the standard architecture: **Mapping ->  Pairing -> Object Formation -> Binning.**

### 1. Genomic Mapping (`fanc map`)
Align both forward and reverse paired-end bisulfite-free reads independently to the pre-indexed `hg19_chr18_19` reference genome using FAN-C's alignment wrapper.

```bash
# Map Read 1 (Forward)
fanc map SRR4271982_chr18_19_1.fastq.gzip index/hg19_chr18_19 out/SRR4271982_chr18_19_1.sam

# Map Read 2 (Reverse)
fanc map SRR4271982_chr18_19_2.fastq.gzip index/hg19_chr18_19 out/SRR4271982_chr18_19_2.sam
```

### 2. Format Conversion & Name-Sorting (`Samtools`)
Convert text alignment files to compressed binary maps, and sort them strictly by read name (`-n`) instead of coordinate location. Name-sorting is mandatory for downstream mate pairing.

```bash
# Convert SAM to BAM
samtools view -bS out/SRR4271982_chr18_19_1.sam > out/SRR4271982_chr18_19_1.bam
samtools view -bS out/SRR4271982_chr18_19_2.sam > out/SRR4271982_chr18_19_2.bam

# Name-Sort alignments
samtools sort -n out/SRR4271982_chr18_19_1.bam -o out/SRR4271982_chr18_19_1.sorted.bam
samtools sort -n out/SRR4271982_chr18_19_2.bam -o out/SRR4271982_chr18_19_2.sorted.bam
```

### 3. Mate Pairing (`fanc pairs`)
Combine the separated name-sorted files back into physical structural ligated pairs. The algorithm filters unligated fragments using a restriction enzyme digest footprint dictionary (`hg19_chr18_19_re_fragments.bed`).

```bash
fanc pairs out/SRR4271982_chr18_19_1.sorted.bam out/SRR4271982_chr18_19_2.sorted.bam out/SRR4271982_chr18_19.pairs -g hg19_chr18_19_re_fragments.bed
```

### 4. Hi-C Object Initialization (`fanc hic`)
Compile the raw textual coordinates file into a highly efficient compiled `.hic` genomic container matrix.

```bash
fanc hic out/SRR4271982_chr18_19.pairs out/fragment.hic
```

### 5. Multi-Resolution Genomic Binning
Bin the multi-contact object arrays into standard fixed intervals to analyze structures at different scales (e.g., larger compartments vs. local loops).

```bash
# Coarse resolution matrix for global compartments (1 Megabase)
fanc hic out/fragment.hic out/fragment_1mb.hic -b 1mb

# Fine resolution matrix for local domain structures (50 Kilobases)
fanc hic out/fragment.hic out/fragment_50kb.hic -b 50kb
```

---

## Matrix Visualizations & Downstream Analyses

### 1. Triangular Matrix Proximity Plots
Generate local interactive proximity triangle plots to evaluate distance-dependent polymer decay and map insulation boundaries.

```bash
fancplot -o 50kb_plot.png chr18:6mb-10mb -p triangular -vmax 30 out/fragment_50kb.hic
```
*   **Structural Insight:** The heatmaps show the strongest signal intensity concentrated tightly along the major central diagonal axis, confirming that linear proximity governs high baseline collision frequency. The distinct clustering of high-intensity isolated triangles along the path indicates the position of stable Topologically Associating Domains (TADs) bounded by structural insulator elements (e.g., CTCF binding sites).

### 2. A/B Eigenvector Compartment Profiling
Perform matrix eigenvector decomposition analysis to calculate the checkboard-like open/closed structural identities across chromosome domains.

```bash
# Calculate eigenvector compartment matrices
fanc compartments out/fragment_1mb.hic out/1mb.ab

# Render square whole-chromosome contact maps
fancplot -o 1mb.png chr18 -p square out/1mb.ab -vmin -0.75 -vmax 0.75 -c RdBu_r

# Annotate genomic domains with GC content correlations
fanc compartments -g hg19_chr18_19.fa -d out/1mb_gc.domains.bed out/1mb.ab
```

---

## Functional Genomics Interpretations

### Whole-Chromosome Square Heatmaps
*   **Checkerboard Matrix Patterns:** The full square contacts display a classic grid landscape mapping interaction preferences across megabase distances. High interaction frequencies occur between domains assigned to the same spatial context, regardless of linear separation distance.
*   **Anti-Correlated Signals:** Blue zones mark strong spatial avoidance between distinct structural designations, confirming that chromatin partitions into distinct, non-mixing micro-environments.

### BED File Track Segmentations (`1mb_gc.domains.bed`)
The compiled segmentation tracks split the sequence landscape into an alternating, organized multi-domain sequence:
*   **Compartment A (Positive Eigenvector Values):** Corresponds to highly transcribed, active, gene-dense euchromatin domains localized away from the nuclear periphery.
*   **Compartment B (Negative Eigenvector Values):** Maps transcriptionally silent, gene-poor heterochromatic regions or nuclear-lamina-associated domains (LADs).
