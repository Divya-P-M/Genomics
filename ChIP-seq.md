# ChIP-seq Analysis Pipeline (hg38)

An end-to-end bioinformatics workflow for processing, aligning, calling peaks, annotating, and identifying motifs from Chromatin Immunoprecipitation Sequencing (ChIP-seq) data.

## Project Objective
To study transcription factor binding sites and chromatin interactions by analyzing raw single-end sequencing reads against the human reference genome (`hg38`).

## Resources & Dependencies
*   **Alignment:** [Bowtie2](https://sourceforge.net)
*   **File Processing:** [Samtools](https://htslib.org)
*   **Peak Calling:** [MACS2](https://github.com)
*   **Annotation:** [ChIPseeker](https://bioconductor.org) (R package)
*   **Visualization & Extraction:** [UCSC Genome Browser](https://genome.ucsc.edu/)
*   **Motif Discovery:** [MEME-ChIP](https://meme-suite.org)

---

## Step-by-Step Workflow Execution

### Data Preparation
All input data and intermediate outputs are managed inside the tracking folder:
```bash
cd AG_Chipseq/
```

### 1. Alignment to Reference Genome (`Bowtie2`)
Align raw single-end fastq files to the `hg38` reference index.
```bash
# Treatment sample alignment
bowtie2 -q -p 4 -k 1 --no-unal -x hg38 -U 6h_IFNa.fastq.gz -S 6h_IFNa.sam

# Control/Input sample alignment
bowtie2 -q -p 4 -k 1 --no-unal -x hg38 -U Input_6h.fastq.gz -S Input_6h.sam
```

### 2. Format Conversion (`SAM` to `BAM`)
Convert text-based alignment data to binary format for speed and storage compression.
```bash
samtools view -bS 6h_IFNa.sam > 6h_IFNa.bam
samtools view -bS Input_6h.sam > Input_6h.bam
```

### 3. Coordinate Sorting
Sort BAM tracking alignments by genomic position.
```bash
samtools sort 6h_IFNa.bam -o 6h_IFNa_sort.bam
samtools sort Input_6h.bam -o Input_6h_sort.bam
```

### 4. Peak Calling (`MACS2`)
Identify statistically significant enriched genomic regions (peaks).
```bash
macs2 callpeak -t 6h_IFNa_sort.bam -c Input_6h_sort.bam -n 6h_IFNa_output -g hs --bdg -q 0.05 -f BAM
```

#### Key Generated Result Files:
*   `*_peaks.narrowPeak`: BED6+4 format containing peak locations, summits, p-values, and q-values.
*   `*_peaks.xls`: Tabular file containing fold enrichment metrics.
*   `*_summits.bed`: Precise point mutation coordinates (Recommended for down-stream motif analysis).
*   `*_treat_pileup.bdg` & `*_control_lambda.bdg`: Signal track files for genome browser track display.

### 5. Coordinate Extraction
Isolate core coordinate metrics for functional genomic annotation.
```bash
cut -f1,2,3,4 6h_IFNa_output_peaks.narrowPeak > peaks.bed
```

---

### 6. Peak Annotation (`ChIPseeker` in R)
Execute the downstream annotation using R to resolve gene assignment and distance mapping to the transcription start site (TSS).

```r
# 1. Install foundational dependencies
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("ChIPseeker", "GO.db", "HDO.db", "clusterProfiler", "org.Hs.eg.db", "TxDb.Hsapiens.UCSC.hg38.knownGene"))

# 2. Load core libraries
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GenomicRanges)

# 3. Read genomic configurations
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
setwd("path/to/AG_Chipseq")

peakfile <- read.table("peaks.bed", sep = "\t", header = FALSE)
gr <- makeGRangesFromDataFrame(peakfile, seqnames.field = "V1", start.field = "V2", end.field = "V3")

# 4. Annotate peaks relative to TSS
peakAnno <- annotatePeak(gr, tssRegion = c(-3000, 3000), TxDb = txdb, annoDb = "org.Hs.eg.db")

# 5. Convert and export data matrices
peakAnno_DF <- as.data.frame(peakAnno)
write.table(peakAnno_DF, file = "annotated_peaks_output.txt", sep = "\t", header = TRUE, row.names = FALSE)

# 6. Generate Profile Plots
plotAnnoBar(peakAnno)
plotAnnoPie(peakAnno)
plotDistToTSS(peakAnno, title = "Distribution of Transcription Factor-Binding Loci \n Relative to TSS")
```

---

### 7. Sequence Extraction (`UCSC Genome Browser`)
To extract DNA strings matching your computed binding targets:
1. Navigate to the [UCSC Genome Browser](https://genome.ucsc.edu/).
2. Select **Genome Browser** -> **Add Custom Tracks** -> Upload your `peaks.bed` file.
3. Launch **Tools** -> **Table Browser** -> Choose `hg38` assembly -> Output Format: `sequence`.

### 8. Motif Analysis (`MEME-ChIP`)
Upload the exported multi-FASTA file to the online [MEME-ChIP portal](https://meme-suite.org) to calculate consensus regulatory bindings.

## Summary Findings
* **Identified Consensus Motif Sequence:** `GAATGGAA`
* **Statistical Significance:** `E-value: 5.3e-258`
* **Biological Context:** Strong biochemical proof indicating a fixed targeted recruitment locus for the tracking transcription factor.

