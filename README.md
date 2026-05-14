# Multi-Omics, Variant Discovery, & Structural Genomics Pipelines

A comprehensive suite of production-ready bioinformatics workflows spanning epigenomics, 3D chromatin conformation, structural variant discovery, and hybrid de novo genome assembly.

---

## Repository Pipeline Index

*   **ChIP-seq Analysis** -> Transcription Factor & Histone Profiling
*   **Bismark BS-Seq** -> DNA Methylation & Epigenetics
*   **Hi-C Conformation** -> 3D Spatial Chromatin Architecture
*   **Hybrid Assembly** -> De Novo Bacterial Genome Reconstruction
*   **Variant Calling** -> Single Nucleotide Polymorphisms & Indels

---

## 1. ChIP-seq Analysis Pipeline

### Objective
To study transcription factor binding loci and chromatin interactions by analyzing raw single-end sequencing reads against the human reference genome (`hg38`).

### Execution Steps
```bash
cd AG_Chipseq/

# 1. Alignment & Formatting
bowtie2 -q -p 4 -k 1 --no-unal -x hg38 -U 6h_IFNa.fastq.gz -S 6h_IFNa.sam
samtools view -bS 6h_IFNa.sam > 6h_IFNa.bam
samtools sort 6h_IFNa.bam -o 6h_IFNa_sort.bam

# (Repeat identical alignment steps for control file: Input_6h.fastq.gz)

# 2. Peak Calling
macs2 callpeak -t 6h_IFNa_sort.bam -c Input_6h_sort.bam -n 6h_IFNa_output -g hs --bdg -q 0.05 -f BAM

# 3. Extract coordinates for annotation
cut -f1,2,3,4 6h_IFNa_output_peaks.narrowPeak > peaks.bed
```

### Functional Genomic Annotation (R)
```r
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

peakfile <- read.table("peaks.bed", sep = "\t", header = FALSE)
gr <- makeGRangesFromDataFrame(peakfile, seqnames.field = "V1", start.field = "V2", end.field = "V3")

peakAnno <- annotatePeak(gr, tssRegion = c(-3000, 3000), TxDb = txdb, annoDb = "org.Hs.eg.db")
plotAnnoBar(peakAnno)
plotDistToTSS(peakAnno, title = "Distribution of TF-Binding Loci Relative to TSS")
```
*   **Consensus Motif Sequence:** `GAATGGAA` (E-value: `5.3e-258` via MEME-ChIP)

---

## 2. Bismark Bisulfite-Seq Pipeline

### Objective
To map bisulfite-treated single-end sequencing reads and perform quantitative cytosine methylation calling against the human reference genome (`GRCh37`) using sample data from human embryonic stem (hES) cells.

### Execution Steps
```bash
# 1. Bisulfite Genome Preparation
bismark_genome_preparation human_GRCh37/

# 2. Alignment & Mapping
bismark human_GRCh37/ test_data.fastq

# 3. Methylation Extraction (Outputs CpG, CHG, and CHH contexts)
bismark_methylation_extractor test_data_bismark_bt2.bam

# 4. Graphical HTML Report Generation
bismark2report
```

---

## 3. Hi-C Chromatin Contact Mapping Pipeline

###  Objective
To analyze paired-end human adrenal tissue sequence reads (Sample `GSM2322539`, chromosomes 18/19) to resolve local 3D architecture, insulation boundaries, and active/inactive structural compartments.

### Execution Steps
```bash
# 1. Genomic Mapping
fanc map SRR4271982_chr18_19_1.fastq.gzip index/hg19_chr18_19 out/SRR4271982_chr18_19_1.sam
fanc map SRR4271982_chr18_19_2.fastq.gzip index/hg19_chr18_19 out/SRR4271982_chr18_19_2.sam

# 2. Format Conversion & Name-Sorting
samtools view -bS out/SRR4271982_chr18_19_1.sam > out/SRR4271982_chr18_19_1.bam
samtools view -bS out/SRR4271982_chr18_19_2.sam > out/SRR4271982_chr18_19_2.bam
samtools sort -n out/SRR4271982_chr18_19_1.bam -o out/SRR4271982_chr18_19_1.sorted.bam
samtools sort -n out/SRR4271982_chr18_19_2.bam -o out/SRR4271982_chr18_19_2.sorted.bam

# 3. Mate Pairing & Object Formation
fanc pairs out/SRR4271982_chr18_19_1.sorted.bam out/SRR4271982_chr18_19_2.sorted.bam out/SRR4271982_chr18_19.pairs -g hg19_chr18_19_re_fragments.bed
fanc hic out/SRR4271982_chr18_19.pairs out/fragment.hic

# 4. Multi-Resolution Binning
fanc hic out/fragment.hic out/fragment_1mb.hic -b 1mb
fanc hic out/fragment.hic out/fragment_50kb.hic -b 50kb

# 5. Matrix Proximity Plots & A/B Compartments
fancplot -o 50kb_plot.png chr18:6mb-10mb -p triangular -vmax 30 out/fragment_50kb.hic
fanc compartments out/fragment_1mb.hic out/1mb.ab
fancplot -o 1mb.png chr18 -p square out/1mb.ab -vmin -0.75 -vmax 0.75 -c RdBu_r
```

---

## 4. Hybrid Genome Assembly Pipeline

### Objective
To construct a complete, unfragmented *Streptococcus thermophilus* genome by pairing accurate Illumina short reads with structural Nanopore long reads.

### Execution Steps
```bash
# 1. Execute Hybrid Co-Assembly
spades.py -1 Short_1.fastq.gz -2 Short_2.fastq.gz --nanopore ont.fastq -o spades_out

# 2. Evaluate Assembly Quality metrics
quast.py spades_out/scaffolds.fasta -o quast_hybrid_out
```

### Assembly Metric Comparison
*   **Contiguity Increase:** Hybrid assembly increased the **N50 from 64,398 bp to 196,927 bp** compared to short-read only.
*   **Fragmentation Drop:** Total genomic contigs fell from **55 down to 14**, allowing long reads to seamlessly bridge repeat regions.

---

## 5. Variant Calling Pipeline (SNPs & Indels)

### Objective
To map paired-end high-quality short reads against the *Mycobacterium tuberculosis* H37Rv reference genome to isolate true genomic variations.

### Execution Steps
```bash
# 1. Reference Indexing
bowtie2-build -f mtb_h37rv_NC_000962.fna h37rv

# 2. Alignment & Optimization
bowtie2 -x ref_genomes/h37rv -1 2287_R1.fastq -2 2287_R2.fastq -S 2287_aln.sam &
samtools view -bS -o 2287_aln.bam 2287_aln.sam &
samtools sort 2287_aln.bam -o 2287_aln_sorted.bam

# 3. Pileup & Quality Filtering (Max Depth = 100)
samtools mpileup -g -u -f ref_genomes/mtb_h37rv_NC_000962.fna 2287_aln_sorted.bam > mtb.bcf &
bcftools view mtb.bcf | vcfutils.pl varFilter -D100 > mtb_varflt.vcf
```

---

