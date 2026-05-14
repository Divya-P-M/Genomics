#  Variant Calling Pipeline (SNPs & Indels)

A robust short-read variant calling workflow optimized for identifying Single Nucleotide Polymorphisms (SNPs) and Insertions/Deletions (Indels) in resequencing datasets using the *Mycobacterium tuberculosis* H37Rv reference genome.

## Project Objective
To efficiently map paired-end high-quality short reads against a microbial reference genome, filter artifacts, and call true genomic variations. Variant outputs are compiled into a standard, annotated Variant Call Format (VCF) database.

## Resources & System Dependencies
*   **Core Aligner:** [Bowtie2](https://sourceforge.net) (Optimized for fast, short-read alignment with minimal reported mappings per read)
*   **Alignment Engineering:** [Samtools](https://htslib.org)
*   **Variant Caller:** [Bcftools & VCFutils](https://github.io)

### Input Requirements
*   **Reference Genome:** `mtb_h37rv_NC_000962.fna` (*M. tuberculosis* H37Rv FASTA file)
*   **Sequencing Data:** `2287_R1.fastq` & `2287_R2.fastq` (High-quality paired-end Illumina reads)

---

## Step-by-Step Workflow Execution

### 1. Reference Genome Indexing (`bowtie2-build`)
Construct a multi-indexed database footprint from the raw reference sequence file to facilitate accelerated downstream coordinate traversals.

```bash
# Usage: bowtie2-build -f <reference_fasta> <index_base_name>
bowtie2-build -f mtb_h37rv_NC_000962.fna h37rv
```
*   **Note:** Keep index files organized within a dedicated subdirectory (e.g., `ref_genomes/`).

### 2. Alignment & Coordinate Optimization

#### A. Paired-End Mapping (`bowtie2`)
Align the paired-end forward and reverse short reads to the generated reference index, outputting a Sequence Alignment Map (SAM).
```bash
bowtie2 -x ref_genomes/h37rv -1 2287_R1.fastq -2 2287_R2.fastq -S 2287_aln.sam &
```
*(The trailing `&` routes the active sequence computation into the terminal background)*

#### B. Compressed Format Conversion (`samtools view`)
Convert the large text-based SAM file into a compressed binary format (BAM) to minimize storage foot-printing.
```bash
samtools view -bS -o 2287_aln.bam 2287_aln.sam &
```

#### C. Coordinate Sorting (`samtools sort`)
Sort the raw alignment coordinates by genomic position. Genomic sort order is an explicit operational requirement for downstream variant calling tools.
```bash
samtools sort 2287_aln.bam -o 2287_aln_sorted.bam
```

### 3. Genotype Likelihood Generation & Variant Filtering

#### A. Pileup Compilation (`samtools mpileup`)
Compute genotype likelihoods across all covered base positions by inspecting the alignment characteristics relative to the physical reference sequence.
```bash
samtools mpileup -g -u -f ref_genomes/mtb_h37rv_NC_000962.fna 2287_aln_sorted.bam > mtb.bcf &
```
*   **Flag Interpretations:** 
    *   `-g`: Instructs the engine to compute genotype likelihoods and store outputs in a Binary Call Format (BCF).
    *   `-u`: Generates uncompressed BCF output to pipeline smoothly into storage matrices.

#### B. Variant Calling & Quality Filtering (`bcftools` & `vcfutils.pl`)
Isolate true polymorphisms from baseline noise by processing the BCF stream through a strict sequencing depth depth-filter ($D \le 100$).
```bash
bcftools view mtb.bcf | vcfutils.pl varFilter -D100 > mtb_varflt.vcf
```

---

## Understanding the Output Format (VCF)

The resulting `mtb_varflt.vcf` file is a structured data sheet divided into three primary components:

1.  **Metadata Lines (`##` prefix):** Rich key-value pairs defining the exact software configurations, pipeline commands, reference paths, and information tags used.
2.  **Header Line (`#` prefix):** A fixed 8-column tracking structure defining the core variant attributes, ending with optional context metrics:
    ```tsv
    #CHROM   POS   ID   REF   ALT   QUAL   FILTER   INFO   [FORMAT]   [SAMPLE]
    ```
3.  **Data Lines:** Row-by-row genomic records containing the precise chromosome position, wild-type allele (**REF**), mutated target marker (**ALT**), Phred-scaled confidence score (**QUAL**), and conditional annotation tags tracking polymorphism categories (SNPs vs. Indels).
