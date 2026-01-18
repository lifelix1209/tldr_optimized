[![DOI](https://zenodo.org/badge/205993555.svg)](https://zenodo.org/badge/latestdoi/205993555)

## tldr

*Transposons from Long DNA Reads*

# Installation

tldr requires python > 3.6 and has the following dependencies:
- HTSLIB/Samtools
- minimap2
- MAFFT
- Exonerate
- some python dependecies in the background. 

## One-step Conda environment setup 
There is a pre-baked Conda (or [mamba](https://anaconda.org/conda-forge/mamba)) environment file provided (tldr.yml) that can be used to create a tldr Conda environment with all of the necessary dependencies. 

```
git clone https://github.com/adamewing/tldr.git
cd tldr
conda env create -f tldr.yml
conda activate tldr
pip install -e $PWD
tldr -h
```
If you use the above method, make sure to activate the Conda environment first with `conda activate tldr` whenever using tldr.

## Installing dependencies seperately
## HTSLIB / SAMtools
Easiest method is via conda:

```
conda install -c bioconda tabix
conda install -c bioconda samtools
```
Manual installation:
```
git clone https://github.com/samtools/htslib.git
git clone https://github.com/samtools/samtools.git

make -C htslib && sudo make install -C htslib
make -C samtools && sudo make install -C samtools
```

## minimap2
Via conda:
```
conda install -c bioconda minimap2
```
For manual installation see [minimap2 github](https://github.com/lh3/minimap2)


## MAFFT
Via conda:
```
conda install -c bioconda mafft
```

For manual installation see the [mafft website](https://mafft.cbrc.jp/alignment/software/linux.html)

Note: different versions of MAFFT may yield different results from tldr. We currently recommend MAFFT v7.480.

## Exonerate
Via conda:

```
conda install -c bioconda exonerate
```
For manual installation see the [exonerate website](https://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate)

# Install

Install tldr package + python dependencies:

```
python setup.py install
```

# Running tldr

## Basic Analysis

Synopsis (minimal input requirements), assuming reads aligned to hg38 using minimap2:
```
tldr -b aligned_reads.bam -e /path/to/tldr/ref/teref.ont.human.fa -r /path/to/minimap2-indexed/reference/genome.fasta --color_consensus
```

## De Novo (Trio) Analysis

For family trio analysis to identify de novo TE insertions:

```
# Step 1: Run on child sample with --denovo flag
tldr -b child.bam -e ref/teref.ont.human.fa -r genome.fasta --denovo -p 8 -o child_output

# Step 2: Analyze parent support
python scripts/parent_support.py \
    --mom_bam mother.bam \
    --dad_bam father.bam \
    --candidates child_output.table.txt \
    --output parent_analysis.json \
    --denovo
```

## Command-line Options

### -b/--bams
Multiple .bam files can provided in a comma-delimited list.

### -e/--elts
Reference elements in .fasta format. The header for each should be formatted as`>Superfamily:Subfamily` e.g. `>ALU:AluYb9`.
If `none` is specifed instead of a filename, tldr will run without a reference TE collection. This is useful for genomes where active mobile element content is not well understood or for unbiased identification of inserted sequenced and is also useful for identifying viral intregration and gene retrocopy insertions.

### -r/--ref
Reference genome .fasta, expects a samtools index i.e. `samtools faidx`.

### -p/--procs
Spread work over _p_ processes. Uses python multiprocessing.

### -m/--minreads
Minimum supporting read count to trigger a consensus / insertion call (default = 3)

### --embed_minreads
Minimum number of reads completely embeddeding the insertion (default = 1, requires at least 1).

### -o/--outbase
Specify a base name for output files. The default is to use the name of the input bam(s) without the .bam extension and joined with "_" if > 1 .bam file given

### -c/--chroms
Specify a text file of chromosome names (one per line) and tldr will focus only on these.

### --max_te_len
Maximum insertion size (default = 10000)

### --min_te_len
Minimum insertion size (default = 200)

### --min_alt_frac
Parameter for allowing base changes in consensus cleanup (default = 0.5)

### --min_alt_depth
Parameter for allowing base changes in consensus cleanup (default = 3)

### --min_total_depth_frac
Parameter for allowing base changes in consensus cleanup (default = 0.25)

### --max_cluster_size
Limit cluster size and downsample clusters larger than the cutoff (default = no limit). Downsampling is biased such that reads completely embedding the inserted sequence are preferred.

### --wiggle
Allows for sloppy breakpoints in initial breakpoint search (default = 200)

### --denovo
Enable de novo candidate generation mode for trio analysis. Outputs lightweight results suitable for parent support analysis. Adds additional columns: `bp_left`, `bp_right`, `wiggle`, `TE_family`, `child_support`. Designed for child-side only candidate discovery in trio (Mom/Dad/Child) workflows.

### --flanksize
Trim reads to contain at most `--flanksize` bases on either side of the insertion. Setting too large makes consensus building slower and more error-prone.

### -n/--nonref
Annotate insertion with known non-reference insertion sites (examples provided in `/path/to/tldr/ref`

### --color_consensus
This will annotate the consensus sequence with ANSI escape characters that yield coloured text on the command-line:
red = TSD, blue = TE insertion sequence, yellow = non-TE insertion sequence
While this looks nice on the command line (try `less -R`) and is helpful for evaluating insertion calls, the output may not translate well to other applications as the escape sequences for the ANSI colours will be embedded in the sequence.

### --detail_output
Creates a directory (name is the output base name) with extended consensus sequences, per-insertion read mapping information and per-insertion .bam files. Required for mCpG analysis.

### --extend_consensus
If --detail_output option is enabled, extend output per-sample consensus by n bases (default 0). This is useful in the analysis of CpG methylation to add context on either end of the insertion.

### --trdcol
Adds 5' and 3' transduction columns needed by the `call_transductions.py` script, if you're into that kind of thing.

### --keep_pickles
Saves pickles for later.

### --use_pickles <folder>
Search specified folder for .pickle files and use them instead of clustering reads. Faster for re-running with different options, requires `--keep_pickles`.

## Output

Some fields in the output table (basename.table.txt) may not be self-explainatory:

### StartTE / EndTE
Start / end position relative to TE consensus provided via `-e/--elts`

### LengthIns
Length of actual inserted sequence. _Not necessarily the same as EndTE-StartTE_

### Inversion
Internal inversion detected in TE

### UnmapCover
Fraction of inserted sequence covered by TE sequence

### MedianMapQ
Median mapping quality score from input .bam(s)

### TEMatch
Overall mean identity to TE in reference library (`-e/--elts`)

### UsedReads
Number of reads used in consensus generation

### SpanReads
Number of reads which completely embed the insertion

### NumSamples
Number of samples (.bam files) in which the insertion was detected

### SampleReads
Per-sample accounting of supporting reads

### EmptyReads
Number of reads spanning both TSDs +/- `--wiggle` parameter with no evidence for insertion, useful for inferring genotype

### EmptyPhase
Phase/Haplotype info for empty site reads (requires PS/HP tags to be present)

### NonRef
If `-n/--nonref` given, annotate whether insertion is a known non-reference insertion ("NA" otherwise)

### TSD
Target site duplication (based on reference genome)

### Consensus
Upper case bases = reference genome sequence, lower case bases = insertion sequence. If `--color_consensus` given TSD will be red, TE will be blue, other inserted sequence (e.g. transduction) will be yellow using ANSI terminal colours (may be affected by specific terminal config)

### Phasing
Phase/Haplotype info for filled site reads (requires PS/HP tags to be present)

### Filter
Annotate whether an insertion call is problematic; "PASS" otherwise (similar to VCF filter column).

### De Novo Output Columns
When running with `--denovo` flag, additional columns are added to the output:

| Column | Description |
|--------|-------------|
| `bp_left` | Left breakpoint coordinate |
| `bp_right` | Right breakpoint coordinate |
| `wiggle` | Wiggle parameter used for clustering |
| `TE_family` | TE family/subfamily classification |
| `child_support` | Read support summary (e.g., "useable:10,embedded:5") |

## De Novo Analysis

The de novo analysis pipeline enables trio analysis to identify potential de novo TE insertions in a child that are not present in either parent.

### Workflow

```
┌─────────────────────────────────────────────────────────────────────┐
│  1. Child Analysis (tldr --denovo)                                  │
│     tldr -b child.bam -e ref.fa -r genome.fa --denovo -o child      │
└─────────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────────┐
│  2. Parent Support Analysis (parent_support.py)                     │
│     python scripts/parent_support.py \                              │
│         --mom_bam mom.bam \                                         │
│         --dad_bam dad.bam \                                         │
│         --candidates child.table.txt \                              │
│         --output parent_analysis.json \                             │
│         --denovo                                                    │
└─────────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────────┐
│  3. Results Evaluation                                              │
│                                                                     │
│     Categories:                                                     │
│     - PASS_DENOVO: Meets all criteria for de novo                   │
│     - UNCERTAIN: Parent coverage insufficient                       │
│     - FAIL: Parent shows alt support                                │
└─────────────────────────────────────────────────────────────────────┘
```

### Child Analysis

```bash
# Run tldr on child sample with --denovo flag
tldr -b child.bam \
     -e ref/teref.ont.human.fa \
     -r genome.fasta \
     --denovo \
     -p 8 \
     -o child_output
```

### Parent Support Analysis

```bash
# Analyze parent BAM files for supporting evidence
python scripts/parent_support.py \
    --mom_bam mother.bam \
    --dad_bam father.bam \
    --candidates child_output.table.txt \
    --output parent_analysis.json \
    --denovo

# With optional TE-enhanced soft-clip alignment
python scripts/parent_support.py \
    --mom_bam mother.bam \
    --dad_bam father.bam \
    --candidates child_output.table.txt \
    --output parent_analysis.json \
    --denovo \
    --te_enhance
```

### Evaluation Criteria

The de novo evaluation uses the following thresholds:

- **Minimum parent depth**: `max(8, 0.25 * median_parent_depth)`
  - Ensures adequate coverage in both parents
  - Floor of 8 reads, scaled by 25% of median depth
- **Maximum parent alt reads**: 0 (strict)
  - No soft-clipped or insertion CIGAR reads allowed near breakpoints
- **Breakpoint precision**: Not explicitly validated by parent_support.py
  - Child-side precision controlled by `--wiggle` parameter

### Output Fields

Parent analysis results include:

| Field | Description |
|-------|-------------|
| `uuid` | Candidate identifier |
| `chrom` | Chromosome |
| `bp_left`, `bp_right` | Breakpoint coordinates |
| `te_family` | TE family |
| `mom_depth_mean` | Median depth in mother |
| `dad_depth_mean` | Median depth in father |
| `mom_softclip` | Soft-clip reads in mother |
| `dad_softclip` | Soft-clip reads in father |
| `mom_insertion` | Insertion CIGAR reads in mother |
| `dad_insertion` | Insertion CIGAR reads in father |
| `mom_total_alt` | Total alt support in mother |
| `dad_total_alt` | Total alt support in father |
| `evaluation` | PASS_DENOVO / UNCERTAIN / FAIL |
| `reasons` | Detailed evaluation reasons |

## Methylation

Non-reference methylation can be assessed through the use of scripts located in the `scripts/` directory:

| script                | description |
|-----------------------|-------------|
| tldr_callmeth.sh      | Must be run from within the diretory where `nanopolish index` was run to index a .fastq file against a set of ONT .fast5 files. Takes as input a .fastq (indexed via `nanopolish index`), an output directory created via the `--detail_output` option, a UUID and a sample name. Creates a tabix indexed table from the output of nanopolish call-methylation on the sample+uuid combination. Can be automated via xargs or GNU parallel. |
| tablemeth_nonref.py   | Creates a table with per-element mCpG summary data given a tldr output table and the directory created by `--detail_output`. Only considers element + sample combinations from the tldr table where `tldr_callmeth.sh` has been run. Requires pysam, pandas, numpy, and scipy. |
| plotmeth_nonref.py    | Makes a plot of a TE (requires running `tldr_callmeth.sh` first) plus the surrounding region if `--extend_consensus` is specified. Tracks include translation to CpG space, raw log-likelihood, and smoothed methylation fraction. Requires pysam, pandas, numpy, scipy, matplotlib, and seaborn. |
| parent_support.py     | Parent support analysis for de novo TE detection. Analyzes Mom/Dad BAM files to compute depth and alt support for child candidate sites identified with `tldr --denovo`. Outputs JSON with PASS_DENOVO/UNCERTAIN/FAIL classification. Requires pysam and numpy. |

## Reference TEs

See https://github.com/adamewing/te-nanopore-tools

## References

Adam D. Ewing, Nathan Smits, Francisco J. Sanchez-Luque, Sandra R. Richardson, Seth W. Cheetham, Geoffrey J. Faulkner. Nanopore Sequencing Enables Comprehensive Transposable Element Epigenomic Profiling. 2020. Molecular Cell, Online ahead of print: https://doi.org/10.1016/j.molcel.2020.10.024

## Getting help

Reporting [issues](https://github.com/adamewing/tldr/issues) and questions through github is preferred versus e-mail.
