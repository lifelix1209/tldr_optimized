#!/usr/bin/env python3
"""
Create mock BAM files for testing de novo TE insertion detection.

This script creates BAM files with soft-clipped reads that tldr can detect.
It simulates a more realistic scenario with:

1. Different coverage depths for child, mom, and dad
2. Various TE insertion types and inheritance patterns
3. Mosaic/germline variations
4. Different read support levels

Scenario (W trio inheritance model):
- Reference genome: chr1 with 50kb sequence
- Multiple TE insertion sites with different inheritance patterns:

| Site | Position | Type | Mom | Dad | Child | Evaluation | Coverage |
|------|----------|------|-----|-----|-------|------------|----------|
| 1 | 5,000 | L1 | Yes | No | Yes | PASS (inherited) | High |
| 2 | 15,000 | L1 | No | Yes | Yes | PASS (inherited) | Medium |
| 3 | 25,000 | L1 | No | No | Yes | PASS_DENOVO | Low |
| 4 | 35,000 | Alu | Yes | No | No | FAIL (private) | High |
| 5 | 40,000 | Alu | No | Yes | No | FAIL (private) | Medium |
| 6 | 10,000 | L1 | Mosaic | No | Mosaic | UNCERTAIN | Low |
| 7 | 20,000 | Alu | No | Mosaic | Mosaic | UNCERTAIN | Low |

TE Insertion Detection by tldr:
- tldr looks for CIGAR operation 4 (soft-clip) with length > min_te_len (~100bp)
- Soft-clipped portion represents the TE sequence that couldn't align to reference
- Pattern: "X S" means X bp match + S bp soft-clipped (TE sequence)
- Or: "S X" means S bp soft-clipped + X bp match
"""

import pysam
import random
import json

# Set random seed for reproducibility
random.seed(42)

# Directory paths
REF_DIR = "ref"
BAM_DIR = "bams"

# TE detection parameters (matching tldr defaults)
MIN_TE_LEN = 100  # Minimum TE length for tldr detection
TE_CLIP_LEN = 150  # Length of soft-clipped TE sequence

# Define TE sequences for different families
TE_SEQUENCES = {
    "L1HS": "GGGGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCCCTGAGATATATACTTCCTTCGAAATGAGGGGACCTTCCCTGCAAGAAGAGAGGTACAGTAAGCCAGAGCCCTTGGATAAGTTTACGCACGTGCTATTCTGTCTCAAATCGCGTACTGTTCCCTTTCTAATGCGAAGCGCGGGCTCAAATACCGTTATTTCCCACAGCCGGGGTTGCTATGCGGGGACCCTCCGAGATCGATTGACCAGGACCGTTGAACCGCTGTCCTTGGATAGCCTCGAAAGCCTGATAGCCAGAGGGTAAGATGTCTACAGGGGTACGTACCCACGCAGCGGGGAGGTCGTTCGCGTCAGGGCAGGTAAGATTTCCGGAACTCGCAGCGGAACCGTTGCGTGCTCGTCGCGAGTTCTCATGTGCGTTGCACCACTACATTTCACGATGGTATTTGAGCGGATGGACAGGCGCTCACCCAACAGGACCAGCAAGCCACGAACAGGTGGCTGGGTTCATTCCGGCAATGTCGGTTCACCCATGAATTCGTGCCGGTCCTGTCTGAGATGCGCGTATAATTCCGATTTGAAAGCGAGCGTCTCTCAATCACCGAGTCTTAGCCCCGGCACTTATGAAAGCCTAGCCTGCAGGCCAGGTAGCGCTGGACGAGTTTCATACGTAGTCACTCGCGCTTCATGGGACAATTGCAGTCACTAAAGACTCGCTTCAATACTTACGAGCGTAACCGTGCGCATATCGTGATGTCTTCCAGCTGGGACGAGGTAAACACAAATGTTGGATAGTGTTGCCTGTCCTTCAATCCATCCCTGACATAGTGAAGCCTAGTGCGTGGTCCTTTGATTCGTCGCGTAGTACTGAACAGTGCAGGAAGGCAGCGTCTAATAATGCG",
    "AluY": "GGGGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCCCTGAGATATATACTTCCTTCGAAATGAGGGGACCTTCCCTGCAAGAAGAGAGGTACAGTAAGCCAGAGCCCTTGGATAAGTTTACGCACGTGCTATTCTGTCTCAAATCGCGTACTGTTCCCTTTCTAATGCGAAGCGCGGGCTCAAATACCGTTATTTCCCACAGCCGGGGTTGCTATGCGGGGACCCTCCGAGATCGATTGACCAGGACCGTTGAACCGCTGTCCTTGGATAGCCTCGAAAGCCTGATAGCCAGAGGGTAAGATGTCTACAGGGGTACGTACCCACGCAGCGGGGAGGTCGTTCGCGTCAGGGCAGGTAAGATTTCCGGAACTCGCAGCGGAACCGTTGCGTGCTCGTCGCGAGTTCTCATGTGCGTTGCACCACTACATTTCACGATGGTATTTGAGCGGATGGACAGGCGCTCACCCAACAGGACCAGCAAGCCACGAACAGGTGGCTGGGTTCATTCCGGCAATGTCGGTTCACCCATGAATTCGTGCCGGTCCTGTCTGAGATGCGCGTATAATTCCGATTTGAAAGCGAGCGTCTCTCAATCACCGAGTCTTAGCCCCGGCACTTATGAAAGCCTAGCCTGCAGGCCAGGTAGCGCTGGACGAGTTTCATACGTAGTCACTCGCGCTTCATGGGACAATTGCAGTCACTAAAGACTCGCTTCAATACTTACGAGCGTAACCGTGCGCATATCGTGATGTCTTCCAGCTGGGACGAGGTAAACACAAATGTTGGATAGTGTTGCCTGTCCTTCAATCCATCCCTGACATAGTGAAGCCTAGTGCGTGGTCCTTTGATTCGTCGCGTAGTACTGAACAGTGCAGGAAGGCAGCGTCTAATAATGCG",
    "SVA": "GGGGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCCCTGAGATATATACTTCCTTCGAAATGAGGGGACCTTCCCTGCAAGAAGAGAGGTACAGTAAGCCAGAGCCCTTGGATAAGTTTACGCACGTGCTATTCTGTCTCAAATCGCGTACTGTTCCCTTTCTAATGCGAAGCGCGGGCTCAAATACCGTTATTTCCCACAGCCGGGGTTGCTATGCGGGGACCCTCCGAGATCGATTGACCAGGACCGTTGAACCGCTGTCCTTGGATAGCCTCGAAAGCCTGATAGCCAGAGGGTAAGATGTCTACAGGGGTACGTACCCACGCAGCGGGGAGGTCGTTCGCGTCAGGGCAGGTAAGATTTCCGGAACTCGCAGCGGAACCGTTGCGTGCTCGTCGCGAGTTCTCATGTGCGTTGCACCACTACATTTCACGATGGTATTTGAGCGGATGGACAGGCGCTCACCCAACAGGACCAGCAAGCCACGAACAGGTGGCTGGGTTCATTCCGGCAATGTCGGTTCACCCATGAATTCGTGCCGGTCCTGTCTGAGATGCGCGTATAATTCCGATTTGAAAGCGAGCGTCTCTCAATCACCGAGTCTTAGCCCCGGCACTTATGAAAGCCTAGCCTGCAGGCCAGGTAGCGCTGGACGAGTTTCATACGTAGTCACTCGCGCTTCATGGGACAATTGCAGTCACTAAAGACTCGCTTCAATACTTACGAGCGTAACCGTGCGCATATCGTGATGTCTTCCAGCTGGGACGAGGTAAACACAAATGTTGGATAGTGTTGCCTGTCCTTCAATCCATCCCTGACATAGTGAAGCCTAGTGCGTGGTCCTTTGATTCGTCGCGTAGTACTGAACAGTGCAGGAAGGCAGCGTCTAATAATGCG"
}

# Define insertion sites with inheritance patterns
# Each site: (position, length, te_family, mom_status, dad_status, child_status, coverage)
# status: 0=no, 1=yes, 2=mosaic
# coverage: "high" (5 reads), "medium" (3 reads), "low" (1-2 reads)
INSERTION_SITES = [
    # Standard inheritance
    {"pos": 5000, "len": 300, "te_family": "L1HS", "mom_has": 1, "dad_has": 0, "child_has": 1, "coverage": "high", "name": "inherited_from_mom"},
    {"pos": 15000, "len": 300, "te_family": "L1HS", "mom_has": 0, "dad_has": 1, "child_has": 1, "coverage": "medium", "name": "inherited_from_dad"},

    # De novo insertion
    {"pos": 25000, "len": 300, "te_family": "L1HS", "mom_has": 0, "dad_has": 0, "child_has": 1, "coverage": "low", "name": "denovo"},

    # Private to parent (not inherited by child)
    {"pos": 35000, "len": 200, "te_family": "AluY", "mom_has": 1, "dad_has": 0, "child_has": 0, "coverage": "high", "name": "mom_private"},
    {"pos": 40000, "len": 200, "te_family": "AluY", "mom_has": 0, "dad_has": 1, "child_has": 0, "coverage": "medium", "name": "dad_private"},

    # Mosaic cases
    {"pos": 10000, "len": 250, "te_family": "L1HS", "mom_has": 2, "dad_has": 0, "child_has": 2, "coverage": "low", "name": "mosaic_mom_to_child"},
    {"pos": 20000, "len": 200, "te_family": "AluY", "mom_has": 0, "dad_has": 2, "child_has": 2, "coverage": "low", "name": "mosaic_dad_to_child"},

    # Low support de novo (near detection limit)
    {"pos": 30000, "len": 250, "te_family": "L1HS", "mom_has": 0, "dad_has": 0, "child_has": 1, "coverage": "very_low", "name": "denovo_low_support"},
]

# Coverage settings per sample type
SAMPLE_COVERAGE = {
    "child": {"softclip": 5, "normal": 10, "span": 3},
    "mom": {"softclip": 4, "normal": 8, "span": 2},
    "dad": {"softclip": 4, "normal": 8, "span": 2},
}

# Mapping of coverage level to read count
COVERAGE_LEVELS = {
    "very_low": 1,
    "low": 2,
    "medium": 3,
    "high": 5,
}


def create_reference():
    """Create mock reference genome."""
    ref_seq = "N" * 50000  # 50kb of N's
    # Add some recognizable sequences
    ref_seq = ref_seq[:1000] + "TGCATGCATGCATGCATGCATGCATGCATGCATGCAT" + ref_seq[1038:]
    ref_seq = ref_seq[:5000] + "TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCAT" + ref_seq[5100:]
    ref_seq = ref_seq[:15000] + "TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCAT" + ref_seq[15300:]
    ref_seq = ref_seq[:25000] + "TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCAT" + ref_seq[25300:]

    ref_file = f"{REF_DIR}/mock_ref.fa"
    with open(ref_file, 'w') as f:
        f.write(">chr1\n")
        for i in range(0, len(ref_seq), 80):
            f.write(ref_seq[i:i+80] + "\n")

    # Create index
    pysam.faidx(ref_file)
    print(f"Created reference: {ref_file}")
    return ref_seq


def create_te_library():
    """Create mock TE reference library."""
    te_file = f"{REF_DIR}/mock_te.fa"
    with open(te_file, 'w') as f:
        for family, seq in TE_SEQUENCES.items():
            f.write(f">{family}:{family}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + "\n")
    print(f"Created TE library: {te_file}")
    return TE_SEQUENCES


def create_read_no_insertion(ref_seq, read_start, read_end, sample_name, mapq=60):
    """Create a normal read without TE insertion."""
    a = pysam.AlignedSegment()
    actual_len = read_end - read_start
    a.query_name = f"{sample_name}_normal_{read_start:05d}"
    a.flag = 0
    a.reference_id = 0
    a.mapping_quality = mapq
    a.reference_start = read_start
    a.cigartuples = [(0, actual_len)]
    a.query_sequence = ref_seq[read_start:read_end]
    a.query_qualities = pysam.qualitystring_to_array("I" * len(a.query_sequence))
    return a


def create_softclip_read_left(ref_seq, te_seq, read_start, read_end, sample_name, clip_len=None):
    """Create a soft-clipped read with TE sequence at the LEFT end (S M pattern)."""
    a = pysam.AlignedSegment()
    clip_len = clip_len or TE_CLIP_LEN
    a.query_name = f"{sample_name}_L{clip_len}_{read_start:05d}"
    a.flag = 0
    a.reference_id = 0
    a.mapping_quality = 60
    a.reference_start = read_start
    query_seq = te_seq[:clip_len] + ref_seq[read_start:read_end]
    a.query_sequence = query_seq
    a.query_qualities = pysam.qualitystring_to_array("I" * len(query_seq))
    match_len = read_end - read_start
    a.cigartuples = [(4, clip_len), (0, match_len)]  # S, M
    return a


def create_softclip_read_right(ref_seq, te_seq, read_start, read_end, sample_name, clip_len=None):
    """Create a soft-clipped read with TE sequence at the RIGHT end (M S pattern)."""
    a = pysam.AlignedSegment()
    clip_len = clip_len or TE_CLIP_LEN
    a.query_name = f"{sample_name}_R{clip_len}_{read_start:05d}"
    a.flag = 0
    a.reference_id = 0
    a.mapping_quality = 60
    a.reference_start = read_start
    query_seq = ref_seq[read_start:read_end] + te_seq[:clip_len]
    a.query_sequence = query_seq
    a.query_qualities = pysam.qualitystring_to_array("I" * len(query_seq))
    match_len = read_end - read_start
    a.cigartuples = [(0, match_len), (4, clip_len)]  # M, S
    return a


def create_span_read(ref_seq, te_seq, read_start, read_end, sample_name, ins_pos_offset=0):
    """Create a spanning read that crosses the insertion point (M I M pattern)."""
    a = pysam.AlignedSegment()
    a.query_name = f"{sample_name}_span_{read_start:05d}"
    a.flag = 0
    a.reference_id = 0
    a.mapping_quality = 60
    a.reference_start = read_start
    ins_len = 300
    ins_in_read = 200 - ins_pos_offset
    pre_genomic = ref_seq[read_start:read_start + ins_in_read]
    post_genomic = ref_seq[read_start + ins_in_read:read_end]
    query_seq = pre_genomic + te_seq[:ins_len] + post_genomic
    a.query_sequence = query_seq
    a.query_qualities = pysam.qualitystring_to_array("I" * len(query_seq))
    pre_ins = ins_in_read
    post_ins = len(post_genomic)
    a.cigartuples = [(0, pre_ins), (1, ins_len), (0, post_ins)]  # M, I, M
    return a


def get_sample_coverage(sample_name, coverage_type="softclip"):
    """Get the number of reads to generate for a sample based on coverage type."""
    base_coverage = SAMPLE_COVERAGE.get(sample_name, SAMPLE_COVERAGE["child"])
    return base_coverage.get(coverage_type, 3)


def create_bam_file(output_bam, ref_seq, te_sequences, sample_name, has_insertions_at, mosaic_prob=0.5):
    """
    Create a BAM file for a sample with soft-clipped reads for TE detection.

    Args:
        output_bam: output BAM path
        ref_seq: reference sequence
        te_sequences: dict of TE family -> sequence
        sample_name: 'child', 'mom', or 'dad'
        has_insertions_at: list of site names this sample has insertions at
        mosaic_prob: probability of mosaic site being present in sample
    """
    header = {
        'HD': {'VN': '1.0'},
        'SQ': [{'LN': len(ref_seq), 'SN': 'chr1'}]
    }

    reads = []
    clip_reads_count = 0
    span_reads_count = 0
    normal_reads_count = 0

    for site in INSERTION_SITES:
        pos = site["pos"]
        site_name = site["name"]
        te_family = site["te_family"]
        te_seq = te_sequences.get(te_family, list(te_sequences.values())[0])
        has_te = site_name in has_insertions_at

        # For mosaic sites, determine if this sample has the insertion
        is_mosaic = site.get("mom_has", 0) == 2 or site.get("dad_has", 0) == 2 or site.get("child_has", 0) == 2
        if is_mosaic and has_te:
            has_te = random.random() < mosaic_prob

        # Get coverage multiplier based on site coverage level
        coverage_level = site.get("coverage", "medium")
        coverage_mult = COVERAGE_LEVELS.get(coverage_level, 1)

        # Get sample-specific coverage multiplier
        sample_mult = get_sample_coverage(sample_name, "softclip") / 3.0  # normalize to medium

        # Calculate number of reads to generate
        num_softclip_reads = max(1, int(3 * coverage_mult * sample_mult))
        num_span_reads = max(1, int(2 * coverage_mult * sample_mult))

        # Generate soft-clipped and span reads for this site
        for i in range(num_softclip_reads):
            # Alternate between left and right soft-clip
            if i % 2 == 0:
                read_start = pos - 200
                read_end = read_start + 750
                if has_te:
                    read = create_softclip_read_left(ref_seq, te_seq, read_start, read_end, sample_name)
                    clip_reads_count += 1
                else:
                    read = create_read_no_insertion(ref_seq, read_start, read_end, sample_name)
                    normal_reads_count += 1
            else:
                read_start = pos
                read_end = read_start + 750
                if has_te:
                    read = create_softclip_read_right(ref_seq, te_seq, read_start, read_end, sample_name)
                    clip_reads_count += 1
                else:
                    read = create_read_no_insertion(ref_seq, read_start, read_end, sample_name)
                    normal_reads_count += 1
            reads.append(read)

        # Generate span reads
        for i in range(num_span_reads):
            read_start = pos - 200
            read_end = read_start + 500
            if has_te:
                read = create_span_read(ref_seq, te_seq, read_start, read_end, sample_name)
                span_reads_count += 1
            else:
                read = create_read_no_insertion(ref_seq, read_start, read_end, sample_name)
                normal_reads_count += 1
            reads.append(read)

        # Generate normal coverage reads around the site
        for i in range(2):
            if i == 0:
                read_start = max(0, pos - 1500)
            else:
                read_start = pos + 100
            read_end = min(read_start + 400, len(ref_seq))
            if read_end > read_start:
                read = create_read_no_insertion(ref_seq, read_start, read_end, sample_name)
                reads.append(read)
                normal_reads_count += 1

    # Add background coverage reads
    num_bg_reads = get_sample_coverage(sample_name, "normal")
    for i in range(num_bg_reads):
        read_start = random.randint(0, len(ref_seq) - 500)
        read_end = min(read_start + random.randint(200, 500), len(ref_seq))
        read = create_read_no_insertion(ref_seq, read_start, read_end, sample_name)
        reads.append(read)
        normal_reads_count += 1

    # Sort reads by position
    reads.sort(key=lambda x: x.reference_start)

    # Write to BAM
    with pysam.AlignmentFile(output_bam, "wb", header=header) as outf:
        for read in reads:
            outf.write(read)

    print(f"Created BAM file: {output_bam}")
    print(f"  Total reads: {len(reads)}")
    print(f"  Soft-clip reads: {clip_reads_count}")
    print(f"  Span reads: {span_reads_count}")
    print(f"  Normal reads: {normal_reads_count}")


def index_bam(bam_file):
    """Index BAM file."""
    pysam.index(bam_file)


def generate_tldr_output():
    """Generate tldr-compatible output table for child insertions."""
    output_file = f"{BAM_DIR}/child.table.txt"

    header = [
        "UUID", "Chrom", "Start", "End", "Strand", "Family", "Subfamily",
        "StartTE", "EndTE", "LengthIns", "Inversion", "UnmapCover",
        "MedianMapQ", "TEMatch", "UsedReads", "SpanReads", "NumSamples",
        "SampleReads", "EmptyReads", "EmptyPhase", "NonRef", "TSD",
        "Consensus", "Phasing", "Remappable", "Filter",
        "bp_left", "bp_right", "wiggle", "TE_family", "child_support"
    ]

    insertions = []
    for site in INSERTION_SITES:
        if site["child_has"]:
            # Determine expected evaluation
            if site["mom_has"] == 1 and site["dad_has"] == 0:
                expected = "Inherited from mom"
                evaluation = "PASS"
            elif site["mom_has"] == 0 and site["dad_has"] == 1:
                expected = "Inherited from dad"
                evaluation = "PASS"
            elif site["child_has"] == 2:
                expected = "Mosaic in child"
                evaluation = "UNCERTAIN"
            elif site["mom_has"] == 0 and site["dad_has"] == 0:
                expected = "De novo"
                evaluation = "PASS_DENOVO"
            else:
                expected = "Unknown"
                evaluation = "UNCERTAIN"

            ins = {
                "UUID": f"ins_{site['name']}",
                "Chrom": "chr1",
                "Start": site["pos"] - 100,
                "End": site["pos"] + site["len"] + 100,
                "Strand": "+",
                "Family": "LINE" if "L1" in site["te_family"] else ("SINE" if "Alu" in site["te_family"] else "SVA"),
                "Subfamily": site["te_family"],
                "StartTE": "1",
                "EndTE": str(site["len"]),
                "LengthIns": str(site["len"]),
                "Inversion": "N",
                "UnmapCover": "0.95",
                "MedianMapQ": "60",
                "TEMatch": "0.98",
                "UsedReads": "5",
                "SpanReads": "5",
                "NumSamples": "1",
                "SampleReads": "child|5",
                "EmptyReads": "0",
                "EmptyPhase": "NA",
                "NonRef": "NA",
                "TSD": "TGAC",
                "Consensus": "TGAC" + "N" * (site["len"] - 8) + "TGAC",
                "Phasing": "NA",
                "Remappable": "True",
                "Filter": "PASS",
                "bp_left": site["pos"],
                "bp_right": site["pos"] + site["len"],
                "wiggle": "50",
                "TE_family": site["te_family"],
                "child_support": f"useable:5,embedded:5 ({expected})"
            }
            insertions.append(ins)

    with open(output_file, 'w') as f:
        f.write('\t'.join(header) + '\n')
        for ins in insertions:
            row = [str(ins[col]) for col in header]
            f.write('\t'.join(row) + '\n')

    print(f"\nGenerated tldr output: {output_file}")
    print(f"  - {len(insertions)} candidate insertions")
    return output_file


def generate_candidates_json():
    """Generate JSON output for parent_analysis script."""
    candidates = []

    for site in INSERTION_SITES:
        if site["child_has"]:
            # Determine expected evaluation
            if site["mom_has"] == 1 and site["dad_has"] == 0:
                evaluation = "PASS"
                support = "Child and Mom have supporting reads, Dad does not"
            elif site["mom_has"] == 0 and site["dad_has"] == 1:
                evaluation = "PASS"
                support = "Child and Dad have supporting reads, Mom does not"
            elif site["child_has"] == 2:
                evaluation = "UNCERTAIN"
                support = "Mosaic insertion in child"
            elif site["mom_has"] == 0 and site["dad_has"] == 0:
                evaluation = "PASS_DENOVO"
                support = "Only Child has supporting reads, both parents lack insertion"
            else:
                evaluation = "UNCERTAIN"
                support = "Unclear inheritance pattern"

            candidate = {
                "uuid": f"ins_{site['name']}",
                "Chrom": "chr1",
                "bp_left": site["pos"],
                "bp_right": site["pos"] + site["len"],
                "TE_family": site["te_family"],
                "Strand": "+",
                "Family": "LINE" if "L1" in site["te_family"] else ("SINE" if "Alu" in site["te_family"] else "SVA"),
                "Subfamily": site["te_family"],
                "evaluation": evaluation,
                "child_support": support,
                "coverage": site.get("coverage", "medium"),
                "inheritance": "inherited_from_mom" if site["mom_has"] == 1 and site["dad_has"] == 0 else
                             ("inherited_from_dad" if site["mom_has"] == 0 and site["dad_has"] == 1 else
                             ("de_novo" if site["mom_has"] == 0 and site["dad_has"] == 0 and site["child_has"] == 1 else
                             ("mosaic" if site["child_has"] == 2 else "unknown")))
            }
            candidates.append(candidate)

    return candidates


def main():
    """Main function to generate mock data."""
    print("=" * 70)
    print("Creating mock BAM files for TE insertion detection")
    print("=" * 70)

    # Create directories
    os.makedirs(REF_DIR, exist_ok=True)
    os.makedirs(BAM_DIR, exist_ok=True)

    # Create reference and TE library
    ref_seq = create_reference()
    te_sequences = create_te_library()

    # Determine which insertions each sample has
    child_insertions = [site["name"] for site in INSERTION_SITES if site["child_has"]]
    mom_insertions = [site["name"] for site in INSERTION_SITES if site["mom_has"] in (1, 2)]
    dad_insertions = [site["name"] for site in INSERTION_SITES if site["dad_has"] in (1, 2)]

    print(f"\nSample configurations:")
    print(f"  Child has insertions at: {child_insertions}")
    print(f"  Mom has insertions at: {mom_insertions}")
    print(f"  Dad has insertions at: {dad_insertions}")

    # Create BAM files for each sample
    print("\n" + "-" * 70)
    print("Generating BAM files:")
    print("-" * 70)

    child_bam = f"{BAM_DIR}/child.bam"
    create_bam_file(child_bam, ref_seq, te_sequences, "child", child_insertions)
    index_bam(child_bam)

    mom_bam = f"{BAM_DIR}/mom.bam"
    create_bam_file(mom_bam, ref_seq, te_sequences, "mom", mom_insertions)
    index_bam(mom_bam)

    dad_bam = f"{BAM_DIR}/dad.bam"
    create_bam_file(dad_bam, ref_seq, te_sequences, "dad", dad_insertions)
    index_bam(dad_bam)

    # Generate tldr-compatible output
    generate_tldr_output()

    # Generate JSON output for parent_analysis
    candidates = generate_candidates_json()
    candidates_file = f"{BAM_DIR}/expected_candidates.json"
    with open(candidates_file, 'w') as f:
        json.dump(candidates, f, indent=2)
    print(f"\nGenerated expected candidates: {candidates_file}")

    print("\n" + "=" * 70)
    print("Mock data creation complete!")
    print("=" * 70)
    print("\nNext steps:")
    print("  1. Run tldr on child sample:")
    print(f"     python3 tldr/tldr -b {child_bam} -e {REF_DIR}/mock_te.fa -r {REF_DIR}/mock_ref.fa \\")
    print(f"       --denovo -m 2 --min_te_len 100 -o {BAM_DIR}/tldr_output")
    print(f"\n  2. Run parent support analysis:")
    print(f"     python scripts/parent_support.py \\")
    print(f"       --mom_bam {mom_bam} \\")
    print(f"       --dad_bam {dad_bam} \\")
    print(f"       --child_bam {child_bam} \\")
    print(f"       --candidates {BAM_DIR}/tldr_output.table.txt \\")
    print(f"       --output {BAM_DIR}/parent_analysis.json \\")
    print(f"       --denovo --visualize")


if __name__ == "__main__":
    import os
    main()
