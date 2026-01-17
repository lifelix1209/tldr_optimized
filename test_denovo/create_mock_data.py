#!/usr/bin/env python3
"""
Create mock BAM files for testing de novo TE insertion detection.

Scenario:
- Reference genome: chr1 with 50kb sequence
- TE insertion: L1 element (300bp) inserted at position 25000 in child only
- Child BAM: 20 reads covering the insertion (10 fully embedding, 10 partially)
- Mom BAM: 20 reads covering the region without insertion
- Dad BAM: 20 reads covering the region without insertion
"""

import pysam
import random
import os

# Set random seed for reproducibility
random.seed(42)

# Directory paths
REF_DIR = "ref"
BAM_DIR = "bams"

# Create reference genome
def create_reference_genome():
    """Create a simple reference genome."""
    ref_file = f"{REF_DIR}/mock_ref.fa"

    # Create 50kb random sequence
    bases = ['A', 'C', 'G', 'T']
    seq = ''.join(random.choices(bases, k=50000))

    with open(ref_file, 'w') as f:
        f.write(">chr1\n")
        # Write in 80 character lines
        for i in range(0, len(seq), 80):
            f.write(seq[i:i+80] + "\n")

    print(f"Created reference genome: {ref_file}")
    return ref_file, seq

# Create TE reference library
def create_te_reference():
    """Create a mock TE reference library."""
    te_file = f"{REF_DIR}/mock_te.fa"

    # Create a 300bp L1-like sequence
    bases = ['A', 'C', 'G', 'T']
    # Add some GC-rich regions typical of L1
    te_seq = 'G' * 20 + 'C' * 20  # GC-rich start
    te_seq += ''.join(random.choices(bases, weights=[0.3, 0.2, 0.2, 0.3], k=260))

    with open(te_file, 'w') as f:
        f.write(">LINE:L1HS\n")
        for i in range(0, len(te_seq), 80):
            f.write(te_seq[i:i+80] + "\n")

    print(f"Created TE reference: {te_file}")
    return te_file, te_seq

# Helper function to create reads
def generate_read_sequence(ref_seq, start, end, with_insertion=False, te_seq=None, ins_pos=None):
    """Generate a read sequence from reference, optionally with TE insertion."""
    if not with_insertion:
        return ref_seq[start:end]

    # Calculate where insertion occurs within the read
    if start <= ins_pos <= end:
        # Insertion is within this read
        left_flank = ref_seq[start:ins_pos]
        right_flank = ref_seq[ins_pos:end]
        return left_flank + te_seq + right_flank
    else:
        return ref_seq[start:end]

def create_cigar_with_insertion(read_len, ins_pos_in_read, ins_len):
    """Create CIGAR string with insertion."""
    if ins_pos_in_read <= 0:
        return [(1, ins_len), (0, read_len)]  # Insertion at start, then match
    elif ins_pos_in_read >= read_len:
        return [(0, read_len), (1, ins_len)]  # Match, then insertion at end
    else:
        return [(0, ins_pos_in_read), (1, ins_len), (0, read_len - ins_pos_in_read)]

# Create BAM file
def create_bam_file(output_bam, ref_seq, te_seq, insertion_pos, has_insertion=False, num_reads=20):
    """Create a BAM file with or without TE insertion."""

    # Create header
    header = {
        'HD': {'VN': '1.0'},
        'SQ': [{'LN': len(ref_seq), 'SN': 'chr1'}]
    }

    # Store reads in a list first, then sort by position
    reads = []

    # Generate reads covering the insertion region
    # Region: insertion_pos - 5000 to insertion_pos + 5000
    start_region = max(0, insertion_pos - 5000)
    end_region = min(len(ref_seq), insertion_pos + 5000)

    for i in range(num_reads):
            # Create read
            a = pysam.AlignedSegment()
            a.query_name = f"read_{i:04d}"
            a.flag = 0  # Forward strand, primary alignment
            a.reference_id = 0  # chr1
            a.mapping_quality = 60

            if has_insertion:
                # Half the reads fully embed the insertion
                # Half partially cover it
                if i < num_reads // 2:
                    # Fully embedding reads
                    read_start = random.randint(insertion_pos - 1000, insertion_pos - 100)
                    read_len = 3000  # Long enough to cover insertion

                    # Calculate position of insertion in read
                    ins_pos_in_read = insertion_pos - read_start

                    # Create sequence with insertion
                    left_flank = ref_seq[read_start:insertion_pos]
                    right_flank = ref_seq[insertion_pos:read_start + read_len]
                    query_seq = left_flank + te_seq + right_flank

                    # Create CIGAR: Match left, Insertion, Match right
                    a.cigartuples = [
                        (0, len(left_flank)),      # Match
                        (1, len(te_seq)),           # Insertion
                        (0, len(right_flank))       # Match
                    ]

                    a.reference_start = read_start
                    a.query_sequence = query_seq
                    a.query_qualities = pysam.qualitystring_to_array("I" * len(query_seq))

                else:
                    # Partially covering reads (soft-clipped at insertion site)
                    if random.random() < 0.5:
                        # Left-side read with right soft-clip
                        read_start = random.randint(insertion_pos - 1000, insertion_pos - 100)
                        match_len = insertion_pos - read_start
                        clip_len = random.randint(150, 250)  # Part of TE sequence

                        query_seq = ref_seq[read_start:insertion_pos] + te_seq[:clip_len]

                        a.cigartuples = [
                            (0, match_len),    # Match
                            (4, clip_len)      # Soft clip (TE sequence)
                        ]

                        a.reference_start = read_start
                    else:
                        # Right-side read with left soft-clip
                        match_start = insertion_pos
                        match_len = random.randint(1000, 1500)
                        clip_len = random.randint(150, 250)  # Part of TE sequence

                        query_seq = te_seq[-clip_len:] + ref_seq[match_start:match_start + match_len]

                        a.cigartuples = [
                            (4, clip_len),     # Soft clip (TE sequence)
                            (0, match_len)     # Match
                        ]

                        a.reference_start = match_start

                    a.query_sequence = query_seq
                    a.query_qualities = pysam.qualitystring_to_array("I" * len(query_seq))

            else:
                # No insertion - normal reads spanning the region
                read_start = random.randint(start_region, insertion_pos - 500)
                read_len = random.randint(2000, 4000)
                read_end = min(read_start + read_len, len(ref_seq))
                actual_len = read_end - read_start

                query_seq = ref_seq[read_start:read_end]

                a.cigartuples = [(0, actual_len)]  # All match
                a.reference_start = read_start
                a.query_sequence = query_seq
                a.query_qualities = pysam.qualitystring_to_array("I" * len(query_seq))

            reads.append(a)

    # Sort reads by reference position
    reads.sort(key=lambda x: x.reference_start)

    # Write sorted reads to BAM
    with pysam.AlignmentFile(output_bam, "wb", header=header) as outf:
        for read in reads:
            outf.write(read)

    print(f"Created BAM file: {output_bam}")

def index_bam(bam_file):
    """Index BAM file."""
    pysam.index(bam_file)
    print(f"Indexed: {bam_file}")

def index_fasta(fasta_file):
    """Index FASTA file."""
    pysam.faidx(fasta_file)
    print(f"Indexed: {fasta_file}")

def main():
    print("Creating mock data for de novo TE insertion testing...")
    print("=" * 60)

    # Create reference genome and TE sequence
    ref_file, ref_seq = create_reference_genome()
    te_file, te_seq = create_te_reference()

    # Index reference files
    index_fasta(ref_file)

    # Insertion position in reference (middle of chr1)
    INSERTION_POS = 25000

    print(f"\nInsertion scenario:")
    print(f"  - Position: chr1:{INSERTION_POS}")
    print(f"  - TE type: L1HS")
    print(f"  - TE length: {len(te_seq)}bp")
    print(f"  - Present in: Child only (de novo)")
    print()

    # Create child BAM with insertion
    child_bam = f"{BAM_DIR}/child.bam"
    create_bam_file(child_bam, ref_seq, te_seq, INSERTION_POS, has_insertion=True, num_reads=20)
    index_bam(child_bam)

    # Create mom BAM without insertion
    mom_bam = f"{BAM_DIR}/mom.bam"
    create_bam_file(mom_bam, ref_seq, te_seq, INSERTION_POS, has_insertion=False, num_reads=20)
    index_bam(mom_bam)

    # Create dad BAM without insertion
    dad_bam = f"{BAM_DIR}/dad.bam"
    create_bam_file(dad_bam, ref_seq, te_seq, INSERTION_POS, has_insertion=False, num_reads=20)
    index_bam(dad_bam)

    print("\n" + "=" * 60)
    print("Mock data creation complete!")
    print("\nCreated files:")
    print(f"  Reference: {ref_file}")
    print(f"  TE library: {te_file}")
    print(f"  Child BAM: {child_bam}")
    print(f"  Mom BAM: {mom_bam}")
    print(f"  Dad BAM: {dad_bam}")
    print("\nNext steps:")
    print(f"  1. Run: tldr -b {child_bam} -e {te_file} -r {ref_file} --denovo -o output/child")
    print(f"  2. Run: python scripts/parent_support.py --mom_bam {mom_bam} --dad_bam {dad_bam} \\")
    print(f"          --candidates output/child.table.txt --output output/parent_analysis.json --denovo")

if __name__ == "__main__":
    main()
