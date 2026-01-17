#!/usr/bin/env python3
"""
Create a mock tldr output table for testing parent_support.py
"""

import os

def create_mock_tldr_table():
    """Create a mock tldr output table with de novo fields."""

    output_file = "output/child.table.txt"
    os.makedirs("output", exist_ok=True)

    # Header with denovo fields
    header = [
        "Chromosome", "Start", "End", "UUID", "Superfamily", "Subfamily",
        "StartTE", "EndTE", "Strand", "LengthIns", "Inversion", "UnmapCover",
        "MedianMapQ", "TEMatch", "UsedReads", "SpanReads", "NumSamples",
        "SampleReads", "EmptyReads", "EmptyPhase", "NonRef", "TSD",
        "Consensus", "Phasing", "Remappable", "Filter",
        "bp_left", "bp_right", "wiggle", "TE_family", "child_support"
    ]

    # Create three test insertions:
    # 1. De novo insertion at chr1:25000 (should be PASS_DENOVO)
    # 2. Insertion at chr1:30000 (coordinate 0 test)
    # 3. Insertion at chr1:35000 with edge case

    insertions = [
        {
            "Chromosome": "chr1",
            "Start": "24900",
            "End": "25100",
            "UUID": "ins_001_denovo",
            "Superfamily": "LINE",
            "Subfamily": "L1HS",
            "StartTE": "1",
            "EndTE": "300",
            "Strand": "+",
            "LengthIns": "300",
            "Inversion": "N",
            "UnmapCover": "0.95",
            "MedianMapQ": "60",
            "TEMatch": "0.98",
            "UsedReads": "10",
            "SpanReads": "10",
            "NumSamples": "1",
            "SampleReads": "child|10",
            "EmptyReads": "0",
            "EmptyPhase": "NA",
            "NonRef": "NA",
            "TSD": "TGAC",
            "Consensus": "TGACatcgatcgatcgatcgatcgTGAC",
            "Phasing": "NA",
            "Remappable": "True",
            "Filter": "PASS",
            "bp_left": "25000",
            "bp_right": "25300",
            "wiggle": "50",
            "TE_family": "L1HS",
            "child_support": "useable:10,embedded:10"
        },
        {
            "Chromosome": "chr1",
            "Start": "0",  # Test coordinate 0
            "End": "200",
            "UUID": "ins_002_coord_zero",
            "Superfamily": "SINE",
            "Subfamily": "Alu",
            "StartTE": "1",
            "EndTE": "200",
            "Strand": "+",
            "LengthIns": "200",
            "Inversion": "N",
            "UnmapCover": "0.92",
            "MedianMapQ": "55",
            "TEMatch": "0.96",
            "UsedReads": "8",
            "SpanReads": "8",
            "NumSamples": "1",
            "SampleReads": "child|8",
            "EmptyReads": "0",
            "EmptyPhase": "NA",
            "NonRef": "NA",
            "TSD": "AGCT",
            "Consensus": "AGCTgatcgatcgatcgatcgAGCT",
            "Phasing": "NA",
            "Remappable": "True",
            "Filter": "PASS",
            "bp_left": "0",  # Should not be treated as falsy!
            "bp_right": "200",
            "wiggle": "50",
            "TE_family": "Alu",
            "child_support": "useable:8,embedded:8"
        },
        {
            "Chromosome": "chr1",
            "Start": "34900",
            "End": "35100",
            "UUID": "ins_003_inherited",
            "Superfamily": "LINE",
            "Subfamily": "L1PA2",
            "StartTE": "1",
            "EndTE": "250",
            "Strand": "-",
            "LengthIns": "250",
            "Inversion": "N",
            "UnmapCover": "0.90",
            "MedianMapQ": "58",
            "TEMatch": "0.94",
            "UsedReads": "12",
            "SpanReads": "12",
            "NumSamples": "1",
            "SampleReads": "child|12",
            "EmptyReads": "0",
            "EmptyPhase": "NA",
            "NonRef": "NA",
            "TSD": "CGAT",
            "Consensus": "CGATtcgatcgatcgatcgatCGAT",
            "Phasing": "NA",
            "Remappable": "True",
            "Filter": "PASS",
            "bp_left": "35000",
            "bp_right": "35250",
            "wiggle": "50",
            "TE_family": "L1PA2",
            "child_support": "useable:12,embedded:12"
        }
    ]

    with open(output_file, 'w') as f:
        # Write header
        f.write('\t'.join(header) + '\n')

        # Write insertions
        for ins in insertions:
            row = [str(ins[col]) for col in header]
            f.write('\t'.join(row) + '\n')

    print(f"Created mock tldr table: {output_file}")
    print(f"  - {len(insertions)} candidate insertions")
    print(f"  - ins_001: Should be PASS_DENOVO (parents have no support)")
    print(f"  - ins_002: Coordinate 0 test (tests falsy value fix)")
    print(f"  - ins_003: Should be FAIL (parents have support)")

    return output_file

if __name__ == "__main__":
    create_mock_tldr_table()
