#!/usr/bin/env python3
"""
IGV-style visualization for de novo TE insertions.

Creates a multi-track view showing child, mom, and dad BAM files
with highlighted insertion regions.
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
import pysam
import numpy as np
from collections import defaultdict
import argparse
import json
import logging

logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)


class IGVStyleVisualizer:
    """IGV-style read alignment visualizer for de novo TE insertions."""

    def __init__(self, figsize=(18, 12)):
        self.figsize = figsize
        self.colors = {
            'match': '#A0A0A0',           # Gray for aligned reads (IGV default)
            'insertion': '#8B008B',        # Dark magenta for insertion markers
            'child': '#4169E1',            # Royal blue for child label
            'mom': '#DC143C',              # Crimson for mom label
            'dad': '#228B22',              # Forest green for dad label
            'highlight': '#808080',        # Gray for insertion region highlight
            'highlight_alpha': 0.25,       # Semi-transparent for IGV-style highlighting
            'coverage': '#90EE90',          # Light green for coverage
            'softclip': '#FF4500',         # Orange-red for soft-clipped reads
        }

    def extract_reads(self, bam_file, chrom, start, end, max_reads=50):
        """
        Extract reads from BAM file in the specified region.

        Returns:
            list of dicts with read information
        """
        reads = []
        try:
            bam = pysam.AlignmentFile(bam_file, 'rb')

            for i, read in enumerate(bam.fetch(chrom, start, end)):
                if i >= max_reads:
                    break

                if read.is_secondary or read.is_supplementary:
                    continue

                # Parse CIGAR to find insertions
                insertions = []
                ref_pos = read.reference_start
                query_pos = 0

                if read.cigartuples:
                    for op, length in read.cigartuples:
                        if op == 1:  # Insertion
                            insertions.append({
                                'ref_pos': ref_pos,
                                'length': length,
                                'query_pos': query_pos
                            })

                        # Update positions
                        if op in (0, 7, 8):  # M, =, X (match/mismatch)
                            ref_pos += length
                            query_pos += length
                        elif op == 1:  # I (insertion)
                            query_pos += length
                        elif op in (2, 3):  # D, N (deletion)
                            ref_pos += length
                        elif op == 4:  # S (soft clip)
                            query_pos += length

                reads.append({
                    'name': read.query_name,
                    'start': read.reference_start,
                    'end': read.reference_end,
                    'insertions': insertions,
                    'cigar': read.cigartuples,
                    'mapq': read.mapping_quality
                })

            bam.close()

        except Exception as e:
            logger.warning(f"Error reading {bam_file}: {e}")

        return reads

    def calculate_coverage(self, reads, start, end, bin_size=10):
        """Calculate coverage across the region."""
        coverage = defaultdict(int)

        for read in reads:
            for pos in range(max(start, read['start']), min(end, read['end'])):
                coverage[pos] += 1

        # Bin the coverage
        positions = []
        depths = []
        for i in range(start, end, bin_size):
            bin_depth = sum(coverage[p] for p in range(i, min(i + bin_size, end)))
            positions.append(i)
            depths.append(bin_depth / bin_size)

        return positions, depths

    def plot_track(self, ax, reads, track_name, track_color, start, end,
                   insertion_region=None, y_offset=0, max_rows=20,
                   is_child_track=False):
        """
        Plot a single BAM track (IGV-style).

        Args:
            ax: matplotlib axes
            reads: list of read dicts
            track_name: name for the track
            track_color: color for the track label
            start, end: genomic coordinates
            insertion_region: (start, end) tuple to highlight
            y_offset: vertical offset for this track
            max_rows: maximum number of read rows to display
            is_child_track: whether this is the child track (for highlighting)
        """
        # Sort reads by start position
        reads = sorted(reads, key=lambda r: r['start'])

        # Assign reads to rows (pileup style - IGV-style packing)
        rows = []
        for read in reads[:max_rows * 3]:
            placed = False
            for row_idx, row in enumerate(rows):
                if len(row) == 0 or row[-1]['end'] < read['start']:
                    row.append(read)
                    placed = True
                    break

            if not placed and len(rows) < max_rows:
                rows.append([read])

        # Draw semi-transparent gray highlight region for child track (IGV-style)
        if is_child_track and insertion_region:
            ins_start, ins_end = insertion_region
            highlight = mpatches.Rectangle(
                (ins_start, y_offset - max_rows * 0.8),
                ins_end - ins_start,
                max_rows * 0.8,
                facecolor=self.colors['highlight'],
                edgecolor='none',
                linewidth=0,
                alpha=self.colors['highlight_alpha'],
                zorder=0
            )
            ax.add_patch(highlight)

            # Add thin border to highlight region for visibility
            highlight_border = mpatches.Rectangle(
                (ins_start, y_offset - max_rows * 0.8),
                ins_end - ins_start,
                max_rows * 0.8,
                facecolor='none',
                edgecolor=self.colors['highlight'],
                linewidth=1,
                alpha=0.5,
                zorder=1
            )
            ax.add_patch(highlight_border)

        # Draw reads (IGV-style: gray rectangles with optional strand arrows)
        read_height = 0.7
        read_spacing = 0.8

        for row_idx, row in enumerate(rows):
            y_pos = y_offset - row_idx * read_spacing

            for read in row:
                read_start = read['start']
                read_end = read['end']
                read_width = read_end - read_start

                # Draw main read body (IGV-style gray rectangle)
                rect = mpatches.Rectangle(
                    (read_start, y_pos - read_height/2),
                    read_width, read_height,
                    facecolor=self.colors['match'],
                    edgecolor='#404040',
                    linewidth=0.4,
                    alpha=0.85,
                    zorder=2
                )
                ax.add_patch(rect)

                # Draw insertions as vertical bars (dark magenta)
                for ins in read['insertions']:
                    ins_x = ins['ref_pos']
                    ins_marker = mpatches.Rectangle(
                        (ins_x - 1, y_pos - read_height/2),
                        2, read_height,
                        facecolor=self.colors['insertion'],
                        edgecolor=self.colors['insertion'],
                        linewidth=0,
                        alpha=0.95,
                        zorder=3
                    )
                    ax.add_patch(ins_marker)

        # Draw track label on the left side (IGV-style)
        label_x = start - (end - start) * 0.02
        ax.text(label_x, y_offset - max_rows * read_spacing / 2,
                track_name,
                fontsize=12, fontweight='bold',
                verticalalignment='center',
                horizontalalignment='right',
                color=track_color)

        return len(rows)

    def visualize_denovo_insertion(self, child_bam, mom_bam, dad_bam,
                                   chrom, bp_left, bp_right,
                                   candidate_info=None,
                                   output_file=None,
                                   window_extend=500):
        """
        Create IGV-style visualization of de novo insertion.

        Displays three BAM tracks (Child, Mom, Dad) with the insertion region
        highlighted in semi-transparent gray on the child track.

        Args:
            child_bam, mom_bam, dad_bam: BAM file paths
            chrom: chromosome name
            bp_left, bp_right: insertion breakpoints
            candidate_info: dict with insertion metadata
            output_file: output PNG path
            window_extend: bases to show on each side
        """
        # Define viewing region
        view_start = max(0, bp_left - window_extend)
        view_end = bp_right + window_extend

        logger.info(f"Visualizing {chrom}:{view_start}-{view_end}")

        # Extract reads from all three BAM files
        child_reads = self.extract_reads(child_bam, chrom, view_start, view_end, max_reads=60)
        mom_reads = self.extract_reads(mom_bam, chrom, view_start, view_end, max_reads=60)
        dad_reads = self.extract_reads(dad_bam, chrom, view_start, view_end, max_reads=60)

        logger.info(f"  Child: {len(child_reads)} reads")
        logger.info(f"  Mom:   {len(mom_reads)} reads")
        logger.info(f"  Dad:   {len(dad_reads)} reads")

        # Create figure with IGV-style layout
        fig, axes = plt.subplots(4, 1, figsize=self.figsize,
                                gridspec_kw={'height_ratios': [0.8, 3, 3, 3]})
        fig.suptitle(f'De Novo Insertion: {chrom}:{bp_left:,}-{bp_right:,}',
                    fontsize=16, fontweight='bold', y=0.98)

        # Plot 1: Genomic coordinates ruler (IGV-style)
        ax_coord = axes[0]
        ax_coord.set_xlim(view_start, view_end)
        ax_coord.set_ylim(-0.5, 0.5)
        ax_coord.axis('off')

        # Draw coordinate ruler line
        coord_y = 0
        ax_coord.plot([view_start, view_end], [coord_y, coord_y], 'k-', linewidth=2)

        # Add tick marks and labels
        num_ticks = 5
        for i in range(num_ticks + 1):
            tick_x = view_start + (view_end - view_start) * i / num_ticks
            ax_coord.plot([tick_x, tick_x], [coord_y - 0.1, coord_y + 0.1], 'k-', linewidth=1)
            ax_coord.text(tick_x, coord_y - 0.25, f'{int(tick_x):,}',
                         ha='center', va='top', fontsize=9)

        # Highlight insertion region on coordinate track
        ax_coord.axvspan(bp_left, bp_right, alpha=0.4, color=self.colors['highlight'])
        ax_coord.text((bp_left + bp_right) / 2, 0.4, 'De Novo Insertion',
                     ha='center', va='bottom', fontsize=10, fontweight='bold',
                     color='darkred')

        # Plot 2-4: Child, Mom, Dad tracks (IGV-style)
        tracks = [
            (axes[1], child_reads, 'Child', self.colors['child'], True),
            (axes[2], mom_reads, 'Mom', self.colors['mom'], False),
            (axes[3], dad_reads, 'Dad', self.colors['dad'], False)
        ]

        for ax, reads, name, color, is_child in tracks:
            ax.set_xlim(view_start, view_end)
            ax.set_ylim(-21, 1)
            ax.set_yticks([])
            ax.spines['top'].set_visible(False)
            ax.spines['left'].set_visible(False)
            ax.spines['right'].set_visible(False)

            # Only show x-axis on bottom track
            if ax != axes[-1]:
                ax.set_xticks([])
                ax.spines['bottom'].set_visible(False)
            else:
                ax.set_xlabel('Genomic Position (bp)', fontsize=11)
                ax.tick_params(axis='x', labelsize=9)

            # Plot the track with highlight on child track
            insertion_region = (bp_left, bp_right) if is_child else None
            num_rows = self.plot_track(ax, reads, name, color,
                                      view_start, view_end,
                                      insertion_region=insertion_region,
                                      y_offset=0, max_rows=20,
                                      is_child_track=is_child)

        # Add metadata panel if provided
        if candidate_info:
            info_text = f"UUID: {candidate_info.get('uuid', 'N/A')}\n"
            info_text += f"TE Family: {candidate_info.get('te_family', 'N/A')}\n"
            info_text += f"Evaluation: {candidate_info.get('evaluation', 'N/A')}"

            if candidate_info.get('child_support'):
                info_text += f"\n{candidate_info['child_support']}"

            fig.text(0.98, 0.98, info_text,
                    transform=fig.transFigure,
                    fontsize=9,
                    verticalalignment='top',
                    horizontalalignment='right',
                    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

        # Add legend (IGV-style)
        legend_elements = [
            mpatches.Patch(facecolor=self.colors['match'], edgecolor='#404040',
                          label='Aligned reads', alpha=0.85),
            mpatches.Patch(facecolor=self.colors['insertion'],
                          label='Insertion (CIGAR)', alpha=0.95),
            mpatches.Patch(facecolor=self.colors['highlight'],
                          edgecolor=self.colors['highlight'],
                          label='De Novo region (gray highlight)', alpha=self.colors['highlight_alpha'])
        ]
        fig.legend(handles=legend_elements, loc='upper left',
                  bbox_to_anchor=(0.01, 0.97), fontsize=9,
                  frameon=True, facecolor='white', edgecolor='gray')

        plt.tight_layout()
        plt.subplots_adjust(top=0.93)

        # Save or show
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
            logger.info(f"Saved visualization: {output_file}")
        else:
            plt.show()

        plt.close()


def visualize_from_parent_analysis(parent_analysis_json, child_bam, mom_bam, dad_bam,
                                   output_dir='visualizations',
                                   filter_evaluation=None):
    """
    Create visualizations for all candidates from parent analysis.

    Args:
        parent_analysis_json: path to parent_analysis.json
        child_bam, mom_bam, dad_bam: BAM file paths
        output_dir: directory for output images
        filter_evaluation: only visualize candidates with this evaluation (e.g., 'PASS_DENOVO')
    """
    import os
    os.makedirs(output_dir, exist_ok=True)

    # Load analysis results
    with open(parent_analysis_json, 'r') as f:
        results = json.load(f)

    visualizer = IGVStyleVisualizer()

    logger.info(f"\nGenerating IGV-style visualizations...")
    logger.info(f"=" * 70)

    for result in results:
        uuid = result.get('uuid', 'unknown')
        evaluation = result.get('evaluation', 'UNKNOWN')

        # Filter by evaluation if specified
        if filter_evaluation and evaluation != filter_evaluation:
            continue

        chrom = result.get('chrom') or result.get('Chrom')
        bp_left = result.get('bp_left')
        bp_right = result.get('bp_right')
        uuid = result.get('uuid') or result.get('UUID')
        te_family = result.get('te_family') or result.get('TE_family')

        if not all([chrom, bp_left is not None, bp_right is not None]):
            logger.warning(f"Skipping {uuid}: missing coordinates")
            continue

        logger.info(f"\n{uuid} ({evaluation})")

        output_file = os.path.join(output_dir, f"{uuid}.png")

        try:
            visualizer.visualize_denovo_insertion(
                child_bam, mom_bam, dad_bam,
                chrom, bp_left, bp_right,
                candidate_info=result,
                output_file=output_file
            )
        except Exception as e:
            logger.error(f"  Error: {e}")
            continue

    logger.info(f"\n" + "=" * 70)
    logger.info(f"Visualizations saved to: {output_dir}/")


def main():
    parser = argparse.ArgumentParser(
        description='Create IGV-style visualizations for de novo TE insertions'
    )
    parser.add_argument('--parent_analysis', required=True,
                       help='parent_analysis.json from parent_support.py')
    parser.add_argument('--child_bam', required=True,
                       help='Child BAM file')
    parser.add_argument('--mom_bam', required=True,
                       help='Maternal BAM file')
    parser.add_argument('--dad_bam', required=True,
                       help='Paternal BAM file')
    parser.add_argument('--output_dir', default='visualizations',
                       help='Output directory for PNG files')
    parser.add_argument('--filter', choices=['PASS_DENOVO', 'UNCERTAIN', 'FAIL'],
                       help='Only visualize candidates with this evaluation')
    parser.add_argument('--window', type=int, default=500,
                       help='Bases to show on each side of insertion (default: 500)')

    args = parser.parse_args()

    visualize_from_parent_analysis(
        args.parent_analysis,
        args.child_bam,
        args.mom_bam,
        args.dad_bam,
        output_dir=args.output_dir,
        filter_evaluation=args.filter
    )


if __name__ == '__main__':
    main()
