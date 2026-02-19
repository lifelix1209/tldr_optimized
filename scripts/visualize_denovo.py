#!/usr/bin/env python3
"""
Enhanced IGV-style visualization for de novo TE insertions.

Features:
- Soft-clip visualization: Mark clipped read ends with conspicuous colors
- Coverage track: Area graph showing depth variations above reads
- Read direction: Arrow-shaped reads distinguishing forward/reverse strands
- Synchronized highlighting across all tracks
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.path import Path
from matplotlib.patches import FancyArrowPatch, Rectangle
import pysam
import numpy as np
from collections import defaultdict
import argparse
import json
import logging

logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)


class EnhancedIGVVisualizer:
    """
    Enhanced IGV-style visualizer with:
    - Soft-clip visualization for TE breakpoint detection
    - Coverage depth area track
    - Strand-specific arrow-shaped reads
    """

    def __init__(self, figsize=(16, 14)):
        self.figsize = figsize

        # Professional color palette (colorblind-friendly)
        self.colors = {
            # Read colors by strand
            'forward': '#2E86AB',           # Blue for forward strand (5'->3')
            'forward_dark': '#1B4965',      # Darker blue for arrow body
            'reverse': '#E94F37',           # Red for reverse strand (3'->5')
            'reverse_dark': '#A82E1E',      # Darker red for arrow body

            # Feature colors
            'softclip_left': '#FF6B35',     # Orange-red for left soft-clips (5' end)
            'softclip_right': '#FF0000',    # Bright red for right soft-clips (3' end)
            'insertion': '#8B008B',         # Magenta for insertion markers

            # Coverage colors
            'coverage_fill': '#90EE90',     # Light green fill
            'coverage_line': '#228B22',     # Forest green line
            'coverage_highlight': '#7CFC00', # Lawn green for high depth

            # Track colors
            'child': '#2E86AB',             # Child track color
            'mom': '#A23B72',               # Mom track color
            'dad': '#F18F01',               # Dad track color

            # Highlight
            'highlight': '#FFE066',         # Yellow highlight for insertion region
            'highlight_border': '#E67E22',  # Orange border

            # Other
            'bg': '#FFFFFF',                # White background
        }

    def extract_reads_with_features(self, bam_file, chrom, start, end, max_reads=80):
        """
        Extract reads with detailed feature extraction for visualization.

        Returns:
            dict with 'reads', 'coverage_x', 'coverage_y'
        """
        reads = []
        coverage = defaultdict(int)

        try:
            bam = pysam.AlignmentFile(bam_file, 'rb')

            for read in bam.fetch(chrom, start, end):
                if read.is_secondary or read.is_supplementary:
                    continue

                # Calculate coverage
                for pos in range(max(start, read.reference_start), min(end, read.reference_end)):
                    coverage[pos] += 1

                # Parse CIGAR for features
                features = {
                    'insertions': [],        # CIGAR: 1 (insertion in read)
                    'softclips': [],         # CIGAR: 4 (soft clip)
                    'deletions': [],         # CIGAR: 2 (deletion in ref)
                    'read_length': len(read.query_sequence) if read.query_sequence else 0
                }

                ref_pos = read.reference_start
                query_pos = 0

                if read.cigartuples:
                    for op, length in read.cigartuples:
                        if op == 1:  # Insertion
                            features['insertions'].append({
                                'ref_pos': ref_pos,
                                'length': length,
                                'query_pos': query_pos
                            })
                            query_pos += length
                        elif op == 4:  # Soft clip
                            # Determine if clip is at start or end of read
                            clip_end = 'left' if query_pos == 0 else 'right'
                            features['softclips'].append({
                                'ref_pos': ref_pos,
                                'length': length,
                                'end': clip_end
                            })
                            query_pos += length
                        elif op == 2:  # Deletion
                            features['deletions'].append({
                                'start': ref_pos,
                                'end': ref_pos + length
                            })
                            ref_pos += length
                        elif op in (0, 7, 8):  # Match/Mismatch
                            ref_pos += length
                            query_pos += length
                        elif op == 3:  # Skip
                            ref_pos += length

                reads.append({
                    'name': read.query_name,
                    'start': read.reference_start,
                    'end': read.reference_end,
                    'is_reverse': read.is_reverse,
                    'mapq': read.mapping_quality,
                    'features': features,
                    'query_seq': read.query_sequence,
                    'clip_len': features['softclips'][0]['length'] if features['softclips'] else 0
                })

                if len(reads) >= max_reads:
                    break

            bam.close()

        except Exception as e:
            logger.warning(f"Error reading {bam_file}: {e}")
            return {'reads': [], 'coverage_x': np.arange(start, end), 'coverage_y': np.zeros(end - start)}

        # Prepare coverage data
        cov_x = np.arange(start, end)
        cov_y = np.array([coverage.get(p, 0) for p in cov_x])

        # Smooth coverage
        if len(cov_y) > 5:
            kernel_size = min(5, len(cov_y) // 2 * 2 + 1)
            if kernel_size >= 3:
                kernel = np.ones(kernel_size) / kernel_size
                cov_y = np.convolve(cov_y, kernel, mode='same')

        return {'reads': reads, 'coverage_x': cov_x, 'coverage_y': cov_y}

    def _create_strand_arrow(self, x_start, x_end, y_pos, is_reverse, height=0.6):
        """
        Create a rectangular read with arrow head at the end indicating strand direction.

        Args:
            x_start, x_end: genomic coordinates
            y_pos: vertical position
            is_reverse: True if reverse strand
            height: arrow height
        """
        arrow_len = min((x_end - x_start) * 0.1, 12)  # Arrow head length
        body_len = x_end - x_start - arrow_len

        # Use light gray for read body (as requested)
        color = '#B0B0B0'  # Light gray

        if is_reverse:
            # Reverse strand: arrow points left
            # Rectangle body from x_start+arrow_len to x_end
            # Small triangle at left end pointing left
            vertices = [
                (x_start + arrow_len, y_pos - height/2),      # Bottom left of body
                (x_end, y_pos - height/2),                    # Bottom right
                (x_end, y_pos + height/2),                    # Top right
                (x_start + arrow_len, y_pos + height/2),      # Top left of body
                (x_start + arrow_len, y_pos + height/4),      # Arrow notch top
                (x_start, y_pos),                             # Arrow tip
                (x_start + arrow_len, y_pos - height/4),      # Arrow notch bottom
            ]
            codes = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO,
                     Path.LINETO, Path.LINETO, Path.CLOSEPOLY]
        else:
            # Forward strand: arrow points right
            # Rectangle body from x_start to x_end-arrow_len
            # Small triangle at right end pointing right
            vertices = [
                (x_start, y_pos - height/2),                  # Bottom left
                (x_end - arrow_len, y_pos - height/2),        # Bottom left of arrow
                (x_end - arrow_len, y_pos - height/4),        # Arrow notch bottom
                (x_end, y_pos),                               # Arrow tip
                (x_end - arrow_len, y_pos + height/4),        # Arrow notch top
                (x_end - arrow_len, y_pos + height/2),        # Top left of arrow
                (x_start, y_pos + height/2),                  # Top left
            ]
            codes = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO,
                     Path.LINETO, Path.LINETO, Path.LINETO]

        path = Path(vertices, codes)
        return path, color

    def _draw_softclip_marker(self, ax, read, y_pos, height=0.6, clip_length_threshold=10):
        """
        Draw conspicuous soft-clip markers at read ends.

        Key feature for TE insertion breakpoint detection:
        - Left soft-clip (5' end): orange-red
        - Right soft-clip (3' end): bright red
        """
        for clip in read['features']['softclips']:
            clip_len = clip['length']
            clip_end = clip['end']

            if clip_len < 3:  # Skip very small clips
                continue

            # Make longer clips more visible
            opacity = min(1.0, 0.5 + clip_len / 50)

            if clip_end == 'left':
                # Left soft-clip (5' end of read) - indicates TE start breakpoint
                clip_x = read['start']
                clip_width = min(clip_len * 1.5, 40)  # Visual width

                # Draw triangular marker at 5' end
                marker = Rectangle(
                    (clip_x - clip_width, y_pos - height/2),
                    clip_width, height,
                    facecolor=self.colors['softclip_left'],
                    edgecolor=self.colors['softclip_left'],
                    linewidth=1.5,
                    alpha=opacity,
                    zorder=5
                )
                ax.add_patch(marker)

                # Add small "S" indicator for soft-clip
                ax.text(clip_x - clip_width/2, y_pos, 'S',
                       fontsize=6, ha='center', va='center',
                       color='white', fontweight='bold', zorder=6)

            else:
                # Right soft-clip (3' end of read) - indicates TE end breakpoint
                clip_x = read['end'] - min(clip_len * 1.5, 40)
                clip_width = min(clip_len * 1.5, 40)

                marker = Rectangle(
                    (clip_x, y_pos - height/2),
                    clip_width, height,
                    facecolor=self.colors['softclip_right'],
                    edgecolor=self.colors['softclip_right'],
                    linewidth=1.5,
                    alpha=opacity,
                    zorder=5
                )
                ax.add_patch(marker)

                # Add small "S" indicator
                ax.text(clip_x + clip_width/2, y_pos, 'S',
                       fontsize=6, ha='center', va='center',
                       color='white', fontweight='bold', zorder=6)

    def _draw_insertion_marker(self, ax, read, y_pos, height=0.6):
        """Draw vertical markers for CIGAR insertions."""
        for ins in read['features']['insertions']:
            ins_x = ins['ref_pos']
            ins_len = ins['length']

            if ins_len < 3:  # Skip small insertions
                continue

            marker = Rectangle(
                (ins_x - 1.5, y_pos - height * 0.7),
                3, height * 1.4,
                facecolor=self.colors['insertion'],
                edgecolor='none',
                alpha=0.9,
                zorder=4
            )
            ax.add_patch(marker)

    def _plot_coverage_track(self, ax, cov_x, cov_y, view_start, view_end,
                             track_height=2, base_y=0):
        """
        Plot coverage as an area graph above the reads.

        This visualizes depth variations that may indicate:
        - Coverage drops (possible deletion/insertion)
        - Coverage spikes (possible duplication)
        """
        # Clip to view
        mask = (cov_x >= view_start) & (cov_x <= view_end)
        plot_x = cov_x[mask]
        plot_y = cov_y[mask]

        if len(plot_y) == 0:
            return

        # Normalize to track height
        max_depth = max(np.max(plot_y), 1)
        scale = track_height / max_depth
        normalized_y = plot_y * scale

        # Fill area under curve
        ax.fill_between(plot_x, base_y, base_y + normalized_y,
                       color=self.colors['coverage_fill'], alpha=0.4, zorder=1)

        # Draw outline
        ax.plot(plot_x, base_y + normalized_y,
               color=self.colors['coverage_line'], linewidth=1.5, alpha=0.8, zorder=2)

        # Add depth annotations
        if max_depth > 0:
            ax.text(view_end - (view_end - view_start) * 0.02,
                   base_y + track_height * 0.9,
                   f'{int(max_depth)}x',
                   fontsize=8, ha='right', va='top',
                   bbox=dict(boxstyle="round,pad=0.2", facecolor='white',
                            edgecolor='none', alpha=0.8))

        # Mark depth drops (possible insertion breakpoints)
        depth_threshold = max_depth * 0.5
        drop_positions = []
        for i in range(1, len(plot_y)):
            if plot_y[i] < depth_threshold and plot_y[i-1] >= depth_threshold:
                drop_positions.append(plot_x[i])

        for pos in drop_positions:
            ax.axvline(x=pos, color='red', linestyle='--', alpha=0.3, linewidth=1)

    def _layout_reads(self, reads, max_rows=15, min_gap=5):
        """
        Pack reads into rows to minimize overlap (IGV-style).

        Returns:
            list of (read, row_idx) tuples
        """
        reads = sorted(reads, key=lambda r: r['start'])
        rows = []

        for read in reads:
            placed = False
            for row_idx, row_reads in enumerate(rows):
                # Check if read fits in this row (non-overlapping)
                fits = True
                for r in row_reads:
                    if not (read['end'] + min_gap < r['start'] or
                            r['end'] + min_gap < read['start']):
                        fits = False
                        break
                if fits:
                    row_reads.append(read)
                    placed = True
                    break

            if not placed and len(rows) < max_rows:
                rows.append([read])

        # Return reads with assigned row indices
        result = []
        for row_idx, row_reads in enumerate(rows):
            for read in row_reads:
                result.append((read, row_idx))
        return result

    def plot_track(self, ax, data, track_name, track_color, view_start, view_end,
                   insertion_region=None, max_rows=15):
        """
        Plot a single BAM track with coverage, reads, and features.

        Args:
            ax: matplotlib axes
            data: dict with 'reads', 'coverage_x', 'coverage_y'
            track_name: name for the track
            track_color: color for the track label
            view_start, end: genomic coordinates
            insertion_region: (start, end) tuple to highlight
        """
        reads = data['reads']
        cov_x = data['coverage_x']
        cov_y = data['coverage_y']

        ax.set_xlim(view_start, view_end)

        # Calculate layout
        layout = self._layout_reads(reads, max_rows=max_rows)
        track_height = max(6, len(layout) * 0.8 + 3)

        # Plot coverage track at top
        self._plot_coverage_track(ax, cov_x, cov_y, view_start, view_end,
                                 track_height=2, base_y=track_height - 2)

        # Highlight insertion region
        if insertion_region:
            ins_start, ins_end = insertion_region
            highlight = Rectangle(
                (ins_start, 0),
                ins_end - ins_start,
                track_height,
                facecolor=self.colors['highlight'],
                edgecolor=self.colors['highlight_border'],
                linewidth=2,
                linestyle='--',
                alpha=0.3,
                zorder=0
            )
            ax.add_patch(highlight)

        # Draw reads with direction arrows
        read_height = 0.6
        for read, row_idx in layout:
            y_pos = track_height - 3 - row_idx * 0.8

            # Only draw if in view
            if read['end'] < view_start or read['start'] > view_end:
                continue

            # Draw strand arrow
            path, color = self._create_strand_arrow(
                max(view_start, read['start']),
                min(view_end, read['end']),
                y_pos,
                read['is_reverse'],
                height=read_height
            )
            patch = mpatches.PathPatch(path, facecolor=color,
                                       edgecolor='black', linewidth=0.5, zorder=2)
            ax.add_patch(patch)

            # Draw soft-clip markers
            self._draw_softclip_marker(ax, read, y_pos, height=read_height)

            # Draw insertion markers
            self._draw_insertion_marker(ax, read, y_pos, height=read_height)

        # Track label
        ax.text(view_start - (view_end - view_start) * 0.01,
               track_height - 1,
               track_name,
               fontsize=11, fontweight='bold',
               verticalalignment='center',
               horizontalalignment='right',
               color=track_color)

        ax.set_ylim(0, track_height)
        ax.set_yticks([])

        # Remove spines
        for spine in ax.spines.values():
            spine.set_visible(False)

        return track_height

    def visualize_denovo_insertion(self, child_bam, mom_bam, dad_bam,
                                   chrom, bp_left, bp_right,
                                   candidate_info=None,
                                   output_file=None,
                                   window_extend=500):
        """
        Create enhanced IGV-style visualization of de novo insertion.

        Features:
        - Three tracks: Child, Mom, Dad
        - Coverage area track above each read track
        - Arrow-shaped reads showing strand direction
        - Conspicuous soft-clip markers for breakpoint detection
        - Synchronized highlighting of insertion region
        """
        # Define viewing region
        view_start = max(0, bp_left - window_extend)
        view_end = bp_right + window_extend

        logger.info(f"Visualizing {chrom}:{view_start:,}-{view_end:,}")

        # Extract reads and coverage from all three BAM files
        child_data = self.extract_reads_with_features(child_bam, chrom, view_start, view_end, max_reads=80)
        mom_data = self.extract_reads_with_features(mom_bam, chrom, view_start, view_end, max_reads=80)
        dad_data = self.extract_reads_with_features(dad_bam, chrom, view_start, view_end, max_reads=80)

        logger.info(f"  Child: {len(child_data['reads'])} reads, max coverage: {int(np.max(child_data['coverage_y']))}x")
        logger.info(f"  Mom:   {len(mom_data['reads'])} reads, max coverage: {int(np.max(mom_data['coverage_y']))}x")
        logger.info(f"  Dad:   {len(dad_data['reads'])} reads, max coverage: {int(np.max(dad_data['coverage_y']))}x")

        # Create figure
        fig, axes = plt.subplots(4, 1, figsize=self.figsize,
                                gridspec_kw={'height_ratios': [0.3, 1, 1, 1]})

        # Title
        uuid = candidate_info.get('uuid', 'N/A') if candidate_info else 'N/A'
        eval_status = candidate_info.get('evaluation', 'UNKNOWN') if candidate_info else 'UNKNOWN'
        fig.suptitle(f'De Novo TE Insertion: {chrom}:{bp_left:,}-{bp_right:,}  |  {uuid}  |  {eval_status}',
                    fontsize=14, fontweight='bold', y=0.98)

        # Track 1: Ruler/coordinates
        ax_ruler = axes[0]
        ax_ruler.set_xlim(view_start, view_end)
        ax_ruler.set_ylim(-0.5, 0.5)
        ax_ruler.axis('off')

        # Draw ruler line
        ax_ruler.plot([view_start, view_end], [0, 0], 'k-', linewidth=2)

        # Add ticks and labels
        num_ticks = 6
        for i in range(num_ticks + 1):
            tick_x = view_start + (view_end - view_start) * i / num_ticks
            ax_ruler.plot([tick_x, tick_x], [-0.1, 0.1], 'k-', linewidth=1)
            ax_ruler.text(tick_x, -0.25, f'{int(tick_x):,}',
                         ha='center', va='top', fontsize=9)

        # Highlight insertion on ruler
        ax_ruler.axvspan(bp_left, bp_right, alpha=0.5, color=self.colors['highlight'])
        ax_ruler.text((bp_left + bp_right) / 2, 0.5, 'De Novo Region',
                     ha='center', va='bottom', fontsize=10, fontweight='bold',
                     color=self.colors['highlight_border'])

        # Tracks 2-4: Child, Mom, Dad
        tracks = [
            (axes[1], child_data, 'Child', self.colors['child']),
            (axes[2], mom_data, 'Mom', self.colors['mom']),
            (axes[3], dad_data, 'Dad', self.colors['dad'])
        ]

        insertion_region = (bp_left, bp_right)

        for ax, data, name, color in tracks:
            ax.set_xlim(view_start, view_end)
            ax.spines['top'].set_visible(False)
            ax.spines['left'].set_visible(False)
            ax.spines['right'].set_visible(False)

            # Only show x-axis on bottom track
            if ax != axes[-1]:
                ax.set_xticks([])
            else:
                ax.set_xlabel('Genomic Position (bp)', fontsize=11)
                ax.tick_params(axis='x', labelsize=9)

            self.plot_track(ax, data, name, color, view_start, view_end,
                          insertion_region=insertion_region, max_rows=12)

        # Add metadata panel
        if candidate_info:
            info_parts = [f"UUID: {candidate_info.get('uuid', 'N/A')}"]
            info_parts.append(f"TE Family: {candidate_info.get('te_family', 'N/A')}")
            info_parts.append(f"Evaluation: {candidate_info.get('evaluation', 'N/A')}")

            if 'mom_depth' in candidate_info:
                info_parts.append(f"Mom Depth: {candidate_info.get('mom_depth', 'N/A')}x")
            if 'dad_depth' in candidate_info:
                info_parts.append(f"Dad Depth: {candidate_info.get('dad_depth', 'N/A')}x")

            if candidate_info.get('child_support'):
                info_parts.append(f"Child Support: {candidate_info['child_support']}")

            info_text = '\n'.join(info_parts)

            fig.text(0.99, 0.98, info_text,
                    transform=fig.transFigure,
                    fontsize=9,
                    verticalalignment='top',
                    horizontalalignment='right',
                    bbox=dict(boxstyle='round', facecolor='lightyellow',
                             alpha=0.9, edgecolor='gray'))

        # Add legend
        legend_elements = [
            mpatches.Patch(facecolor=self.colors['forward'], edgecolor='black',
                          label='Forward strand (5\'→3\')'),
            mpatches.Patch(facecolor=self.colors['reverse'], edgecolor='black',
                          label='Reverse strand (3\'→5\')'),
            mpatches.Patch(facecolor=self.colors['softclip_left'], edgecolor='black',
                          label='Soft-clip (5\' end)'),
            mpatches.Patch(facecolor=self.colors['softclip_right'], edgecolor='black',
                          label='Soft-clip (3\' end)'),
            mpatches.Patch(facecolor=self.colors['insertion'], edgecolor='black',
                          label='Insertion (CIGAR=1)'),
            mpatches.Patch(facecolor=self.colors['coverage_fill'], edgecolor=self.colors['coverage_line'],
                          label='Coverage depth', alpha=0.5),
            mpatches.Patch(facecolor=self.colors['highlight'], edgecolor=self.colors['highlight_border'],
                          linestyle='--', label='De Novo region'),
        ]

        fig.legend(handles=legend_elements, loc='upper left',
                  bbox_to_anchor=(0.01, 0.97), fontsize=9,
                  frameon=True, facecolor='white', edgecolor='gray',
                  ncol=2)

        plt.tight_layout()
        plt.subplots_adjust(top=0.92)

        # Save or show
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
            logger.info(f"Saved: {output_file}")
        else:
            plt.show()

        plt.close()


def visualize_from_tldr_output(tldr_table, child_bam, mom_bam, dad_bam,
                                output_dir='enhanced_visualizations',
                                filter_evaluation=None):
    """
    Create enhanced visualizations for all candidates from tldr output table.

    Uses the actual detected positions (bp_left, bp_right) from tldr
    instead of expected positions from mock data.
    """
    import os
    import csv

    os.makedirs(output_dir, exist_ok=True)

    # Read tldr output table
    candidates = []
    with open(tldr_table, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            candidates.append(row)

    visualizer = EnhancedIGVVisualizer()

    logger.info(f"\nGenerating enhanced IGV-style visualizations from tldr output...")
    logger.info(f"{'=' * 70}")

    for row in candidates:
        uuid = row.get('UUID', 'unknown')
        filter_status = row.get('Filter', 'UNKNOWN')

        # Use bp_left and bp_right from tldr output (actual detected positions)
        try:
            bp_left = int(float(row.get('bp_left', 0)))
            bp_right = int(float(row.get('bp_right', 0)))
            chrom = row.get('Chrom', 'chr1')

            if bp_left is None or bp_right is None or bp_left >= bp_right:
                logger.warning(f"Skipping {uuid}: invalid coordinates")
                continue
        except (ValueError, TypeError) as e:
            logger.warning(f"Skipping {uuid}: coordinate error: {e}")
            continue

        # Build candidate info from tldr output
        candidate_info = {
            'uuid': uuid,
            'chrom': chrom,
            'bp_left': bp_left,
            'bp_right': bp_right,
            'te_family': row.get('TE_family', 'NA'),
            'evaluation': 'PASS_DENOVO' if filter_status == 'PASS' else filter_status,
            'child_support': row.get('child_support', 'NA'),
            'mom_depth': row.get('mom_depth', 'NA'),
            'dad_depth': row.get('dad_depth', 'NA'),
            'mom_alt': row.get('mom_alt', '0'),
            'dad_alt': row.get('dad_alt', '0'),
        }

        # Filter by evaluation if requested
        if filter_evaluation and candidate_info['evaluation'] != filter_evaluation:
            continue

        logger.info(f"\n{uuid} ({filter_status})")
        logger.info(f"  Detected region: {chrom}:{bp_left:,}-{bp_right:,}")
        logger.info(f"  Expected mock data: ~{bp_left + 249}-{bp_right + 249}")  # Approximate offset

        output_file = os.path.join(output_dir, f"{uuid}_{filter_status}.png")

        try:
            visualizer.visualize_denovo_insertion(
                child_bam, mom_bam, dad_bam,
                chrom, bp_left, bp_right,
                candidate_info=candidate_info,
                output_file=output_file,
                window_extend=500
            )
        except Exception as e:
            logger.error(f"  Error: {e}")
            import traceback
            traceback.print_exc()
            continue

    logger.info(f"\n{'=' * 70}")
    logger.info(f"Visualizations saved to: {output_dir}/")


def visualize_from_parent_analysis(parent_analysis_json, child_bam, mom_bam, dad_bam,
                                   output_dir='enhanced_visualizations',
                                   filter_evaluation=None):
    """
    Create enhanced visualizations for all candidates from parent analysis.
    NOTE: This uses expected positions from parent_analysis.json.
    For actual detected positions, use visualize_from_tldr_output() instead.
    """
    import os
    os.makedirs(output_dir, exist_ok=True)

    with open(parent_analysis_json, 'r') as f:
        results = json.load(f)

    visualizer = EnhancedIGVVisualizer()

    logger.info(f"\nGenerating enhanced IGV-style visualizations...")
    logger.info(f"{'=' * 70}")

    for result in results:
        uuid = result.get('uuid', 'unknown')
        evaluation = result.get('evaluation', 'UNKNOWN')

        if filter_evaluation and evaluation != filter_evaluation:
            continue

        chrom = result.get('chrom') or result.get('Chrom')
        bp_left = result.get('bp_left')
        bp_right = result.get('bp_right')

        try:
            if bp_left is None or bp_right is None:
                logger.warning(f"Skipping {uuid}: missing coordinates")
                continue
            bp_left = int(float(bp_left))
            bp_right = int(float(bp_right))
            if bp_left >= bp_right or bp_left < 0:
                logger.warning(f"Skipping {uuid}: invalid coordinates")
                continue
        except (ValueError, TypeError) as e:
            logger.warning(f"Skipping {uuid}: coordinate error: {e}")
            continue

        logger.info(f"\n{uuid} ({evaluation})")
        logger.info(f"  Region: {chrom}:{bp_left:,}-{bp_right:,}")

        output_file = os.path.join(output_dir, f"{uuid}_{evaluation}.png")

        try:
            visualizer.visualize_denovo_insertion(
                child_bam, mom_bam, dad_bam,
                chrom, bp_left, bp_right,
                candidate_info=result,
                output_file=output_file,
                window_extend=500
            )
        except Exception as e:
            logger.error(f"  Error: {e}")
            import traceback
            traceback.print_exc()
            continue

    logger.info(f"\n{'=' * 70}")
    logger.info(f"Visualizations saved to: {output_dir}/")


def main():
    parser = argparse.ArgumentParser(
        description='Enhanced IGV-style visualization for de novo TE insertions'
    )
    parser.add_argument('--tldr_table',
                       help='tldr output table (uses actual detected positions)')
    parser.add_argument('--parent_analysis',
                       help='parent_analysis.json from parent_support.py (uses expected positions)')
    parser.add_argument('--child_bam', required=True,
                       help='Child BAM file')
    parser.add_argument('--mom_bam', required=True,
                       help='Maternal BAM file')
    parser.add_argument('--dad_bam', required=True,
                       help='Paternal BAM file')
    parser.add_argument('--output_dir', default='enhanced_visualizations',
                       help='Output directory for PNG files')
    parser.add_argument('--filter', choices=['PASS_DENOVO', 'UNCERTAIN', 'FAIL'],
                       help='Only visualize candidates with this evaluation')
    parser.add_argument('--window', type=int, default=500,
                       help='Bases to show on each side of insertion (default: 500)')

    args = parser.parse_args()

    # Prefer tldr_table if provided (uses actual detected positions)
    if args.tldr_table:
        visualize_from_tldr_output(
            args.tldr_table,
            args.child_bam,
            args.mom_bam,
            args.dad_bam,
            output_dir=args.output_dir,
            filter_evaluation=args.filter
        )
    elif args.parent_analysis:
        visualize_from_parent_analysis(
            args.parent_analysis,
            args.child_bam,
            args.mom_bam,
            args.dad_bam,
            output_dir=args.output_dir,
            filter_evaluation=args.filter
        )
    else:
        parser.error("Either --tldr_table or --parent_analysis is required")


if __name__ == '__main__':
    main()
