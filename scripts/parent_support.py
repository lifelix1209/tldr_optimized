#!/usr/bin/env python
"""
Parent support analysis module for tldr.
Analyzes Mom/Dad BAM files to compute depth and alt support for child candidate sites.

Input format: tldr --denovo output (dict or table row)
Expected fields:
    - UUID: candidate identifier
    - Chrom: chromosome
    - bp_left: left breakpoint (or Start)
    - bp_right: right breakpoint (or End)
    - Strand: '+' or '-'
    - TE_family: TE family name (or Family)
    - child_support: support info (optional)
"""

import pysam
import numpy as np
from collections import defaultdict
import logging
import re

logger = logging.getLogger(__name__)

NA_STRINGS = {'', 'NA', 'N/A', 'NONE', 'NULL'}


def _as_positive_int_or_none(value):
    """Parse a positive integer; return None when unavailable/invalid."""
    if value is None:
        return None
    if isinstance(value, str) and value.strip().upper() in NA_STRINGS:
        return None
    try:
        parsed = int(value)
    except (ValueError, TypeError):
        return None
    return parsed if parsed > 0 else None


def _safe_fetch_interval(bp_left, bp_right, pad):
    """
    Build a fetch interval that is valid even for zero-length events.
    """
    start = max(0, min(bp_left, bp_right) - pad)
    end = max(bp_left, bp_right) + pad + 1
    if end <= start:
        end = start + 1
    return start, end


def _canonical_te_family_labels(te_family):
    """
    Generate robust lookup labels for TE family/subfamily names.
    """
    if not te_family:
        return []

    labels = []

    def _add(label):
        if label and label not in labels:
            labels.append(label)

    raw = str(te_family).strip()
    if raw.upper() in NA_STRINGS:
        return []

    raw_upper = raw.upper()
    _add(raw)
    _add(raw_upper)

    # token splits such as ALU:AluYb9
    if ':' in raw:
        left, right = raw.split(':', 1)
        _add(left)
        _add(right)
        _add(left.upper())
        _add(right.upper())

    # common family aliases
    if raw_upper.startswith('ALU'):
        _add('ALU')
    if raw_upper.startswith('L1'):
        _add('L1')
        _add('LINE')
    if raw_upper.startswith('SVA'):
        _add('SVA')

    return labels


def resolve_te_sequence(te_family, te_family_seqs):
    """Resolve TE sequence by exact/family aliases."""
    if not te_family or not te_family_seqs:
        return None

    for key in _canonical_te_family_labels(te_family):
        seq = te_family_seqs.get(key)
        if seq:
            return seq

    return None


def calculate_window_depth(bam, chrom, bp_left, bp_right, window_size=100):
    """
    Calculate median depth around breakpoints.

    Args:
        bam: pysam AlignmentFile
        chrom: chromosome name
        bp_left: left breakpoint
        bp_right: right breakpoint
        window_size: window size around each breakpoint

    Returns:
        dict with 'depth_left', 'depth_right', 'depth_mean'
    """
    result = {
        'depth_left': 0,
        'depth_right': 0,
        'depth_mean': 0
    }

    try:
        left_start = max(0, bp_left - window_size)
        left_end = bp_left + window_size
        right_start = bp_right - window_size
        right_end = bp_right + window_size

        left_depths = [pc.n for pc in bam.pileup(chrom, left_start, left_end, truncate=True)]
        right_depths = [pc.n for pc in bam.pileup(chrom, right_start, right_end, truncate=True)]

        if left_depths:
            result['depth_left'] = float(np.median(left_depths))
        if right_depths:
            result['depth_right'] = float(np.median(right_depths))

        all_depths = left_depths + right_depths
        if all_depths:
            result['depth_mean'] = float(np.median(all_depths))

    except Exception as e:
        logger.warning(f'Error calculating depth for {chrom}:{bp_left}-{bp_right}: {e}')

    return result


def parse_child_support_count(child_support):
    """
    Parse child support count from tldr --denovo child_support field.

    Supported formats:
        - "useable:10,embedded:5"
        - "usable:10,embedded:5"
        - numeric values

    Returns:
        int or None if not parseable
    """
    if child_support is None:
        return None

    if isinstance(child_support, (int, float)):
        return int(child_support)

    if isinstance(child_support, dict):
        for key in ('useable', 'usable', 'support', 'total'):
            if key in child_support:
                try:
                    return int(child_support[key])
                except (ValueError, TypeError):
                    pass
        return None

    if not isinstance(child_support, str):
        return None

    m = re.search(r'(?:useable|usable)\s*:\s*(\d+)', child_support, flags=re.IGNORECASE)
    if m:
        return int(m.group(1))

    # Fallback: first integer in string
    m = re.search(r'(\d+)', child_support)
    if m:
        return int(m.group(1))

    return None


def find_softclip_reads(bam, chrom, bp_left, bp_right, wiggle=30, min_softclip_len=30, min_mapq=20):
    """
    Find reads with soft-clips near breakpoints.

    Returns:
        list of dicts with read info including clip_end, ref_pos, bp_dist
    """
    softclip_reads = []
    search_start, search_end = _safe_fetch_interval(bp_left, bp_right, wiggle)

    try:
        for read in bam.fetch(chrom, search_start, search_end):
            if read.is_secondary or read.is_supplementary or not read.cigartuples:
                continue
            if read.mapping_quality < min_mapq:
                continue

            ref_pos = read.reference_start

            for i, cig in enumerate(read.cigartuples):
                if cig[0] == 4:  # Soft clip
                    clip_length = cig[1]
                    if clip_length < min_softclip_len:
                        continue
                    clip_end = 'left' if i == 0 else 'right'

                    # Calculate position in query sequence correctly
                    # CIGAR ops that consume query: M(0), I(1), S(4), =(7), X(8)
                    if i == 0:  # Left clip
                        clip_pos_in_read = 0
                    else:  # Right clip
                        # Sum all query-consuming operations before this clip
                        clip_pos_in_read = sum(c[1] for c in read.cigartuples[:i] if c[0] in [0, 1, 4, 7, 8])

                    if clip_end == 'left':
                        anchor_pos = ref_pos
                    else:
                        anchor_pos = read.reference_end if read.reference_end is not None else ref_pos

                    left_bp_dist = abs(anchor_pos - bp_left)
                    right_bp_dist = abs(anchor_pos - bp_right)
                    bp_dist = min(left_bp_dist, right_bp_dist)
                    bp_side = 'left' if left_bp_dist <= right_bp_dist else 'right'

                    if bp_dist <= wiggle:
                        # Extract the full clipped sequence (preferred for TE filtering).
                        clip_seq = ''
                        if read.query_sequence:
                            if clip_end == 'left':
                                start_pos = 0
                                end_pos = min(clip_length, len(read.query_sequence))
                            else:
                                start_pos = max(0, clip_pos_in_read)
                                end_pos = min(clip_pos_in_read + clip_length, len(read.query_sequence))
                            clip_seq = read.query_sequence[start_pos:end_pos].upper()

                        softclip_reads.append({
                            'read_name': read.query_name,
                            'clip_end': clip_end,
                            'clip_length': clip_length,
                            'ref_pos': ref_pos,
                            'bp_pos_obs': anchor_pos,
                            'bp_side': bp_side,
                            'bp_dist': bp_dist,
                            'seq': clip_seq
                        })

    except Exception as e:
        logger.warning(f'Error finding soft-clips for {chrom}:{bp_left}-{bp_right}: {e}')

    return softclip_reads


def find_insertion_cigar_reads(bam, chrom, bp_left, bp_right, wiggle=30, min_insertion_len=50, min_mapq=20):
    """
    Find reads with insertion CIGAR operations near breakpoints.

    Returns:
        list of dicts with ins_pos, bp_dist
    """
    insertion_reads = []
    search_start, search_end = _safe_fetch_interval(bp_left, bp_right, wiggle)

    try:
        for read in bam.fetch(chrom, search_start, search_end):
            if read.is_secondary or read.is_supplementary or not read.cigartuples:
                continue
            if read.mapping_quality < min_mapq:
                continue

            ref_pos = read.reference_start

            for cig in read.cigartuples:
                op, length = cig

                if op == 1:  # Insertion
                    if length < min_insertion_len:
                        continue
                    # Check if insertion is near either breakpoint
                    near_left = bp_left - wiggle <= ref_pos <= bp_left + wiggle
                    near_right = bp_right - wiggle <= ref_pos <= bp_right + wiggle

                    if near_left or near_right:
                        if near_left and near_right:
                            left_dist = abs(ref_pos - bp_left)
                            right_dist = abs(ref_pos - bp_right)
                            bp_side = 'left' if left_dist <= right_dist else 'right'
                        else:
                            bp_side = 'left' if near_left else 'right'

                        bp_target = bp_left if bp_side == 'left' else bp_right
                        bp_dist = abs(ref_pos - bp_target)

                        insertion_reads.append({
                            'read_name': read.query_name,
                            'ins_pos': ref_pos,
                            'bp_pos_obs': ref_pos,
                            'bp_side': bp_side,
                            'ins_length': length,
                            'bp_dist': bp_dist
                        })
                    # Insertion does not consume reference, so ref_pos stays the same

                elif op in (0, 2, 3, 7, 8):  # Consumes reference: M, D, N, =, X
                    ref_pos += length
                # Note: I(1), S(4), H(5), P(6) do not consume reference

    except Exception as e:
        logger.warning(f'Error finding insertions for {chrom}:{bp_left}-{bp_right}: {e}')

    return insertion_reads


def lightweight_softclip_te_align(softclip_seq, te_seq, min_match_frac=0.6):
    """
    Lightweight alignment of softclip sequence to TE using k-mer matching.
    """
    if not softclip_seq or not te_seq:
        return {'aligned': False, 'match_frac': 0}

    k = min(10, len(softclip_seq), len(te_seq))
    total_positions = len(softclip_seq) - k + 1
    if total_positions <= 0:
        return {'aligned': False, 'match_frac': 0}

    matches = sum(1 for i in range(total_positions) if softclip_seq[i:i+k] in te_seq)
    match_frac = matches / total_positions

    return {
        'aligned': match_frac >= min_match_frac,
        'match_frac': match_frac
    }


def calculate_alt_support_parent(
    bam, chrom, bp_left, bp_right, wiggle=30, use_te_enhance=False, te_seq=None,
    min_softclip_len=30, min_insertion_len=50, min_mapq=20, dedup_by_read=True
):
    """
    Calculate alternative support from parent BAM.

    Returns:
        dict with parent alt support counts and breakpoint consistency metrics
    """
    result = {
        'softclip_count': 0,
        'insertion_count': 0,
        'softclip_left': 0,
        'softclip_right': 0,
        'insertion_left': 0,
        'insertion_right': 0,
        'total_left': 0,
        'total_right': 0,
        'has_both_sides': False,
        'bp_position_std': 0.0,
        'bp_position_iqr': 0.0,
        'total_alt': 0
    }
    bp_positions = []

    def _init_event():
        return {
            'types': set(),
            'sides': set(),
            'softclip_sides': set(),
            'insertion_sides': set(),
            'positions': []
        }

    event_by_read = defaultdict(_init_event)

    def _record_event(read_name, evidence_type, side, bp_pos):
        evt = event_by_read[read_name]
        evt['types'].add(evidence_type)
        if side in ('left', 'right'):
            evt['sides'].add(side)
            if evidence_type == 'softclip':
                evt['softclip_sides'].add(side)
            elif evidence_type == 'insertion':
                evt['insertion_sides'].add(side)
        if bp_pos is not None:
            evt['positions'].append(bp_pos)

    # Soft-clip support
    softclip_reads = find_softclip_reads(
        bam, chrom, bp_left, bp_right, wiggle,
        min_softclip_len=min_softclip_len, min_mapq=min_mapq
    )
    for sc in softclip_reads:
        is_supporting = True
        if use_te_enhance and te_seq and sc.get('seq'):
            align_result = lightweight_softclip_te_align(sc['seq'], te_seq)
            is_supporting = align_result['aligned']
        if is_supporting:
            side = sc.get('bp_side')
            bp_pos = sc.get('bp_pos_obs')
            if dedup_by_read:
                _record_event(sc['read_name'], 'softclip', side, bp_pos)
            else:
                result['softclip_count'] += 1
                if side == 'left':
                    result['softclip_left'] += 1
                elif side == 'right':
                    result['softclip_right'] += 1
                if bp_pos is not None:
                    bp_positions.append(bp_pos)

    # Insertion CIGAR support
    insertion_reads = find_insertion_cigar_reads(
        bam, chrom, bp_left, bp_right, wiggle,
        min_insertion_len=min_insertion_len, min_mapq=min_mapq
    )
    if not dedup_by_read:
        result['insertion_count'] = len(insertion_reads)
    for ins in insertion_reads:
        side = ins.get('bp_side')
        bp_pos = ins.get('bp_pos_obs')
        if dedup_by_read:
            _record_event(ins['read_name'], 'insertion', side, bp_pos)
        else:
            if side == 'left':
                result['insertion_left'] += 1
            elif side == 'right':
                result['insertion_right'] += 1
            if bp_pos is not None:
                bp_positions.append(bp_pos)

    if dedup_by_read:
        for evt in event_by_read.values():
            if 'softclip' in evt['types']:
                result['softclip_count'] += 1
            if 'insertion' in evt['types']:
                result['insertion_count'] += 1
            if 'left' in evt['softclip_sides']:
                result['softclip_left'] += 1
            if 'right' in evt['softclip_sides']:
                result['softclip_right'] += 1
            if 'left' in evt['insertion_sides']:
                result['insertion_left'] += 1
            if 'right' in evt['insertion_sides']:
                result['insertion_right'] += 1
            if evt['positions']:
                bp_positions.append(int(np.median(evt['positions'])))

        result['total_alt'] = len(event_by_read)
        result['total_left'] = sum(1 for evt in event_by_read.values() if 'left' in evt['sides'])
        result['total_right'] = sum(1 for evt in event_by_read.values() if 'right' in evt['sides'])
    else:
        result['total_alt'] = result['softclip_count'] + result['insertion_count']
        result['total_left'] = result['softclip_left'] + result['insertion_left']
        result['total_right'] = result['softclip_right'] + result['insertion_right']

    result['has_both_sides'] = result['total_left'] > 0 and result['total_right'] > 0

    if len(bp_positions) >= 2:
        result['bp_position_std'] = float(np.std(bp_positions))
        result['bp_position_iqr'] = float(calculate_iqr(bp_positions))

    return result


def _normalize_candidate(cand):
    """
    Normalize tldr output dict to internal format with type conversion.

    Args:
        cand: dict from tldr --denovo output

    Returns:
        dict with normalized fields, or None if invalid
    """
    bp_left_val = cand.get('bp_left')
    bp_right_val = cand.get('bp_right')
    te_family_val = cand.get('TE_family')
    candidate_wiggle_val = cand.get('wiggle')
    raw_filter = cand.get('Filter')

    # Determine breakpoint values with fallback
    bp_left = bp_left_val if bp_left_val is not None else cand.get('Start')
    bp_right = bp_right_val if bp_right_val is not None else cand.get('End')

    # Convert to int with validation
    try:
        # Handle 'NA' and empty strings
        if bp_left == 'NA' or bp_left == '' or bp_left is None:
            logger.warning(f"Invalid bp_left for candidate {cand.get('UUID')}: {bp_left}")
            return None
        if bp_right == 'NA' or bp_right == '' or bp_right is None:
            logger.warning(f"Invalid bp_right for candidate {cand.get('UUID')}: {bp_right}")
            return None

        bp_left = int(bp_left)
        bp_right = int(bp_right)

        # Sanity check
        if bp_left < 0 or bp_right < 0:
            logger.warning(f"Negative coordinates for {cand.get('UUID')}: {bp_left}, {bp_right}")
            return None
        if bp_right < bp_left:
            logger.warning(f"Invalid coordinate order for {cand.get('UUID')}: {bp_left} > {bp_right}")
            return None

    except (ValueError, TypeError) as e:
        logger.warning(f"Failed to convert breakpoints for {cand.get('UUID')}: {e}")
        return None

    candidate_wiggle = _as_positive_int_or_none(candidate_wiggle_val)
    te_family = te_family_val if te_family_val is not None else cand.get('Family')
    if isinstance(te_family, str):
        te_family = te_family.strip()

    chrom = cand.get('Chrom') or cand.get('Chromosome')
    if not chrom:
        logger.warning(f"Missing chromosome for candidate {cand.get('UUID')}")
        return None

    return {
        'uuid': cand.get('UUID'),
        'chrom': chrom,
        'bp_left': bp_left,
        'bp_right': bp_right,
        'strand': cand.get('Strand', '+'),
        'te_family': te_family,
        'child_support': cand.get('child_support'),
        'wiggle': candidate_wiggle,
        'raw_filter': raw_filter
    }


def _candidate_is_eligible(cand, allow_na_te_family=False, require_filter_pass=False):
    """
    Decide if normalized candidate should enter parent-side analysis.
    """
    te_family = (cand.get('te_family') or '').strip().upper()
    raw_filter = (cand.get('raw_filter') or '').strip().upper()

    if not allow_na_te_family and te_family in NA_STRINGS:
        return False, 'TE_family is NA'

    if require_filter_pass and raw_filter and raw_filter != 'PASS':
        return False, f'Filter={cand.get("raw_filter")}'

    return True, None


def analyze_parent_support(
    mom_bam_path, dad_bam_path, child_candidates,
    wiggle=80, window_size=100, use_te_enhance=False, te_family_seqs=None,
    min_softclip_len=30, min_insertion_len=50, min_mapq=20, dedup_by_read=True,
    allow_na_te_family=False, require_filter_pass=False
):
    """
    Analyze parent support for child candidate sites.

    Args:
        mom_bam_path: path to maternal BAM file
        dad_bam_path: path to paternal BAM file
        child_candidates: list of dicts (tldr --denovo output format)
        wiggle: default window around breakpoint for alt support search
        window_size: window size for depth calculation
        use_te_enhance: whether to use TE alignment enhancement
        te_family_seqs: dict mapping TE family names to sequences (optional)
        min_softclip_len: minimum soft-clip size to count as alt
        min_insertion_len: minimum CIGAR insertion size to count as alt
        min_mapq: minimum MAPQ for parent-side evidence
        dedup_by_read: if True, count each read at most once in total alt support
        allow_na_te_family: if True, keep candidates with TE_family=NA
        require_filter_pass: if True, keep only candidates with Filter=PASS

    Returns:
        list of dicts with parent support analysis results
    """
    results = []

    # Normalize candidates and filter out invalid ones
    normalized = [_normalize_candidate(c) for c in child_candidates]
    normalized = [n for n in normalized if n is not None]

    filtered = []
    dropped = defaultdict(int)
    for cand in normalized:
        keep, reason = _candidate_is_eligible(
            cand,
            allow_na_te_family=allow_na_te_family,
            require_filter_pass=require_filter_pass
        )
        if keep:
            filtered.append(cand)
        else:
            dropped[reason] += 1

    if dropped:
        summary = ', '.join(f'{k}: {v}' for k, v in sorted(dropped.items()))
        logger.info(f'Dropped {sum(dropped.values())} candidates before parent analysis ({summary})')

    if not filtered:
        logger.error('No valid candidates after normalization/filtering')
        return results

    mom_bam = None
    dad_bam = None

    try:
        mom_bam = pysam.AlignmentFile(mom_bam_path, 'rb')
        dad_bam = pysam.AlignmentFile(dad_bam_path, 'rb')

        logger.info(f'Analyzing {len(filtered)} valid candidates from Mom: {mom_bam_path}, Dad: {dad_bam_path}')

        for i, cand in enumerate(filtered):
            if (i + 1) % 100 == 0:
                logger.info(f'Processing candidate {i + 1}/{len(filtered)}')

            chrom = cand['chrom']
            bp_left = cand['bp_left']
            bp_right = cand['bp_right']
            te_family = cand['te_family']
            candidate_wiggle = cand.get('wiggle')
            analysis_wiggle = candidate_wiggle if candidate_wiggle is not None else wiggle

            # Get TE sequence for enhancement
            te_seq = None
            if use_te_enhance and te_family and te_family_seqs:
                te_seq = resolve_te_sequence(te_family, te_family_seqs)

            result = {
                'uuid': cand['uuid'],
                'chrom': chrom,
                'bp_left': bp_left,
                'bp_right': bp_right,
                'te_family': te_family,
                'child_support': cand.get('child_support'),
                'analysis_wiggle': analysis_wiggle
            }

            # Mom analysis
            mom_depth = calculate_window_depth(mom_bam, chrom, bp_left, bp_right, window_size)
            mom_alt = calculate_alt_support_parent(mom_bam, chrom, bp_left, bp_right,
                                                    wiggle=analysis_wiggle, use_te_enhance=use_te_enhance, te_seq=te_seq,
                                                    min_softclip_len=min_softclip_len, min_insertion_len=min_insertion_len,
                                                    min_mapq=min_mapq, dedup_by_read=dedup_by_read)

            result['mom_depth_mean'] = mom_depth['depth_mean']
            result['mom_softclip'] = mom_alt['softclip_count']
            result['mom_insertion'] = mom_alt['insertion_count']
            result['mom_total_alt'] = mom_alt['total_alt']
            result['mom_total_left'] = mom_alt['total_left']
            result['mom_total_right'] = mom_alt['total_right']
            result['mom_has_both_sides'] = mom_alt['has_both_sides']
            result['mom_bp_position_std'] = mom_alt['bp_position_std']
            result['mom_bp_position_iqr'] = mom_alt['bp_position_iqr']

            # Dad analysis
            dad_depth = calculate_window_depth(dad_bam, chrom, bp_left, bp_right, window_size)
            dad_alt = calculate_alt_support_parent(dad_bam, chrom, bp_left, bp_right,
                                                    wiggle=analysis_wiggle, use_te_enhance=use_te_enhance, te_seq=te_seq,
                                                    min_softclip_len=min_softclip_len, min_insertion_len=min_insertion_len,
                                                    min_mapq=min_mapq, dedup_by_read=dedup_by_read)

            result['dad_depth_mean'] = dad_depth['depth_mean']
            result['dad_softclip'] = dad_alt['softclip_count']
            result['dad_insertion'] = dad_alt['insertion_count']
            result['dad_total_alt'] = dad_alt['total_alt']
            result['dad_total_left'] = dad_alt['total_left']
            result['dad_total_right'] = dad_alt['total_right']
            result['dad_has_both_sides'] = dad_alt['has_both_sides']
            result['dad_bp_position_std'] = dad_alt['bp_position_std']
            result['dad_bp_position_iqr'] = dad_alt['bp_position_iqr']

            results.append(result)

        logger.info(f'Completed analysis for {len(results)} candidates')

    finally:
        if mom_bam:
            mom_bam.close()
        if dad_bam:
            dad_bam.close()

    return results


def load_candidates_from_tldr_table(tldr_table_file):
    """
    Load child candidates from tldr output table.

    Args:
        tldr_table_file: path to tldr .table.txt file

    Returns:
        list of dicts suitable for analyze_parent_support()
    """
    candidates = []

    with open(tldr_table_file, 'r') as f:
        header = f.readline().strip().split('\t')

        for line in f:
            cols = line.strip().split('\t')
            rec = dict(zip(header, cols))
            candidates.append(rec)

    logger.info(f'Loaded {len(candidates)} candidates from {tldr_table_file}')
    return candidates


def load_te_family_sequences(te_fasta):
    """
    Load TE family sequences from FASTA.

    FASTA header first token is used as key, e.g. '>L1HS' -> 'L1HS'.
    """
    te_seqs = {}
    current_name = None
    current_seq = []

    with open(te_fasta, 'r') as f:
        for raw_line in f:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_name and current_seq:
                    seq = ''.join(current_seq).upper()
                    for label in _canonical_te_family_labels(current_name):
                        te_seqs.setdefault(label, seq)
                current_name = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)

    if current_name and current_seq:
        seq = ''.join(current_seq).upper()
        for label in _canonical_te_family_labels(current_name):
            te_seqs.setdefault(label, seq)

    logger.info(f'Loaded {len(te_seqs)} TE family sequences from {te_fasta}')
    return te_seqs


def calculate_iqr(values):
    """Calculate interquartile range."""
    if len(values) < 2:
        return float('inf')
    sorted_vals = sorted(values)
    q1 = sorted_vals[int(len(sorted_vals) * 0.25)]
    q3 = sorted_vals[int(len(sorted_vals) * 0.75)]
    return q3 - q1


def evaluate_denovo(parent_results, m_child=1,
                    min_parent_depth_default=8, min_parent_depth_frac=0.25,
                    max_parent_alt=0, require_both_sides=False,
                    bp_position_std_threshold=15, bp_position_iqr_threshold=20,
                    adaptive_bp_thresholds=True,
                    bp_std_wiggle_factor=0.25, bp_iqr_wiggle_factor=0.5):
    """
    Evaluate de novo status for candidate insertions.

    Thresholds:
        min_parent_depth = max(8, 0.25 * median_parent_depth)
        max_parent_alt = 0 (strict)

    Categories:
        PASS_DENOVO: Meets all criteria for de novo
        UNCERTAIN: Insufficient parent depth / child support / diffuse parent alt
        FAIL: Parent shows coherent alt support above threshold

    Args:
        parent_results: list from analyze_parent_support()
        m_child: minimum child support reads
        min_parent_depth_default: minimum depth floor (default 8)
        min_parent_depth_frac: fraction of median parent depth (default 0.25)
        max_parent_alt: maximum allowed alt reads in parent (default 0)
        require_both_sides: require parent alt evidence from both left/right breakpoints to call FAIL
        bp_position_std_threshold: max std of parent alt breakpoint positions to call FAIL
        bp_position_iqr_threshold: max IQR of parent alt breakpoint positions to call FAIL
        adaptive_bp_thresholds: adapt position thresholds to candidate wiggle
        bp_std_wiggle_factor: std threshold lower bound as fraction of wiggle
        bp_iqr_wiggle_factor: iqr threshold lower bound as fraction of wiggle

    Returns:
        list of dicts with evaluation results
    """
    # Calculate adaptive depth threshold
    all_depths = []
    for p in parent_results:
        if p.get('mom_depth_mean', 0) > 0:
            all_depths.append(p['mom_depth_mean'])
        if p.get('dad_depth_mean', 0) > 0:
            all_depths.append(p['dad_depth_mean'])

    median_depth = np.median(all_depths) if all_depths else 32
    min_parent_depth = max(min_parent_depth_default, min_parent_depth_frac * median_depth)

    logger.info(f'De novo evaluation: min_parent_depth={min_parent_depth:.1f} (median={median_depth:.1f})')

    results = []

    for p in parent_results:
        result = {
            'uuid': p.get('uuid'),
            'chrom': p.get('chrom'),
            'bp_left': p.get('bp_left'),
            'bp_right': p.get('bp_right'),
            'te_family': p.get('te_family'),
            'child_support': p.get('child_support'),
            'analysis_wiggle': p.get('analysis_wiggle')
        }

        mom_depth = p.get('mom_depth_mean', 0)
        dad_depth = p.get('dad_depth_mean', 0)
        mom_alt = p.get('mom_total_alt', 0)
        dad_alt = p.get('dad_total_alt', 0)
        child_support_count = parse_child_support_count(p.get('child_support'))

        mom_total_left = p.get('mom_total_left', 0)
        mom_total_right = p.get('mom_total_right', 0)
        dad_total_left = p.get('dad_total_left', 0)
        dad_total_right = p.get('dad_total_right', 0)

        mom_has_both_sides = bool(p.get('mom_has_both_sides', mom_total_left > 0 and mom_total_right > 0))
        dad_has_both_sides = bool(p.get('dad_has_both_sides', dad_total_left > 0 and dad_total_right > 0))

        mom_bp_std = float(p.get('mom_bp_position_std', 0.0) or 0.0)
        dad_bp_std = float(p.get('dad_bp_position_std', 0.0) or 0.0)
        mom_bp_iqr = float(p.get('mom_bp_position_iqr', 0.0) or 0.0)
        dad_bp_iqr = float(p.get('dad_bp_position_iqr', 0.0) or 0.0)
        analysis_wiggle = _as_positive_int_or_none(p.get('analysis_wiggle'))

        std_threshold = bp_position_std_threshold
        iqr_threshold = bp_position_iqr_threshold
        if adaptive_bp_thresholds and analysis_wiggle is not None:
            std_threshold = max(std_threshold, float(analysis_wiggle) * float(bp_std_wiggle_factor))
            iqr_threshold = max(iqr_threshold, float(analysis_wiggle) * float(bp_iqr_wiggle_factor))

        result['min_parent_depth'] = min_parent_depth
        result['mom_depth'] = mom_depth
        result['dad_depth'] = dad_depth
        result['mom_alt'] = mom_alt
        result['dad_alt'] = dad_alt
        result['child_support_count'] = child_support_count
        result['mom_total_left'] = mom_total_left
        result['mom_total_right'] = mom_total_right
        result['dad_total_left'] = dad_total_left
        result['dad_total_right'] = dad_total_right
        result['mom_has_both_sides'] = mom_has_both_sides
        result['dad_has_both_sides'] = dad_has_both_sides
        result['mom_bp_position_std'] = mom_bp_std
        result['dad_bp_position_std'] = dad_bp_std
        result['mom_bp_position_iqr'] = mom_bp_iqr
        result['dad_bp_position_iqr'] = dad_bp_iqr
        result['bp_std_threshold'] = std_threshold
        result['bp_iqr_threshold'] = iqr_threshold

        mom_covered = mom_depth >= min_parent_depth
        dad_covered = dad_depth >= min_parent_depth

        evaluation = 'FAIL'
        fail_reasons = []
        uncertain_reasons = []

        # Check parent alt support with optional side/consistency constraints
        if mom_alt > max_parent_alt:
            mom_side_ok = (not require_both_sides) or mom_has_both_sides
            mom_pos_ok = (mom_bp_std <= std_threshold) and (mom_bp_iqr <= iqr_threshold)
            if mom_side_ok and mom_pos_ok:
                fail_reasons.append(f'Mom alt too high ({mom_alt})')
            else:
                if require_both_sides and not mom_has_both_sides:
                    uncertain_reasons.append('Mom alt present but not supported on both breakpoints')
                if not mom_pos_ok:
                    uncertain_reasons.append(
                        f'Mom alt is position-diffuse (std={mom_bp_std:.1f}, iqr={mom_bp_iqr:.1f}, '
                        f'thr={std_threshold:.1f}/{iqr_threshold:.1f})'
                    )

        if dad_alt > max_parent_alt:
            dad_side_ok = (not require_both_sides) or dad_has_both_sides
            dad_pos_ok = (dad_bp_std <= std_threshold) and (dad_bp_iqr <= iqr_threshold)
            if dad_side_ok and dad_pos_ok:
                fail_reasons.append(f'Dad alt too high ({dad_alt})')
            else:
                if require_both_sides and not dad_has_both_sides:
                    uncertain_reasons.append('Dad alt present but not supported on both breakpoints')
                if not dad_pos_ok:
                    uncertain_reasons.append(
                        f'Dad alt is position-diffuse (std={dad_bp_std:.1f}, iqr={dad_bp_iqr:.1f}, '
                        f'thr={std_threshold:.1f}/{iqr_threshold:.1f})'
                    )

        if fail_reasons:
            evaluation = 'FAIL'
            reasons = fail_reasons
        else:
            # Check depth coverage
            if not mom_covered and not dad_covered:
                uncertain_reasons.append('Both parents have insufficient depth')
            elif not mom_covered:
                uncertain_reasons.append(f'Mom depth too low ({mom_depth:.1f} < {min_parent_depth:.1f})')
            elif not dad_covered:
                uncertain_reasons.append(f'Dad depth too low ({dad_depth:.1f} < {min_parent_depth:.1f})')

            # Child support threshold check (if parseable from child_support field)
            if m_child > 0 and child_support_count is not None and child_support_count < m_child:
                uncertain_reasons.append(f'Child support too low ({child_support_count} < {m_child})')

            if uncertain_reasons:
                evaluation = 'UNCERTAIN'
                reasons = uncertain_reasons
            else:
                evaluation = 'PASS_DENOVO'
                reasons = ['Passed all criteria']

        result['evaluation'] = evaluation
        result['reasons'] = reasons
        results.append(result)

    # Summary
    eval_counts = {}
    for r in results:
        eval_counts[r['evaluation']] = eval_counts.get(r['evaluation'], 0) + 1
    logger.info(f'De novo evaluation summary: {eval_counts}')

    return results


def filter_candidates(denovo_results, evaluation='PASS_DENOVO'):
    """Filter candidates by evaluation result."""
    return [r for r in denovo_results if r.get('evaluation') == evaluation]


if __name__ == '__main__':
    import argparse
    import json

    parser = argparse.ArgumentParser(description='Parent support analysis for tldr candidates')
    parser.add_argument('--mom_bam', required=True, help='Maternal BAM file')
    parser.add_argument('--dad_bam', required=True, help='Paternal BAM file')
    parser.add_argument('--candidates', required=True, help='tldr table.txt or JSON file')
    parser.add_argument('--output', required=True, help='Output JSON file')
    parser.add_argument('--wiggle', type=int, default=80, help='Default window around breakpoint (candidate wiggle has priority)')
    parser.add_argument('--depth_window', type=int, default=100, help='Window for depth')
    parser.add_argument('--min_softclip_len', type=int, default=30,
                       help='Minimum soft-clip length to count as parent alt support')
    parser.add_argument('--min_insertion_len', type=int, default=50,
                       help='Minimum CIGAR insertion length to count as parent alt support')
    parser.add_argument('--min_mapq', type=int, default=20,
                       help='Minimum MAPQ to count parent alt support')
    parser.add_argument('--no_dedup_by_read', action='store_true',
                       help='Disable read-level de-duplication for parent alt support counting')
    parser.add_argument('--allow_na_te_family', action='store_true',
                       help='Keep candidates where TE_family is NA')
    parser.add_argument('--require_filter_pass', action='store_true',
                       help='Keep only candidates with Filter=PASS from tldr output')
    parser.add_argument('--te_enhance', action='store_true',
                       help='Use TE sequence-aware soft-clip filtering')
    parser.add_argument('--te_fasta',
                       help='FASTA with TE family sequences (header should match TE_family)')
    parser.add_argument('--denovo', action='store_true', help='Run de novo evaluation')
    parser.add_argument('--m_child', type=int, default=1,
                       help='Minimum child support reads (parsed from child_support)')
    parser.add_argument('--min_parent_depth_default', type=float, default=8,
                       help='Minimum parent depth floor')
    parser.add_argument('--min_parent_depth_frac', type=float, default=0.25,
                       help='Fraction of median parent depth for adaptive threshold')
    parser.add_argument('--max_parent_alt', type=int, default=0,
                       help='Maximum allowed parent alt reads')
    parser.add_argument('--require_both_sides', action='store_true',
                       help='Require parent alt support on both left/right breakpoints to call FAIL')
    parser.add_argument('--bp_position_std_threshold', type=float, default=15,
                       help='Max breakpoint position standard deviation for coherent parent alt support')
    parser.add_argument('--bp_position_iqr_threshold', type=float, default=20,
                       help='Max breakpoint position IQR for coherent parent alt support')
    parser.add_argument('--no_adaptive_bp_thresholds', action='store_true',
                       help='Use fixed bp std/iqr thresholds instead of wiggle-adaptive thresholds')
    parser.add_argument('--bp_std_wiggle_factor', type=float, default=0.25,
                       help='Adaptive std threshold lower bound = wiggle * factor')
    parser.add_argument('--bp_iqr_wiggle_factor', type=float, default=0.5,
                       help='Adaptive iqr threshold lower bound = wiggle * factor')
    parser.add_argument('--child_bam', help='Child BAM file (required for visualization)')
    parser.add_argument('--visualize', action='store_true',
                       help='Generate IGV-style visualizations for candidates')
    parser.add_argument('--visualize_dir', default='denovo_visualizations',
                       help='Directory for visualization output')
    parser.add_argument('--visualize_window', type=int, default=500,
                       help='Base pairs to show on each side of insertion')
    parser.add_argument('--visualize_filter', choices=['PASS_DENOVO', 'UNCERTAIN', 'FAIL'],
                       help='Only visualize candidates with this evaluation')

    args = parser.parse_args()

    # Validate visualization arguments
    if args.visualize and not args.child_bam:
        parser.error('--visualize requires --child_bam')

    # Load candidates
    if args.candidates.endswith('.txt'):
        candidates = load_candidates_from_tldr_table(args.candidates)
    else:
        with open(args.candidates) as f:
            candidates = json.load(f)

    te_family_seqs = None
    if args.te_enhance:
        if args.te_fasta:
            te_family_seqs = load_te_family_sequences(args.te_fasta)
        else:
            logger.warning('--te_enhance enabled without --te_fasta; TE sequence filter will be skipped')

    # Run analysis
    results = analyze_parent_support(
        args.mom_bam, args.dad_bam, candidates,
        wiggle=args.wiggle,
        window_size=args.depth_window,
        use_te_enhance=args.te_enhance,
        te_family_seqs=te_family_seqs,
        min_softclip_len=args.min_softclip_len,
        min_insertion_len=args.min_insertion_len,
        min_mapq=args.min_mapq,
        dedup_by_read=(not args.no_dedup_by_read),
        allow_na_te_family=args.allow_na_te_family,
        require_filter_pass=args.require_filter_pass
    )

    if args.denovo:
        results = evaluate_denovo(
            results,
            m_child=args.m_child,
            min_parent_depth_default=args.min_parent_depth_default,
            min_parent_depth_frac=args.min_parent_depth_frac,
            max_parent_alt=args.max_parent_alt,
            require_both_sides=args.require_both_sides,
            bp_position_std_threshold=args.bp_position_std_threshold,
            bp_position_iqr_threshold=args.bp_position_iqr_threshold,
            adaptive_bp_thresholds=(not args.no_adaptive_bp_thresholds),
            bp_std_wiggle_factor=args.bp_std_wiggle_factor,
            bp_iqr_wiggle_factor=args.bp_iqr_wiggle_factor
        )

    with open(args.output, 'w') as f:
        json.dump(results, f, indent=2)

    print(f'Saved {len(results)} results to {args.output}')

    # Generate visualizations if requested
    if args.visualize:
        try:
            from visualize_denovo import visualize_from_parent_analysis

            print(f'\nGenerating IGV-style visualizations...')
            visualize_from_parent_analysis(
                args.output,
                args.child_bam,
                args.mom_bam,
                args.dad_bam,
                output_dir=args.visualize_dir,
                filter_evaluation=args.visualize_filter
            )

            # Print summary
            filter_desc = f" (filter: {args.visualize_filter})" if args.visualize_filter else ""
            print(f'\nVisualization complete! Output saved to: {args.visualize_dir}/{filter_desc}')

        except ImportError as e:
            print(f'Warning: Could not import visualization module: {e}')
            print('To enable visualization, ensure visualize_denovo.py is in the same directory.')
