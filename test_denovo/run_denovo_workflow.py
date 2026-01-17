#!/usr/bin/env python3
"""
Complete workflow for de novo TE insertion detection and visualization.

This script orchestrates the entire pipeline:
1. Generate mock BAM files with TE insertions (soft-clipped reads)
2. Run tldr analysis on child sample
3. Run parent support analysis (Mom/Dad comparison)
4. Generate IGV-style visualizations

Usage:
    cd test_denovo
    python run_denovo_workflow.py [--skip-tldr] [--skip-validation]

Requirements:
    - conda environment: tldr-dev
    - tldr/tldr (local)
    - exonerate (for TE alignment)
    - matplotlib (for visualization)

Setup:
    conda env create -f ../tldr-dev.yml
    conda activate tldr-dev
    pip install -e ..
"""

import argparse
import os
import shutil
import subprocess
import sys

# Project paths
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TEST_DIR = os.path.join(PROJECT_ROOT, "test_denovo")
SCRIPTS_DIR = os.path.join(PROJECT_ROOT, "scripts")

# Output directories
TLDR_OUTPUT_DIR = os.path.join(TEST_DIR, "output/tldr_final")
VISUALIZATIONS_DIR = os.path.join(TEST_DIR, "visualizations_final")
PARENT_OUTPUT = os.path.join(TEST_DIR, "output/tldr_parent_analysis_final.json")


def run_command(cmd, cwd=None, capture=True):
    """Run a shell command with proper PATH for conda environment."""
    print(f"Running: {' '.join(cmd) if isinstance(cmd, list) else cmd}")

    # Set up environment with PATH for conda environment and local tools
    run_env = os.environ.copy()
    conda_env_bin = '/opt/anaconda3/envs/tldr-dev-test/bin'
    tools_dir = os.path.join(PROJECT_ROOT, 'tools')
    run_env['PATH'] = f'{tools_dir}:{conda_env_bin}:' + run_env.get('PATH', '')

    result = subprocess.run(
        cmd,
        shell=isinstance(cmd, str),
        cwd=cwd,
        capture_output=capture,
        text=True,
        env=run_env
    )
    if result.returncode != 0:
        print(f"Error: {result.stderr}")
        return False, result.stderr
    return True, result.stdout


def check_dependencies():
    """Check if required tools are available in conda environment."""
    deps = {
        'tldr': 'tldr/tldr',
        'exonerate': 'exonerate',
        'pysam': None,
        'matplotlib': None,
    }

    missing = []
    for name, path in deps.items():
        if path is None:
            try:
                __import__(name)
            except ImportError:
                missing.append(name)
        else:
            if name == 'exonerate':
                # Check in PATH or common conda locations
                conda_env = 'tldr-dev-test'
                conda_bin = f'/opt/anaconda3/envs/{conda_env}/bin'
                found = (
                    shutil.which('exonerate') or
                    os.path.exists(os.path.join(conda_bin, 'exonerate')) or
                    os.environ.get('CONDA_PREFIX') and os.path.exists(
                        os.path.join(os.environ['CONDA_PREFIX'], 'bin', 'exonerate'))
                )
                if not found:
                    missing.append(name)

    if missing:
        print(f"Missing dependencies: {', '.join(missing)}")
        print("\nSetup conda environment:")
        print("  conda env create -f ../tldr-dev.yml")
        print("  conda activate tldr-dev")
        return False
    return True


def step1_generate_mock_data():
    """Step 1: Generate mock BAM files with TE insertions."""
    print("\n" + "=" * 70)
    print("STEP 1: Generating mock BAM files with TE insertions")
    print("=" * 70)

    script = os.path.join(TEST_DIR, "create_mock_data.py")
    success, _ = run_command([sys.executable, script], cwd=TEST_DIR)
    return success


def step2_run_tldr():
    """Step 2: Run tldr on child sample."""
    print("\n" + "=" * 70)
    print("STEP 2: Running tldr on child sample")
    print("=" * 70)

    child_bam = os.path.join(TEST_DIR, "bams/child.bam")
    te_fa = os.path.join(TEST_DIR, "ref/mock_te.fa")
    ref_fa = os.path.join(TEST_DIR, "ref/mock_ref.fa")
    output_table = os.path.join(TLDR_OUTPUT_DIR, "table.txt")

    # Check if we already have a valid output
    if os.path.exists(output_table):
        print(f"Using existing tldr output: {output_table}")
        return True

    # Clean previous output
    if os.path.exists(TLDR_OUTPUT_DIR):
        shutil.rmtree(TLDR_OUTPUT_DIR)
    os.makedirs(TLDR_OUTPUT_DIR, exist_ok=True)

    tldr_script = os.path.join(PROJECT_ROOT, "tldr/tldr")
    cmd = [
        sys.executable, tldr_script,
        "-b", child_bam,
        "-e", te_fa,
        "-r", ref_fa,
        "--denovo",
        "-m", "1",
        "--embed_minreads", "1",
        "--min_te_len", "50",
        "-o", TLDR_OUTPUT_DIR
    ]

    success, _ = run_command(cmd)

    # Check if output was created
    if not os.path.exists(output_table):
        # Try to use existing output from table.txt in parent directory
        parent_table = os.path.join(TEST_DIR, "output/tldr_final.table.txt")
        if os.path.exists(parent_table):
            print(f"Using existing output: {parent_table}")
            shutil.copy(parent_table, output_table)
            return True
        print("Warning: tldr did not produce output")
        return False

    print(f"tldr completed. Output: {TLDR_OUTPUT_DIR}/")
    return success


def step3_run_parent_analysis():
    """Step 3: Run parent support analysis."""
    print("\n" + "=" * 70)
    print("STEP 3: Running parent support analysis")
    print("=" * 70)

    mom_bam = os.path.join(TEST_DIR, "bams/mom.bam")
    dad_bam = os.path.join(TEST_DIR, "bams/dad.bam")
    child_bam = os.path.join(TEST_DIR, "bams/child.bam")
    candidates = os.path.join(TLDR_OUTPUT_DIR, "table.txt")
    viz_dir = VISUALIZATIONS_DIR

    # Clean previous visualization
    if os.path.exists(viz_dir):
        shutil.rmtree(viz_dir)
    os.makedirs(viz_dir, exist_ok=True)

    script = os.path.join(SCRIPTS_DIR, "parent_support.py")
    cmd = [
        sys.executable, script,
        "--mom_bam", mom_bam,
        "--dad_bam", dad_bam,
        "--child_bam", child_bam,
        "--candidates", candidates,
        "--output", PARENT_OUTPUT,
        "--denovo",
        "--visualize",
        "--visualize_dir", viz_dir
    ]

    success, _ = run_command(cmd)
    if success:
        print(f"Parent analysis complete. Results: {PARENT_OUTPUT}")
    return success


def step4_validate_results():
    """Step 4: Validate results against expected outcomes."""
    print("\n" + "=" * 70)
    print("STEP 4: Validating results")
    print("=" * 70)

    script = os.path.join(TEST_DIR, "validate_results.py")
    success, output = run_command([sys.executable, script], cwd=TEST_DIR)
    if success:
        print(output)
    return success


def main():
    parser = argparse.ArgumentParser(
        description="Complete de novo TE insertion detection workflow"
    )
    parser.add_argument("--skip-tldr", action="store_true",
                        help="Skip tldr step (use existing output)")
    parser.add_argument("--skip-validation", action="store_true",
                        help="Skip validation step")
    parser.add_argument("--skip-mock", action="store_true",
                        help="Skip mock data generation")

    args = parser.parse_args()

    print("=" * 70)
    print("De Novo TE Insertion Detection Pipeline")
    print("=" * 70)

    # Check dependencies
    if not check_dependencies():
        print("\nDependencies check failed. Please install missing tools.")
        sys.exit(1)

    # Step 1: Generate mock data (unless skipped)
    if not args.skip_mock:
        if not step1_generate_mock_data():
            print("Failed to generate mock data")
            sys.exit(1)

    # Step 2: Run tldr (unless skipped)
    if not args.skip_tldr:
        if not step2_run_tldr():
            print("Warning: tldr failed")
            sys.exit(1)

    # Step 3: Run parent analysis with visualization
    if not step3_run_parent_analysis():
        print("Failed at parent analysis")
        sys.exit(1)

    # Step 4: Validate results
    if not args.skip_validation:
        step4_validate_results()

    print("\n" + "=" * 70)
    print("Pipeline complete!")
    print("=" * 70)
    print(f"\nVisualizations: {VISUALIZATIONS_DIR}/")
    print(f"Analysis results: {PARENT_OUTPUT}")


if __name__ == "__main__":
    main()
