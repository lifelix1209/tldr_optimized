# Changelog

All notable changes to this project should be documented in this file.

## [1.3.0] - 2026-03-04

### Added
- De novo parent-support v1 improvements in `scripts/parent_support.py`
  - Candidate wiggle-aware parent analysis
  - Parent evidence quality thresholds and read-level de-duplication
  - Better TE family matching for `--te_enhance`
  - Adaptive breakpoint consistency thresholds
- `doc/denovo_evaluation.md` updated to reflect current de novo logic
- Repository governance docs:
  - `CONTRIBUTING.md`
  - `SECURITY.md`
  - `RELEASE_CHECKLIST.md`
  - GitHub issue/PR templates

### Changed
- Packaging metadata/build config aligned for release readiness (`pyproject.toml`, `setup.py`)
- Development environment made reproducible (`tldr-dev.yml`)
- README installation/dev sections refreshed for this repository
