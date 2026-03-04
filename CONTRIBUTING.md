# Contributing

Thanks for contributing to `tldr_optimized`.

## Development Setup

```bash
git clone https://github.com/lifelix1209/tldr_optimized.git
cd tldr_optimized
conda env create -f tldr-dev.yml
conda activate tldr-dev
pip install -e .
```

## Local Checks

Run the same syntax checks used in CI:

```bash
python -m py_compile tldr/tldr
find scripts test test_denovo -name "*.py" -print0 | xargs -0 python -m py_compile
```

Optional test scripts:

```bash
bash test/run_test.sh
bash test/run_test_targets.sh
python test_denovo/run_denovo_workflow.py --skip-validation
```

## Pull Requests

1. Create a branch from `main`.
2. Keep PRs focused and small when possible.
3. Update documentation for user-facing behavior changes.
4. Include commands you ran for validation in the PR description.
5. Ensure CI passes before requesting review.

## Style

- Keep changes compatible with Python >= 3.9 unless discussed otherwise.
- Prefer clear names over compact code.
- Avoid introducing heavyweight dependencies without justification.
