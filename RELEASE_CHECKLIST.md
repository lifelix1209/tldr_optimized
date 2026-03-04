# Release Checklist

Use this checklist before creating a GitHub release.

## 1. Code and Docs

- [ ] `README.md` reflects current installation and CLI behavior
- [ ] `doc/denovo_evaluation.md` updated for de novo logic changes
- [ ] `CHANGELOG`/release notes draft prepared (GitHub release notes is acceptable)
- [ ] No local-only files are tracked

## 2. Versioning

- [ ] `tldr/tldr` `__version__` updated
- [ ] `pyproject.toml` version updated

## 3. Validation

- [ ] CI passes on `main`
- [ ] Local syntax checks pass:
  - `python -m py_compile tldr/tldr`
  - `find scripts test test_denovo -name "*.py" -print0 | xargs -0 python -m py_compile`
- [ ] Core smoke test command executed on representative data

## 4. Tag and Publish

- [ ] Create annotated tag: `git tag -a vX.Y.Z -m "Release vX.Y.Z"`
- [ ] Push commits and tag: `git push origin main --tags`
- [ ] Create GitHub Release from the tag
- [ ] Attach release notes with major changes and any breaking behavior
