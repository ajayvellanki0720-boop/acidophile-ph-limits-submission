#!/bin/sh
set -eu

mode=${1:-verify}

case "$mode" in
  verify)
    sha256sum -c CHECKSUMS.sha256
    python -m pytest -q
    python scripts/69_document_consistency.py
    python scripts/77_audit.py
    ;;
  figures)
    python scripts/100_figures.py
    ;;
  manuscript)
    python scripts/105_submission_docx.py
    ;;
  *)
    echo "usage: $0 {verify|figures|manuscript}" >&2
    exit 2
    ;;
esac
