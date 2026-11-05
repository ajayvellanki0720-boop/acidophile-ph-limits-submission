#!/bin/sh
set -eu

mode=${1:-verify}

case "$mode" in
  verify)
    sha256sum -c CHECKSUMS.sha256
    python3 -m pytest -q
    python3 scripts/69_document_consistency.py
    python3 scripts/77_audit.py
    ;;
  figures)
    python3 scripts/100_figures.py
    ;;
  manuscript)
    python3 scripts/105_submission_docx.py
    python3 scripts/108_peerj_submission.py
    ;;
  *)
    echo "usage: $0 {verify|figures|manuscript}" >&2
    exit 2
    ;;
esac
