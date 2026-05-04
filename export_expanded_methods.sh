#!/usr/bin/env bash
# Exports three Word documents from the supplementary source files:
#   1. expanded_methods.docx  — Expanded Methods text with inline citations resolved
#   2. supplementary_table3.docx — ST3 caption + table
#   3. references.docx        — Full reference list (all entries in the bib file)
# Run from the repo root: bash export_expanded_methods.sh

set -e

METHODS_TEX="Supplementary/Components/Sections/methods.tex"
ST3_CAPTION_TEX="Supplementary/Components/Sections/tables.tex"
ST3_TABLE_TEX="Supplementary/Components/Tables/ST3.tex"
BIB="Supplementary/Components/References/supplementary.bib"
CSL="Supplementary/Components/References/thyroid.csl"
OUT_DIR="/Users/JoshsMacbook2015/Library/CloudStorage/OneDrive-EmoryUniversity/Research/Manuscripts and Projects/Active Projects/Thyroid Exposomics Paper/2023_Smith_Thyroid Cancer/Env International Submission_Smith"

OUT_METHODS="$OUT_DIR/expanded_methods.docx"
OUT_ST3="$OUT_DIR/supplementary_table3.docx"
OUT_REFS="$OUT_DIR/references.docx"

#--- 1: Expanded Methods ---
TMP_METHODS=$(mktemp /tmp/exp_methods_XXXX.tex)
cat "$METHODS_TEX" > "$TMP_METHODS"
pandoc "$TMP_METHODS" \
  --from=latex \
  --to=docx \
  --citeproc \
  --bibliography="$BIB" \
  --csl="$CSL" \
  --output="$OUT_METHODS"
rm -f "$TMP_METHODS"
echo "Done: $OUT_METHODS"

#--- 2: Supplementary Table 3 (caption + table) ---
TMP_ST3=$(mktemp /tmp/exp_st3_XXXX.tex)
# Extract the ST3 caption block (lines 31–46), dropping the \input line
sed -n '31,46p' "$ST3_CAPTION_TEX" | grep -v '\\input{Tables/ST3' > "$TMP_ST3"
# Append the actual table content inline
cat "$ST3_TABLE_TEX" >> "$TMP_ST3"
pandoc "$TMP_ST3" \
  --from=latex \
  --to=docx \
  --citeproc \
  --bibliography="$BIB" \
  --csl="$CSL" \
  --output="$OUT_ST3"
rm -f "$TMP_ST3"
echo "Done: $OUT_ST3"

#--- 3: References ---
TMP_REFS=$(mktemp /tmp/exp_refs_XXXX.md)
cat > "$TMP_REFS" << 'EOF'
---
nocite: '@*'
---
EOF
pandoc "$TMP_REFS" \
  --from=markdown \
  --to=docx \
  --citeproc \
  --bibliography="$BIB" \
  --csl="$CSL" \
  --output="$OUT_REFS"
rm -f "$TMP_REFS"
echo "Done: $OUT_REFS"
