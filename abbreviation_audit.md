# Chemical Name Capitalization Audit

**Date:** 2026-07-14
**Scope:** Consistency of *displayed* chemical-name capitalization across the manuscript
text, manuscript tables, supplementary tables (ST1–ST3), figure captions, and
figure-label code.
**Status:** Audit complete. **No edits applied** — this is a proposal/reference document.

---

## TL;DR

Capitalization is **substantially more consistent than feared**. The manuscript already
follows the correct scientific convention (Title-Case in tables/figure labels, lowercase
in running prose). There is **one genuine inconsistency**: **`o-aminoazotoluene`** is the
lone lowercase-root name in otherwise Title-Case tables and disagrees with the supplement.

---

## Method & coverage

| Source | How inspected | Coverage |
|---|---|---|
| Manuscript text + tables | R1 docx → pandoc markdown (`--track-changes=accept`) | full text + tables |
| Supplementary tables ST1–ST3 | `.tex` sources | ST1 sampled (710 rows; sampled rows uniform), ST2/ST3 full |
| Figure captions | `Supplementary/Components/Sections/figures.Rmd` | full |
| Figure-label code | `12_IARC_vis.R`, `R/Utilities/Visualization/*` | full (code, not baked pixels) |

**Not fully verified (see Coverage gaps):** pixels baked into rendered figure grobs; every
one of ST1's 710 rows; a full reverse-direction prose scan (mid-sentence Title-Case errors).

---

## The correct convention (canonical rule)

- **Locant / stereo prefixes** (`o-`, `m-`, `p-`, `n-`, `N-`, `O-`, `sec-`, `tert-`,
  `cis-`, `trans-`, `α-`, `β-`, `γ-`): **lowercase** (italic when typeset).
- **Parent (root) name:**
  - **Title-Case** in tables, figure/axis labels, and at the **start of a sentence**.
  - **lowercase** in running prose **mid-sentence** (chemical names are common nouns).
- **Full names** used consistently; **abbreviations** (4-ABP, MOCA, DEET) only where
  space-constrained and defined.

> Example — all correct: table cell `o-Toluidine`; prose "…the carcinogen *o-toluidine*…";
> sentence start "*o*-Toluidine was detected…".

---

## ⚠️ Genuine inconsistency (proposed fix)

### `o-aminoazotoluene` → `o-Aminoazotoluene`

The only lowercase-root chemical name appearing inside Title-Case tables. It disagrees with
the supplement, which already has it correct.

| Location | Current | Proposed | Who fixes |
|---|---|---|---|
| Manuscript table (the `^†‡¶^` row) | `o-aminoazotoluene` ❌ | `o-Aminoazotoluene` | user (docx, find/replace) |
| Validation spectrum figure label (short-name in `figure_order`) | `o-aminoazotoluene` ❌ | `o-Aminoazotoluene` | needs `figure_order` edit **+ grob re-render** |
| Supplement ST tables | `o-Aminoazotoluene` ✅ | (no change) | — |

Every sibling (`o-Toluidine`, `o-Cresol`, etc.) is Title-Case in the same tables, so this is
an oversight, not a deliberate style choice.

---

## Looks inconsistent but is CORRECT — do **not** change

- **Table Title-Case vs. prose lowercase** for the same chemical
  (`o-Toluidine` / `o-toluidine`, `4-Aminobiphenyl` / `4-aminobiphenyl`,
  `Fluoranthene` / `fluoranthene`, `Methoxychlor` / `methoxychlor`,
  `2,4-Dimethylphenol` / `2,4-dimethylphenol`, `o-Cresol` / `o-cresol`). This is proper
  scientific style. Forcing Title-Case into prose would be an error.
- **Abbreviation vs. full name** (`4-ABP` vs `4-Aminobiphenyl`, `MOCA`, `2-ABP`, `DEET`):
  different purposes (space-constrained figure short-names vs. table names), not miscasing.
- **Different chemicals sharing a root / class plurals**: `Benzo(b)fluoranthene`,
  `1-Naphthylamine` vs `2-Naphthylamine`, `Anthracenes`, `Pyrenes`, `Carbofuran` vs
  `3-Hydroxycarbofuran` — distinct entities, all internally consistent.

---

## Cross-source extraction (reference data)

Per root fragment: distinct case-forms found in each source. `***` = a genuine same-name
casing difference (all others are convention-correct or different chemicals).

```
toluidine ***    MS: o-Toluidine(table) / o-toluidine(prose) | 5-Nitro-o-toluidine
                 ST: o-Toluidine | 3-Chloro-o-toluidine | 4-Chloro-o-toluidine | 5-Nitro-o-toluidine
                 -> convention-correct (table vs prose)

aminobiphenyl *** MS: 4-Aminobiphenyl(table) / 4-aminobiphenyl(prose)
                 ST: 4-Aminobiphenyl | 2-Aminobiphenyl
                 -> convention-correct (table vs prose)

aminoazotoluene *** MS: o-aminoazotoluene(table)   <-- FIX to o-Aminoazotoluene
                 ST: o-Aminoazotoluene (correct)

cresol ***       MS: o-Cresol(table) / o-cresol(prose)
                 ST: m-Cresol | o-Cresol | p-Cresol      -> convention-correct

pyrene ***       MS: Pyrene(table) / pyrene(only in "Benzo[a]pyrene", a reference) -> fine
fluoranthene *** MS: Fluoranthene(table) / fluoranthene(prose)                     -> fine
methoxychlor *** MS: Methoxychlor(table) / methoxychlor(prose) | 2,4'-Methoxychlor -> fine
dimethylphenol***MS: 2,4-Dimethylphenol(table) / 2,4-dimethylphenol(prose)         -> fine

naphthylamine    1-Naphthylamine, 2-Naphthylamine (distinct)                       -> fine
anthracene       Anthracene, Benz(a)anthracene, Dibenz(a,h)anthracene (distinct)   -> fine
phthalate        class label "Phthalate"/"phthalate" in fig usage-class            -> verify (class category, not a chemical name)
```

Greek-letter chemicals (`gamma-BHC`, etc.) appear only in validation spectra (supplement
figures), **not** in the manuscript tables/text — consistent within their scope, but the
symbol-vs-spelled-out form (`γ-` vs `gamma-`) was not cross-checked in the rendered grobs.

---

## Coverage gaps / follow-ups

1. **Rendered figure grobs** — the `o-Aminoazotoluene` spectrum label (and any other baked-in
   figure text) was not pixel-inspected. Verify on the actual rendered validation figure.
2. **ST1 full 710 rows** — sampled only; a full pass would confirm no stragglers.
3. **Reverse prose check** — a chemical wrongly Title-Cased mid-sentence can't be fully
   auto-detected without a chemical dictionary. Frequently-used names spot-checked clean.
4. **Greek forms** in rendered spectra (`γ-` vs `gamma-`) — cross-check for consistency.

---

## Proposed action items

- [ ] **Manuscript docx:** find/replace `o-aminoazotoluene` → `o-Aminoazotoluene` (table row).
- [ ] **`validation.xlsx` → `figure_order`:** change short-name `o-aminoazotoluene` →
      `o-Aminoazotoluene`, then re-render that spectrum so the figure label updates.
- [ ] **Supplement ST tables:** no change (already correct).
- [ ] (Optional) Exhaustive ST1 (710-row) pass + full reverse prose check.
- [ ] (Optional) Verify `γ-`/`gamma-` and other locants in rendered spectra grobs.
