# CoxeterDimBound

Exact-arithmetic programs, finite input catalogues, and reproducibility data
for local Coxeter-diagram weights and dimension-bound comparisons.
The repository is self-contained for these computations. It does not require
a manuscript, TeX installation, or network connection after dependencies are installed.

## Requirements

The pinned environment was checked with Python 3.9.12. Install its dependencies:

```sh
python -m pip install -r requirements.txt
```

The dependencies are NetworkX 2.7.1 and SymPy 1.14.0. The default stored-result
validation was also checked with Python 3.12.14. The certificate generator
and independent replayer use the Python standard library only. Computations
use exact rational or symbolic arithmetic; stored timing fields are diagnostic.

## Verification

From the repository root:

```sh
python verification/verify_all.py
```

This checks the integrity manifest, validates the supplied result summaries,
and recomputes the final rational comparisons. It does **not** rerun the large
enumerations or independently replay every tree certificate.

To add the short product-link, coordinate-example, and certificate-arithmetic audits:

```sh
python verification/verify_all.py --quick
```

To independently reconstruct and replay all supplied tree certificates:

```sh
python verification/verify_all.py --replay-tree
```

To regenerate the complete finite computations, including the large cusp
enumeration and the tree records:

```sh
python verification/verify_all.py --regenerate
```

Full regeneration is substantially more expensive than checking stored
results. Do not run multiple regeneration commands against the same checkout
at once: the programs write their result files in place. Use ordinary Python,
without `-O`, because the programs use assertions for verification.

Commands may be launched from another directory by giving the absolute path
to the driver. Paths to inputs and outputs are resolved from each script.
The driver writes a report to `verification/verification_summary.json` and
local execution logs under `verification/run_logs/`. It refreshes
`SHA256SUMS.json` after a successful run to account for regenerated metadata.
The manifest is a consistency check, not an independently authenticated proof.

## Contents

| Location | Purpose |
|---|---|
| `supplementary/edge/verification/` | Reduced cores, three-face markings, extension checks, exact metric and cusp calculations |
| `supplementary/edge/tables/` | Generated catalogues, marking tables, and exact local results |
| `supplementary/edge/logs/` | Structured results for the edge computations and short algebra checks |
| `supplementary/vertex/verification/` | Planar-angle budgets, triangle and quadrilateral checks, cycle checks, cusp enumeration, final comparisons |
| `supplementary/tree/` | Tree-support generator, 79,113 records, and an independent certificate replayer |
| `supplementary/shared_moment_audit.py` | Product-link polynomial and moment identities in a finite test range |
| `supplementary/vertex_tree_certificate_audit.py` | Separate exact checks of certificate inequalities and preserved short paths |
| `verification/verify_all.py` | Unified verification driver |

The historical filename `chapter4_readability_checks.py` denotes exact checks
of illustrative examples and local-to-global algebra; it does not read or
edit a document. `terminal_audit.py` reconstructs the eight final comparison
rows from the stored finite-vertex and cusp budgets and checks the stated
dimension-81 candidate separately.

Required catalogue inputs include `edge/verification/reduced_cores.json`,
the two copies of `high_cores.json`, and
`vertex/verification/expected_tail_positions.json`. Supplied result files
allow inexpensive validation without rerunning the complete census.

## Scope

The stored computations include 61 reduced core entries, 5,270 three-face
markings, 79,113 tree certificates, and 2,770,038 reduced cusp lists.
The edge and planar extension checks impose different conditions and have
different outcomes: respectively `92 = 76 + 3 + 13` and
`92 = 3 + 65 + 3 + 21`.

These programs verify finite arithmetic under their stated classification
and geometric hypotheses. They do not establish the classification's
completeness, geometric realizability of every candidate, or the
arbitrary-length and arbitrary-rank reductions by computation alone.
They are not a proof-assistant formalization. See `supplementary/tree/README.txt`
for the precise certificate inputs, checks, and coverage assumptions.

