"""Run the available exact checks for the two-face dimension manuscript.

The tree-support certificates are supplied separately in ../../tree.
This driver verifies the other exact components; the top-level driver
also regenerates and independently replays the tree certificates.
"""
from pathlib import Path
import argparse
import hashlib
import json
import subprocess
import sys

HERE = Path(__file__).absolute().parent
CHECKS = (
    ('budget_audit.py', 'budget_results.json'),
    ('tail_gluing_audit.py', 'tail_gluing_results.json'),
    ('triangle_audit.py', 'triangle_results.json'),
    ('cycle_audit.py', 'cycle_results.json'),
    ('envelope_audit.py', 'envelope_results.json'),
    ('cusp_enumeration.py', 'cusp_results.json'),
)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cached', action='store_true',
                        help='Validate already generated exact summaries without rerunning enumeration.')
    args = parser.parse_args()
    entries = []
    for script, result in CHECKS:
        if not args.cached:
            subprocess.run([sys.executable, str(HERE / script)], cwd=HERE, check=True)
        path = HERE / result
        data = json.loads(path.read_text())
        if script == 'budget_audit.py':
            assert data['stable_formula_checks'] == 156
            assert len(data['finite_DP']) == 8
        elif script == 'triangle_audit.py':
            assert len(data['records']) == 1265
        elif script == 'tail_gluing_audit.py':
            assert len(data['tails']) == 21 and data['marked_records'] == 16
            assert data['candidates'] == 92
            assert data['counts'] == {
                'S_not_PD': 65, 'angle_not_PSD': 3,
                'conflicting_labels': 3, 'survivor': 21,
            }
            assert data['minimum_score'] == '1'
        elif script == 'cycle_audit.py':
            assert len(data['records']) == 192
            assert sum(row['negative_witnesses'] for row in data['summaries']) == 185
        elif script == 'envelope_audit.py':
            assert data['records'] == 341
            assert data['max'] == '254166/24565'
        elif script == 'cusp_enumeration.py':
            assert data['enumerated_reduced_partitions'] == 2770038
            assert len(data['cusp_table']) == 8
        entries.append({
            'script': script, 'summary': result,
            'summary_sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
        })
    report = {
        'mode': 'validate_stored_exact_summaries' if args.cached else 'regenerate_exact_checks',
        'completed_checks': entries,
        'tree_support_certificates_79113': 'supplied_in_supplementary_tree',
        'geometric_arbitrary_length_and_concavity_coverage': 'mathematical_proof_in_manuscript',
    }
    (HERE / 'verification_summary.json').write_text(json.dumps(report, indent=2) + '\n')
    print('Available exact checks passed:', len(entries), flush=True)
    print('Tree certificates are supplied in supplementary/tree and checked by the unified driver.', flush=True)

if __name__ == '__main__':
    main()
