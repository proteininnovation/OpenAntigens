"""Build small, offline interaction fixtures; coordinates are synthetic."""
import json
import sys
from pathlib import Path

from Bio.SeqUtils import seq3
from agdesign2.portal import build_portal
from tests.test_portal_db import _write_fixture
from tests.test_construct_tabs import tab_constructs
from tests.test_report_navigation import with_optional_sections


def structure(sequence):
    return '\n'.join(
        f'ATOM  {i:5d}  CA  {seq3(aa).upper()} A{i:4d}    {float(i):8.3f}{13.207:8.3f}{2.100:8.3f}  1.00 91.00           C'
        for i, aa in enumerate(sequence, 1)
    ) + '\nEND\n'


def build(root):
    root.mkdir(parents=True, exist_ok=True)
    batch = root / 'batch'
    batch.mkdir(exist_ok=True)
    summary, payload = _write_fixture(batch, with_structure=True)
    sequence = 'MMCM' + 'M' * 6
    payload['target']['sequence'] = sequence
    payload['construct_details'][0]['sequence'] = sequence
    payload['construct_details'][0]['cysteine_analysis'] = [{'position': 3, 'paired_with': None}]
    payload['cysteine_analysis'] = [{'position': 3, 'paired_with': None}]
    payload['ptms'] = [{'type': 'MOD_RES', 'start': 3, 'end': 3, 'description': '<img id="injected-annotation" src="bad" onerror="window.annotationExecuted=true">'}]
    (batch / 'test_human_report.json').write_text(json.dumps(payload))
    (batch / 'alphafold/P00001.pdb').write_text(structure(sequence))
    (batch / 'alphafold/P00001.pae.json').write_text(json.dumps({'pae': [[1] * 10 for _ in range(10)]}))
    build_portal(summary, output_dir=root / 'site', refresh_stale_reports=False, fetch_literature=False, bundle_vendor_assets=True)

    tab_report = {**payload, 'construct_details': tab_constructs(payload['construct_details'][0])}
    (batch / 'test_human_report.json').write_text(json.dumps(tab_report))
    build_portal(summary, output_dir=root / 'tabs-site', refresh_stale_reports=False, fetch_literature=False, bundle_vendor_assets=True)
    (batch / 'test_human_report.json').write_text(json.dumps(with_optional_sections(tab_report)))
    build_portal(summary, output_dir=root / 'navigation-site', refresh_stale_reports=False, fetch_literature=False, bundle_vendor_assets=True)
    (batch / 'test_human_report.json').write_text(json.dumps(payload))

    fixture = json.loads((Path(__file__).parent / 'fixtures/ptprc_mapping.json').read_text())
    human = 'M' * (fixture['human_start'] - 1) + fixture['human_sequence']
    mouse = fixture['construct']
    full_homolog = {**mouse, 'start': fixture['mouse_start'], 'end': fixture['mouse_start'] + len(fixture['mouse_sequence']) - 1, 'sequence': fixture['mouse_sequence']}
    ptprc = {
        'target': {'accession': 'P08575', 'entry_name': 'PTPRC_HUMAN', 'gene_symbol': 'PTPRC', 'sequence': human},
        'ectodomain': {'start': fixture['human_start'], 'end': len(human), 'label': 'ectodomain'},
        'construct_details': [
            {'name': 'full_ectodomain', 'start': fixture['human_start'], 'end': len(human), 'sequence': fixture['human_sequence'], 'homologs': [full_homolog]},
            {'name': 'strict_domain_2', 'start': 390, 'end': 481, 'sequence': human[389:481], 'homologs': [mouse]},
        ],
    }
    (batch / 'ptprc_human_report.json').write_text(json.dumps(ptprc))
    (batch / 'alphafold/P08575.pdb').write_text(structure(human))
    summary.write_text(json.dumps([{'query': 'PTPRC_HUMAN', 'status': 'ok', 'json_report': str(batch / 'ptprc_human_report.json')}]))
    build_portal(summary, output_dir=root / 'ptprc-site', refresh_stale_reports=False, fetch_literature=False, bundle_vendor_assets=True)
    # Exercise the independent no-structure card renderer as well.
    (batch / 'alphafold/P00001.pdb').unlink()
    summary.write_text(json.dumps([{'query': 'TEST_HUMAN', 'status': 'ok', 'json_report': str(batch / 'test_human_report.json')}]))
    build_portal(summary, output_dir=root / 'no-structure-site', refresh_stale_reports=False, fetch_literature=False)


if __name__ == '__main__':
    build(Path(sys.argv[1]).resolve())
