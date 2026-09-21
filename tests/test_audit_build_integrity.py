import json
import sqlite3
import tempfile
import unittest
from contextlib import ExitStack
from pathlib import Path
from unittest.mock import patch

from fastapi.testclient import TestClient
from agdesign2 import deployment, portal
from agdesign2.dynamic_portal import create_app
from agdesign2.module_refresh import _module_refresh_config
from agdesign2.portal_db import build_portal_database, build_portal_from_database, create_db_portal_app
from tests.test_portal_db import _write_fixture


class AuditBuildIntegrityTests(unittest.TestCase):
    def test_download_links_and_manifest_follow_available_files(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            summary, _ = _write_fixture(root)
            for available in (False, True, False):
                source = root / 'open_targets_disease_associations.tsv'
                if available:
                    source.write_text('targetId\n')
                else:
                    source.unlink(missing_ok=True)
                site = root / 'site'
                portal.build_portal(summary, output_dir=site, refresh_stale_reports=False, fetch_literature=False)
                db = build_portal_database(summary, root / 'portal.sqlite')
                db_site = root / 'db-site'
                build_portal_from_database(db, db_site)
                documents = [(p / 'downloads.html').read_text() for p in (site, db_site)]
                manifests = [json.loads((p / 'downloads/download_manifest.json').read_text()) for p in (site, db_site)]
                for client in (TestClient(create_app(summary)), TestClient(create_db_portal_app(db))):
                    documents.append(client.get('/downloads.html').text)
                    manifests.append(client.get('/downloads/download_manifest.json').json())
                    self.assertEqual(client.get('/open_targets_disease_associations.tsv').status_code, 200 if available else 404)
                for p in (site, db_site):
                    self.assertEqual((p / source.name).is_file(), available)
                for document, manifest in zip(documents, manifests, strict=True):
                    with self.subTest(available=available):
                        self.assertEqual('href="open_targets_disease_associations.tsv"' in document, available)
                        self.assertNotIn('href="open_targets_disease_associations.json"', document)
                        paths = {item['path'] for item in manifest['files']}
                        self.assertEqual('open_targets_disease_associations.tsv' in paths, available)
                        self.assertNotIn('open_targets_disease_associations.json', paths)

    def test_live_summary_refresh(self):
        with tempfile.TemporaryDirectory() as tmp:
            summary, _ = _write_fixture(Path(tmp))
            client = TestClient(create_app(summary))
            self.assertEqual(len(client.get('/api/genes').json()), 1)
            rows = json.loads(summary.read_text())
            rows.append({'query': 'SECOND_HUMAN', 'status': 'pending'})
            summary.write_text(json.dumps(rows))
            self.assertEqual(len(client.get('/api/genes').json()), 2)

    def test_db_assets_and_indexed_case_insensitive_lookup(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            summary, _ = _write_fixture(root)
            db = build_portal_database(summary, root / 'portal.sqlite')
            site = root / 'site'
            build_portal_from_database(db, site)
            self.assertTrue((site / 'portal-disease-index.js').is_file())
            self.assertEqual(TestClient(create_db_portal_app(db)).get('/portal-disease-index.js').status_code, 200)
            with sqlite3.connect(db) as conn:
                plan = conn.execute('EXPLAIN QUERY PLAN SELECT report_json FROM reports WHERE UPPER(entry_name)=UPPER(?) OR LOWER(detail_page)=LOWER(?)', ('test_human', 'test_human')).fetchall()
            self.assertFalse(any('SCAN reports' in row[3] for row in plan), plan)

    def test_removed_targets_are_pruned_from_rebuilt_portal(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            summary, _ = _write_fixture(root, with_structure=True)
            site = root / 'site'
            portal.build_portal(summary, output_dir=site, refresh_stale_reports=False, fetch_literature=False)
            self.assertTrue((site / 'reports/test_human.html').exists())
            summary.write_text('[]')
            portal.build_portal(summary, output_dir=site, refresh_stale_reports=False, fetch_literature=False)
            self.assertFalse((site / 'reports/test_human.html').exists())
            self.assertFalse((site / 'report_scripts/test_human.js').exists())
            self.assertEqual(list((site / 'structures').glob('*.pdb')), [])

    def test_incomplete_snapshot_resume_enters_analysis_then_blocks_packaging(self):
        with tempfile.TemporaryDirectory() as tmp, ExitStack() as stack:
            root = Path(tmp)
            snapshot = root / 'partial'
            batch = snapshot / 'outputs/surfy_batch'
            batch.mkdir(parents=True)
            (snapshot / 'inputs').mkdir()
            (snapshot / 'inputs/accessible_targets.tsv').write_text('uniprot_name\nAAA_HUMAN\nBBB_HUMAN\n')
            (batch / 'batch_summary.json').write_text('[{"query":"AAA_HUMAN","status":"ok","source_column":"uniprot_name"}]')
            mocks = {}
            for name in ['AntigenAnalyzer', 'prefetch_targets_from_tsv', 'build_ortholog_table_from_tsv', 'build_family_alignments_from_tsv', 'build_paralog_reference_from_tsv', 'run_batch_from_tsv', '_package_public_site', 'validate_release_data']:
                mocks[name] = stack.enter_context(patch.object(deployment, name))
            stack.enter_context(patch.object(deployment, '_resolve_default_af3_catalog', return_value=None))
            with self.assertRaisesRegex(ValueError, 'every requested target'):
                deployment.build_fresh_snapshot(snapshot_root=root, snapshot_name='partial', resume=True, include_mouse=False)
            mocks['run_batch_from_tsv'].assert_called_once()
            mocks['_package_public_site'].assert_not_called()

    def test_reanalysis_rebuilds_existing_mouse_reports_and_portal(self):
        from types import SimpleNamespace

        for reanalyze in (False, True):
            with self.subTest(reanalyze=reanalyze), tempfile.TemporaryDirectory() as tmp, ExitStack() as stack:
                root = Path(tmp)
                snapshot = root / 'existing'
                (snapshot / 'inputs').mkdir(parents=True)
                (snapshot / 'inputs/accessible_targets.tsv').write_text('uniprot_name\nAAA_HUMAN\n')
                source = snapshot / 'mouse_data/portal'
                source.mkdir(parents=True)
                (source / 'index.html').write_text('rebuilt mouse')
                published = snapshot / 'public_site/mouse'
                published.mkdir(parents=True)
                (published / 'index.html').write_text('old mouse')
                mocks = {}
                for name in ['AntigenAnalyzer', 'prefetch_targets_from_tsv', 'build_ortholog_table_from_tsv', 'build_family_alignments_from_tsv', 'build_paralog_reference_from_tsv', 'run_batch_from_tsv', '_validate_batch_completion', 'refresh_report_modules', 'generate_assets_for_summary', 'build_open_targets_associations_from_bulk_downloads', 'build_portal', '_package_public_site', '_write_snapshot_manifest', 'validate_release_data']:
                    mocks[name] = stack.enter_context(patch.object(deployment, name))
                mouse = stack.enter_context(patch.object(deployment, 'build_mouse_portal_from_human_orthologs', return_value=SimpleNamespace(portal_index_path=source / 'index.html')))
                deployment.build_fresh_snapshot(snapshot_root=root, snapshot_name='existing', resume=True, reanalyze_reports=reanalyze)
                self.assertEqual(mocks['run_batch_from_tsv'].call_args.kwargs['resume'], not reanalyze)
                if reanalyze:
                    mouse.assert_called_once()
                    self.assertFalse(mouse.call_args.kwargs['resume'])
                    self.assertEqual((published / 'index.html').read_text(), 'rebuilt mouse')
                else:
                    mouse.assert_not_called()
                    self.assertEqual((published / 'index.html').read_text(), 'old mouse')

    def test_refresh_paths_stay_in_snapshot_before_directories_exist(self):
        with tempfile.TemporaryDirectory() as tmp:
            snapshot = Path(tmp) / 'snapshot'
            config = _module_refresh_config(batch_dir=snapshot / 'outputs/surfy_batch', verbose=False, enable_complex_portal=False)
            for path in [config.cache_dir, config.data_dir, config.blast_db_dir, config.ortholog_fasta_dir, config.ortholog_blast_db_dir]:
                self.assertTrue(path.is_relative_to(snapshot), path)
