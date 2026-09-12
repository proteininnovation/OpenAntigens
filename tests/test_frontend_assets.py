import tempfile
import unittest
from pathlib import Path

from agdesign2.portal import _version_portal_assets
from scripts.build_public_release import validate_frontend_assets


class FrontendAssetTests(unittest.TestCase):
    def test_versions_change_with_served_bytes_and_validate_references(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            page = root / 'index.html'
            asset = root / 'data.js'
            page.write_text('<script src="data.js?v=old"></script>')
            asset.write_text('window.rows = [1];')
            _version_portal_assets(root)
            first = page.read_text()
            _version_portal_assets(root)
            self.assertEqual(page.read_text(), first)
            asset.write_text('window.rows = [1, 2];')
            _version_portal_assets(root)
            self.assertNotEqual(page.read_text(), first)
            validate_frontend_assets(root)
            asset.unlink()
            with self.assertRaisesRegex(ValueError, 'Missing local frontend asset'):
                validate_frontend_assets(root)
