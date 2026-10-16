import importlib.util
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path

from agdesign2.portal import _portal_brand_header, _portal_footer

ROOT = Path(__file__).resolve().parents[1]
BRANDING = ROOT / "assets" / "branding"


class BrandAssetTests(unittest.TestCase):
    def test_header_pairs_receptor_mark_with_ipi_logo_in_both_portals(self):
        for title, prefix in (("OpenAntigens", ""), ("OpenAntigens Mouse", "../")):
            with self.subTest(title=title):
                header = _portal_brand_header(prefix=prefix, portal_title=title)
                self.assertIn(f'src="{prefix}openantigens-mark-reversed.svg"', header)
                self.assertIn(f'src="{prefix}ipi-logo-light-800.png"', header)
                self.assertIn('alt="Institute for Protein Innovation"', header)
                self.assertIn('href="https://proteininnovation.org/"', header)
        self.assertNotIn("brand-species", _portal_brand_header(portal_title="OpenAntigens"))
        self.assertIn('<span class="brand-species">Mouse</span>', _portal_brand_header(portal_title="OpenAntigens Mouse"))

    def test_footer_uses_species_mark_and_ipi_logo(self):
        human = _portal_footer(portal_title="OpenAntigens")
        mouse = _portal_footer(prefix="../", portal_title="OpenAntigens Mouse")
        self.assertIn('src="openantigens-mark.svg"', human)
        self.assertIn('src="../openantigens-mark-mouse.svg"', mouse)
        for footer in (human, mouse):
            self.assertIn("ipi-logo-dark-800.png", footer)

    def test_brand_files_exist_and_svgs_define_their_masks(self):
        for name in ("openantigens-mark.svg", "openantigens-mark-reversed.svg", "openantigens-mark-mouse.svg",
                     "favicon.svg", "manrope-latin.woff2", "MANROPE-OFL.txt"):
            self.assertTrue((BRANDING / name).is_file(), name)
        self.assertEqual((BRANDING / "manrope-latin.woff2").read_bytes()[:4], b"wOF2")
        spec = importlib.util.spec_from_file_location("build_brand_assets", ROOT / "scripts" / "build_brand_assets.py")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        for palette in ("light", "reversed", "mouse"):
            for size in ("full", "tiny"):
                text = module.svg(palette, size)
                root = ET.fromstring(text)
                ns = "{http://www.w3.org/2000/svg}"
                masks = {m.get("id") for m in root.iter(f"{ns}mask")}
                used = {el.get("mask")[5:-1] for el in root.iter() if el.get("mask")}
                self.assertTrue(used and used <= masks, (palette, size))
        self.assertEqual(module.svg("light"), (BRANDING / "openantigens-mark.svg").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
