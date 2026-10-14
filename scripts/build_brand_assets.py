"""Generate the OpenAntigens receptor-mark files in assets/branding/.

One geometry (a 64-unit grid) drives both outputs, so the SVG logos and the
raster favicons/app icons cannot drift apart:

- openantigens-mark.svg / -reversed.svg / -mouse.svg: full mark for light,
  navy and mouse-portal light grounds; knockout gaps are masked (transparent).
- favicon.svg, favicon.ico (16/32/48), favicon-32.png: the small-size mark.
- apple-touch-icon.png (full-bleed navy), icon-192.png, icon-512.png (rounded tile).

Run: python scripts/build_brand_assets.py
"""

from __future__ import annotations

from pathlib import Path

from PIL import Image, ImageDraw

OUT = Path(__file__).resolve().parents[1] / "assets" / "branding"
NAVY = "#0d1f3c"
GAP = 1.6  # the one knockout gap between touching shapes
DOMAINS = [(31, 40.3), (38, 26.9), (31, 13.5)]  # membrane-proximal -> distal
EPITOPE = (37.4, 7.9, 3.6)
HEAD_ROWS = (52.4, 58.6)
HEAD_X = (12.4, 18.6, 24.8, 37.2, 43.4, 49.6, 55.8)

PALETTES = {
    "light": dict(prox="#0d1f3c", mid="#0e8c88", dist="#4ec3dc", epi="#e8684a", lip="#a9bccd", stalk="#0d1f3c"),
    "reversed": dict(prox="#ffffff", mid="#2bb0a8", dist="#92efe5", epi="#e8684a", lip="#3d5878", stalk="#ffffff"),
    "mouse": dict(prox="#33262b", mid="#bd4d61", dist="#e8909c", epi="#0e8c88", lip="#d9bcc3", stalk="#33262b"),
}


def shapes(size: str, tile: bool) -> dict:
    """Geometry for 'full' (>=48px) or 'tiny' (<=32px) renditions."""
    tiny = size == "tiny"
    r = 9.0 if tiny else 8.2
    circles = [(x, y, r) for x, y in DOMAINS] + ([] if tiny else [EPITOPE])
    if tiny:
        membrane = [("rect", 4, 51, 56, 7)]
    else:
        xs = ((6.2,) if tile else ()) + HEAD_X
        membrane = [("head", x, y, 2.3) for y in HEAD_ROWS for x in xs]
    return dict(
        circles=circles,
        gap=1.5 if tiny else GAP,
        membrane=membrane,
        stalk=(46, 64 if tile else (62 if tiny else 61), 6 if tiny else 4.2),
    )


def svg(palette: str, size: str = "full") -> str:
    p, g = PALETTES[palette], shapes(size, tile=False)
    cols = [p["prox"], p["mid"], p["dist"], p["epi"]]
    defs, body = [], []

    def mask(mid: str, cutters: list[tuple[float, float, float]]) -> str:
        holes = "".join(f'<circle cx="{x}" cy="{y}" r="{r}" fill="#000"/>' for x, y, r in cutters)
        return (f'<mask id="{mid}" maskUnits="userSpaceOnUse" x="-4" y="-4" width="72" height="72">'
                f'<rect x="-4" y="-4" width="72" height="72" fill="#fff"/>{holes}</mask>')

    for kind, *a in g["membrane"]:
        if kind == "head":
            body.append(f'<circle cx="{a[0]}" cy="{a[1]}" r="{a[2]}" fill="{p["lip"]}"/>')
        else:
            body.append(f'<rect x="{a[0]}" y="{a[1]}" width="{a[2]}" height="{a[3]}" rx="{a[3] / 2}" fill="{p["lip"]}"/>')
    x0, y0, r0 = g["circles"][0]
    y_top, y_bottom, width = g["stalk"]
    defs.append(mask("k0", [(x0, y0, r0 + g["gap"])]))
    body.append(f'<line x1="31" y1="{y_top}" x2="31" y2="{y_bottom}" stroke="{p["stalk"]}" '
                f'stroke-width="{width}" stroke-linecap="round" mask="url(#k0)"/>')
    for i, (x, y, r) in enumerate(g["circles"]):
        later = [(a, b, c + g["gap"]) for a, b, c in g["circles"][i + 1:]]
        attr = ""
        if later:
            defs.append(mask(f"k{i + 1}", later))
            attr = f' mask="url(#k{i + 1})"'
        body.append(f'<circle cx="{x}" cy="{y}" r="{r}" fill="{cols[i]}"{attr}/>')
    return ('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 64 64" width="64" height="64" '
            'role="img" aria-label="OpenAntigens"><title>OpenAntigens</title>'
            f'<defs>{"".join(defs)}</defs>{"".join(body)}</svg>\n')


def raster(px: int, palette: str, size: str, tile: str | None = None) -> Image.Image:
    """Draw the mark at px x px; tile is None, 'rounded' or 'bleed'. Supersampled 8x."""
    ss = 8
    scale = px * ss / 64
    p, g = PALETTES[palette], shapes(size, tile=tile is not None)
    ground = (0, 0, 0, 0)
    img = Image.new("RGBA", (px * ss, px * ss), ground)
    d = ImageDraw.Draw(img)
    if tile:
        ground = NAVY
        if tile == "bleed":
            d.rectangle([0, 0, px * ss, px * ss], fill=NAVY)
        else:
            d.rounded_rectangle([0, 0, px * ss - 1, px * ss - 1], radius=14 * scale, fill=NAVY)

    def disc(x: float, y: float, r: float, fill) -> None:
        d.ellipse([(x - r) * scale, (y - r) * scale, (x + r) * scale, (y + r) * scale], fill=fill)

    for kind, *a in g["membrane"]:
        if kind == "head":
            disc(a[0], a[1], a[2], p["lip"])
        else:
            d.rounded_rectangle([a[0] * scale, a[1] * scale, (a[0] + a[2]) * scale, (a[1] + a[3]) * scale],
                                radius=a[3] / 2 * scale, fill=p["lip"])
    y_top, y_bottom, width = g["stalk"]
    d.line([31 * scale, y_top * scale, 31 * scale, y_bottom * scale], fill=p["stalk"], width=round(width * scale))
    disc(31, y_top, width / 2, p["stalk"])
    disc(31, y_bottom, width / 2, p["stalk"])
    cols = [p["prox"], p["mid"], p["dist"], p["epi"]]
    for i, (x, y, r) in enumerate(g["circles"]):
        disc(x, y, r + g["gap"], ground)  # cut the gap out of everything already drawn
        disc(x, y, r, cols[i])
    if tile == "rounded":  # restore transparent corners the gap cuts may have painted navy
        corner = Image.new("L", img.size, 0)
        ImageDraw.Draw(corner).rounded_rectangle([0, 0, px * ss - 1, px * ss - 1], radius=14 * scale, fill=255)
        img.putalpha(Image.composite(img.getchannel("A"), Image.new("L", img.size, 0), corner))
    return img.resize((px, px), Image.LANCZOS)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for name, palette in (("openantigens-mark", "light"), ("openantigens-mark-reversed", "reversed"),
                          ("openantigens-mark-mouse", "mouse")):
        (OUT / f"{name}.svg").write_text(svg(palette), encoding="utf-8")
    (OUT / "favicon.svg").write_text(svg("light", "tiny"), encoding="utf-8")
    raster(32, "light", "tiny").save(OUT / "favicon-32.png")
    raster(48, "light", "tiny").save(OUT / "favicon.ico", sizes=[(16, 16), (32, 32), (48, 48)])
    raster(180, "reversed", "full", tile="bleed").convert("RGB").save(OUT / "apple-touch-icon.png")
    raster(192, "reversed", "full", tile="rounded").save(OUT / "icon-192.png")
    raster(512, "reversed", "full", tile="rounded").save(OUT / "icon-512.png")
    for path in sorted(OUT.glob("*")):
        if path.name.startswith(("openantigens-mark", "favicon", "apple-touch", "icon-")):
            print(f"{path.name:34s} {path.stat().st_size:>7,d} bytes")


if __name__ == "__main__":
    main()
