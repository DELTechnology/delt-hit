from __future__ import annotations

import argparse
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont, ImageOps


BG = "white"
INK = "#111111"
MUTED = "#4f5b66"
BORDER = "#d9dee3"
SECTION_BG = "#fafbfd"


def parse_args() -> argparse.Namespace:
    script_dir = Path(__file__).resolve().parent
    default_visualization_root = (
        script_dir.parent
        / "supporting_material"
        / "experiments"
        / "example-single-stranded"
        / "campaign"
        / "library"
        / "visualization"
    )
    default_output_dir = script_dir / "figures"
    parser = argparse.ArgumentParser(
        description=(
            "Generate the composite manuscript figures for the single-stranded "
            "example workflow."
        )
    )
    parser.add_argument(
        "--visualization-root",
        type=Path,
        default=default_visualization_root,
        help="Root directory containing the visualization PNG outputs.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=default_output_dir,
        help="Directory where the composite figures should be written.",
    )
    return parser.parse_args()


def _get_font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont | ImageFont.ImageFont:
    candidates: list[str] = []
    if bold:
        candidates.extend(
            [
                "/System/Library/Fonts/Supplemental/Arial Bold.ttf",
                "/System/Library/Fonts/Supplemental/Helvetica Bold.ttf",
            ]
        )
    candidates.extend(
        [
            "/System/Library/Fonts/Supplemental/Arial Unicode.ttf",
            "/System/Library/Fonts/Supplemental/Arial.ttf",
            "/System/Library/Fonts/Supplemental/Helvetica.ttc",
            "/System/Library/Fonts/SFNS.ttf",
        ]
    )
    for candidate in candidates:
        path = Path(candidate)
        if not path.exists():
            continue
        try:
            return ImageFont.truetype(str(path), size=size)
        except OSError:
            continue
    return ImageFont.load_default()


TITLE = _get_font(64, bold=True)
SUB = _get_font(40, bold=True)
LABEL = _get_font(30, bold=True)
SMALL = _get_font(24)


def contain(image: Image.Image, size: tuple[int, int]) -> Image.Image:
    return ImageOps.contain(image.convert("RGB"), size, method=Image.Resampling.LANCZOS)


def draw_header(
    draw: ImageDraw.ImageDraw,
    x: int,
    y: int,
    title: str,
    subtitle: str,
) -> None:
    draw.text((x, y), title, font=TITLE, fill=INK)
    draw.text((x, y + 78), subtitle, font=SMALL, fill=MUTED)


def panel(
    image_path: Path,
    size: tuple[int, int],
    *,
    title: str | None = None,
    subtitle: str | None = None,
) -> Image.Image:
    image = Image.open(image_path)
    canvas = Image.new("RGB", size, BG)
    draw = ImageDraw.Draw(canvas)
    draw.rounded_rectangle(
        (2, 2, size[0] - 3, size[1] - 3),
        22,
        outline=BORDER,
        width=3,
        fill=BG,
    )
    inner_w = size[0] - 40
    inner_h = size[1] - 90
    top_offset = 0
    if title or subtitle:
        inner_h = size[1] - 120
        top_offset = 30
    body = contain(image, (inner_w, inner_h))
    x = (size[0] - body.width) // 2
    y = top_offset + (inner_h - body.height) // 2 + (30 if title or subtitle else 20)
    canvas.paste(body, (x, y))
    if title:
        draw.text((22, 18), title, font=LABEL, fill=INK)
    if subtitle:
        draw.text((22, 56), subtitle, font=SMALL, fill=MUTED)
    return canvas


def build_enumeration_composite(visualization_root: Path, output_path: Path) -> None:
    reaction_graph = visualization_root / "reaction_graph.png"
    reaction_abf1 = visualization_root / "reactions" / "ABF1.png"
    reaction_sr = visualization_root / "reactions" / "SR.png"

    b0_examples = [visualization_root / "building_blocks" / "B0" / f"{n}.png" for n in [0, 1, 10]]
    b1_examples = [visualization_root / "building_blocks" / "B1" / f"{n}.png" for n in [0, 1, 10]]
    compound_examples = [
        visualization_root / "libraries" / "AG24_4_top_hits" / name
        for name in [
            "B0=11-B1=274.png",
            "B0=23-B1=314.png",
            "B0=219-B1=314.png",
        ]
    ]

    width, height = 3600, 3000
    figure = Image.new("RGB", (width, height), BG)
    draw = ImageDraw.Draw(figure)
    draw_header(
        draw,
        120,
        70,
        "Enumeration workflow outputs",
        "Reaction graph, representative reactions, building blocks, and enumerated top-hit structures",
    )

    graph_panel = panel(
        reaction_graph,
        (2050, 1320),
        title="Library reaction graph",
        subtitle="Visualized from the example-single-stranded workflow",
    )
    figure.paste(graph_panel, (100, 220))

    figure.paste(panel(reaction_abf1, (1250, 390), title="Reaction template ABF1"), (2230, 220))
    figure.paste(panel(reaction_sr, (1250, 390), title="Reaction template SR"), (2230, 640))

    compound_y = 1080
    compound_w = 390
    for i, compound_path in enumerate(compound_examples, start=1):
        compound_panel = panel(
            compound_path,
            (compound_w, 520),
            title=f"Top compound {i}",
            subtitle=compound_path.stem,
        )
        figure.paste(compound_panel, (2230 + (i - 1) * (compound_w + 20), compound_y))

    section_y = 1680
    for x, title in [
        (100, "Representative B0 building blocks"),
        (1840, "Representative B1 building blocks"),
    ]:
        draw.rounded_rectangle(
            (x, section_y, x + 1660, section_y + 1040),
            28,
            fill=SECTION_BG,
            outline=BORDER,
            width=3,
        )
        draw.text((x + 30, section_y + 24), title, font=SUB, fill=INK)

    for i, building_block_path in enumerate(b0_examples):
        block_panel = panel(building_block_path, (500, 860), title=f"B0 code {building_block_path.stem}")
        figure.paste(block_panel, (140 + i * 520, 1820))

    for i, building_block_path in enumerate(b1_examples):
        block_panel = panel(building_block_path, (500, 860), title=f"B1 code {building_block_path.stem}")
        figure.paste(block_panel, (1880 + i * 520, 1820))

    figure.save(output_path, quality=95)


def build_library_examples_composite(visualization_root: Path, output_path: Path) -> None:
    library_examples = [
        visualization_root / "libraries" / "AG24_4_top_hits" / name
        for name in [
            "B0=11-B1=274.png",
            "B0=23-B1=314.png",
            "B0=219-B1=314.png",
            "B0=290-B1=314.png",
            "B0=295-B1=556.png",
            "B0=222-B1=314.png",
        ]
    ]

    width, height = 3300, 2600
    figure = Image.new("RGB", (width, height), BG)
    draw = ImageDraw.Draw(figure)
    draw_header(
        draw,
        110,
        70,
        "Library visualization examples",
        "Six representative enumerated top-hit structures from AG24_4_top_hits",
    )

    positions = [
        (110, 260),
        (1160, 260),
        (2210, 260),
        (110, 1360),
        (1160, 1360),
        (2210, 1360),
    ]
    size = (980, 980)
    for i, (image_path, (x, y)) in enumerate(zip(library_examples, positions), start=1):
        example_panel = panel(image_path, size, title=f"Example {i}", subtitle=image_path.stem)
        figure.paste(example_panel, (x, y))

    figure.save(output_path, quality=95)


def main() -> None:
    args = parse_args()
    visualization_root = args.visualization_root
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    build_enumeration_composite(visualization_root, output_dir / "image2.png")
    build_library_examples_composite(visualization_root, output_dir / "image3.png")

    print(f"Wrote {output_dir / 'image2.png'}")
    print(f"Wrote {output_dir / 'image3.png'}")


if __name__ == "__main__":
    main()
