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
    default_properties_root = (
        script_dir.parent
        / "supporting_material"
        / "experiments"
        / "example-single-stranded"
        / "campaign"
        / "library"
        / "properties"
        / "AG24_4_top_hits"
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
    parser.add_argument(
        "--properties-root",
        type=Path,
        default=default_properties_root,
        help="Directory containing molecular property distribution PNG outputs.",
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


SUB = _get_font(66, bold=True)
SECTION = _get_font(64, bold=True)
LABEL = _get_font(50, bold=True)
LABEL_PLAIN = _get_font(50, bold=False)
SMALL = _get_font(44)


def contain(image: Image.Image, size: tuple[int, int]) -> Image.Image:
    return ImageOps.contain(image.convert("RGB"), size, method=Image.Resampling.LANCZOS)


def require_paths(paths: list[Path]) -> list[Path]:
    missing = [path for path in paths if not path.exists()]
    if missing:
        joined = ", ".join(str(path) for path in missing)
        raise FileNotFoundError(f"Missing expected PNG files: {joined}")
    return paths


def first_building_blocks(directory: Path, count: int) -> list[Path]:
    indexed_paths = sorted(directory.glob("*.png"), key=lambda path: int(path.stem))
    if len(indexed_paths) < count:
        raise FileNotFoundError(f"Expected at least {count} PNG files in {directory}, found {len(indexed_paths)}")
    return indexed_paths[:count]


def top_hit_paths(visualization_root: Path) -> list[Path]:
    library_dir = visualization_root / "libraries" / "AG24_4_top_hits"
    hits = [
        library_dir / "B0=11-B1=274.png",
        library_dir / "B0=23-B1=314.png",
        library_dir / "B0=219-B1=314.png",
        library_dir / "B0=290-B1=314.png",
        library_dir / "B0=295-B1=556.png",
        library_dir / "B0=222-B1=314.png",
    ]
    return require_paths(hits)


def panel(
    image_path: Path,
    size: tuple[int, int],
    *,
    title: str | None = None,
    subtitle: str | None = None,
    title_font: ImageFont.FreeTypeFont | ImageFont.ImageFont = LABEL,
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
        inner_h = size[1] - 210
        top_offset = 58
    body = contain(image, (inner_w, inner_h))
    x = (size[0] - body.width) // 2
    y = top_offset + (inner_h - body.height) // 2 + (68 if title or subtitle else 20)
    canvas.paste(body, (x, y))
    if title:
        draw.text((28, 24), title, font=title_font, fill=INK)
    if subtitle:
        draw.text((28, 96), subtitle, font=SMALL, fill=MUTED)
    return canvas


def build_enumeration_composite(visualization_root: Path, output_path: Path) -> None:
    reaction_graph = visualization_root / "reaction_graph.png"
    reaction_abf1 = visualization_root / "reactions" / "ABF1.png"
    reaction_sr = visualization_root / "reactions" / "SR.png"
    reaction_aba2 = visualization_root / "reactions" / "ABF2.png"

    b0_examples = first_building_blocks(visualization_root / "building_blocks" / "B0", 3)
    b1_examples = first_building_blocks(visualization_root / "building_blocks" / "B1", 3)
    require_paths([reaction_graph, reaction_abf1, reaction_sr, reaction_aba2])

    width, height = 3600, 3000
    figure = Image.new("RGB", (width, height), BG)
    draw = ImageDraw.Draw(figure)
    graph_panel = panel(
        reaction_graph,
        (2050, 1420),
        title="Library reaction graph",
    )
    figure.paste(graph_panel, (100, 100))

    figure.paste(panel(reaction_abf1, (1250, 400), title="Reaction template ABF1"), (2230, 100))
    figure.paste(panel(reaction_sr, (1250, 400), title="Reaction template SR"), (2230, 520))
    figure.paste(panel(reaction_aba2, (1250, 400), title="Reaction template ABF2"), (2230, 940))

    section_y = 1560
    for x, title in [
        (100, "Representative B0 building blocks"),
        (1840, "Representative B1 building blocks"),
    ]:
        draw.rounded_rectangle(
            (x, section_y, x + 1660, section_y + 720),
            28,
            fill=SECTION_BG,
            outline=BORDER,
            width=3,
        )
        draw.text((x + 30, section_y + 30), title, font=SECTION, fill=INK)

    for i, building_block_path in enumerate(b0_examples):
        block_panel = panel(
            building_block_path,
            (500, 540),
            title=f"Code {building_block_path.stem}",
            title_font=LABEL_PLAIN,
        )
        figure.paste(block_panel, (140 + i * 520, 1680))

    for i, building_block_path in enumerate(b1_examples):
        block_panel = panel(
            building_block_path,
            (500, 540),
            title=f"Code {building_block_path.stem}",
            title_font=LABEL_PLAIN,
        )
        figure.paste(block_panel, (1880 + i * 520, 1680))

    figure.save(output_path, quality=95)


def build_library_examples_composite(visualization_root: Path, output_path: Path) -> None:
    library_examples = top_hit_paths(visualization_root)

    width, height = 3300, 2600
    figure = Image.new("RGB", (width, height), BG)

    positions = [
        (110, 110),
        (1160, 110),
        (2210, 110),
        (110, 1260),
        (1160, 1260),
        (2210, 1260),
    ]
    size = (980, 980)
    for i, (image_path, (x, y)) in enumerate(zip(library_examples, positions), start=1):
        example_panel = panel(image_path, size, title=f"Example {i}", subtitle=image_path.stem)
        figure.paste(example_panel, (x, y))

    figure.save(output_path, quality=95)


def build_properties_and_examples_composite(
    visualization_root: Path,
    properties_root: Path,
    output_path: Path,
) -> None:
    property_plot = properties_root / "prop_mw.png"
    require_paths([property_plot])
    library_examples = top_hit_paths(visualization_root)[:5]

    width, height = 3600, 3000
    figure = Image.new("RGB", (width, height), BG)
    draw = ImageDraw.Draw(figure)

    property_panel = panel(
        property_plot,
        (3380, 1280),
        title="Top-hit molecular weight distribution",
        subtitle="supporting_material/experiments/example-single-stranded/campaign/library/properties/AG24_4_top_hits/prop_mw.png",
    )
    figure.paste(property_panel, (110, 100))

    draw.rounded_rectangle(
        (110, 1460, 3490, 2820),
        28,
        fill=SECTION_BG,
        outline=BORDER,
        width=3,
    )
    draw.text((140, 1500), "Representative top-hit structures", font=SECTION, fill=INK)
    draw.text((140, 1585), "Examples from AG24_4_top_hits visualized with delt-hit visualize library", font=SMALL, fill=MUTED)

    positions = [
        (140, 1695),
        (810, 1695),
        (1480, 1695),
        (2150, 1695),
        (2820, 1695),
    ]
    size = (620, 980)
    for i, (image_path, (x, y)) in enumerate(zip(library_examples, positions), start=1):
        example_panel = panel(image_path, size, title=f"Example {i}", subtitle=image_path.stem)
        figure.paste(example_panel, (x, y))

    figure.save(output_path, quality=95)


def main() -> None:
    args = parse_args()
    visualization_root = args.visualization_root
    properties_root = args.properties_root
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    build_enumeration_composite(visualization_root, output_dir / "enumeration-summary.png")
    build_library_examples_composite(visualization_root, output_dir / "top-hit-structures.png")
    build_properties_and_examples_composite(
        visualization_root,
        properties_root,
        output_dir / "top-hit-properties-and-structures.png",
    )

    print(f"Wrote {output_dir / 'enumeration-summary.png'}")
    print(f"Wrote {output_dir / 'top-hit-structures.png'}")
    print(f"Wrote {output_dir / 'top-hit-properties-and-structures.png'}")


if __name__ == "__main__":
    main()
