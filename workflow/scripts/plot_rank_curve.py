#!/usr/bin/env python3
"""Descending rank curve with a linear rank axis and logarithmic size/support axis."""
import argparse
import csv
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter, MaxNLocator

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", required=True, type=Path)
    p.add_argument("--id-column", required=True)
    p.add_argument("--value-column", required=True)
    p.add_argument("--title", required=True)
    p.add_argument("--x-label", required=True)
    p.add_argument("--y-label", required=True)
    p.add_argument("--data-output", required=True, type=Path)
    p.add_argument("--png-output", required=True, type=Path)
    p.add_argument("--svg-output", required=True, type=Path)
    p.add_argument("--pdf-output", required=True, type=Path)
    p.add_argument("--dpi", type=int, default=300)
    p.add_argument("--color", default="#2C7FB8")
    a = p.parse_args()
    with a.input.open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        if not {a.id_column, a.value_column}.issubset(reader.fieldnames or []):
            raise ValueError("Missing plot identifier or value column")
        data = [(r[a.id_column], int(r[a.value_column])) for r in reader]
    if not data or any(v < 1 for _, v in data):
        raise ValueError("Rank curves require positive integer values and nonempty input")
    data.sort(key=lambda x: (-x[1], x[0]))
    a.data_output.parent.mkdir(parents=True, exist_ok=True)
    with a.data_output.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["rank", a.id_column, a.value_column])
        w.writerows((i, identifier, value) for i, (identifier, value) in enumerate(data, 1))
    plt.rcParams.update({"font.family": "DejaVu Sans", "svg.fonttype": "none", "pdf.fonttype": 42})
    fig, ax = plt.subplots(figsize=(9.23, 6.15))
    values = [v for _, v in data]
    ax.plot(range(1, len(data) + 1), values, color=a.color, linewidth=0.85)
    ax.set_xlim(0, max(2, len(data) * 1.05))
    ax.set_ylim(0.7, max(values) * 1.55)
    ax.set_yscale("log")
    ax.yaxis.set_major_locator(LogLocator(base=10))
    ax.yaxis.set_minor_locator(LogLocator(base=10, subs=range(2, 10)))
    ax.yaxis.set_minor_formatter(NullFormatter())
    formatter = FuncFormatter(lambda x, _: f"{x:,.0f}" if x >= 1 else "")
    ax.yaxis.set_major_formatter(formatter)
    ax.xaxis.set_major_formatter(formatter)
    ax.xaxis.set_major_locator(MaxNLocator(nbins=7, integer=True))
    if max(values) < 10:
        ax.set_yticks([value for value in [1, 2, 5] if value <= max(values) * 1.55])
    ax.set_xlabel(a.x_label, fontsize=14)
    ax.set_ylabel(a.y_label, fontsize=14)
    ax.set_title(a.title + "\n" + f"n = {len(values):,}; maximum = {max(values):,}; singletons = {sum(v == 1 for v in values):,}", fontsize=15, pad=14)
    ax.grid(which="major", color="#D9D9D9", linestyle=":", linewidth=0.6)
    ax.grid(which="minor", axis="y", color="#ECECEC", linestyle=":", linewidth=0.4)
    ax.tick_params(labelsize=11)
    fig.tight_layout()
    for path in (a.png_output, a.svg_output, a.pdf_output):
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=a.dpi, facecolor="white")
    plt.close(fig)
    print(f"Rank curve: observations={len(data)}; column={a.value_column}")

if __name__ == "__main__":
    main()
