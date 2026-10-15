#!/usr/bin/env python3
"""Plot the saved histogram tables with the step-line style of the reference.

Requires numpy and matplotlib. Run from any directory:
    python3 /path/to/data-for-delphes/plot_histograms.py

Signal/background PNGs are saved alongside their input tables. Overlays and 2-by-2
panels are saved in histograms/comparisons. Only PNG files are replaced;
the C++ analysis and histogram tables are never modified.

Density normalization uses the sum of the saved bins and their widths.
Underflow/overflow bins are absent from the tables and cannot be included.
The angle definitions from the analysis are preserved: only the b-jet receives
the reference histo.cpp transformation (marked with a star). Background tables
are read from histograms/m<MASS>/background/<observable>_m<MASS>.txt.
Names ending in _background.txt are also accepted. Each signal requires a
background with the same observable and mass. Backgrounds may differ across
mass directories: comparisons then use separate panels, each with its own
signal and background. Identical backgrounds share one comparison curve.
"""

import argparse
from dataclasses import dataclass
from pathlib import Path
import re
from typing import Optional

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MaxNLocator, ScalarFormatter
import numpy as np


LABELS = {
    "cost_bbar": r"$\cos\theta^{*}(j_b)$",
    "dphi_j1_bbar": r"$\Delta\phi(j_1,j_b^{*})$",
    "dphi_a1_bbar": r"$\Delta\phi(\gamma_1,j_b^{*})$",
    "dphi_a2_bbar": r"$\Delta\phi(\gamma_2,j_b^{*})$",
    "dphi_a1_j1": r"$\Delta\phi(\gamma_1,j_1)$",
    "dphi_j1_j2": r"$\Delta\phi(j_1,j_2)$",
    "dphi_a1_j2": r"$\Delta\phi(\gamma_1,j_2)$",
}
COLORS = {10: "#0000ff", 50: "#ff0000", 100: "#ff9933"}


@dataclass
class Histogram:
    source: Path
    observable: str
    mass: Optional[int]
    edges: np.ndarray
    contents: np.ndarray
    is_background: bool = False

    def heights(self, normalization):
        if normalization == "none" or not np.any(self.contents):
            return self.contents.copy()
        total = self.contents.sum()
        if total <= 0:
            raise ValueError(f"Cannot normalize a non-positive bin sum: {self.source}")
        return self.contents / (total * np.diff(self.edges))


def load_histograms(hist_dir):
    histograms = []
    sources = list(hist_dir.glob("m*/*.txt")) + list(hist_dir.glob("m*/background/*.txt"))
    for source in sorted(sources):
        is_background = source.parent.name == "background"
        pattern = r"(.+?)(?:_m(\d+))?(?:_background)?" if is_background else r"(.+)_m(\d+)"
        match = re.fullmatch(pattern, source.stem)
        if not match or match[1] not in LABELS:
            continue
        observable = match[1]
        if is_background:
            directory_mass = re.fullmatch(r"m(\d+)", source.parent.parent.name)
            if not directory_mass:
                continue
            mass = int(directory_mass[1])
            if match[2] is not None and int(match[2]) != mass:
                raise ValueError(f"Background filename mass differs from its directory: {source}")
        else:
            mass = int(match[2])
        data = np.loadtxt(source, skiprows=1, ndmin=2)
        if data.shape[1] != 2 or len(data) < 2 or not np.isfinite(data).all():
            raise ValueError(f"Expected finite bin-center/content pairs: {source}")

        # These are the original analysis ranges, not a new binning choice.
        limits = (-1.0, 1.0) if observable == "cost_bbar" else (-3.14159, 3.14159)
        edges = np.linspace(*limits, len(data) + 1)
        if not np.allclose(data[:, 0], (edges[:-1] + edges[1:]) / 2, atol=1e-5):
            raise ValueError(f"Bin centers do not match the analysis range: {source}")
        histograms.append(Histogram(source, observable, mass, edges, data[:, 1], is_background))
    if not histograms:
        raise ValueError(f"No recognized histogram tables in {hist_dir}")
    seen = set()
    for hist in histograms:
        key = (hist.observable, hist.mass, hist.is_background)
        if key in seen:
            raise ValueError(f"Duplicate histogram for {key}: {hist.source}")
        seen.add(key)
    return sorted(histograms, key=lambda hist: (hist.is_background, hist.mass or 0, hist.observable))


def comparison_groups(histograms):
    """Collapse identical backgrounds; retain distinct ones for separate panels."""
    groups = {}
    for name in LABELS:
        group = [hist for hist in histograms if hist.observable == name and not hist.is_background]
        backgrounds = [hist for hist in histograms if hist.observable == name and hist.is_background]
        if backgrounds:
            first = backgrounds[0]
            identical = all(
                np.array_equal(first.edges, background.edges) and
                np.array_equal(first.contents, background.contents)
                for background in backgrounds[1:]
            )
            group.extend([first] if identical else backgrounds)
        if group:
            groups[name] = group
    return groups


def comparison_panels(group):
    """Each panel contains at most one background distribution."""
    backgrounds = {hist.mass: hist for hist in group if hist.is_background}
    if len(backgrounds) <= 1:
        return [group]
    return [[signal, backgrounds[signal.mass]] for signal in group if not signal.is_background]


def configure_style():
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "mathtext.fontset": "stixsans",
        "axes.labelsize": 12,
        "axes.labelpad": 3,
        "axes.linewidth": 0.8,
        "axes.grid": False,
        "axes.spines.top": True,
        "axes.spines.right": True,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "xtick.major.size": 6,
        "ytick.major.size": 6,
        "xtick.minor.size": 3,
        "ytick.minor.size": 3,
        "xtick.major.width": 0.7,
        "ytick.major.width": 0.7,
        "xtick.minor.width": 0.6,
        "ytick.minor.width": 0.6,
        "legend.frameon": False,
        "legend.fontsize": 9,
        "legend.handlelength": 2.1,
        "legend.handletextpad": 0.5,
        "legend.labelspacing": 1.2,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
    })


def draw_histograms(ax, histograms, normalization, annotation=""):
    observable = histograms[0].observable
    background_count = sum(hist.is_background for hist in histograms)
    if background_count > 1:
        raise ValueError(f"Expected one shared background for {observable}")
    minimum, maximum = 0.0, 0.0
    for hist in sorted(histograms, key=lambda item: (item.is_background, item.mass or 0)):
        heights = hist.heights(normalization)
        minimum = min(minimum, float(heights.min()))
        maximum = max(maximum, float(heights.max()))
        if hist.is_background:
            label = "Background"
            color, linestyle = "black", ":"
        else:
            label = rf"$m_a = {hist.mass}\ \mathrm{{GeV}}$"
            color, linestyle = COLORS.get(hist.mass, "#00aa00"), "-"
        if not np.any(hist.contents):
            label += " (no events)"
        ax.stairs(
            heights, hist.edges, baseline=None, fill=False,
            linewidth=1.0 if hist.is_background else 0.9,
            color=color, linestyle=linestyle, label=label,
        )

    ax.set_xlim(histograms[0].edges[[0, -1]])
    headroom = 1.5 if len(histograms) > 3 else 1.35
    ax.set_ylim(minimum * 1.15, maximum * headroom if maximum > 0 else 1.0)
    ax.set_xlabel(LABELS[observable], loc="right")
    ylabel = r"$(1/\sigma)\,\mathrm{d}\sigma/\mathrm{d}\xi$" if normalization == "density" else "Weighted bin content"
    ax.set_ylabel(ylabel, loc="top")
    ax.xaxis.set_major_locator(MaxNLocator(nbins=6))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=7))
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax.tick_params(which="both", top=True, right=True)
    formatter = ScalarFormatter(useMathText=True)
    formatter.set_powerlimits((-3, 3))
    formatter.set_useOffset(False)
    ax.yaxis.set_major_formatter(formatter)
    ax.legend(loc="upper right", borderaxespad=0.9,
              fontsize=8 if len(histograms) > 3 else 9,
              labelspacing=0.65 if len(histograms) > 3 else 1.2)
    if annotation:
        ax.text(0.025, 0.96, annotation, transform=ax.transAxes,
                va="top", fontsize=8, fontweight="bold")
    if maximum == 0 and minimum == 0:
        ax.text(0.5, 0.5, "No selected events", transform=ax.transAxes, ha="center")


def save_figure(fig, destination, dpi, normalization):
    normalization_note = "Each distribution normalized independently" if normalization == "density" else "Original weighted bin contents"
    fig.savefig(destination, dpi=dpi, metadata={
        "Description": f"Saved signal/background histogram bins; normalization={normalization}. "
                       f"{normalization_note}; underflow/overflow excluded. "
                       "Signals are paired with backgrounds from the matching mass directory. "
                       "Only the b-jet has the reference histo.cpp transformation."
    })
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--hist-dir", type=Path, default=Path(__file__).resolve().parent / "histograms")
    parser.add_argument("--normalization", choices=("density", "none"), default="density",
                        help="density: unit area over saved bins (default); none: original weighted contents")
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--annotation", default="",
                        help="Optional verified experiment label; no beam energy or luminosity is assumed")
    args = parser.parse_args()
    if args.dpi <= 0:
        parser.error("--dpi must be positive")
    try:
        histograms = load_histograms(args.hist_dir)
        groups = comparison_groups(histograms)
    except (ValueError, OSError) as error:
        parser.error(str(error))
    configure_style()

    signals = [hist for hist in histograms if not hist.is_background]
    backgrounds = {(hist.observable, hist.mass): hist for hist in histograms if hist.is_background}
    if not signals:
        parser.error("No signal histograms were found")
    missing = [hist.source.parent / "background" / hist.source.name
               for hist in signals if (hist.observable, hist.mass) not in backgrounds]
    if missing:
        parser.error("Missing matching background tables:\n" + "\n".join(map(str, missing)))
    for hist in signals:
        group = [hist, backgrounds[(hist.observable, hist.mass)]]
        fig, ax = plt.subplots(figsize=(5.4, 4.6))
        fig.subplots_adjust(left=0.15, right=0.965, bottom=0.15, top=0.95)
        draw_histograms(ax, group, args.normalization, args.annotation)
        save_figure(fig, hist.source.with_suffix(".png"), args.dpi, args.normalization)

    for hist in backgrounds.values():
        fig, ax = plt.subplots(figsize=(5.4, 4.6))
        fig.subplots_adjust(left=0.15, right=0.965, bottom=0.15, top=0.95)
        draw_histograms(ax, [hist], args.normalization, args.annotation)
        save_figure(fig, hist.source.with_suffix(".png"), args.dpi, args.normalization)

    comparison_dir = args.hist_dir / "comparisons"
    comparison_dir.mkdir(parents=True, exist_ok=True)
    for name, group in groups.items():
        panels = comparison_panels(group)
        if len(panels) == 1:
            fig, axes = plt.subplots(figsize=(5.4, 4.6))
            fig.subplots_adjust(left=0.15, right=0.965, bottom=0.15, top=0.95)
            axes = [axes]
        else:
            fig, axes = plt.subplots(1, len(panels), figsize=(5.4 * len(panels), 4.6))
            fig.subplots_adjust(left=0.06, right=0.985, bottom=0.15, top=0.95, wspace=0.28)
        for ax, panel in zip(axes, panels):
            draw_histograms(ax, panel, args.normalization, args.annotation)
        save_figure(fig, comparison_dir / f"{name}.png", args.dpi, args.normalization)

    panels = [panel for group in groups.values() for panel in comparison_panels(group)]
    for start in range(0, len(panels), 4):
        fig, axes = plt.subplots(2, 2, figsize=(10.4, 8.7))
        fig.subplots_adjust(left=0.085, right=0.975, bottom=0.08, top=0.97, wspace=0.24, hspace=0.23)
        batch = panels[start:start + 4]
        for ax, group in zip(axes.flat, batch):
            draw_histograms(ax, group, args.normalization, args.annotation)
        for ax in list(axes.flat)[len(batch):]:
            ax.set_visible(False)
        destination = comparison_dir / f"comparison_panel_{start // 4 + 1}.png"
        save_figure(fig, destination, args.dpi, args.normalization)

    matched = sum((hist.observable, hist.mass) in backgrounds for hist in signals)
    print(f"Updated {len(signals)} signal PNGs; {matched} include background from their mass directory.")
    print(f"Saved {len(backgrounds)} background PNGs alongside their tables in {args.hist_dir}/m*/background.")
    print(f"Saved {len(groups)} mass comparisons and {(len(panels) + 3) // 4} panel sheets in {comparison_dir}.")
    print(f"Normalization: {args.normalization}. Analysis code and input tables were not modified.")


if __name__ == "__main__":
    main()
