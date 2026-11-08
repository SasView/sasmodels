import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch
from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter

METRICS = dict(wall="Wall clock time", overhead="Overhead (wall clock - on device)", device="On-device time")
COLORS = dict(float64="#2a78d6", float32="#eb6834")
INK, SECONDARY, MUTED, GRID, AXIS, SURFACE, BAND = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7", "#fcfcfb", "#f2f1ec"
NAMES = dict(dll_1core="DLL 1 core", dll_omp="DLL OpenMP", opencl_cpu="OpenCL CPU", opencl_gpu="OpenCL GPU",
             numba="numba", torch_cpu="torch CPU", torch_gpu="torch GPU")
EPSILON = dict(float64=2.2e-16, float32=1.2e-7)
DIMENSIONS = dict(sphere="1d", ellipsoid="2d", triaxial="3d")
BAR, FAMILY_GAP, CASE_GAP, ROW_INCH, LABELED_ROW_INCH = 0.38, 0.4, 0.9, 0.17, 0.25


def per_set(row: dict) -> dict:
    n, wall, device = row["sets"], row["wall"], row.get("device")
    out = dict(wall=wall / n, wall_std=np.sqrt(row["wall_var"]) / n)
    if device is not None:
        out.update(device=device / n, device_std=np.sqrt(row["device_var"]) / n,
                   overhead=max(wall - device, 0.0) / n,
                   overhead_std=np.sqrt(row["wall_var"] + row["device_var"]) / n)
    return out


def device_kind(backend: str) -> str:
    return "gpu" if backend.endswith("gpu") else "cpu"


def references(rows: list[dict]) -> tuple[dict, dict]:
    done = [r for r in rows if "wall" in r]
    dll = {r["case"]: per_set(r)["wall"] for r in done if r["backend"] == "dll_1core"}
    ocl = {(r["case"], r["precision"], device_kind(r["backend"])): per_set(r)["wall"]
           for r in done if r["backend"].startswith("opencl")}
    return dll, ocl


def speedups(row: dict, dll: dict, ocl: dict) -> dict:
    wall, kind = per_set(row)["wall"], device_kind(row["backend"])
    out = {}
    if row["backend"] != "dll_1core" and row["case"] in dll:
        out["DLL"] = dll[row["case"]] / wall
    if not row["backend"].startswith("opencl") and (row["case"], row["precision"], kind) in ocl:
        out["OCL"] = ocl[row["case"], row["precision"], kind] / wall
    return out


def _duration(ms: float) -> str:
    for unit, scale in (("s", 1e3), ("ms", 1.0), ("µs", 1e-3)):
        if ms >= scale:
            return f"{ms / scale:.3g} {unit}"
    return f"{ms / 1e-6:.3g} ns"


def _times(x: float) -> str:
    return f"x{x:.0f}" if x >= 100 else f"x{x:.3g}"


def family(backend: str) -> str:
    return backend if backend.startswith("torch") else backend.split("_")[0]


def run_label(backend: str, mode: str, by_precision: dict) -> str:
    sets = list(dict.fromkeys(r["sets"] for r, *_ in by_precision.values()))
    count = " / ".join(map(str, sets)) + (" set" if sets == [1] else " sets")
    return NAMES.get(backend, backend) + ("" if mode == "single" else f", {mode} ({count})")


def _log_axis(ax: plt.Axes, values: list[float], fmt: str) -> None:
    finite = [v for v in values if np.isfinite(v) and v > 0]
    ax.set_xscale("log")
    if finite:
        ax.set_xlim(10 ** np.floor(np.log10(min(finite))), 10 ** np.ceil(np.log10(max(finite) * 1.5)))
    ax.xaxis.set_major_locator(LogLocator(base=10, numticks=12))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: format(v, fmt)))
    ax.xaxis.set_minor_formatter(NullFormatter())


def _style(ax: plt.Axes, title: str, label: str, grid: str = "y") -> None:
    ax.set_facecolor(SURFACE)
    ax.set_title(title, color=INK, fontsize=10, loc="left")
    (ax.set_xlabel if grid == "x" else ax.set_ylabel)(label, color=SECONDARY, fontsize=8)
    ax.tick_params(colors=SECONDARY, labelsize=7, length=0)
    ax.grid(axis=grid, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(AXIS)


def _groups(rows: list[dict], value) -> dict:
    groups = {}
    for r in rows:
        got = value(r)
        if got is not None:
            groups.setdefault(r["case"], {}).setdefault((r["backend"], r["mode"]), {})[r["precision"]] = (r, *got)
    order = list(NAMES)
    return {case: dict(sorted(runs.items(), key=lambda kv: order.index(kv[0][0]) if kv[0][0] in order else len(order)))
            for case, runs in groups.items()}


def _layout(groups: dict) -> tuple[list, list, list, float]:
    items, headers, bands, y = [], [], [], 0.0
    for case, runs in groups.items():
        start, previous = y, None
        headers.append((y, case))
        y += 1
        for (backend, mode), by_precision in runs.items():
            if previous is not None and family(backend) != previous:
                y += FAMILY_GAP
            previous = family(backend)
            items.append((y, run_label(backend, mode, by_precision), by_precision))
            y += 1
        bands.append((start - 0.5, y - 0.5))
        y += CASE_GAP
    return items, headers, bands, y - CASE_GAP


def plot_bars(rows: list[dict], value, title: str, xlabel: str, path: Path, fmt: str = "g",
              epsilon: bool = False, note=None) -> None:
    groups = _groups(rows, value)
    if not groups:
        return
    items, headers, bands, height = _layout(groups)
    precisions = [p for p in COLORS if any(p in by for *_, by in items)]
    row_inch = ROW_INCH if note is None else LABELED_ROW_INCH
    fig, ax = plt.subplots(figsize=(9.0 if note is None else 10.5, max(3.0, height * row_inch + 1.4)))
    fig.patch.set_facecolor(SURFACE)
    for k, (lo, hi) in enumerate(bands):
        if k % 2:
            ax.axhspan(lo, hi, color=BAND, linewidth=0, zorder=0)
    values = []
    floor = min(EPSILON[p] for p in precisions) / 1e3 if epsilon else None
    for y, _, by_precision in items:
        present = [p for p in precisions if p in by_precision]
        for k, p in enumerate(present):
            r, v, std = by_precision[p]
            at = y + (k - (len(present) - 1) / 2) * BAR
            if epsilon and v == 0:
                ax.text(floor, at, " exact", color=SECONDARY, fontsize=6, va="center")
                continue
            v = max(v, floor) if epsilon else v
            xerr = None if std is None else [[min(std, 0.9 * v)], [std]]
            ax.barh(at, v, BAR * 0.92, xerr=xerr, color=COLORS[p], zorder=2,
                    error_kw=dict(ecolor=SECONDARY, elinewidth=0.8, capsize=1.5, zorder=3))
            values.append(v)
            if note is not None:
                ax.annotate(note(r, v), (v + (std or 0.0), at), xytext=(3, 0), textcoords="offset points",
                            va="center", fontsize=5.5, color=SECONDARY, zorder=4)
    handles = [Patch(color=COLORS[p], label=p) for p in precisions]
    if epsilon:
        for p in precisions:
            handles.append(ax.axvline(EPSILON[p], color=COLORS[p], linewidth=0.9, linestyle="--", zorder=1,
                                      label=f"{p} epsilon"))
        values += [floor, max(values + [EPSILON[p] for p in precisions]) * 100]
    _log_axis(ax, values, fmt)
    if note is not None:
        lo, hi = ax.get_xlim()
        ax.set_xlim(lo, hi * 10 ** np.ceil(np.log10(hi / lo) / 4))
    ticks = [y for y, _ in headers] + [y for y, *_ in items]
    labels = [case for _, case in headers] + [text for _, text, _ in items]
    ax.set_yticks(ticks, labels)
    for tick, bold in zip(ax.get_yticklabels(), [True] * len(headers) + [False] * len(items)):
        tick.set(color=INK if bold else SECONDARY, fontweight="bold" if bold else "normal", fontsize=8 if bold else 7)
    ax.set_ylim(height + 0.5, -0.5)
    _style(ax, "", xlabel, grid="x")
    ax.tick_params(axis="x", labeltop=True)
    ax.spines["bottom"].set_visible(False)
    top = 1 - 0.62 / fig.get_figheight()
    fig.text(0.01, 1 - 0.18 / fig.get_figheight(), title, color=INK, fontsize=11, va="top")
    fig.legend(handles=handles, frameon=False, fontsize=8, labelcolor=SECONDARY, ncol=len(handles),
               loc="upper left", bbox_to_anchor=(0.005, 1 - 0.42 / fig.get_figheight()))
    fig.tight_layout(rect=(0, 0, 1, top))
    fig.savefig(path, dpi=150, facecolor=SURFACE)
    plt.close(fig)


def _metric(metric: str):
    def value(row: dict):
        p = per_set(row)
        if metric not in p or not p[metric] > 0:
            return None
        return p[metric] * 1e3, p[metric + "_std"] * 1e3
    return value


def _error(row: dict):
    if row["backend"] == "dll_1core" or "error_rel" not in row:
        return None
    return row["error_rel"], None


def plot_gpu(rows: list[dict], path: Path) -> None:
    gpu = [r for r in rows if "timeline" in r]
    if not gpu:
        return
    cases = list(dict.fromkeys(r["case"] for r in gpu))
    modes = list(dict.fromkeys(r["mode"] for r in gpu))
    fig, axes = plt.subplots(len(cases), len(modes), figsize=(4.2 * len(modes), 1.9 * len(cases)), squeeze=False)
    fig.patch.set_facecolor(SURFACE)
    for i, case in enumerate(cases):
        for j, mode in enumerate(modes):
            ax = axes[i, j]
            runs = [r for r in gpu if r["case"] == case and r["mode"] == mode]
            for r in runs:
                ax.plot(r["timeline"]["time_ms"], np.asarray(r["timeline"]["busy"]) * 100, color=COLORS[r["precision"]],
                        linewidth=1.5, label=r["precision"])
            ax.set_ylim(0, 105)
            sets = " / ".join(dict.fromkeys(str(r["sets"]) for r in runs))
            _style(ax, f"{case}, {mode} ({sets} sets)", "GPU busy [%]")
            ax.set_xlabel("time in call [ms]", color=SECONDARY, fontsize=7)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, frameon=False, fontsize=8, labelcolor=SECONDARY, loc="upper center", ncol=len(labels))
    fig.tight_layout(rect=(0, 0, 1, 1 - 0.3 / fig.get_figheight()))
    fig.savefig(path, dpi=150, facecolor=SURFACE)
    plt.close(fig)


def _ms(seconds: float | None, std: float | None = None) -> str:
    if seconds is None or np.isnan(seconds):
        return "-"
    return f"{seconds * 1e3:.4g}" + (f" ± {std * 1e3:.2g}" if std is not None else "")


def table(results: dict) -> str:
    rows, meta = results["rows"], results["meta"]
    dll, ocl = references(rows)
    head = ("case", "prec", "backend", "mode", "sets", "wall/set [ms]", "overhead/set [ms]", "device/set [ms]",
            "startup [s]", "warmup [ms]", "x 1-core DLL", "x OpenCL (same device)", "error norm", "rel. error")
    lines = [f"speed summary {meta['date']}: python {meta['python']}, numpy {meta['numpy']}, numba {meta['numba']} "
             f"({meta['threading_layer']}), torch {meta['torch']} ({meta['torch_compile']}), GPU {meta['gpu']}, "
             f"{meta['threads']} threads, mean ± std of {meta['repeats']} calls after {meta['warmup']} warm-up calls",
             "", "| " + " | ".join(head) + " |", "|" + "---|" * len(head)]
    for r in rows:
        if "wall" not in r:
            lines.append(f"| {r['case']} | {r['precision']} | {r['backend']} | {r['mode']} | "
                         f"{r.get('skipped') or r.get('error')} |" + " |" * (len(head) - 5))
            continue
        p = per_set(r)
        speedup = dll.get(r["case"], np.nan) / p["wall"]
        speedup_ocl = ocl.get((r["case"], r["precision"], device_kind(r["backend"])), np.nan) / p["wall"]
        cells = (r["case"], r["precision"], r["backend"], r["mode"], str(r["sets"]), _ms(p["wall"], p["wall_std"]),
                 _ms(p.get("overhead"), p.get("overhead_std")), _ms(p.get("device"), p.get("device_std")),
                 f"{r['startup']:.3g}", _ms(r.get("warmup")), f"{speedup:.3g}", f"{speedup_ocl:.3g}",
                 f"{r.get('error_norm', np.nan):.2e}", f"{r.get('error_rel', np.nan):.2e}")
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


def plot(results: dict | str | Path, out_dir: str | Path) -> None:
    if not isinstance(results, dict):
        results = json.loads(Path(results).read_text())
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    rows = [r for r in results["rows"] if "wall" in r]
    dll, ocl = references(rows)

    def wall_note(row: dict, v: float) -> str:
        fast = speedups(row, dll, ocl)
        return " | ".join([_duration(v)] + [f"{_times(fast[name])} {name}" for name in ("OCL", "DLL") if name in fast])

    split = {}
    for r in rows:
        prefix = r["case"].split("_")[0]
        split.setdefault((DIMENSIONS.get(prefix, prefix), prefix), []).append(r)
    for (dim, prefix), subset in split.items():
        plot_bars(subset, _metric("wall"), f"{METRICS['wall']} per parameter set, {dim} ({prefix})",
                  "ms per parameter set", out_dir / f"timing_wall_{dim}.png", note=wall_note)
        plot_bars(subset, _error, f"Relative error vs 1-core DLL, {dim} ({prefix})", "relative L2 error",
                  out_dir / f"errors_{dim}.png", ".0e", epsilon=True)
    for metric in ("overhead", "device"):
        plot_bars(rows, _metric(metric), f"{METRICS[metric]} per parameter set", "ms per parameter set",
                  out_dir / f"timing_{metric}.png", note=lambda row, v: _duration(v))
    plot_gpu(rows, out_dir / "gpu_usage.png")
    text = table(results)
    (out_dir / "table.md").write_text(text + "\n")
    print(text)


if __name__ == "__main__":
    import sys

    plot(sys.argv[1] if len(sys.argv) > 1 else "results.json", sys.argv[2] if len(sys.argv) > 2 else "plots")
