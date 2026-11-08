import gc
import json
import os
import platform
import tempfile
import time
from collections.abc import Callable
from contextlib import nullcontext
from pathlib import Path

import numpy as np

import utils


def measure(build: Callable[[], Callable[[], object]], repeats: int, warmup: int, cuda: bool) -> dict:
    import torch

    sync = torch.cuda.synchronize if cuda else (lambda: None)
    start = time.perf_counter()
    fn = build()
    out = fn()
    sync()
    startup = time.perf_counter() - start
    warm, wall, device = [], [], []
    for i in range(warmup + repeats):
        if cuda:
            begin, end = torch.cuda.Event(enable_timing=True), torch.cuda.Event(enable_timing=True)
            begin.record()
        start = time.perf_counter()
        out = fn()
        if cuda:
            end.record()
        sync()
        (warm if i < warmup else wall).append(time.perf_counter() - start)
        if cuda and i >= warmup:
            device.append(begin.elapsed_time(end) / 1e3)
    return dict(startup=startup, warmup=float(np.mean(warm)) if warm else None, wall=float(np.mean(wall)),
                wall_var=float(np.var(wall)), device=float(np.mean(device)) if cuda else None,
                device_var=float(np.var(device)) if cuda else None, fn=fn, out=out)


def gpu_timeline(fn: Callable[[], object], bins: int) -> dict:
    import torch
    from torch.profiler import ProfilerActivity, profile, record_function

    with profile(activities=[ProfilerActivity.CPU, ProfilerActivity.CUDA]) as prof:
        with record_function("call"):
            fn()
            torch.cuda.synchronize()
    with tempfile.TemporaryDirectory() as tmp:
        prof.export_chrome_trace(str(Path(tmp) / "trace.json"))
        events = json.loads((Path(tmp) / "trace.json").read_text())["traceEvents"]
    call = next(e for e in events if e.get("name") == "call" and e.get("cat") == "user_annotation")
    kernels = [(e["ts"] - call["ts"], e["ts"] - call["ts"] + e["dur"]) for e in events if e.get("cat") == "kernel"]
    edges = np.linspace(0.0, call["dur"], bins + 1)
    busy = [sum(max(0.0, min(b, e) - max(a, s)) for s, e in kernels) / (b - a) for a, b in zip(edges[:-1], edges[1:])]
    return dict(time_ms=(edges[:-1] / 1e3).tolist(), busy=np.minimum(busy, 1.0).tolist())


def builders(case: dict, precision: str, s: dict, controls: bool) -> list[tuple]:
    import torch

    import numba_kernel as nk
    import torch_kernel as tk

    model, pars, kind = utils.MODELS[case["model"]], case["pars"], case["kind"]
    q, q_bins = (getattr(np, k)(a, b, n) for k, a, b, n in (s["q_1d"], s["q_bins_2d"]))
    qx, qy = utils.image_grid(q_bins)

    def control(platform: str, dtype: str, device: str | None = None) -> Callable:
        def build() -> Callable[[], object]:
            name = utils.opencl_device(device) if device else None
            call, index, info = utils.sasmodels_control(model, kind, q, q_bins, s["cutoff"], dtype, platform)
            fn = lambda: call(pars)  # noqa: E731
            fn.index, fn.info = index, dict(info, device=name)
            return fn
        return build

    def numba() -> Callable[[], object]:
        m = nk.load_model(case["model"], precision)
        grid = dict(q=q.astype(precision)) if kind == "Fq" else dict(qx=qx.astype(precision), qy=qy.astype(precision))
        return lambda: nk.calculate(m, pars, s["cutoff"], **grid)

    def torch_run(device: torch.device, factor: float | None) -> Callable:
        def build() -> Callable[[], object]:
            m, dtype = tk.load_model(case["model"]), getattr(torch, precision)
            arrays = (q,) if kind == "Fq" else (qx, qy)
            grid = {k: torch.as_tensor(a, dtype=dtype, device=device)
                    for k, a in zip(("q",) if kind == "Fq" else ("qx", "qy"), arrays)}
            calc = tk.Calculator(m, pars, s["cutoff"], **grid)
            if factor is None:
                return lambda: calc(torch.from_numpy(utils.pack(m, [pars])))
            sets = tk.batch_size(m, pars, s["cutoff"], factor, **grid)
            pool = [torch.from_numpy(utils.pack(m, utils.batch_sets(m, pars, sets, s["spread"], seed)))
                    for seed in range(s["warmup"] + s["repeats"] + 2)]
            pool = [p.pin_memory() for p in pool] if device.type == "cuda" else pool
            calls = [0]

            def call() -> object:
                calls[0] += 1
                return calc(pool[calls[0] % len(pool)])
            call.sets = sets
            return call
        return build

    cl_dtype = "double" if precision == "float64" else "single"
    runs = [("dll_1core", "single", control("dll", "double"), utils.single_core),
            ("dll_omp", "single", control("dll", "double"), utils.openmp_build)] if controls else []
    runs += [("opencl_cpu", "single", control("ocl", cl_dtype, "cpu"), None),
             ("opencl_gpu", "single", control("ocl", cl_dtype, "gpu"), None),
             ("numba", "single", numba, None)]
    for device in [torch.device("cpu")] + ([torch.device("cuda", 0)] if torch.cuda.is_available() else []):
        label = "torch_gpu" if device.type == "cuda" else "torch_cpu"
        runs.append((label, "single", torch_run(device, None), None))
        for factor in s["batch_factors"] if device.type == "cuda" else s["batch_factors"][:1]:
            runs.append((label, "batched" if factor == 1 else f"batched x{factor:g}", torch_run(device, factor), None))
    return runs


def environment(threads: int) -> None:
    path = os.pathsep.join(p for p in os.environ.get("PATH", "").split(os.pathsep) if not os.path.isfile(p))
    os.environ.update(PATH=path, OMP_NUM_THREADS=str(threads), NUMBA_NUM_THREADS=str(threads),
                      NUMBA_THREADING_LAYER="omp", OMP_PROC_BIND="close", OMP_PLACES="cores",
                      POCL_MAX_PTHREAD_COUNT=str(threads))
    import numba
    import torch

    numba.set_num_threads(threads)
    torch.set_num_threads(threads)


def run(s: dict) -> dict:
    environment(s["threads"])
    import numba
    import torch

    import torch_kernel as tk

    tk.WORK_PER_THREAD, tk.MEMORY_FRACTION = s["work_per_thread"], s["memory_fraction"]
    gpu = torch.device("cuda", 0) if torch.cuda.is_available() else None
    compiled = tk.compile_engine(gpu or torch.device("cpu"))
    rows, references = [], {}
    for name, case in s["cases"].items():
        for precision in s["precisions"]:
            if gpu:
                tk.tune(gpu, getattr(torch, precision))
            for backend, mode, build, context in builders(case, precision, s, name not in references):
                row = dict(case=name, precision="float64" if backend.startswith("dll") else precision,
                           backend=backend, mode=mode)
                print(f"{name} {row['precision']} {backend} {mode}", flush=True)
                try:
                    with context() if context else nullcontext():
                        result = measure(build, s["repeats"], s["warmup"], backend == "torch_gpu")
                    fn, out = result.pop("fn"), result.pop("out")
                    out = out.detach()[0].double().cpu().numpy() if hasattr(out, "detach") else np.asarray(out)
                    if backend == "dll_1core":
                        references[name] = (out.reshape(-1), fn.index)
                    row.update(result, sets=getattr(fn, "sets", 1), info=getattr(fn, "info", {}))
                    if name in references:
                        row["error_norm"], row["error_rel"] = utils.errors(out, *references[name])
                    if backend == "torch_gpu":
                        row["timeline"] = gpu_timeline(fn, s["timeline_bins"])
                    del fn, out
                except utils.Skip as reason:
                    row["skipped"] = str(reason)
                except Exception as error:
                    row["error"] = f"{type(error).__name__}: {str(error).strip()[:300]}"
                gc.collect()
                if gpu:
                    torch.cuda.empty_cache()
                rows.append(row)
    meta = dict(date=time.strftime("%Y-%m-%d %H:%M"), python=platform.python_version(), numpy=np.__version__,
                numba=numba.__version__, threading_layer=numba.threading_layer(), torch=torch.__version__,
                torch_compile=compiled, gpu=torch.cuda.get_device_name(gpu) if gpu else None,
                threads=s["threads"], repeats=s["repeats"], warmup=s["warmup"], q_1d=s["q_1d"],
                q_bins_2d=s["q_bins_2d"])
    results = dict(meta=meta, rows=rows)
    Path(s["results_file"]).write_text(json.dumps(results, indent=1, default=str))
    return results
