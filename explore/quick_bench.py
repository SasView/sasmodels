from pathlib import Path

import numpy as np

import timing
import utils


def entry(case: str, backend: str, path: str, sets: int, result: dict, reference: dict | None, error: float | None
          ) -> dict:
    wall, std = result["wall"] / sets, np.sqrt(result["wall_var"]) / sets
    speedup = reference["wall"] / result["wall"] * sets if reference else None
    return dict(case=case, backend=backend, path=path, sets=sets, wall=wall, std=std, speedup=speedup, error=error)


def table(rows: list[dict]) -> str:
    head = "| case | backend | path | sets | wall/set [ms] | x 1-core DLL | rel. error |\n|---|---|---|---|---|---|---|"
    lines = [f"| {r['case']} | {r['backend']} | {r['path']} | {r['sets']} | {r['wall'] * 1e3:.4g} ± {r['std'] * 1e3:.2g} "
             f"| {'-' if r['speedup'] is None else format(r['speedup'], '.3g')} | {'-' if r['error'] is None else format(r['error'], '.1e')} |" for r in rows]
    return "\n".join([head, *lines])


if __name__ == "__main__":
    THREADS = 12
    REPEATS = 5
    WARMUP = 2
    PRECISION = "float64"
    PD_N, PD_NSIGMA = 35, 3
    BATCH = dict(cpu=16, cuda=1024)
    SPREAD = 0.2
    Q_1D = np.logspace(-3.0, -1.0, 200)
    Q_BINS_2D = np.linspace(-0.1, 0.1, 65)
    KINDS = ("Fq", "image")
    RESULTS_FILE = Path(__file__).resolve().parent / "quick.md"
    PD = dict(pd=0.1, pd_n=PD_N, pd_nsigma=PD_NSIGMA)
    CASES = dict(
        sphere=dict(sld=4.0, sld_solvent=1.0, background=0.0, radius=400.0,
                    **{f"radius_{k}": v for k, v in PD.items()}),
        ellipsoid=dict(sld=4.0, sld_solvent=1.0, background=0.0, radius_polar=400.0, radius_equatorial=20.0,
                       theta=60.0, phi=60.0, **{f"radius_polar_{k}": v for k, v in PD.items()}),
        triaxial_ellipsoid=dict(sld=4.0, sld_solvent=1.0, background=0.0, radius_equat_minor=20.0,
                                radius_equat_major=400.0, radius_polar=10.0, theta=60.0, phi=60.0, psi=60.0,
                                **{f"radius_equat_major_{k}": v for k, v in PD.items()},
                                **{f"radius_polar_{k}": v for k, v in PD.items()}))

    timing.environment(THREADS)
    import torch

    import numba_kernel as nk
    import torch_kernel as tk

    devices = ["cpu"] + (["cuda"] if torch.cuda.is_available() else [])
    rows = []
    for name, pars in CASES.items():
        model = utils.MODELS[name]
        for kind in KINDS:
            case = f"{name} {'1D' if kind == 'Fq' else '2D'}"
            print(case, flush=True)
            with utils.single_core():
                call, index, _ = utils.sasmodels_control(model, kind, Q_1D, Q_BINS_2D, 0.0, "double", "dll")
                reference = timing.measure(lambda: lambda: call(pars), REPEATS, WARMUP, False)
            expected = reference["out"]
            rows.append(entry(case, "DLL 1 core", "DirectModel", 1, reference, reference, None))
            backends = [("numba", nk.KernelModel(name, PRECISION))]
            backends += [(f"torch {d}", tk.KernelModel(name, PRECISION, d)) for d in devices]
            for backend, kernel_model in backends:
                call, _ = utils.direct(kernel_model, kind, Q_1D, Q_BINS_2D)
                result = timing.measure(lambda: lambda: call(pars), REPEATS, WARMUP, False)
                rows.append(entry(case, backend, "DirectModel", 1, result, reference,
                                  utils.errors(result["out"], expected, index)[1]))
                if not backend.startswith("torch"):
                    continue
                device = kernel_model.device.type
                sets = utils.batch_sets(model, pars, BATCH[device], SPREAD, 0)
                kernel = kernel_model.make_kernel(call.q_vectors)
                args = timing.measure(lambda: lambda: utils.batch_args(kernel, sets), REPEATS, WARMUP, False)
                rows.append(entry(case, "sasmodels", "get_mesh + make_kernel_args only", len(sets), args, None,
                                  None))
                details, values, magnetic = args["out"]
                result = timing.measure(lambda: lambda: kernel.Iq(details, values, 0.0, magnetic), REPEATS, WARMUP,
                                        False)
                rows.append(entry(case, backend, "kernel batch", len(sets), result, reference,
                                  utils.errors(result["out"][0], expected, index)[1]))
                m, dtype = kernel_model.module, getattr(torch, PRECISION)
                arrays = dict(q=Q_1D) if kind == "Fq" else dict(zip(("qx", "qy"), utils.image_grid(Q_BINS_2D)))
                grid = {k: torch.as_tensor(a, dtype=dtype, device=device) for k, a in arrays.items()}
                calc = tk.Calculator(m, pars, 0.0, **grid)
                packed = torch.from_numpy(utils.pack(m, sets))
                packed = packed.pin_memory() if device == "cuda" else packed
                result = timing.measure(lambda: lambda: calc(packed), REPEATS, WARMUP, device == "cuda")
                out = result["out"][0].double().cpu().numpy()
                rows.append(entry(case, backend, "own batch (Calculator)", len(sets), result, reference,
                                  utils.errors(out, expected, index)[1]))
    text = table(rows)
    RESULTS_FILE.write_text(text + "\n")
    print(text)
