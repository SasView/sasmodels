import math
import sys
import types
from collections.abc import Callable
from contextlib import nullcontext

import numpy as np
import torch
from torch import Tensor

import utils

ADAPTIVE_MAX_N = 100000
ADAPTIVE_MAX_OUTER = ADAPTIVE_MAX_N // 500
N_POINTS = (20, 76, 500, 5000)
CHUNK_ELEMS = 1 << 24
MAX_BATCH = 1 << 24
WORK_PER_THREAD = 4096
WAVES_PER_LAUNCH = 128
MEMORY_FRACTION = 0.5
POINT_ARRAYS = 24
OUTPUT_COPIES = 4
CPU_LANES = 8
MAX_SETS = 4096
COMPILE_MIN_CAPABILITY = (7, 0)
IMAGE_MODES = dict(Iq=0, Iqac=1, Iqabc=2)
GRAPH_WARMUP = 2
COUNT = dict(evaluations=0, integrations=0)
COUNT_ONLY = False
_TABLES: dict[tuple, tuple[Tensor, Tensor]] = {}


def roots_legendre(n: int) -> tuple[np.ndarray, np.ndarray]:
    i = np.arange(1, n + 1, dtype=np.float64)
    x = (1 - 1 / (8 * n**2) + 1 / (8 * n**3)) * np.cos(
        np.pi * (4 * i - 1) / (4 * n + 2)
    )
    for _ in range(20):
        p0, p1 = np.ones_like(x), x
        for k in range(2, n + 1):
            p0, p1 = p1, ((2 * k - 1) * x * p1 - (k - 1) * p0) / k
        dp = n * (x * p1 - p0) / (x * x - 1)
        dx = p1 / dp
        x = x - dx
        if np.abs(dx).max() < 1e-15:
            break
    return x[::-1].copy(), (2 / ((1 - x * x) * dp * dp))[::-1].copy()


def square(x: Tensor) -> Tensor:
    return x * x


def cube(x: Tensor) -> Tensor:
    return x * x * x


def cbrt(x: Tensor) -> Tensor:
    return torch.sign(x) * torch.abs(x) ** (1.0 / 3.0)


def sas_3j1x_x(x: Tensor) -> Tensor:
    small = torch.abs(x) < (0.1 if x.dtype == torch.float64 else 0.7)
    x2 = x * x
    series = 1.0 + x2 * (-3.0 / 30.0 + x2 * (3.0 / 840.0 + x2 * (-3.0 / 45360.0)))
    s = torch.where(small, 1.0, x)
    return torch.where(small, series, 3.0 * (torch.sin(s) / s - torch.cos(s)) / (s * s))


def gauss_table(n: int, lo: float, hi: float, like: Tensor) -> tuple[Tensor, Tensor]:
    key = (n, float(lo), float(hi), like.dtype, like.device)
    if key not in _TABLES:
        z, w = roots_legendre(n)
        half, center = 0.5 * (hi - lo), 0.5 * (hi + lo)
        _TABLES[key] = (
            torch.as_tensor(half * z + center, dtype=like.dtype, device=like.device),
            torch.as_tensor(half * w, dtype=like.dtype, device=like.device),
        )
    return _TABLES[key]


def gauss_weights(qr: Tensor, n_outer: int | Tensor) -> Tensor:
    n = torch.as_tensor(n_outer, device=qr.device)
    rule = (qr >= 10.0).long() + (qr >= 100.0).long() + (qr >= 750.0).long()
    cap = (
        3
        - (n > ADAPTIVE_MAX_N // 5000).long()
        - (n > ADAPTIVE_MAX_N // 500).long()
        - (n > ADAPTIVE_MAX_N // 76).long()
    )
    return torch.minimum(rule, cap)


def _sum_1d(f: Callable, x: Tensor, w: Tensor, *args: Tensor) -> tuple[Tensor, Tensor]:
    F1, F2 = f(x, *[t[:, None] for t in args])
    return F1 @ w, F2 @ w


def _sum_2d(
    f: Callable,
    x_outer: Tensor,
    w_outer: Tensor,
    x_inner: Tensor,
    w_inner: Tensor,
    *args: Tensor,
) -> tuple[Tensor, Tensor]:
    F1, F2 = f(x_outer[:, None], x_inner[None, :], *[t[:, None, None] for t in args])
    return (F1 @ w_inner) @ w_outer, (F2 @ w_inner) @ w_outer


SUM = {1: _sum_1d, 2: _sum_2d}


def _integrate(
    f: Callable, rules: tuple[Tensor, ...], ranges: tuple, args: tuple[Tensor, ...]
) -> tuple[Tensor, Tensor]:
    shape = torch.broadcast_shapes(*[t.shape for t in (*rules, *args)])
    n_rules = len(N_POINTS)
    key = (
        sum(r * n_rules ** (len(rules) - 1 - d) for d, r in enumerate(rules))
        .expand(shape)
        .reshape(-1)
    )
    args = [t.expand(shape).reshape(-1) for t in args]
    values, order = torch.sort(key, stable=True)
    bounds = torch.searchsorted(
        values,
        torch.arange(n_rules ** len(rules) + 1, device=key.device) - 1,
        right=True,
    ).tolist()
    F1, F2 = torch.zeros_like(args[0]), torch.zeros_like(args[0])
    COUNT["integrations"] += 1
    for g in range(len(bounds) - 1):
        start, stop = bounds[g], bounds[g + 1]
        if start == stop:
            continue
        counts = [
            N_POINTS[(g // n_rules ** (len(rules) - 1 - d)) % n_rules]
            for d in range(len(rules))
        ]
        COUNT["evaluations"] += (stop - start) * math.prod(counts)
        if COUNT_ONLY:
            continue
        tables = [
            x
            for n, (lo, hi) in zip(counts, ranges)
            for x in gauss_table(n, lo, hi, args[0])
        ]
        step = max(1, CHUNK_ELEMS // math.prod(counts))
        for s in range(start, stop, step):
            sel = order[s : min(s + step, stop)]
            sub = [t[sel] for t in args]
            for t in sub:
                torch._dynamo.maybe_mark_dynamic(t, 0)
            F1[sel], F2[sel] = SUM[len(rules)](f, *tables, *sub)
    return F1.reshape(shape), F2.reshape(shape)


@torch._dynamo.disable
def integrate(
    f: Callable, qr: Tensor, lo: float, hi: float, *args: Tensor
) -> tuple[Tensor, Tensor]:
    return _integrate(f, (gauss_weights(qr, 1),), ((lo, hi),), args)


@torch._dynamo.disable
def integrate_fixed(
    f: Callable, n: int, lo: float, hi: float, *args: Tensor
) -> tuple[Tensor, Tensor]:
    return _integrate(
        f, (torch.tensor(N_POINTS.index(n), device=args[0].device),), ((lo, hi),), args
    )


@torch._dynamo.disable
def integrate_2d(
    f: Callable,
    qr_outer: Tensor,
    lo_outer: float,
    hi_outer: float,
    qr_inner: Tensor,
    lo_inner: float,
    hi_inner: float,
    *args: Tensor,
) -> tuple[Tensor, Tensor]:
    outer = gauss_weights(qr_outer, ADAPTIVE_MAX_OUTER)
    inner = gauss_weights(qr_inner, torch.tensor(N_POINTS, device=outer.device)[outer])
    return _integrate(
        f, (outer, inner), ((lo_outer, hi_outer), (lo_inner, hi_inner)), args
    )


API = dict(
    kernel=lambda f: f,
    M_4PI_3=4.0 * math.pi / 3.0,
    M_PI_2=math.pi / 2.0,
    sqrt=torch.sqrt,
    sin=torch.sin,
    cos=torch.cos,
    log=torch.log,
    asin=torch.asin,
    cbrt=cbrt,
    square=square,
    cube=cube,
    fmin=torch.minimum,
    fmax=torch.maximum,
    where=torch.where,
    sas_3j1x_x=sas_3j1x_x,
    integrate=integrate,
    integrate_fixed=integrate_fixed,
    integrate_2d=integrate_2d,
)


def load_model(name: str) -> types.ModuleType:
    model = utils.MODELS[name]
    m = sys.modules[f"torch_{name}"] = types.ModuleType(f"torch_{name}")
    m.__dict__.update(API, **utils.model_metadata(model))
    for fn, source in utils.model_functions(model):
        exec(compile(source, f"{name}.{fn}", "exec"), m.__dict__)
    sizes = [m.PARAMETERS.index(n) for n in m.VOLUME_PARAMETERS]
    m.kernel_Fq = lambda q, p: m.Fq(q, *p)
    m.kernel_Iqxy = dict(
        Iq=lambda qa, qb, qc, p: m.Fq(qa, *p)[1],
        Iqac=lambda qa, qb, qc, p: m.Iqac(qa, qb, *p),
        Iqabc=lambda qa, qb, qc, p: m.Iqabc(qa, qb, qc, *p),
    )[m.IMAGE]
    m.kernel_volume = lambda p: (m.form_volume(*[p[i] for i in sizes]),) * 2
    m.kernel_radius_effective = lambda mode, p: m.radius_effective(
        mode, *[p[i] for i in sizes]
    )
    m.kernel_valid = lambda p: m.valid(*[p[i] for i in sizes])
    m.mode = IMAGE_MODES[m.IMAGE]
    return m


def pd_tables(layout: tuple, like: Tensor) -> tuple[tuple[Tensor, Tensor], ...]:
    return tuple(
        tuple(
            torch.as_tensor(a, dtype=like.dtype, device=like.device)
            for a in utils.pd_table(n, nsigma)
        )
        for n, nsigma in layout
    )


def pd_inputs(
    model: types.ModuleType, packed: Tensor, tables: tuple
) -> tuple[tuple[Tensor, Tensor], ...]:
    npars, inputs = len(model.PARAMETERS), []
    for k, (z, w) in enumerate(tables):
        center, width = packed[:, k : k + 1], packed[:, npars + k : npars + k + 1]
        x = center + width * center * z
        lb, ub = model.LIMITS.get(model.PARAMETERS[k], (-math.inf, math.inf))
        ok = (x >= lb) & (x <= ub)
        inputs.append(
            (torch.where(ok, x, center), torch.where(ok, w, 0.0).expand(x.shape))
        )
    return tuple(inputs)


def pd_points(
    volume: Callable,
    radius_effective: Callable,
    valid: Callable,
    inputs: tuple,
    radius_effective_mode: int,
    cutoff: float,
) -> tuple[list[Tensor], Tensor, tuple[Tensor, ...]]:
    sets, npars = inputs[0][0].shape[0], len(inputs)
    shape = [sets, *(x.shape[1] for x, _ in inputs)]
    grid, weight = [], 1.0
    for k, (x, w) in enumerate(inputs):
        view = [sets] + [1] * npars
        view[k + 1] = x.shape[1]
        grid.append(x.reshape(view).expand(shape).reshape(sets, -1))
        weight = weight * w.reshape(view)
    weight = weight.expand(shape).reshape(sets, -1)
    center = weight.argmax(-1, keepdim=True)
    ok = valid(grid) & (weight > cutoff)
    grid = [torch.where(ok, g, g.gather(1, center)) for g in grid]
    weight = torch.where(ok, weight, 0.0)
    form, shell = volume(grid)
    radius = (
        radius_effective(radius_effective_mode, grid)
        if radius_effective_mode > 0
        else torch.zeros_like(weight)
    )
    return (
        grid,
        weight,
        (
            weight.sum(-1),
            (weight * form).sum(-1),
            (weight * shell).sum(-1),
            (weight * radius).sum(-1),
        ),
    )


def _Fq_sum(
    Fq: Callable, q: Tensor, p: list[Tensor], w: Tensor
) -> tuple[Tensor, Tensor]:
    F1, F2 = Fq(q, p)
    return (F1 * w).sum(-1), (F2 * w).sum(-1)


def orient(
    mode: int, qx: Tensor, qy: Tensor, angles: Tensor
) -> tuple[Tensor, Tensor, Tensor]:
    theta, phi, psi = (
        torch.deg2rad(a)[:, None]
        for a in torch.nn.functional.pad(angles, (0, 3 - angles.shape[1])).T
    )
    x, y = qx.reshape(1, -1), qy.reshape(1, -1)
    if mode == 0:
        qa = torch.sqrt(x * x + y * y).expand(len(theta), -1)
        return qa, torch.zeros_like(qa), torch.zeros_like(qa)
    st, ct, sp, cp = torch.sin(theta), torch.cos(theta), torch.sin(phi), torch.cos(phi)
    if mode == 1:
        qc = st * (cp * x + sp * y)
        return (
            torch.sqrt(torch.clamp(x * x + y * y - qc * qc, min=0.0)),
            qc,
            torch.zeros_like(qc),
        )
    ss, cs = torch.sin(psi), torch.cos(psi)
    return (
        (-sp * ss + cp * cs * ct) * x + (sp * cs * ct + ss * cp) * y,
        (-sp * cs - ss * cp * ct) * x + (-sp * ss * ct + cp * cs) * y,
        (st * cp) * x + (sp * st) * y,
    )


def _Iqxy_sum(
    Iqxy: Callable,
    mode: int,
    qx: Tensor,
    qy: Tensor,
    angles: Tensor,
    p: list[Tensor],
    w: Tensor,
) -> Tensor:
    qa, qb, qc = orient(mode, qx, qy, angles)
    return (Iqxy(qa[:, :, None], qb[:, :, None], qc[:, :, None], p) * w).sum(-1)


CHUNKS = dict(inputs=pd_inputs, points=pd_points, Fq=_Fq_sum, Iqxy=_Iqxy_sum)


def pd_Fq(
    Fq: Callable,
    volume: Callable,
    radius_effective: Callable,
    valid: Callable,
    q: Tensor,
    inputs: tuple,
    radius_effective_mode: int = 0,
    cutoff: float = 0.0,
) -> tuple[Tensor, Tensor, tuple[Tensor, ...]]:
    grid, weight, norms = CHUNKS["points"](
        volume, radius_effective, valid, inputs, radius_effective_mode, cutoff
    )
    F1, F2 = CHUNKS["Fq"](
        Fq, q[None, :, None], [g[:, None, :] for g in grid], weight[:, None, :]
    )
    return F1, F2, norms


def pd_Iqxy(
    Iqxy: Callable,
    volume: Callable,
    radius_effective: Callable,
    valid: Callable,
    mode: int,
    qx: Tensor,
    qy: Tensor,
    angles: Tensor,
    inputs: tuple,
    radius_effective_mode: int = 0,
    cutoff: float = 0.0,
) -> tuple[Tensor, tuple[Tensor, ...]]:
    grid, weight, norms = CHUNKS["points"](
        volume, radius_effective, valid, inputs, radius_effective_mode, cutoff
    )
    Iq = CHUNKS["Iqxy"](
        Iqxy, mode, qx, qy, angles, [g[:, None, :] for g in grid], weight[:, None, :]
    )
    return Iq, norms


def calculate(
    model: types.ModuleType,
    packed: Tensor,
    tables: tuple,
    cutoff: float = 0.0,
    q: Tensor | None = None,
    qx: Tensor | None = None,
    qy: Tensor | None = None,
) -> Tensor:
    m, npars = model, len(model.PARAMETERS)
    like = q if q is not None else qx
    step = max(1, MAX_BATCH // (like.numel() * math.prod(len(z) for z, _ in tables)))
    out = []
    for chunk in packed.split(step):
        inputs = CHUNKS["inputs"](m, chunk, tables)
        if q is not None:
            _, Iq, norms = pd_Fq(
                m.kernel_Fq,
                m.kernel_volume,
                m.kernel_radius_effective,
                m.kernel_valid,
                q,
                inputs,
                0,
                cutoff,
            )
        else:
            angles = chunk[:, 2 * npars : 2 * npars + len(m.ORIENTATION)]
            Iq, norms = pd_Iqxy(
                m.kernel_Iqxy,
                m.kernel_volume,
                m.kernel_radius_effective,
                m.kernel_valid,
                m.mode,
                qx,
                qy,
                angles,
                inputs,
                0,
                cutoff,
            )
        out.append(chunk[:, -2:-1] * Iq / norms[2][:, None] + chunk[:, -1:])
    return torch.cat(out).reshape(len(packed), *like.shape)


class Calculator:
    def __init__(
        self,
        model: types.ModuleType,
        pars: dict,
        cutoff: float = 0.0,
        graph: bool = True,
        **grid: Tensor,
    ):
        self.model, self.cutoff, self.grid = model, cutoff, grid
        self.like = grid.get("q", grid.get("qx"))
        utils.check_pars(model, pars, ("gaussian",))
        self.tables = pd_tables(utils.pd_layout(model, pars), self.like)
        self.graph = graph and self.like.device.type == "cuda"
        self._graph = self._static = self._out = None

    def compute(self, packed: Tensor) -> Tensor:
        return calculate(self.model, packed, self.tables, self.cutoff, **self.grid)

    def _capture(self, packed: Tensor) -> Tensor:
        stream = torch.cuda.Stream()
        stream.wait_stream(torch.cuda.current_stream())
        before = COUNT["integrations"]
        with torch.cuda.stream(stream):
            for _ in range(GRAPH_WARMUP):
                out = self.compute(packed)
        torch.cuda.current_stream().wait_stream(stream)
        self.graph = COUNT["integrations"] == before
        if self.graph:
            try:
                self._static, graph = packed.clone(), torch.cuda.CUDAGraph()
                with torch.cuda.graph(
                    graph, stream=stream, capture_error_mode="thread_local"
                ):
                    self._out = self.compute(self._static)
                self._graph = graph
            except Exception:
                self.graph = False
        return out

    def __call__(self, packed: Tensor) -> Tensor:
        packed = packed.to(self.like.device, self.like.dtype, non_blocking=True)
        if self.graph and self._graph is None:
            return self._capture(packed)
        if self._graph is None:
            return self.compute(packed)
        self._static.copy_(packed)
        self._graph.replay()
        return self._out.clone()


def compile_engine(device: torch.device) -> str:
    if (
        device.type == "cuda"
        and torch.cuda.get_device_capability(device) < COMPILE_MIN_CAPABILITY
    ):
        return f"eager ({torch.cuda.get_device_name(device)} cannot run Triton)"
    if SUM[1] is _sum_1d:
        import torch._dynamo.config as config

        for name in ("cache_size_limit", "recompile_limit"):
            if hasattr(config, name):
                setattr(config, name, max(getattr(config, name), 256))
        config.suppress_errors = True
        SUM.update({k: torch.compile(f, dynamic=None) for k, f in SUM.items()})
        CHUNKS.update({k: torch.compile(f, dynamic=False) for k, f in CHUNKS.items()})
    return "compiled"


def parallelism(device: torch.device) -> int:
    if device.type != "cuda":
        return torch.get_num_threads() * CPU_LANES
    props = torch.cuda.get_device_properties(device)
    return props.multi_processor_count * getattr(
        props, "max_threads_per_multi_processor", 2048
    )


def budget(device: torch.device) -> int:
    return (
        int(MEMORY_FRACTION * torch.cuda.mem_get_info(device)[0])
        if device.type == "cuda"
        else 1 << 34
    )


def tune(device: torch.device, dtype: torch.dtype) -> dict:
    global CHUNK_ELEMS, MAX_BATCH
    if device.type == "cuda":
        points = budget(device) // (POINT_ARRAYS * torch.finfo(dtype).bits // 8)
        CHUNK_ELEMS, MAX_BATCH = (
            min(WAVES_PER_LAUNCH * parallelism(device), points),
            points,
        )
    return dict(
        parallel_threads=parallelism(device),
        budget_gb=round(budget(device) / 2**30, 2),
        chunk_elems=CHUNK_ELEMS,
        max_batch=MAX_BATCH,
    )


def batch_size(
    model: types.ModuleType,
    pars: dict,
    cutoff: float = 0.0,
    factor: float = 1.0,
    **grid: Tensor,
) -> int:
    global COUNT_ONLY
    like = grid.get("q", grid.get("qx"))
    calc = Calculator(model, pars, cutoff, graph=False, **grid)
    COUNT["evaluations"], COUNT_ONLY = 0, True
    try:
        with getattr(torch.compiler, "set_stance", lambda _: nullcontext())(
            "force_eager"
        ):
            calc(torch.from_numpy(utils.pack(model, [pars])))
    finally:
        COUNT_ONLY = False
    work = max(
        like.numel() * math.prod(len(z) for z, _ in calc.tables), COUNT["evaluations"]
    )
    output = like.numel() * like.element_size() * OUTPUT_COPIES
    sets = 2 ** max(
        0, round(math.log2(max(1.0, WORK_PER_THREAD * parallelism(like.device) / work)))
    )
    return max(1, min(int(factor * min(MAX_SETS, sets)), budget(like.device) // output))


class Kernel(utils.SasKernel):
    def __init__(self, owner: "KernelModel", q_vectors: list[np.ndarray]):
        super().__init__(owner, q_vectors)
        self.grid = [torch.as_tensor(q, device=owner.device) for q in self.q]

    def sums(
        self,
        call_details: object,
        values: np.ndarray,
        cutoff: float,
        radius_effective_mode: int,
    ) -> tuple[np.ndarray, np.ndarray | None, np.ndarray]:
        m, like = self.owner.module, self.grid[0]
        v = torch.as_tensor(values, dtype=like.dtype).to(like.device, non_blocking=True)
        inputs, angles = utils.mesh_inputs(self.model, self.info, call_details, v)
        inputs = tuple((x.contiguous(), w.contiguous()) for x, w in inputs)
        kernels = (m.kernel_volume, m.kernel_radius_effective, m.kernel_valid)
        step = max(
            1, MAX_BATCH // (like.numel() * math.prod(x.shape[1] for x, _ in inputs))
        )
        F2, F1, norms = [], [], []
        for lo in range(0, len(v), step):
            pd = tuple((x[lo : lo + step], w[lo : lo + step]) for x, w in inputs)
            if self.dim == "1d":
                f1, f2, n = pd_Fq(
                    m.kernel_Fq, *kernels, like, pd, radius_effective_mode, cutoff
                )
                F1.append(f1)
            else:
                f2, n = pd_Iqxy(
                    m.kernel_Iqxy,
                    *kernels,
                    m.mode,
                    *self.grid,
                    angles[lo : lo + step].contiguous(),
                    pd,
                    radius_effective_mode,
                    cutoff,
                )
            F2.append(f2)
            norms.append(torch.stack(n, -1))

        def host(parts: list[Tensor]) -> np.ndarray:
            return torch.cat(parts).double().cpu().numpy()

        return host(F2), host(F1) if F1 else None, host(norms)


class KernelModel(utils.SasKernelModel):
    kernel_class = Kernel

    def __init__(self, name: str, dtype: str = "float64", device: str | None = None):
        super().__init__(name, dtype)
        self.device = torch.device(
            device or ("cuda" if torch.cuda.is_available() else "cpu")
        )
        compile_engine(self.device)
        if self.device.type == "cuda":
            tune(self.device, getattr(torch, self.dtype.name))
        self.module = load_model(name)
