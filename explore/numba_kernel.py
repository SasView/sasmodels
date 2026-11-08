import ast
import inspect
import sys
import types
import warnings

import numba
import numpy as np
from numba import njit, prange

import utils

warnings.filterwarnings("ignore", category=numba.NumbaExperimentalFeatureWarning)

FASTMATH: bool | set = {"nsz", "arcp", "contract", "afn", "reassoc"}
ADAPTIVE_MAX_N = 100000
ADAPTIVE_MAX_OUTER = ADAPTIVE_MAX_N // 500
IMAGE_MODES = dict(Iq=0, Iqac=1, Iqabc=2)


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


def fmin(a: float, b: float) -> float:
    return a if a < b else b


def fmax(a: float, b: float) -> float:
    return a if a > b else b


def where(condition: bool, a: float, b: float) -> float:
    return a if condition else b


def square(x: float) -> float:
    return x * x


def cube(x: float) -> float:
    return x * x * x


def sas_3j1x_x(x: float) -> float:
    small = abs(x) < SMALL_X
    s = 1.0 if small else x
    x2 = x * x
    series = 1.0 + x2 * (-3.0 / 30.0 + x2 * (3.0 / 840.0 + x2 * (-3.0 / 45360.0)))
    full = 3.0 * (np.sin(s) / s - np.cos(s)) / (s * s)
    return series if small else full


def gauss_weights(qr: float, n_outer: int) -> np.ndarray:
    if qr < 10.0 or n_outer > ADAPTIVE_MAX_N // 76:
        return GAUSS_WEIGHTS_20
    elif qr < 100.0 or n_outer > ADAPTIVE_MAX_N // 500:
        return GAUSS_WEIGHTS_76
    elif qr < 750.0 or n_outer > ADAPTIVE_MAX_N // 5000:
        return GAUSS_WEIGHTS_500
    return GAUSS_WEIGHTS_5000


def gauss_fixed(n: int) -> np.ndarray:
    if n == 20:
        return GAUSS_WEIGHTS_20
    elif n == 500:
        return GAUSS_WEIGHTS_500
    elif n == 5000:
        return GAUSS_WEIGHTS_5000
    return GAUSS_WEIGHTS_76


def _integrate(f, z_w: np.ndarray, lo: float, hi: float, *args) -> tuple[float, float]:
    half, center = 0.5 * (hi - lo), 0.5 * (hi + lo)
    total_F1, total_F2 = 0.0, 0.0
    for i in range(z_w.shape[1]):
        F1, F2 = f(half * z_w[0, i] + center, *args)
        total_F1 += z_w[1, i] * F1
        total_F2 += z_w[1, i] * F2
    return half * total_F1, half * total_F2


def integrate(f, qr: float, lo: float, hi: float, *args) -> tuple[float, float]:
    return _integrate(f, gauss_weights(qr, 1), lo, hi, *args)


def integrate_fixed(f, n: int, lo: float, hi: float, *args) -> tuple[float, float]:
    return _integrate(f, gauss_fixed(n), lo, hi, *args)


def integrate_2d(
    f,
    qr_outer: float,
    lo_outer: float,
    hi_outer: float,
    qr_inner: float,
    lo_inner: float,
    hi_inner: float,
    *args,
) -> tuple[float, float]:
    z_w_outer = gauss_weights(qr_outer, ADAPTIVE_MAX_OUTER)
    z_w_inner = gauss_weights(qr_inner, z_w_outer.shape[1])
    half_outer, center_outer = 0.5 * (hi_outer - lo_outer), 0.5 * (hi_outer + lo_outer)
    half_inner, center_inner = 0.5 * (hi_inner - lo_inner), 0.5 * (hi_inner + lo_inner)
    total_F1, total_F2 = 0.0, 0.0
    for i in range(z_w_outer.shape[1]):
        x = half_outer * z_w_outer[0, i] + center_outer
        inner_F1, inner_F2 = 0.0, 0.0
        for j in range(z_w_inner.shape[1]):
            F1, F2 = f(x, half_inner * z_w_inner[0, j] + center_inner, *args)
            inner_F1 += z_w_inner[1, j] * F1
            inner_F2 += z_w_inner[1, j] * F2
        total_F1 += z_w_outer[1, i] * inner_F1
        total_F2 += z_w_outer[1, i] * inner_F2
    scale = half_outer * half_inner
    return scale * total_F1, scale * total_F2


def pd_points(
    volume, radius_effective, valid, inputs, radius_effective_mode: int, cutoff: float
) -> tuple[np.ndarray, np.ndarray, tuple[float, float, float, float]]:
    npars = len(inputs)
    stride = np.empty(npars, dtype=np.int64)
    length = np.empty(npars, dtype=np.int64)
    for k in range(npars):
        length[k] = inputs[k][0].shape[0]
        stride[k] = 1 if k == 0 else stride[k - 1] * length[k - 1]
    num_pd = stride[-1] * length[-1]
    grid = np.empty((num_pd, npars), dtype=DTYPE)
    weight = np.empty(num_pd, dtype=DTYPE)
    weight_norm, weighted_form, weighted_shell, weighted_radius = 0.0, 0.0, 0.0, 0.0
    n = 0
    for pd_index in range(num_pd):
        w = 1.0
        for k in range(npars):
            par_index = (pd_index // stride[k]) % length[k]
            grid[n, k] = inputs[k][0][par_index]
            w *= inputs[k][1][par_index]
        if w > cutoff and valid(grid[n]):
            form, shell = volume(grid[n])
            weight_norm += w
            weighted_form += w * form
            weighted_shell += w * shell
            if radius_effective_mode > 0:
                weighted_radius += w * radius_effective(radius_effective_mode, grid[n])
            weight[n] = w
            n += 1
    return (
        grid[:n],
        weight[:n],
        (weight_norm, weighted_form, weighted_shell, weighted_radius),
    )


def pd_Fq(
    Fq,
    volume,
    radius_effective,
    valid,
    q: np.ndarray,
    inputs,
    radius_effective_mode: int = 0,
    cutoff: float = 0.0,
) -> tuple[np.ndarray, np.ndarray, tuple[float, float, float, float]]:
    grid, weight, norms = pd_points(
        volume, radius_effective, valid, inputs, radius_effective_mode, cutoff
    )
    nq, npd = len(q), len(weight)
    items, threads = nq * npd, numba.get_num_threads()
    block = (items + threads - 1) // threads
    F1_k, F2_k = np.empty(items, dtype=DTYPE), np.empty(items, dtype=DTYPE)
    for k in prange(block * threads):
        item = (k % block) * threads + k // block
        if item < items:
            F1_k[item], F2_k[item] = Fq(q[item // npd], grid[item % npd])
    F1, F2 = np.zeros(nq, dtype=DTYPE), np.zeros(nq, dtype=DTYPE)
    for q_index in range(nq):
        for pd_index in range(npd):
            F1[q_index] += weight[pd_index] * F1_k[q_index * npd + pd_index]
            F2[q_index] += weight[pd_index] * F2_k[q_index * npd + pd_index]
    return F1, F2, norms


def orient(
    mode: int,
    x: float,
    y: float,
    st: float,
    ct: float,
    sp: float,
    cp: float,
    ss: float,
    cs: float,
) -> tuple[float, float, float]:
    if mode == 0:
        return np.sqrt(x * x + y * y), 0.0, 0.0
    elif mode == 1:
        qc = st * (cp * x + sp * y)
        r = x * x + y * y - qc * qc
        return np.sqrt(r) if r > 0.0 else 0.0, qc, 0.0
    return (
        (-sp * ss + cp * cs * ct) * x + (sp * cs * ct + ss * cp) * y,
        (-sp * cs - ss * cp * ct) * x + (-sp * ss * ct + cp * cs) * y,
        (st * cp) * x + (sp * st) * y,
    )


def pd_Iqxy(
    Iqxy,
    volume,
    radius_effective,
    valid,
    mode: int,
    qx: np.ndarray,
    qy: np.ndarray,
    theta: float,
    phi: float,
    psi: float,
    inputs,
    radius_effective_mode: int = 0,
    cutoff: float = 0.0,
) -> tuple[np.ndarray, tuple[float, float, float, float]]:
    grid, weight, norms = pd_points(
        volume, radius_effective, valid, inputs, radius_effective_mode, cutoff
    )
    st, ct = np.sin(np.deg2rad(theta)), np.cos(np.deg2rad(theta))
    sp, cp = np.sin(np.deg2rad(phi)), np.cos(np.deg2rad(phi))
    ss, cs = np.sin(np.deg2rad(psi)), np.cos(np.deg2rad(psi))
    Iq = np.empty(len(qx), dtype=DTYPE)
    for i in prange(len(qx)):
        qa, qb, qc = orient(mode, qx[i], qy[i], st, ct, sp, cp, ss, cs)
        total = 0.0
        for pd_index in range(len(weight)):
            total += weight[pd_index] * Iqxy(qa, qb, qc, grid[pd_index])
        Iq[i] = total
    return Iq, norms


ENGINE = (
    fmin,
    fmax,
    where,
    square,
    cube,
    sas_3j1x_x,
    gauss_weights,
    gauss_fixed,
    _integrate,
    integrate,
    integrate_fixed,
    integrate_2d,
    pd_points,
    pd_Fq,
    orient,
    pd_Iqxy,
)
PARALLEL = ("pd_Fq", "pd_Iqxy")
_ENGINES: dict[str, dict] = {}


class _Single(ast.NodeTransformer):
    def visit_Constant(self, node: ast.Constant) -> ast.AST:
        if isinstance(node.value, float):
            return ast.copy_location(
                ast.Call(ast.Name("float32", ast.Load()), [node], []), node
            )
        return node


def _exec(source: str, filename: str, namespace: dict, single: bool) -> None:
    tree = ast.parse(source)
    if single:
        tree.body = [
            _Single().visit(n) if isinstance(n, ast.FunctionDef) else n
            for n in tree.body
        ]
        ast.fix_missing_locations(tree)
    exec(compile(tree, filename, "exec"), namespace)


def engine(dtype: str) -> dict:
    if dtype not in _ENGINES:
        dt = np.dtype(dtype).type
        single = dt is np.float32
        ns = dict(
            np=np,
            numba=numba,
            prange=prange,
            float32=np.float32,
            DTYPE=dt,
            ADAPTIVE_MAX_N=ADAPTIVE_MAX_N,
            ADAPTIVE_MAX_OUTER=ADAPTIVE_MAX_OUTER,
            SMALL_X=dt(0.7 if single else 0.1),
            M_4PI_3=dt(4.0 * np.pi / 3.0),
            M_PI_2=dt(np.pi / 2.0),
            sqrt=np.sqrt,
            sin=np.sin,
            cos=np.cos,
            log=np.log,
            asin=np.arcsin,
            cbrt=np.cbrt,
        )
        for n in (20, 76, 500, 5000):
            ns[f"GAUSS_WEIGHTS_{n}"] = np.vstack(roots_legendre(n)).astype(dt)
        for fn in ENGINE:
            _exec(inspect.getsource(fn), __file__, ns, single)
            ns[fn.__name__] = njit(
                parallel=fn.__name__ in PARALLEL, fastmath=FASTMATH, error_model="numpy"
            )(ns[fn.__name__])
        ns["kernel"] = njit(fastmath=FASTMATH, error_model="numpy")
        _ENGINES[dtype] = ns
    return _ENGINES[dtype]


def load_model(name: str, dtype: str = "float64") -> types.ModuleType:
    model = utils.MODELS[name]
    module = sys.modules[f"numba_{name}_{dtype}"] = types.ModuleType(
        f"numba_{name}_{dtype}"
    )
    module.__dict__.update(engine(dtype), **utils.model_metadata(model))
    for fn, source in utils.model_functions(model):
        _exec(source, f"{name}.{fn}", module.__dict__, np.dtype(dtype) == np.float32)
        module.__dict__[fn] = module.kernel(module.__dict__[fn])
    p = ", ".join(f"p[{i}]" for i in range(len(module.PARAMETERS)))
    sizes = ", ".join(
        f"p[{module.PARAMETERS.index(n)}]" for n in module.VOLUME_PARAMETERS
    )
    Iqxy = dict(
        Iq=f"Fq(qa, {p})[1]", Iqac=f"Iqac(qa, qb, {p})", Iqabc=f"Iqabc(qa, qb, qc, {p})"
    )[module.IMAGE]
    exec(
        f"def kernel_Fq(q, p):\n    return Fq(q, {p})\n"
        f"def kernel_Iqxy(qa, qb, qc, p):\n    return {Iqxy}\n"
        f"def kernel_volume(p):\n    v = form_volume({sizes})\n    return v, v\n"
        f"def kernel_radius_effective(mode, p):\n    return radius_effective(mode, {sizes})\n"
        f"def kernel_valid(p):\n    return valid({sizes})\n",
        module.__dict__,
    )
    for fn in (
        "kernel_Fq",
        "kernel_Iqxy",
        "kernel_volume",
        "kernel_radius_effective",
        "kernel_valid",
    ):
        module.__dict__[fn] = module.kernel(module.__dict__[fn])
    module.dtype, module.mode = np.dtype(dtype), IMAGE_MODES[module.IMAGE]
    return module


def calculate(
    model: types.ModuleType,
    pars: dict,
    cutoff: float = 0.0,
    q: np.ndarray | None = None,
    qx: np.ndarray | None = None,
    qy: np.ndarray | None = None,
) -> np.ndarray:
    inputs = utils.pd_inputs(model, pars, model.dtype)
    scale, background = (
        pars.get(k, model.DEFAULTS[k]) for k in ("scale", "background")
    )
    kernels = (model.kernel_volume, model.kernel_radius_effective, model.kernel_valid)
    if q is not None:
        F1, F2, norms = model.pd_Fq(model.kernel_Fq, *kernels, q, inputs, 0, cutoff)
        return scale * F2 / norms[2] + background
    angles = [
        model.dtype.type(pars.get(a, model.DEFAULTS[a])) for a in model.ORIENTATION
    ]
    angles += [model.dtype.type(0.0)] * (3 - len(angles))
    Iq, norms = model.pd_Iqxy(
        model.kernel_Iqxy,
        *kernels,
        model.mode,
        qx.reshape(-1),
        qy.reshape(-1),
        *angles,
        inputs,
        0,
        cutoff,
    )
    return (scale * Iq / norms[2] + background).reshape(qx.shape)


class Kernel(utils.SasKernel):
    def sums(
        self,
        call_details: object,
        values: np.ndarray,
        cutoff: float,
        radius_effective_mode: int,
    ) -> tuple[np.ndarray, np.ndarray | None, np.ndarray]:
        m = self.owner.module
        inputs, angles = utils.mesh_inputs(
            self.model, self.info, call_details, values.astype(m.dtype, copy=False)
        )
        kernels = (m.kernel_volume, m.kernel_radius_effective, m.kernel_valid)
        F2, F1, norms = [], [], []
        for s in range(len(values)):
            pd = tuple(
                (np.ascontiguousarray(x[s]), np.ascontiguousarray(w[s]))
                for x, w in inputs
            )
            if self.dim == "1d":
                f1, f2, n = m.pd_Fq(
                    m.kernel_Fq,
                    *kernels,
                    self.q[0],
                    pd,
                    radius_effective_mode,
                    float(cutoff),
                )
                F1.append(f1)
            else:
                a = [m.dtype.type(x) for x in angles[s]] + [m.dtype.type(0.0)] * (
                    3 - angles.shape[1]
                )
                f2, n = m.pd_Iqxy(
                    m.kernel_Iqxy,
                    *kernels,
                    m.mode,
                    self.q[0],
                    self.q[1],
                    *a,
                    pd,
                    radius_effective_mode,
                    float(cutoff),
                )
            F2.append(f2)
            norms.append(n)
        return (
            np.array(F2),
            np.array(F1) if F1 else None,
            np.array(norms, dtype=np.float64),
        )


class KernelModel(utils.SasKernelModel):
    kernel_class = Kernel

    def __init__(self, name: str, dtype: str = "float64"):
        super().__init__(name, dtype)
        self.module = load_model(name, self.dtype.name)
