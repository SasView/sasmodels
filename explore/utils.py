import inspect
import os
import textwrap
import types
from collections.abc import Callable, Iterator
from contextlib import contextmanager
from pathlib import Path

import numpy as np

if not hasattr(np, "NaN"):
    np.NaN = np.nan
if not hasattr(np, "trapz"):
    np.trapz = np.trapezoid

import sasmodels
from sasmodels import core, data, direct_model, weights
from sasmodels.details import make_kernel_args
from sasmodels.kernel import Kernel, KernelModel

HERE = Path(__file__).resolve().parent


class Skip(Exception):
    pass


class Sphere:
    MODEL_NAME = "sphere"

    def form_volume(radius: float) -> float:
        return M_4PI_3 * cube(radius)

    def radius_effective(mode: int, radius: float) -> float:
        return radius

    def Fq(
        q: float, sld: float, sld_solvent: float, radius: float
    ) -> tuple[float, float]:
        bes = sas_3j1x_x(q * radius)
        contrast = sld - sld_solvent
        form = contrast * form_volume(radius) * bes
        return 1.0e-2 * form, 1.0e-4 * form * form

    def valid(radius: float) -> bool:
        return radius >= 0.0


class Ellipsoid:
    MODEL_NAME = "ellipsoid"
    SASMODELS_1_1_POINTS = True

    def form_volume(radius_polar: float, radius_equatorial: float) -> float:
        return M_4PI_3 * radius_polar * radius_equatorial * radius_equatorial

    def radius_from_volume(radius_polar: float, radius_equatorial: float) -> float:
        return cbrt(radius_polar * radius_equatorial * radius_equatorial)

    def radius_from_curvature(radius_polar: float, radius_equatorial: float) -> float:
        equal = radius_polar == radius_equatorial
        zero = radius_polar * radius_equatorial == 0.0
        ratio = where(
            equal | zero,
            0.5,
            where(
                radius_polar < radius_equatorial,
                radius_polar / radius_equatorial,
                radius_equatorial / radius_polar,
            ),
        )
        e1 = sqrt(1.0 - ratio * ratio)
        b1 = 1.0 + asin(e1) / (e1 * ratio)
        bL = (1.0 + e1) / (1.0 - e1)
        b2 = 1.0 + 0.5 * ratio * ratio / e1 * log(bL)
        delta = 0.75 * b1 * b2
        ddd = 2.0 * (delta + 1.0) * radius_polar * radius_equatorial * radius_equatorial
        return where(equal, radius_polar, where(zero, 0.0, 0.5 * cbrt(ddd)))

    def radius_effective(
        mode: int, radius_polar: float, radius_equatorial: float
    ) -> float:
        if mode == 2:
            return radius_from_volume(radius_polar, radius_equatorial)
        if mode == 3:
            return fmin(radius_polar, radius_equatorial)
        if mode == 4:
            return fmax(radius_polar, radius_equatorial)
        return radius_from_curvature(radius_polar, radius_equatorial)

    def _Fq_point(
        u: float, q: float, radius_equatorial: float, v_square_minus_one: float
    ) -> tuple[float, float]:
        r = radius_equatorial * sqrt(1.0 + u * u * v_square_minus_one)
        f = sas_3j1x_x(q * r)
        return f, f * f

    def Fq(
        q: float,
        sld: float,
        sld_solvent: float,
        radius_polar: float,
        radius_equatorial: float,
    ) -> tuple[float, float]:
        v_square_minus_one = square(radius_polar / radius_equatorial) - 1.0
        if SASMODELS_1_1_POINTS:
            qr_max = q * fmax(radius_polar, radius_equatorial)
            total_F1, total_F2 = integrate(
                _Fq_point, qr_max, 0.0, 1.0, q, radius_equatorial, v_square_minus_one
            )
        else:
            total_F1, total_F2 = integrate_fixed(
                _Fq_point, 76, 0.0, 1.0, q, radius_equatorial, v_square_minus_one
            )
        s = (sld - sld_solvent) * form_volume(radius_polar, radius_equatorial)
        return 1.0e-2 * s * total_F1, 1.0e-4 * s * s * total_F2

    def Iqac(
        qab: float,
        qc: float,
        sld: float,
        sld_solvent: float,
        radius_polar: float,
        radius_equatorial: float,
    ) -> float:
        qr = sqrt(square(radius_equatorial * qab) + square(radius_polar * qc))
        f = sas_3j1x_x(qr)
        s = (sld - sld_solvent) * form_volume(radius_polar, radius_equatorial)
        return 1.0e-4 * square(f * s)

    def valid(radius_polar: float, radius_equatorial: float) -> bool:
        return (radius_polar >= 0.0) & (radius_equatorial >= 0.0)


class TriaxialEllipsoid:
    MODEL_NAME = "triaxial_ellipsoid"
    SASMODELS_1_1_POINTS = True

    def form_volume(
        radius_equat_minor: float, radius_equat_major: float, radius_polar: float
    ) -> float:
        return M_4PI_3 * radius_equat_minor * radius_equat_major * radius_polar

    def radius_from_curvature(
        radius_equat_minor: float, radius_equat_major: float, radius_polar: float
    ) -> float:
        lo, hi = fmin(radius_equat_minor, radius_equat_major), fmax(
            radius_equat_minor, radius_equat_major
        )
        radius_1, radius_3 = fmin(lo, radius_polar), fmax(hi, radius_polar)
        radius_2 = fmax(lo, fmin(hi, radius_polar))
        upper = radius_2 - radius_1 > radius_3 - radius_2
        r_equat_equiv = where(
            upper, sqrt(radius_2 * radius_3), sqrt(radius_1 * radius_2)
        )
        r_polar_equiv = where(upper, radius_1, radius_3)
        equal = radius_1 == radius_3
        zero = radius_equat_minor * radius_equat_major * radius_polar == 0.0
        ratio = where(
            equal | zero,
            0.5,
            where(
                r_polar_equiv < r_equat_equiv,
                r_polar_equiv / r_equat_equiv,
                r_equat_equiv / r_polar_equiv,
            ),
        )
        e1 = sqrt(1.0 - ratio * ratio)
        b1 = 1.0 + asin(e1) / (e1 * ratio)
        bL = (1.0 + e1) / (1.0 - e1)
        b2 = 1.0 + 0.5 * ratio * ratio / e1 * log(bL)
        delta = 0.75 * b1 * b2
        ddd = 2.0 * (delta + 1.0) * r_polar_equiv * r_equat_equiv * r_equat_equiv
        return where(equal, radius_polar, where(zero, 0.0, 0.5 * cbrt(ddd)))

    def radius_effective(
        mode: int,
        radius_equat_minor: float,
        radius_equat_major: float,
        radius_polar: float,
    ) -> float:
        if mode == 2:
            return cbrt(radius_equat_minor * radius_equat_major * radius_polar)
        if mode == 3:
            return fmin(fmin(radius_equat_minor, radius_equat_major), radius_polar)
        if mode == 4:
            return fmax(fmax(radius_equat_minor, radius_equat_major), radius_polar)
        return radius_from_curvature(
            radius_equat_minor, radius_equat_major, radius_polar
        )

    def _Fq_point(
        phi: float, u: float, q: float, radius_equat_major: float, pa: float, pc: float
    ) -> tuple[float, float]:
        usq = u * u
        r = radius_equat_major * sqrt(
            pa * square(sin(phi)) * (1.0 - usq) + 1.0 + pc * usq
        )
        fq = sas_3j1x_x(q * r)
        return fq, fq * fq

    def _Fq_point_u_phi(
        u: float, phi: float, q: float, radius_equat_major: float, pa: float, pc: float
    ) -> tuple[float, float]:
        return _Fq_point(phi, u, q, radius_equat_major, pa, pc)

    def Fq(
        q: float, sld: float, sld_solvent: float, Ra: float, Rb: float, Rc: float
    ) -> tuple[float, float]:
        maybe_min, maybe_max = fmin(Ra, Rb), fmax(Ra, Rb)
        maybe_mid = fmax(maybe_min, Rc)
        radius_equat_minor = fmin(maybe_min, Rc)
        radius_equat_major = fmin(maybe_max, maybe_mid)
        radius_polar = fmax(maybe_max, maybe_mid)
        pa = square(radius_equat_minor / radius_equat_major) - 1.0
        pc = square(radius_polar / radius_equat_major) - 1.0
        qr_max_inner = fmax(q * radius_equat_minor, q * radius_equat_major)
        qr_max = fmax(qr_max_inner, q * radius_polar)
        if SASMODELS_1_1_POINTS:
            outer_sum_F1, outer_sum_F2 = integrate_2d(
                _Fq_point,
                qr_max,
                0.0,
                M_PI_2,
                qr_max_inner,
                0.0,
                1.0,
                q,
                radius_equat_major,
                pa,
                pc,
            )
        else:
            outer_sum_F1, outer_sum_F2 = integrate_2d(
                _Fq_point_u_phi,
                qr_max,
                0.0,
                1.0,
                qr_max_inner,
                0.0,
                M_PI_2,
                q,
                radius_equat_major,
                pa,
                pc,
            )
        volume = form_volume(radius_equat_minor, radius_equat_major, radius_polar)
        contrast = sld - sld_solvent
        return (
            1.0e-2 * contrast * volume * outer_sum_F1 / M_PI_2,
            1.0e-4 * square(contrast * volume) * outer_sum_F2 / M_PI_2,
        )

    def Iqabc(
        qa: float,
        qb: float,
        qc: float,
        sld: float,
        sld_solvent: float,
        radius_equat_minor: float,
        radius_equat_major: float,
        radius_polar: float,
    ) -> float:
        qr = sqrt(
            square(radius_equat_minor * qa)
            + square(radius_equat_major * qb)
            + square(radius_polar * qc)
        )
        fq = sas_3j1x_x(qr)
        vol = form_volume(radius_equat_minor, radius_equat_major, radius_polar)
        drho = sld - sld_solvent
        return 1.0e-4 * square(vol * drho * fq)

    def valid(
        radius_equat_minor: float, radius_equat_major: float, radius_polar: float
    ) -> bool:
        return (
            (radius_equat_minor >= 0.0)
            & (radius_equat_major >= 0.0)
            & (radius_polar >= 0.0)
        )


def from_sasmodels(model: type) -> type:
    table = core.load_model_info(model.MODEL_NAME).parameters
    kernel = table.kernel_parameters
    model.PARAMETERS = tuple(p.name for p in kernel if p.type != "orientation")
    model.VOLUME_PARAMETERS = tuple(p.name for p in kernel if p.type == "volume")
    model.ORIENTATION = tuple(p.name for p in kernel if p.type == "orientation")
    model.DEFAULTS = {
        p.name: float(p.default) for p in (*table.common_parameters, *kernel)
    }
    model.LIMITS = {p.name: (float(p.limits[0]), float(p.limits[1])) for p in kernel}
    model.IMAGE = ("Iq", "Iq", "Iqac", "Iqabc")[len(model.ORIENTATION)]
    return model


MODELS = {
    m.MODEL_NAME: from_sasmodels(m) for m in (Sphere, Ellipsoid, TriaxialEllipsoid)
}


def model_metadata(model: type) -> dict:
    return {
        k: v
        for k, v in vars(model).items()
        if not k.startswith("__") and not inspect.isfunction(v)
    }


def model_functions(model: type) -> list[tuple[str, str]]:
    return [
        (k, textwrap.dedent(inspect.getsource(v)))
        for k, v in vars(model).items()
        if inspect.isfunction(v)
    ]


def check_pars(
    model: type, pars: dict, pd_types: tuple[str, ...] | None = None
) -> None:
    pd_types = pd_types or tuple(k for k in weights.DISTRIBUTIONS if k != "array")
    for name in (*model.VOLUME_PARAMETERS, *model.ORIENTATION):
        if not (
            float(pars.get(name + "_pd", 0.0)) and int(pars.get(name + "_pd_n", 0)) > 1
        ):
            continue
        if name in model.ORIENTATION:
            raise ValueError(f"orientation dispersity ({name}_pd) is not supported")
        kind = pars.get(name + "_pd_type", "gaussian")
        if kind not in pd_types:
            raise ValueError(
                f"{name}_pd_type={kind!r} is not supported here, use one of {', '.join(pd_types)}"
            )


def pd_inputs(
    model: types.ModuleType, pars: dict, dtype: str = "float64"
) -> tuple[tuple[np.ndarray, np.ndarray], ...]:
    check_pars(model, pars)
    inputs = []
    for name in model.PARAMETERS:
        width = (
            float(pars.get(name + "_pd", 0.0))
            if name in model.VOLUME_PARAMETERS
            else 0.0
        )
        kind = pars.get(name + "_pd_type", "gaussian") if width else "gaussian"
        x, w = weights.get_weights(
            kind,
            int(pars.get(name + "_pd_n", 0)),
            width,
            float(pars.get(name + "_pd_nsigma", 3.0)),
            float(pars.get(name, model.DEFAULTS[name])),
            model.LIMITS.get(name, (-np.inf, np.inf)),
            True,
        )
        inputs.append(
            (np.ascontiguousarray(x, dtype=dtype), np.ascontiguousarray(w, dtype=dtype))
        )
    return tuple(inputs)


def pd_table(n: int, nsigma: float) -> tuple[np.ndarray, np.ndarray]:
    return weights.get_weights(
        "gaussian", n, 1.0, nsigma, 0.0, (-np.inf, np.inf), False
    )


def pd_layout(model: types.ModuleType, pars: dict) -> tuple[tuple[int, float], ...]:
    out = []
    for name in model.PARAMETERS:
        width = (
            float(pars.get(name + "_pd", 0.0))
            if name in model.VOLUME_PARAMETERS
            else 0.0
        )
        n = int(pars.get(name + "_pd_n", 0)) if width else 1
        out.append(
            (n, float(pars.get(name + "_pd_nsigma", 3.0))) if n > 1 else (1, 0.0)
        )
    return tuple(out)


def pack(
    model: types.ModuleType,
    pars_list: list[dict],
    pd_types: tuple[str, ...] = ("gaussian",),
) -> np.ndarray:
    def row(p: dict) -> list[float]:
        check_pars(model, p, pd_types)
        values = [float(p.get(n, model.DEFAULTS[n])) for n in model.PARAMETERS]
        widths = [
            float(p.get(n + "_pd", 0.0)) if n in model.VOLUME_PARAMETERS else 0.0
            for n in model.PARAMETERS
        ]
        angles = [float(p.get(a, model.DEFAULTS[a])) for a in model.ORIENTATION]
        return (
            values
            + widths
            + angles
            + [float(p.get(k, model.DEFAULTS[k])) for k in ("scale", "background")]
        )

    return np.asarray([row(p) for p in pars_list], dtype=np.float64)


def batch_sets(
    model: types.ModuleType, pars: dict, sets: int, spread: float, seed: int
) -> list[dict]:
    rng = np.random.default_rng(seed)
    names = (*model.VOLUME_PARAMETERS, *model.ORIENTATION)
    return [dict(pars)] + [
        dict(
            pars,
            **{
                n: float(pars.get(n, model.DEFAULTS[n]))
                * rng.uniform(1 - spread, 1 + spread)
                for n in names
            },
        )
        for _ in range(sets - 1)
    ]


def image_grid(q_bins: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    return tuple(np.meshgrid(q_bins, q_bins))


def errors(
    Iq: np.ndarray, reference: np.ndarray, index: np.ndarray | None
) -> tuple[float, float]:
    Iq = np.asarray(Iq, dtype=np.float64).reshape(-1)
    Iq = Iq if index is None or len(Iq) == len(reference) else Iq[index]
    norm = float(np.linalg.norm(Iq - reference))
    return norm, norm / float(np.linalg.norm(reference))


@contextmanager
def single_core() -> Iterator[None]:
    if not hasattr(os, "sched_setaffinity"):
        yield
        return
    cores = os.sched_getaffinity(0)
    os.sched_setaffinity(0, {min(cores)})
    try:
        yield
    finally:
        os.sched_setaffinity(0, cores)


@contextmanager
def openmp_build() -> Iterator[None]:
    from sasmodels import kerneldll

    if not hasattr(kerneldll, "compiler"):
        raise Skip(
            "this sasmodels version has no compiler flag list for an OpenMP DLL build"
        )
    path, compiler = kerneldll.SAS_DLL_PATH, list(kerneldll.compiler)
    kerneldll.SAS_DLL_PATH = os.path.join(path, "openmp")
    os.makedirs(kerneldll.SAS_DLL_PATH, exist_ok=True)
    kerneldll.compiler.append("-fopenmp")
    try:
        yield
    finally:
        kerneldll.SAS_DLL_PATH, kerneldll.compiler[:] = path, compiler


def opencl_device(kind: str) -> str:
    try:
        import pyopencl as cl
    except ImportError as error:
        raise Skip("pyopencl is not installed") from error
    wanted = cl.device_type.GPU if kind == "gpu" else cl.device_type.CPU
    for i, platform in enumerate(cl.get_platforms()):
        for j, device in enumerate(platform.get_devices()):
            if device.type & wanted:
                from sasmodels import kernelcl

                os.environ["SAS_OPENCL"] = f"{i}:{j}"
                kernelcl.reset_environment()
                return device.name.strip()
    raise Skip(f"no OpenCL {kind.upper()} device")


def sasmodels_control(
    model: type,
    kind: str,
    q: np.ndarray,
    q_bins: np.ndarray,
    cutoff: float,
    dtype: str,
    platform: str,
) -> tuple[Callable[[dict], np.ndarray], np.ndarray | None, dict]:
    kernel = core.build_model(
        core.load_model_info(model.MODEL_NAME), dtype=dtype, platform=platform
    )
    if platform == "ocl":
        from sasmodels import kernelcl

        if not isinstance(kernel, kernelcl.GpuModel):
            raise Skip(
                f"sasmodels can not run {dtype} on this OpenCL device and fell back to the DLL"
            )
    call, index = direct(kernel, kind, q, q_bins, cutoff)
    dll = getattr(kernel, "dllpath", None)
    info = dict(
        sasmodels=sasmodels.__version__,
        openmp=(
            b"GOMP_" in Path(dll).read_bytes() if dll and os.path.exists(dll) else None
        ),
    )
    return call, index, info


def direct(
    kernel_model: object,
    kind: str,
    q: np.ndarray,
    q_bins: np.ndarray,
    cutoff: float = 0.0,
) -> tuple[Callable[[dict], np.ndarray], np.ndarray | None]:
    dataset = (
        data.Data1D(x=q, dx=0.0 * q)
        if kind == "Fq"
        else data.empty_data2D(q_bins, q_bins)
    )
    calculator = direct_model.DirectModel(dataset, kernel_model, cutoff=cutoff)
    index = None if kind == "Fq" else np.asarray(calculator.index, dtype=bool)
    q_calc = calculator.resolution.q_calc

    def call(pars: dict) -> np.ndarray:
        return np.asarray(calculator(**pars), dtype=np.float64)

    call.q_vectors = list(q_calc) if isinstance(q_calc, (list, tuple)) else [q_calc]
    return call, index


def mesh_inputs(
    model: type, info: object, call_details: object, values: object
) -> tuple[tuple, object]:
    names = [p.name for p in info.parameters.kernel_parameters]
    base, num_weights = info.parameters.nvalues, int(call_details.num_weights)
    for name in model.ORIENTATION:
        if call_details.length[names.index(name)] > 1:
            raise NotImplementedError(
                f"orientation dispersity ({name}_pd) is not supported"
            )

    def column(k: int, shift: int) -> object:
        lo = base + shift + int(call_details.offset[k])
        return values[:, lo : lo + int(call_details.length[k])]

    inputs = tuple(
        (column(k, 0), column(k, num_weights))
        for k in map(names.index, model.PARAMETERS)
    )
    return inputs, values[:, [2 + names.index(a) for a in model.ORIENTATION]]


class SasKernelModel(KernelModel):
    kernel_class: type

    def __init__(self, name: str, dtype: str = "float64"):
        self.info, self.dtype, self.model = (
            core.load_model_info(name),
            np.dtype(dtype),
            MODELS[name],
        )

    def make_kernel(self, q_vectors: list[np.ndarray]) -> "SasKernel":
        return self.kernel_class(self, q_vectors)

    def release(self) -> None:
        pass


class SasKernel(Kernel):
    def __init__(self, owner: SasKernelModel, q_vectors: list[np.ndarray]):
        self.owner, self.info, self.dtype, self.model = (
            owner,
            owner.info,
            owner.dtype,
            owner.model,
        )
        self.dim = "2d" if len(q_vectors) == 2 else "1d"
        self.q = [np.ascontiguousarray(v, dtype=self.dtype) for v in q_vectors]
        self.q_input = types.SimpleNamespace(nq=len(self.q[0]))
        self.nout = 2 if self.info.have_Fq and self.dim == "1d" else 1
        self.result = np.empty(self.nout * self.q_input.nq + 4)

    def sums(
        self,
        call_details: object,
        values: np.ndarray,
        cutoff: float,
        radius_effective_mode: int,
    ) -> tuple[np.ndarray, np.ndarray | None, np.ndarray]:
        raise NotImplementedError

    def _call_kernel(
        self,
        call_details: object,
        values: np.ndarray,
        cutoff: float,
        magnetic: bool,
        radius_effective_mode: int,
    ) -> None:
        if magnetic:
            raise NotImplementedError("magnetism is not supported")
        F2, F1, norms = self.sums(
            call_details, np.atleast_2d(values), cutoff, radius_effective_mode
        )
        nq, nout = self.q_input.nq, self.nout
        result = np.empty((len(F2), nout * nq + 4))
        result[:, 0 : nout * nq : nout] = F2
        if nout == 2:
            result[:, 1 : nout * nq : nout] = F1
        result[:, nout * nq :] = norms
        self.result = result if np.ndim(values) == 2 else result[0]

    def Fq(
        self,
        call_details: object,
        values: np.ndarray,
        cutoff: float,
        magnetic: bool,
        radius_effective_mode: int = 0,
    ) -> tuple:
        self._call_kernel(call_details, values, cutoff, magnetic, radius_effective_mode)
        r, nq, nout = np.atleast_2d(self.result), self.q_input.nq, self.nout
        total = r[:, nout * nq]
        total = np.where(total == 0.0, 1.0, total)
        form, shell, radius = (r[:, nout * nq + k] / total for k in (1, 2, 3))
        shell = np.where(shell == 0.0, 1.0, shell)
        F1 = r[:, 1 : nout * nq : nout] / total[:, None] if nout == 2 else None
        F2 = r[:, 0 : nout * nq : nout] / total[:, None]
        out = (F1, F2, radius, shell, form / shell)
        return (
            out
            if np.ndim(values) == 2
            else tuple(None if x is None else x[0] for x in out)
        )

    def Iq(
        self, call_details: object, values: np.ndarray, cutoff: float, magnetic: bool
    ) -> np.ndarray:
        _, F2, _, shell, _ = self.Fq(call_details, values, cutoff, magnetic, 0)
        v = np.asarray(values, dtype=np.float64)
        return (v[..., 0] / shell)[..., None] * F2 + v[..., 1][..., None]

    __call__ = Iq


def batch_args(
    kernel: SasKernel, pars_list: list[dict], mono: bool = False
) -> tuple[object, np.ndarray, bool]:
    args = [
        make_kernel_args(
            kernel,
            direct_model.get_mesh(kernel.info, dict(p), dim=kernel.dim, mono=mono),
        )
        for p in pars_list
    ]
    details = args[0][0]
    if any(not np.array_equal(d.length, details.length) for d, _, _ in args):
        raise ValueError(
            "all parameter sets in a batch need the same number of dispersity points"
        )
    return details, np.stack([v for _, v, _ in args]), any(m for _, _, m in args)
