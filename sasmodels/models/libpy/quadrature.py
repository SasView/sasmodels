import numpy as np
import numpy.typing as npt
import typing as tp
import h5py
import os

module_dir = os.path.dirname(__file__)
tuned_quad_model_loc = os.path.join(module_dir, "tuned_quad.h5")

FixedQuad = tp.NamedTuple("FixedQuad", [("points", npt.NDArray[np.float64]), ("weights", npt.NDArray[np.float64])])


def fixed_quad_init(n, kronrod_points):
    if n in kronrod_points:
        x = kronrod_points[n][0][1:-1:2].copy().astype(np.float64)  # numba requires contiguous arrays
        w = kronrod_points[n][2][1:-1:2].copy().astype(np.float64)
        return FixedQuad(x, w)
    else:
        raise ValueError(f"n = {n} not found in Kronrod points.")


def fixed_quad_integrate(
        fixed_quad: FixedQuad,
        func: tp.Callable[
            [npt.NDArray[np.float64], np.float64, np.float64, tp.Dict[str, np.float64]], npt.NDArray[np.float64]],
        a: np.float64,
        b: np.float64,
        params: tp.Dict[str, np.float64] = None) -> tp.Union[np.float64]:
    if params is None:
        params = {}
    y = (b - a) * (fixed_quad.points + 1) / 2.0 + a
    return (b - a) / 2.0 * np.sum(fixed_quad.weights * func(y, params), axis=-1)


TunedQuad = tp.NamedTuple("TunedQuad",
                          [("reg_params", tp.Dict[str, npt.NDArray[np.float64]]), ("tuned_mat", npt.NDArray[np.int32]),
                           ("quad_cache", tp.Dict[int, FixedQuad]), ("dims", np.ndarray)])


def tuned_quad_init(
        reg_params: tp.Dict[str, npt.NDArray[np.float64]],
        tuned_mat: npt.NDArray[np.int32],
        kronrod_points_dict) -> TunedQuad:
    quad_cache = {}
    for n in kronrod_points_dict:
        integrator = fixed_quad_init(n, kronrod_points_dict)
        quad_cache[n] = integrator
    # Get the dimensions of the tuned matrix
    dims = np.empty(len(reg_params), dtype=np.int32)
    for i, value in enumerate(reg_params.values()):
        dims[i] = len(value)
    return TunedQuad(reg_params, tuned_mat, quad_cache, dims)


def tuned_quad_ravel_multi_index(
        tuned_quad: TunedQuad,
        multi_index: npt.NDArray[np.int32]) -> int:
    strides = np.cumprod(np.concatenate((np.array([1], dtype=np.int32), tuned_quad.dims[::-1][:-1].astype(np.int32))))[
              ::-1]
    return (multi_index * strides).sum()


def tuned_quad_get_n_kronrod(
        tuned_quad: TunedQuad,
        params: tp.Dict[str, np.float64]) -> np.int32:
    index = np.empty(len(tuned_quad.reg_params), dtype=np.int32)

    for k, key in enumerate(tuned_quad.reg_params):
        i = np.searchsorted(tuned_quad.reg_params[key], params[key])

        i = min(i, len(tuned_quad.reg_params[key]) - 1)
        index[k] = i

    loc = tuned_quad_ravel_multi_index(tuned_quad, index)

    if loc >= len(tuned_quad.tuned_mat):
        raise ValueError("Index out of bounds")

    return tuned_quad.tuned_mat[loc]


def tuned_quad_integrate(
        tuned_quad: TunedQuad,
        func: tp.Callable[
            [tp.Union[np.float64, npt.NDArray[np.float64]], np.float64, np.float64, tp.Dict[str, np.float64]], tp.Union[
                np.float64, npt.NDArray[np.float64]]],
        a: np.float64,
        b: np.float64,
        params: tp.Dict[str, np.float64]) -> tp.Union[np.float64]:
    tuned_quad_check_params(tuned_quad, params)
    n = tuned_quad_get_n_kronrod(tuned_quad, params)
    integrator = tuned_quad.quad_cache[n]
    return fixed_quad_integrate(integrator, func, a, b, params)


def tuned_quad_check_params(
        tuned_quad: TunedQuad,
        params: tp.Dict[str, np.float64]
):
    for key in params:
        if key not in tuned_quad.reg_params:
            raise ValueError(f"Parameter {key} is not registered")


def load_tuned_quad_h5(
        model_name: str,
        kronrod_points_dict: tp.Dict[int, tp.Tuple[np.ndarray, np.ndarray, np.ndarray]]

) -> TunedQuad:
    with h5py.File(tuned_quad_model_loc, "r") as f:
        if model_name not in f:
            raise FileNotFoundError(f"{model_name} not found.")  # Check if the model exists

        file = f[model_name]
        reg_params: tp.Dict[str, npt.NDArray[np.float64]] = {}
        params_order = file["params_order"][:]
        for k in params_order:
            reg_params[str(k.decode('ascii'))] = np.array(file["reg_params"][k][:])
        tuned_mat = np.array(file["tuned_mat"][:])

        tuned_quad = tuned_quad_init(reg_params, tuned_mat, kronrod_points_dict)

    return tuned_quad
