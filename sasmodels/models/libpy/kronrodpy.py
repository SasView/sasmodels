import numpy as np
import typing as tp
import h5py
import os



# KronrodPointsDictType = (types.int32, types.UniTuple(types.float64[:], 3))
# KronrodPointsType = numba.extending.as_numba_type(types.DictType(*KronrodPointsDictType))
# kronrod_points_dict = Dict.empty(*KronrodPointsDictType)
kronrod_points_dict = {}
# kronrod_loc_h5 = "C:/Users/dnwob/GitHub/SasViewDev/sasmodels/sasmodels/models/libpy/kronrod.h5"
module_dir = os.path.dirname(__file__)
kronrod_loc_h5 = os.path.join(module_dir, "kronrod.h5")
""" Kronrod quadrature points. These ahve been precomputed and stored in `kronrod.h5`.
    If you want to generate more points, there is a cpp file that can be used to generate the points.
    store them in `kronrod.h5` so they are loaded here.
"""


def get_kronrod_h5(n: int) -> tp.Tuple[np.ndarray, np.ndarray, np.ndarray]:
    with h5py.File(kronrod_loc_h5, "r") as file:
        if f"{n}" not in file:
            raise FileNotFoundError(f"{n} not found.")
        x = np.array(file[f"{n}/x"][:])
        w = np.array(file[f"{n}/w"][:])
        wg = np.array(file[f"{n}/wg"][:])

    x = np.concatenate([-x[:-1], x[::-1]])
    w = np.concatenate([w[:-1], w[::-1]])
    wg = np.concatenate([wg[:-1], wg[::-1]])
    return x, w, wg


with h5py.File(kronrod_loc_h5, "r") as f:
    all_kron_p = np.array(sorted([int(k) for k in list(f.keys())]))

for n in all_kron_p:
    x, w, wg = get_kronrod_h5(n)
    kronrod_points_dict[n] = (x, w, wg)
