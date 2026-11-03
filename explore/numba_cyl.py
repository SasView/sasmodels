from numba import njit, prange
import numpy as np
from numpy.typing import NDArray
from scipy.special import j1, p_roots

Vector = NDArray[np.float64]

_ADAPTIVE_MAX_N = 100000
ADAPTIVE_MAX_76 = (_ADAPTIVE_MAX_N / 76)
ADAPTIVE_MAX_OUTER = (_ADAPTIVE_MAX_N / 500)

def _build_weights(n):
    # Stack weights and nodes into a (2, N) array for Numba compatibility
    return np.vstack(p_roots(n))

# Pre‑compute weight tables as separate constants (Numba doesn't like the global dict)
GAUSS_WEIGHTS_20 = _build_weights(20)
GAUSS_WEIGHTS_76 = _build_weights(76)
GAUSS_WEIGHTS_500 = _build_weights(500)
GAUSS_WEIGHTS_5000 = _build_weights(5000)

# Turn off numba
#def njit(*args, **kw): return args[0] if args and callable(args[0]) else lambda f: f
#prange = range

@njit(inline="always", cache=False)
def gauss_weights(qr: float, n_outer: int) -> Vector:
    """
    Adaptive integration uses nested gaussian integration with order determined by the
    length scale qr, which drives the oscillation frequency of the scattering amplitude.
    Using polar coordinates, the oscillation with longitude appears to be independent of
    of latitude, so efficient spherical integration methods that sample roughly uniformly
    in solid angle don't perform as well as expected. The theta-phi uniform grid only
    uses 2x more points than the efficient grids, so the improvements from using
    Gauss-Legendre integration outweigh the increase in the number of integration points.

    The qr for the loop over phi may be different from that over theta. This is particularly
    the case for long rods, where the loop over phi depends on the cross section of the
    rod, so you can use far fewer phi points than the max dimension of the shape would
    suggest. Models such as the triaxial ellipsoid can relabel their axis so that c the
    long axis of the shape, allowing a smaller number of total evaluation points.

    If size of the grid is too large even after relabeling then the evaluation will be
    very slow. To prevent this, we pass the number of elements in the outer loop to
    limit the size of the inner loop, giving n_outer*n_inner < _ADAPTIVE_MAX_N. This works
    well enough for existing models, giving high accuracy results even for 20 micro shapes
    at q=1e-3 (the upper limit of qr for USANS). Even out to q=1e-1 the calculations are
    mostly within 10% of the full grid value (good enough to estimate slit resolution for
    USANS measurements) without being exceedingly slow::

        qr_outer = max dimension
        qr_inner = max cross section
        n_outer = gauss_weights(qr_outer, ADAPTIVE_MAX_OUTER, &w_outer, &z_outer);
        n_inner = gauss_weights(qr_inner, n_outer, &w_inner, &z_inner);

    For a standard integral over a hemisphere the loops will look like::

       const double *w_outer, *z_outer;
       const int n_outer = gauss_weights(qr_outer, 1, &w_outer, &z_outer);
       double total_outer = 0.;
       for (int j=0; j < n_outer; j++) {
           cos_theta = 0.5*z_outer[j] + 0.5; [-1, 1] => [0, 1]
           sin_theta = sqrt(1.0 - cos_theta*cos_theta)
           ...
           const double *w_inner, *z_inner;
           const int n_inner = gauss_weights(qr_inner, n_outer, &w_inner, &z_inner);
           double total_inner = 0.;
           for (int k=0; k < n_inner; k++) {
               phi = M_PI*(z_inner[k] + 1.0);  [-1, 1] => [0, 2 pi]
               ...
               total_inner += w_inner[k]*Fq;
           }
           total_outer += w_outer[k]*total_inner;
        }
        total_outer /= 2.0*M_PI; correct [-1, 1] => [0, 1] for u and [0, 2 pi] for phi

    Some models (barbell, capped cylinder, pringle) are not using spherical integration.
    Instead we find that limiting the outer integral to 76 points leads to okay results::

            n_outer = gauss_weights(qr_outer, ADAPTIVE_MAX_76, &w_outer, &z_outer);
    """
    # TODO: add more stages (150, 1500) so the n^2 model slowdown is a little less brutal
    # TODO: adjust cutoff for n=500
    #    For the pringle model I get npoints = 0.365 qr_max + 43
    #    and empirical values 14 => 20, 120 => 76, 1200 => 500, 14500 => 5000
    #    The error threshhold is abrupt: if n is too low the result is bad, but
    #    there is little benefit for having too large an n.
    #    These results are specific to the pringle inner integral and may not hold in general.
    #    Check that all models have reasonable cutoff.
    # *w = Gauss5000Wt; *z = Gauss5000Z; return 5000; # max precision
    # *w = Gauss76Wt; *z = Gauss76Z; return 76; # default

    if qr < 10 or n_outer > _ADAPTIVE_MAX_N // 76:
        return GAUSS_WEIGHTS_20
    elif qr < 100 or n_outer > _ADAPTIVE_MAX_N // 500:
        return GAUSS_WEIGHTS_76
    elif qr < 750 or n_outer > _ADAPTIVE_MAX_N // 5000:
        return GAUSS_WEIGHTS_500
    else:
        return GAUSS_WEIGHTS_5000

M_PI_4 = np.pi / 4.0

@njit(inline="always", cache=True)
def _sas_2J1x_x(x: float) -> float:
    return 2.0 * j1(x) / x if x != 0.0 else 1.0

@njit(inline="always", cache=True)
def _sas_sinx_x(x: float) -> float:
    return np.sin(x) / x if x != 0.0 else 1.0


# *** cylinder specific ***

@njit(inline="always", cache=True)
def _fq(qab:float , qc: float, radius: float, length: float) -> float:
    return _sas_2J1x_x(qab * radius) * _sas_sinx_x(qc * 0.5 * length)

@njit(inline="always", cache=True)
def _integrate(z_w: Vector, q: float, radius: float, length: float) -> tuple[float, float]:
    total_F1 = 0.0
    total_F2 = 0.0
    for i in range(z_w.shape[1]):
        value = z_w[0, i]
        weight = z_w[1, i]
        theta = value * M_PI_4 + M_PI_4
        sin_theta = np.sin(theta)
        cos_theta = np.cos(theta)
        form = _fq(q * sin_theta, q * cos_theta, radius, length)
        #if i == 10: print("integrate", i, z_w.shape[1], q, radius, length, theta, form)
        total_F1 += weight * form * sin_theta
        total_F2 += weight * form * form * sin_theta
    total_F1 *= M_PI_4
    total_F2 *= M_PI_4
    return total_F1, total_F2


@njit(inline="always", cache=True)
def form_volume(radius: float, length: float) -> float:
    """Cylinder volume – used for the contrast factor.
    ``radius`` and ``length`` are in the same units as the scattering vector.
    """
    return np.pi * radius * radius * length

@njit(cache=False)
def Fq(q: float, sld: float, solvent_sld: float, radius: float, length: float) -> tuple[float, float]:
    """Calculate the form factor amplitudes.

    Parameters
    ----------
    q : float
        Scattering vector magnitude.
    sld : float
        Scattering length density of the cylinder material.
    solvent_sld : float
        Scattering length density of the surrounding solvent.
    radius : float
        Cylinder radius.
    length : float
        Cylinder length.

    Returns
    -------
    (F1, F2) : tuple of floats
        Amplitudes corresponding to the original C implementation.
    """
    qr_max = q * max(radius, 0.5 * length)
    z_w = gauss_weights(qr=qr_max, n_outer=1)
    total_F1, total_F2 = _integrate(z_w, q, radius, length)
    s = (sld - solvent_sld) * form_volume(radius, length)
    F1 = 1e-2 * s * total_F1
    F2 = 1e-4 * s * s * total_F2
    return F1, F2

@njit(inline="always", cache=True)
def radius_from_excluded_volume(radius: float, length: float) -> float:
    """Effective radius based on excluded volume (no shell thickness)."""
    return 0.5 * np.cbrt(
        0.75 * radius * (2.0 * radius * length + (radius + length) * (np.pi * radius + length))
    )

@njit(inline="always", cache=True)
def radius_from_volume(radius: float, length: float) -> float:
    """Effective radius of a sphere having the same volume as the cylinder."""
    volume = np.pi * radius * radius * length
    return np.cbrt(volume / (4.0 * np.pi / 3.0))

@njit(inline="always", cache=True)
def radius_from_diagonal(radius: float, length: float) -> float:
    """Half the outer diagonal of the cylinder (radius, length)."""
    return np.sqrt(radius * radius + 0.25 * length * length)

@njit(cache=True)
def radius_effective(mode: int, radius: float, length: float) -> float:
    # Emulate the C switch‑case logic with Python if/elif.
    if mode == 1 or mode == 0:  # default / excluded‑volume
        return radius_from_excluded_volume(radius, length)
    elif mode == 2:
        return radius_from_volume(radius, length)
    elif mode == 3:
        return radius
    elif mode == 4:
        return 0.5 * length
    elif mode == 5:
        # half of the smaller of radius or length
        return radius if radius < 0.5 * length else 0.5 * length
    elif mode == 6:
        # half of the larger of radius or length
        return radius if radius > 0.5 * length else 0.5 * length
    elif mode == 7:
        return radius_from_diagonal(radius, length)
    else:
        # fallback to default case (same as case 1)
        return radius_from_excluded_volume(radius, length)

@njit(cache=True)
def valid(radius: float, length: float) -> bool:
    return radius >= 0.0 and length >= 0.0


# === cylinder model functions ===
@njit(cache=False)
def cylinder_Fq(q: float, p: Vector) -> tuple[float, float]:
    F1, F2 = Fq(q, p[0], p[1], p[2], p[3])
    return F1, F2

@njit(cache=True)
def cylinder_volume(p: Vector) -> tuple[float, float]:
    # form volume == shell volume
    v = form_volume(p[2], p[3])
    return v, v

@njit(cache=True)
def cylinder_radius_effective(mode: int, p: Vector) -> float:
    return radius_effective(mode, p[2], p[3])

@njit(cache=True)
def cylinder_valid(p: Vector) -> bool:
    return valid(p[2], p[3])

# *** generic polydispersity loop ***

@njit(parallel=True, cache=False)
def pd_Fq(
    Fq, volume, radius_effective, valid,
    q: Vector,
    inputs: tuple[tuple[Vector, Vector]],
    radius_effective_mode: int = 0,
    cutoff: float = 0.,
) -> tuple[Vector, Vector, tuple[float, float, float, float]]:
    npars, nq = len(inputs), len(q)
    F1, F2 = np.zeros(nq, dtype='d'), np.zeros(nq, dtype='d')
    # Compute strides: number of points for each polydisperse dimension
    stride = np.empty(npars, dtype=np.int64)
    length = np.empty(npars, dtype=np.int64)
    for k in range(npars):
        # each input is a tuple (values, weights)
        length[k] = inputs[k][0].shape[0]
        stride[k] = 1 if k == 0 else stride[k-1]*length[k-1]
    num_pd = stride[-1] * length[-1]
    #print("length", length, "stride", stride, "total", num_pd)
    pars = np.empty(npars, dtype='d')
    weight_norm = 0.
    weighted_form = 0.
    weighted_shell = 0.
    weighted_radius = 0.
    # TODO: need to chunk the calculation to allow interrupt. See for example the code in kerneldll.py
    # TODO: add magnetism
    # TODO: add support for 2D data and orientation dispersion
    #step = 100
    for pd_index in range(num_pd):
        weight = 1.
        for k in range(npars):
            par_index = (pd_index//stride[k])%length[k]
            pars[k] = inputs[k][0][par_index]
            weight *= inputs[k][1][par_index]
        if weight > cutoff and valid(pars):
            form, shell = volume(pars)
            weight_norm += weight
            weighted_form += weight * form
            weighted_shell += weight * shell
            if radius_effective_mode > 0:
                weighted_radius += weight * radius_effective(radius_effective_mode, pars)

            for q_index in prange(len(q)):
                F1_k, F2_k = Fq(q[q_index], pars)
                F1[q_index] += weight*F1_k
                F2[q_index] += weight*F2_k
        #else: print("cutoff", pars[2:])
        #if pd_index%step == step-1 or pd_index == num_pd-1:
        #    print("numba", pd_index, F2, [weight_norm, weighted_form, weighted_shell, weighted_radius])

    return F1, F2, (weight_norm, weighted_form, weighted_shell, weighted_radius)

def demo():
    import time
    from sasmodels import core, weights, direct_model, data

    q = np.logspace(-3, -1, 15000)
    dq = q * 0.0
    cutoff = 0

    # Input parameters
    n_length, n_radius = 31, 31
    sld = np.array([1.0]), np.array([1.0])
    sld_solvent = np.array([0.0]), np.array([1.0])
    radius = weights.get_weights("gaussian", n=n_radius, width=0.2, nsigmas=3,
                                 value=150, limits=[0, np.inf], relative=True)
    length = weights.get_weights("gaussian", n=n_length, width=0.2, nsigmas=3,
                                 value=15, limits=[0, np.inf], relative=True)

    # Warm‑up both implementations (first call includes compilation overhead)
    F1, F2, norms = pd_Fq(
        cylinder_Fq, cylinder_volume, cylinder_radius_effective, cylinder_valid,
        q=q, inputs=(sld, sld_solvent, radius, length),
        radius_effective_mode=0, cutoff=cutoff
    )
    # Compute reference result for sanity check
    weight_norm, weighted_form, weighted_shell, weighted_radius = norms
    Iq_numba = (F2 / weight_norm) / (weighted_shell / weight_norm)

    # Existing sasmodels cylinder for comparison (with warm-up).
    info = core.load_model_info("cylinder")
    kernel = core.build_model(info, dtype="double", platform="dll")
    #kernel = core.build_model(info, dtype="single", platform="dll")
    #kernel = core.build_model(info, dtype="single", platform="ocl")
    inp = data.Data1D(x=q, dx=dq)
    calculator = direct_model.DirectModel(inp, kernel, cutoff=cutoff)
    Iq_dll = calculator(radius=150, radius_pd=0.2, radius_pd_n=n_radius, radius_pd_nsigma=3,
                  length=15, length_pd=0.2, length_pd_n=n_length, length_pd_nsigma=3,
                  sld=sld[0][0], sld_solvent=sld_solvent[0][0],
                  scale=1.0, background=0.0)
    if 1: # timing
        # Benchmark loops – repeat a few times for stable timing
        repeats = 3
        # pd_Fq timing
        #time.sleep(5)
        start = time.perf_counter()
        for _ in range(repeats):
            F1, F2, norms = pd_Fq(cylinder_Fq, cylinder_volume,
                                cylinder_radius_effective, cylinder_valid,
                                q=q, inputs=(sld, sld_solvent, radius, length),
                                radius_effective_mode=0, cutoff=cutoff)
            # Compute reference result for sanity check
            weight_norm, weighted_form, weighted_shell, weighted_radius = norms
            Iq_numba = (F2 / weight_norm) / (weighted_shell / weight_norm)
        pd_time = (time.perf_counter() - start) / repeats
        #time.sleep(5)
        # Direct calculator timing (DLL backend)
        start = time.perf_counter()
        for _ in range(repeats):
            Iq_dll = calculator(radius=150, radius_pd=0.2, radius_pd_n=n_radius, radius_pd_nsigma=3,
                                length=15, length_pd=0.2, length_pd_n=n_length, length_pd_nsigma=3,
                                sld=sld[0][0], sld_solvent=sld_solvent[0][0],
                                scale=1.0, background=0.0)
        calc_time = (time.perf_counter() - start) / repeats

        # Report timings and a simple accuracy metric
        print(f"Average pd_Fq time   : {pd_time:.6f} s")
        print(f"Average calculator time: {calc_time:.6f} s")
        print("Speed‑up factor (calculator / pd_Fq):",
            f"{calc_time / pd_time:.2f}" if pd_time != 0 else "inf")

    #print(Iq_numba)
    #print(Iq_dll)
    print("Norm of difference (Numba vs DLL):", np.linalg.norm(Iq_numba - Iq_dll))

if __name__ == "__main__":
    demo()