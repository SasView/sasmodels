"""
Check the speed and accuracy of the 1-D model integration routines,
comparing them to a fixed 76-point gaussian.

You can look at the overall performance by fitting the latex_sphere
sasdata example with a triaxial ellipsoid. This can be done in sasview
for a real world check, but it includes the overhead of the sasview shim.
You can skip over the shim by running a fit directly in bumps:

    $ bumps example/simul_fit.py --fit=dream --pop=-100

Modify simul_fit.py to choose SANS, USANS or both for
triaxial_ellipsoid, ellipsoid or sphere.

Sasmodels compare help:

    $ python -m sasmodels.compare -?

Environment variables are listed in the sasmodels.compare help. For example,
to control use of the GPU in the execution engine:

    SAS_OPENCL=vendor:device|cuda:device|none sets the target GPU device

Check speed and accuracy via explore/check_adaptive.py (this file).

There are multiple control parameters for running this program, but you need
to change the code to set them. Search for "runtime control flags" below.

I've checked a couple of models using big_n=15000 by listing them on the command line:

    $ python explore/check_adaptive.py capped_cylinder barbell

Target sizes are 20 μm diameter for big (USANS) and 20 nm diameter for small (SANS/SAXS).
Small models are tested over the SANS q range, 1e-3/Ang to 1/Ang. USANS only tests to 0.1/Ang

For each model the generic a,b,c dimensions are translated to model specific parameters.

There are other model parameters such as "sld_core" for core-shell models, which is set to
sld_solvent=0, and the thickness/core_size ratio "shell" which defaults to 0.1. Maybe a
future version could take these parameters on the command line. The latter is relevant
to USANS since large objects with thin shells have a 1/q rather than 1/q^4 falloff, and
so slit resolution values will be higher.

Three different aspect ratios are tested: rod-like, disk-like and cube-like

You can check random parameter sets for a model using:

    $ python -m sasmodels.compare background=0 -ngauss=0,1000 -sets=5 -double triaxial_ellipsoid

This compares adaptive (ngauss=0) to a 1000x1000 fixed grid.

If you find a dataset that is anomalous then look for the -random=n output on the console
and add it to the command line. He will reproduce a model with the same parameters.

You may want something like "-random=n -ngauss=0,12000 -nq=5 q=0.001:0.1" to look at
five q values over two decades with 12000x12000 grid

Really small models are faster on the CPU. Large models are faster on the GPU.
"""

from math import sqrt

from sasmodels import core, data, generate
from sasmodels.compare import tic
from sasmodels.direct_model import DirectModel

# == Unnested ==

MODELS = {}
KERNELS = {}
NESTED = {}
UNNESTED = {}
def register(loops):
    def wrapper(par_fn):
        model_name = par_fn.__name__
        if loops == 1:
            UNNESTED[model_name] = par_fn
        else:
            NESTED[model_name] = par_fn
        return par_fn
    return wrapper

@register(1)
def core_shell_bicelle(ab, c, shell=0.1, sld_core=0):
    thick_rim = thick_face = ab*shell/2 if shell < 1 else shell
    length = c - 2*thick_face
    radius = ab/2 - thick_rim
    pars = dict(
        radius=radius, thick_rim=thick_rim, thick_face=thick_face, length=length,
        sld_core=sld_core, sld_face=1, sld_rim=1, sld_solvent=0,
    )
    return pars

@register(1)
def core_shell_cylinder(ab, c, shell=0.1, sld_core=0):
    thickness = ab*shell/2 if shell < 1 else shell
    length = c - 2*thickness
    radius = ab/2 - thickness
    pars = dict(
        radius=radius, thickness=thickness, length=length,
        sld_core=sld_core, sld_shell=1, sld_solvent=0,
    )
    return pars

@register(1)
def core_shell_ellipsoid(ab, c, shell=0.1, sld_core=0):
    thick_shell = ab*shell/2 if shell < 1 else shell
    radius_equat_core = ab/2 - thick_shell
    x_polar_shell = 1
    x_core = (c/2 - thick_shell)/(ab/2 - thick_shell)
    pars = dict(
        radius_equat_core=radius_equat_core, x_core=x_core,
        thick_shell=thick_shell, x_polar_shell=x_polar_shell,
        sld_core=sld_core, sld_shell=1, sld_solvent=0,
    )
    return pars

@register(1)
def cylinder(ab, c):
    pars = dict(sld=1, sld_solvent=0, radius=ab/2, length=c)
    return pars

@register(1)
def ellipsoid(ab, c):
    pars = dict(sld=1, sld_solvent=0, radius_equatorial=ab/2, radius_polar=c/2)
    return pars

@register(1)
def flexible_cylinder(ab, c, n=10):
    pars = dict(
        length=c*n, kuhn_length=c,
        radius=ab,
        sld=1, sld_solvent=0,
    )
    return pars

@register(1)
def hollow_cylinder(ab, c, shell=0.1):
    thickness = ab/2*shell
    pars = dict(
        radius=ab/2 - thickness,
        thickness=thickness,
        length=c,
        sld=1, sld_solvent=0,
    )
    return pars

@register(1)
def stacked_disks(ab, c, shell=0.1, n=10, sld_core=0):
    # Each stack (layer+core+layer) has total thickness length/n
    # The shell is a portion of that total thickness, with the core being (1-portion)
    # Therefore thick_core = length/n*(1-shell) and thick_layer = length/n*shell/2
    jitter = 0.05
    pars = dict(
        thick_core=c/n*(1-shell),
        thick_layer=c/n*shell/2,
        radius=ab/2,
        n_stacking=n,
        sigma_d=jitter*(c/n),
        sld_core=sld_core,
        sld_layer=1,
        sld_solvent=0,
    )
    return pars


# == Nested ==

@register(2)
def barbell(a, b, c):
    radius, radius_bell = a/2, b/2
    h = sqrt(radius_bell**2 - radius**2)
    # c/2 = length/2 + radius_bell + h
    length = max(c - 2*(radius_bell + h), 0)  # "disk" fails because radius_bell > c
    pars = dict(radius=radius, radius_bell=radius_bell, length=length, sld=1, sld_solvent=0)
    return pars

@register(2)
def capped_cylinder(a, b, c):
    radius, radius_cap = a/2, b/2
    h = sqrt(radius_cap**2 - radius**2)
    # c/2 = length/2 + radius_bell + h
    length = max(c - 2*(radius_cap - h), 0)
    pars = dict(radius=radius, radius_cap=radius_cap, length=length, sld=1, sld_solvent=0)
    return pars

@register(2)
def core_shell_bicelle_elliptical_belt_rough(a, b, c, shell=0.1, sld_core=0, roughness=0.02):
    thick_rim = thick_face = a*shell/2 if shell < 1 else shell
    length = c - 2*thick_face
    radius = a/2 - thick_rim
    x_core = (b/2 - thick_rim) / (a/2 - thick_rim)
    pars = dict(
        radius=radius, x_core=x_core, thick_rim=thick_rim, thick_face=thick_face, length=length,
        sld_core=sld_core, sld_face=1, sld_rim=1, sld_solvent=0,
        sigma=thick_rim*roughness,
    )
    return pars

@register(2)
def core_shell_bicelle_elliptical(a, b, c, shell=0.1, sld_core=0):
    thick_rim = thick_face = a*shell/2 if shell < 1 else shell
    length = c - 2*thick_face
    radius = a/2 - thick_rim
    x_core = (b/2 - thick_rim) / (a/2 - thick_rim)
    pars = dict(
        radius=radius, x_core=x_core, thick_rim=thick_rim, thick_face=thick_face, length=length,
        sld_core=sld_core, sld_face=1, sld_rim=1, sld_solvent=0,
    )
    return pars

@register(2)
def core_shell_parallelepiped(a, b, c, shell=0.1, sld_core=0):
    thickness = a*shell/2 if shell < 1 else shell
    length_a, length_b, length_c = a - 2*thickness, b - 2*thickness, c - 2*thickness
    pars = dict(
        sld_core=sld_core, sld_a=1, sld_b=1, sld_c=1, sld_solvent=0,
        length_a=length_a, length_b=length_b, length_c=length_c,
        thick_rim_a=thickness, thick_rim_b=thickness, thick_rim_c=thickness,
    )
    return pars

@register(2)
def elliptical_cylinder(a, b, c):
    pars = dict(radius_minor=a/2, axis_ratio=b/c, length=c, sld=1, sld_solvent=0)
    return pars

@register(2)
def flexible_cylinder_elliptical(a, b, c, n=10):
    pars = dict(
        length=c*n, kuhn_length=c,
        radius=a, axis_ratio=b/a,
        sld=1, sld_solvent=0,
    )
    return pars

@register(2)
def hollow_rectangular_prism_thin_walls(a, b, c):
    pars = dict(
        sld=1, sld_solvent=0,
        length_a=a, b2a_ratio=b/a, c2a_ratio=c/a,
    )
    return pars

@register(2)
def hollow_rectangular_prism(a, b, c, shell=0.1):
    thickness = a/2*shell
    pars = dict(
        sld=1, sld_solvent=0,
        length_a=a, b2a_ratio=b/a, c2a_ratio=c/a,
        thickness=thickness,
    )
    return pars

@register(2)
def parallelepiped(a, b, c):
    pars = dict(
        sld=1, sld_solvent=0,
        length_a=a, length_b=b, length_c=c,
    )
    return pars

@register(2)
def pringle(a, b, c):
    # It is nonsensical to have thickness >> radius, so c ought not be the "long" direction,
    # however, this is how the loops are structured in the pringle model, so go with it.
    # Let a be the radius of the disk under curvature (alpha==beta), and b be the uncurved radius.
    alpha = beta = max(abs(1 - a/b), 0.001)
    pars = dict(radius=b/2, thickness=c, alpha=alpha, beta=beta, sld=1, sld_solvent=0)
    return pars

@register(2)
def rectangular_prism(a, b, c):
    pars = dict(sld=1, sld_solvent=0, length_a=a, b2a_ratio=b/a, c2a_ratio=c/a)
    return pars

@register(2)
def triaxial_ellipsoid(a, b, c):
    pars = dict(
        sld=1, sld_solvent=0,
        radius_equat_minor=a, radius_equat_major=b, radius_polar=c,
    )
    return pars


def check(models=[]):
    import numpy as np

    def longlines(fn):
        def wrapped(*args, **kw):
            with np.printoptions(
                linewidth=1000,
                formatter={'float': '{:.14e}'.format},
            ):
                return fn(*args, **kw)
        return wrapped

    def get_kernel(model_name, n_gauss=0, dtype="double", platform="dll"):
        # print(f"{model_name=} {n_gauss=} {dtype=} {platform=}")
        if model_name not in MODELS:
            model_info = core.load_model_info(model_name)
            MODELS[model_name] = model_info
            # set_integration_size() below overwrites info. Cache the adaptive version immediately
            KERNELS[(model_name, 0, dtype, platform)] = core.build_model(model_info, dtype=dtype, platform=platform)
        if (model_name, n_gauss, dtype, platform) not in KERNELS:
            model_info = MODELS[model_name]
            if n_gauss:
                # ocl adaptive... this had better be the first or second get_kernel since set_integration_size is stateful
                generate.set_integration_size(model_info, n_gauss)
            KERNELS[(model_name, n_gauss, dtype, platform)] = core.build_model(model_info, dtype=dtype, platform=platform)
        return KERNELS[(model_name, n_gauss, dtype, platform)]

    def timer(fn):
        # Clear startup overhead
        toc = tic()
        fn()
        dt = toc()
        if dt > 1:
            return dt
        # loop
        toc = tic()
        fn()
        dt = toc()
        n = int(0.5 // dt)
        for _ in range(n): fn()
        return toc()/(n+1)

    @longlines
    def run_test(model_name, pars, speed_q, test_n_gauss, test_q, tol):

        # Parameter string suitable for use with "python -m sasmodels.compare model par=value..."
        par_str = " ".join(f"{k}={v}" for k, v in pars.items())
        # print(f"{model_name} {par_str} -ngauss={test_n_gauss}")
        k_adaptive = get_kernel(model_name, n_gauss=0, dtype="double", platform="dll")
        speed_data = data.empty_data1D(speed_q)
        test_data = data.empty_data1D(test_q)

        Iq_test_adaptive = DirectModel(test_data, k_adaptive)(**pars) # clear startup overhead
        dt_adaptive = timer(lambda: DirectModel(speed_data, k_adaptive)(**pars))

        labelled = False
        def print_label():
            nonlocal labelled
            if not labelled:
                print(f"{model_name} {par_str}")
                labelled = True

        if dt_adaptive > speed_target:
            print_label()
            print(f"! ** {model_name} is slow: {dt_adaptive:.1f} s for {len(speed_q)} points in [{speed_q.min()}, {speed_q.max()}]")

        if speed_check: # compare speed of adaptive with
            if speed_check == "gpu":
                k_speed = get_kernel(model_name, n_gauss=0, dtype="single", platform="ocl")
            elif speed_check == "gauss-76":
                k_speed = get_kernel(model_name, n_gauss=76)
            else:
                raise ValueError(f"Unknown speed check: {speed_check}. Use gpu or gauss-76")
            dt_speed = timer(lambda: DirectModel(speed_data, k_speed)(**pars))

            delta = dt_speed/dt_adaptive - 1
            if delta >= 1 or delta <= -0.5 or dt_adaptive > 0.5:
                print_label()
                value = (
                    f"{dt_speed/dt_adaptive:.1f}x slower" if delta >= 1
                    else f"{dt_adaptive/dt_speed:.1f}x faster" if delta < -0.5
                    else f"{int(delta*100)}% slower" if delta > 0.0
                    else f"{int(-delta*100)}% faster"
                )
                if delta > 0.2:
                    comment = f"  [{speed_check} {value}]" + (" ***" if speed_check != "gauss-76" else "")
                elif delta < -0.2:
                    comment = f"  [{speed_check} {value}]" + (" ***" if speed_check == "gauss-76" else "")
                elif dt_adaptive > 0.2:
                    comment = "  [slow]"
                else:
                    comment = ""
                print(f"  cpu double: {dt_adaptive*1000:.1f} ms  {speed_check}: {dt_speed*1000:.1f} ms{comment}")

        if speed_only: # speed_only test
            return

        k_accurate = get_kernel(model_name, n_gauss=test_n_gauss, dtype="double", platform="dll")
        Iq_test_accurate = DirectModel(test_data, k_accurate)(**pars)
        if not np.isfinite(Iq_test_accurate).any() or (Iq_test_accurate == 0).any():
            print(f">>> Bad target values for\n    {model_name} {par_str}\n    q   : {test_q}\n    I(q): {Iq_test_accurate}")
        if not np.isfinite(Iq_test_adaptive).any() or (Iq_test_adaptive == 0).any():
            print(f">>> Bad target values for\n    {model_name} {par_str}\n    q   : {test_q}\n    I(q): {Iq_test_adaptive}")
        relerr = abs(Iq_test_adaptive - Iq_test_accurate)/Iq_test_accurate
        index = relerr > tol
        if index.any():
            print_label()
            print("  qvalue", test_q[index])
            print("  target", Iq_test_accurate[index])
            print("  actual", Iq_test_adaptive[index])
            print("  relerr", relerr[index])

        # Compare against 76-point gaussian
        k_76 = get_kernel(model_name, n_gauss=76, dtype="double", platform="dll")
        Iq_test_76 = DirectModel(test_data, k_76)(**pars)
        if not np.isfinite(Iq_test_76).any() or (Iq_test_76 == 0).any():
            print(f">>> Bad target values for\n    {model_name} {par_str}\n    q   : {test_q}\n    I(q): {Iq_test_76}")
        relerr_76 = abs(Iq_test_76 - Iq_test_accurate)/Iq_test_accurate
        index = (relerr_76 > tol) #& (relerr > tol) & (relerr > 10*relerr_76)
        if index.any():
            print_label()
            #print("! ** gauss-76 is better than adaptive for {model_name} at some q values")
            print("! ** gauss-76 is out of tolerance")
            print("  qvalue", test_q[index])
            print("  relerr", relerr[index])
            print("  rel-76", relerr_76[index])


    # === runtime control flags ===
    # speed_target=n warns if adaptive is n times slower than gauss-76
    # speed_only only compares the calculation speed, not the accuracy
    # speed_check is what to compare adaptive against:
    #     "gauss-76" fixed grid,
    #     "gpu" adaptive running on the GPU, or
    #     None for no speed check
    # big_n is the grid size for large shapes. Time scales as n²
    #     for NESTED models, so only a few q values are tested,
    #     ranging over USANS values 5e-4 to 1e-1. The high q is
    #     large so we can look at the cost of computing resolution
    #     effects in slit geometry
    # small_n is the grid size for small shapes.
    # n_per_decade=n is the number of q points per decade unsed for
    #     small shapes. With three decades, n=1 means four q points,
    #     which is fewer than the number used for big_n

    speed_target = 2
    speed_only = False
    speed_check = "gauss-76"  # None | "" | "gpu" | "gauss-76"
    #speed_only = True
    #speed_check = None
    big_n = 5000
    #small_n = 1000
    small_n = big_n
    n_per_decade = 1
    print(f"""\
== Speed and accuracy tests for all adaptive integration models ==
* target evaluation time is {speed_target} s (running on a mac M2 chip)
*     q in [1e-5, 1] with 40 points per decade for 201 points total
* large models tested against {big_n} point gaussian integration
*     q=[5e-4, 1e-3, 2e-3] with tol=1e-5 relative (measured q)
*     q=[0.01, 0.1] with tol=0.2 relative (slit resolution limits)
* small models tested against {small_n} point gaussian integration
*     q in [1e-3, 1] with {n_per_decade} points per decade
!!!! These tests run very slowly --- don't use as part of CI !!!!
""")
    def loop(aspect, a, b, c):
        q = np.logspace(-5, 0, 201) # 5 decades, 40 points per decade
        big = max(a, b, c) > 10000  # 1 μm
        if big:
            n_gauss = big_n
            test_q = np.array([0.0005, 0.001, 0.002, 0.01, 0.1])
            test_tol = np.array([5e-5, 5e-5, 5e-5, 2e-1, 2e-1])
        else:
            n_gauss = small_n
            test_q = np.logspace(-3, 0, 3*n_per_decade + 1) # 3 decades, n points per decade
            test_tol = np.full_like(test_q, 1e-10)

        print(f"\n\n=== {'big' if big else 'small'} {aspect}: {a=} {b=} {c=} ===")
        ab = max(a, b)
        for name, fn in UNNESTED.items():  # 1D integrals
            # skip models not listed
            if models and name not in models:
                continue
            pars = dict(background=0, **fn(ab=ab, c=c))
            run_test(name, pars, q, n_gauss, test_q, test_tol)

        ab = max(a, b)
        for name, fn in NESTED.items():  # 2D integrals
            # skip models not listed
            if models and name not in models:
                continue
            pars = dict(background=0, **fn(a=a, b=b, c=c))
            run_test(name, pars, q, n_gauss, test_q, test_tol)

    # small size (20 nm) tests
    loop("rods", a=20, b=40, c=200)
    loop("disks", a=180, b=200, c=40)
    loop("cubes", a=200, b=200, c=200)

    # large size (20 μm) tests
    loop("rods", a=1_000, b=2000, c=200_000)
    loop("disks", a=180_000, b=200_000, c=1_000)
    loop("cubes", a=200_000, b=200_000, c=200_000)

if __name__ == "__main__":
    import sys
    check(models=sys.argv[1:])
