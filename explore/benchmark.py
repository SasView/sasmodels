from pathlib import Path

import plot
import timing

if __name__ == "__main__":
    PATH = Path(__file__).resolve().parent
    PRECISSONS = ("float64", "float32")
    THREADS = 12
    REPEATS = 5
    WARMUP = 2
    CUTOFF = 0.0
    SPREAD = 0.2
    BATCH_FACTORS = (1, 4)
    WORK_PER_THREAD = 4096
    MEMORY_FRACTION = 0.5
    TIMELINE_BINS = 200
    Q_1D = ("logspace", -3.0, -1.0, 1000)
    Q_BINS_2D = ("linspace", -0.1, 0.1, 129)
    RESULTS_FILE = PATH / "results.json"
    OUTPUT_DIR = PATH / "plots"
    PD_N, PD_NSIGMA = 35, 3  # same gauss in sv gui
    SPHERE = dict(
        sld=4.0,
        sld_solvent=1.0,
        scale=1.0,
        background=0.0,
        radius=400.0,
        radius_pd=0.1,
        radius_pd_n=PD_N,
        radius_pd_nsigma=PD_NSIGMA,
    )
    ELLIPSOID = dict(
        sld=4.0,
        sld_solvent=1.0,
        scale=1.0,
        background=0.0,
        radius_polar=400.0,
        radius_polar_pd=0.1,
        radius_polar_pd_n=PD_N,
        radius_polar_pd_nsigma=PD_NSIGMA,
        radius_equatorial=20.0,
        theta=60.0,
        phi=60.0,
    )
    TRIAXIAL = dict(
        sld=4.0,
        sld_solvent=1.0,
        scale=1.0,
        background=0.0,
        radius_equat_minor=20.0,
        radius_equat_major=400.0,
        radius_equat_major_pd=0.1,
        radius_equat_major_pd_n=PD_N,
        radius_equat_major_pd_nsigma=PD_NSIGMA,
        radius_polar=10.0,
        radius_polar_pd=0.1,
        radius_polar_pd_n=PD_N,
        radius_polar_pd_nsigma=PD_NSIGMA,
        theta=60.0,
        phi=60.0,
        psi=60.0,
    )
    CASES = dict(
        sphere_Fq=dict(model="sphere", kind="Fq", pars=SPHERE),
        sphere_Iq=dict(model="sphere", kind="image", pars=SPHERE),
        ellipsoid_Fq=dict(model="ellipsoid", kind="Fq", pars=ELLIPSOID),
        ellipsoid_Iqac=dict(model="ellipsoid", kind="image", pars=ELLIPSOID),
        triaxial_Fq=dict(model="triaxial_ellipsoid", kind="Fq", pars=TRIAXIAL),
        triaxial_Iqabc=dict(model="triaxial_ellipsoid", kind="image", pars=TRIAXIAL),
    )
    results = timing.run(
        dict(
            cases=CASES,
            precisions=PRECISSONS,
            threads=THREADS,
            repeats=REPEATS,
            warmup=WARMUP,
            cutoff=CUTOFF,
            spread=SPREAD,
            batch_factors=BATCH_FACTORS,
            work_per_thread=WORK_PER_THREAD,
            memory_fraction=MEMORY_FRACTION,
            timeline_bins=TIMELINE_BINS,
            q_1d=Q_1D,
            q_bins_2d=Q_BINS_2D,
            results_file=RESULTS_FILE,
        )
    )
    plot.plot(results, OUTPUT_DIR)
