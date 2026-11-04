"""
Unit tests for multi_shell_vesicle.

The single-shell case (n=1) must be numerically identical to the
vesicle model with matching parameters:

    vesicle                    multi_shell_vesicle
    -------------------------------- ------------------------
    radius                  ->  radius
    thickness               ->  thickness1
    sld                     ->  sld1
    sld_solvent             ->  sld_solvent
    volfraction             ->  volfraction
    (core = solvent by def) ->  (core = solvent by def)

Run with:  python -m pytest sasmodels/models/tests/test_multi_shell_vesicle.py
"""
import numpy as np
import pytest

from sasmodels.core import load_model
from sasmodels.data import empty_data1D
from sasmodels.direct_model import DirectModel, call_Fq


def _run(model_name, pars):
    model = load_model(model_name)
    q = np.linspace(0.001, 1.0, 200)
    data = empty_data1D(q)
    finder = DirectModel(data, model)
    return finder(**pars)


def _vesicle_pars(radius=100.0, thickness=30.0, sld=0.5,
                  sld_solvent=6.36, volfraction=0.05):
    return dict(radius=radius, thickness=thickness, sld=sld,
                sld_solvent=sld_solvent, volfraction=volfraction,
                background=0.0, scale=1.0)


def _multi_shell_vesicle_pars(radius=100.0, thickness=30.0, sld=0.5,
                              sld_solvent=6.36, volfraction=0.05):
    return dict(radius=radius, n=1, thickness1=thickness, sld1=sld,
                sld_solvent=sld_solvent, volfraction=volfraction,
                background=0.0, scale=1.0)


@pytest.mark.parametrize(
    "kwargs",
    [dict(),
     dict(radius=50.0, thickness=80.0, sld=2.0, sld_solvent=1.0,
          volfraction=0.2),
     dict(radius=250.0, thickness=5.0, sld=-0.4, sld_solvent=6.3,
          volfraction=0.01),
     ])
def test_single_shell_matches_vesicle(kwargs):
    """multi_shell_vesicle with n=1 == vesicle with same parameters."""
    ves = _run("vesicle", _vesicle_pars(**kwargs))
    lp = _run("multi_shell_vesicle", _multi_shell_vesicle_pars(**kwargs))
    np.testing.assert_allclose(ves, lp, rtol=1e-10, atol=1e-10)


def test_two_shells_matches_manual_sum():
    """Two thin shells must equal one thick shell when SLDs are equal."""
    ves = _run("vesicle", _vesicle_pars(radius=100.0, thickness=40.0,
                                        sld=0.5))
    lp = _run("multi_shell_vesicle",
              dict(radius=100.0, n=2, thickness1=20.0, thickness2=20.0,
                   sld1=0.5, sld2=0.5, sld_solvent=6.36,
                   volfraction=0.05, background=0.0, scale=1.0))
    np.testing.assert_allclose(ves, lp, rtol=1e-10, atol=1e-10)


def test_zero_contrast_gives_zero_scattering():
    """All SLDs equal to solvent -> only background remains."""
    pars = dict(radius=100.0, n=3, thickness1=30.0, thickness2=15.0,
                thickness3=10.0, sld1=6.36, sld2=6.36, sld3=6.36,
                sld_solvent=6.36, volfraction=0.05,
                background=0.0, scale=1.0)
    result = _run("multi_shell_vesicle", pars)
    np.testing.assert_allclose(result, 0.0, atol=1e-10)


def test_effective_radius_modes():
    """Mode 1 (outer radius) and mode 2 (core radius).

    NB: sasmodels radius_effective modes are 1-based; mode 0 means the
    user-provided radius_effective parameter is used instead.
    """
    model = load_model("multi_shell_vesicle")
    assert model.info.radius_effective_modes == ["outer radius", "core radius"]
    pars = dict(radius=100.0, n=2, thickness1=30.0, thickness2=15.0)
    kernel = model.make_kernel([np.array([1e-3])])
    _, _, r_outer, _, _ = call_Fq(kernel, dict(pars, radius_effective_mode=1))
    _, _, r_core, _, _ = call_Fq(kernel, dict(pars, radius_effective_mode=2))
    assert r_outer == pytest.approx(145.0)
    assert r_core == pytest.approx(100.0)
