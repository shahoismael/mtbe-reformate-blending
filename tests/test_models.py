"""Checks of the property models. Run from the project folder:
    python -m pytest tests -q
"""
import json, os, sys
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
from models import tbp_to_d86, d86_to_tbp, rvp_d323, density_15C, d86_from_composition

P = json.load(open(os.path.join(ROOT, "data/processed/component_props.json")))
S = json.load(open(os.path.join(ROOT, "data/processed/surrogates.json")))


def test_d86_tbp_roundtrip():
    d86 = [43.2, 72.3, 93.6, 114.4, 138.2, 167.7, 192.0]
    assert np.allclose(tbp_to_d86(d86_to_tbp(d86)), d86, atol=1e-6)


def test_rvp_pure_component_equals_vapour_pressure():
    p = P["n-pentane"]
    rvp = rvp_d323([100], [p["rho15_kg_m3"]], [p["MW"]], [p["psat378_Pa"]])
    assert abs(rvp - p["psat378_Pa"] / 1000) < 0.01


def test_density_ideal_mixing():
    assert abs(density_15C([50, 50], [700, 800]) - 750) < 1e-9


def _model(which):
    c = S["surrogate_" + which]; n = [x["component"].lower() for x in c]
    v = [x["final_vol_pct"] for x in c]
    return v, [P[i]["Tb_C"] for i in n], [P[i]["rho15_kg_m3"] for i in n], \
        [P[i]["MW"] for i in n], [P[i]["psat378_Pa"] for i in n]


def test_model_reproduces_hysys_d86_for_original_surrogate_A():
    # HYSYS Boiling Point Curves utility, original surrogate A: T10..T90
    v, nbp, *_ = _model("A")
    d = d86_from_composition(v, nbp)[1:6]
    assert np.allclose(d, [62.7, 94.4, 114.2, 137.0, 169.2], atol=1.0)


def test_model_reproduces_hysys_rvp_for_original_surrogates():
    for which, hysys in (("A", 48.4), ("B", 49.4)):
        v, nbp, rho, MW, ps = _model(which)
        assert abs(rvp_d323(v, rho, MW, ps) - hysys) < 1.0


def test_nrtl_reproduces_measured_binary_mtbe_toluene():
    # Chamorro et al. (2004), 313.15 K: x_MTBE = 0.4992, P = 34.903 kPa; pure Psat 59.92 / 7.897 kPa
    from thermo_models import Mixture, load_props
    mix = Mixture(["mtbe", "toluene"], ["Ether", "Aromatic"], load_props())
    x = np.array([0.4992, 0.5008]); g = mix.gammas(x, 313.15, "nrtl")
    P = x[0] * g[0] * 59.92 + x[1] * g[1] * 7.897
    assert abs(P / 34.903 - 1) < 0.01


def test_nrtl_is_ideal_without_mtbe():
    from thermo_models import Mixture, load_props
    mix = Mixture(["mtbe", "toluene", "n-heptane"], ["Ether", "Aromatic", "n-Paraffin"], load_props())
    g = mix.gammas(np.array([0.0, 0.5, 0.5]), 310.0, "nrtl")
    assert np.allclose(g[1:], 1.0)


def _d86(names, classes, vol, R):
    from thermo_models import Mixture, load_props
    from distillation import D86Model
    P = load_props()
    m = D86Model(Mixture(names, classes, P), P, model="ideal", R=R)
    return m.run(vol, step=1.0, end=99), P


def test_d86_pure_component_is_flat_at_boiling_point():
    (E, T), P = _d86(["toluene", "mtbe"], ["Aromatic", "Ether"], [100, 0], 0.0)
    assert np.allclose(T[:50], P["toluene"]["Tb_C"], atol=0.5)


def test_d86_reflux_sharpens_separation():
    names, cls = ["benzene", "toluene", "mtbe"], ["Aromatic", "Aromatic", "Ether"]
    (E0, T0), _ = _d86(names, cls, [50, 50, 0], 0.0)
    (E5, T5), _ = _d86(names, cls, [50, 50, 0], 5.0)
    # more reflux: lower temperature early, higher late
    assert np.interp(10, E5, T5) < np.interp(10, E0, T0)
    assert np.interp(90, E5, T5) > np.interp(90, E0, T0)
