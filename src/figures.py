"""Part 7: figures for the revised manuscript (600 dpi PNG and TIFF).

  fig1_d86_curves       measured base D86 curve and delta-corrected 10 and 20 vol% MTBE blends
  fig2_volatility       (a) RVP, (b) E70 and E100, (c) DI against MTBE content
  fig4_windows          MTBE windows with octane uncertainty and distillation range
  fig3_validation       model responses at 10 vol% MTBE against published measurements

Main case: surrogate A, NRTL with fitted MTBE parameters, reflux calibrated on measured blends.
Bands: range over the other model settings and surrogate B.
Width: single column 90 mm, double column 180 mm (change COL1/COL2 for the target journal).

Run from the project folder (about 10 s):
    python src/figures.py
Outputs in figures/.
"""
import csv, json, os, sys
import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from thermo_models import Mixture, load_props, ROOT
from distillation import D86Model
from blend_distillation import corrected_curves, measured_base, MEASURED, STEP

FIG = os.path.join(ROOT, "figures")
MM = 1 / 25.4
COL1, COL2 = 90, 180
OI = {"black": "#000000", "orange": "#E69F00", "sky": "#56B4E9", "green": "#009E73",
      "blue": "#0072B2", "verm": "#D55E00", "purple": "#CC79A7", "grey": "#999999"}
mpl.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7,
    "legend.fontsize": 7, "axes.linewidth": 0.6, "lines.linewidth": 1.0, "lines.markersize": 3.5,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6, "axes.spines.top": False, "axes.spines.right": False,
    "legend.frameon": False, "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.bbox": "tight", "savefig.pad_inches": 0.02})


def save(fig, name):
    os.makedirs(FIG, exist_ok=True)
    fig.savefig(os.path.join(FIG, name + ".png"), dpi=600)
    fig.savefig(os.path.join(FIG, name + ".tiff"), dpi=600, pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
    print("saved figures/" + name + ".png and .tiff")


def label(ax, s):
    ax.text(-0.18, 1.04, s, transform=ax.transAxes, fontweight="bold", va="bottom")


def read_csv(path):
    with open(os.path.join(ROOT, path), encoding="utf-8") as f:
        return list(csv.DictReader(f))


def series(rows, key_case, case, col):
    r = [x for x in rows if x[key_case] == case]
    return np.array([float(x["MTBE vol%"]) for x in r]), np.array([float(x[col]) for x in r])


# ------------------------------------------------------------------ figure 1
def fig1():
    which = "A"
    fit = json.load(open(os.path.join(ROOT, "results/surrogate_fit", f"surrogate_{which}_rayleigh_fit.json")))
    R = json.load(open(os.path.join(ROOT, "results/distillation", f"reflux_calibration_{which}.json")))["R_blend"]
    sur = json.load(open(os.path.join(ROOT, "data/processed/surrogates.json")))["surrogate_" + which]
    names = [c["component"] for c in sur]; classes = [c["class"] for c in sur]
    v = np.array([fit["fitted_vol_pct"][n] for n in names]); P = load_props()
    m = D86Model(Mixture(names + ["mtbe"], classes + ["Ether"], P), P, model="nrtl", R=R)
    Es, Ts = [], []
    for x in (0.0, 10.0, 20.0):
        E, T = m.run(np.append(v * (100 - x) / v.sum(), x), step=STEP); Es.append(E); Ts.append(T)
    grid, Ec, _ = corrected_curves(Es, Ts)
    Em, Tm = measured_base()
    fig, ax = plt.subplots(figsize=(COL1 * MM, 70 * MM))
    for (lo, hi), t, txt in (((20, 48), 70, "E70"), ((46, 71), 100, "E100"), ((75, 100), 150, "E150")):
        ax.plot([lo, hi], [t, t], color=OI["grey"], lw=3, alpha=0.35, solid_capstyle="butt")
        ax.text(hi + 1, t, txt, va="center", fontsize=7, color="#555555")
    ax.plot(Em, Tm, "o-", color=OI["black"], label="Base, measured (ASTM D86)")
    ax.plot(Ec[1], grid, "--", color=OI["blue"], label="10 vol% MTBE")
    ax.plot(Ec[2], grid, "-", color=OI["verm"], label="20 vol% MTBE")
    ax.set_xlabel("Volume evaporated (vol%)"); ax.set_ylabel("Temperature (°C)")
    ax.set_xlim(0, 100); ax.set_ylim(40, 200); ax.legend(loc="upper left")
    save(fig, "fig1_d86_curves")


# ------------------------------------------------------------------ figure 2
def fig2():
    vp = []
    for s in ("A", "B"):
        vp += read_csv(f"results/vapor_pressure/vp_{s}_fitted.csv")
    d86 = {s: read_csv(f"results/distillation/d86_{s}.csv") for s in ("A", "B")}
    fig, axs = plt.subplots(1, 3, figsize=(COL2 * MM, 58 * MM))
    # (a) RVP
    ax = axs[0]
    ax.axhspan(44, 60, color=OI["grey"], alpha=0.15, lw=0)
    ax.text(0.5, 59, "Summer limits 44–60 kPa", va="top", fontsize=7, color="#555555")
    styles = {"nrtl": (OI["verm"], "-", "o", "NRTL, fitted"), "unifac": (OI["blue"], "--", "s", "UNIFAC"),
              "ideal": (OI["green"], ":", "^", "Ideal")}
    for mdl, (c, ls, mk, lab) in styles.items():
        x, y = series([r for r in vp if r["surrogate"] == "A_fitted"], "model", mdl, "RVP delta-corrected (kPa)")
        ax.plot(x, y, ls=ls, color=c, marker=mk, markevery=4, label=lab, clip_on=False)
    ax.set_xlabel("MTBE (vol%)"); ax.set_ylabel("RVP (kPa)"); ax.set_ylim(40, 62); ax.set_xlim(0, 20)
    ax.legend(loc="lower right"); label(ax, "(a)")
    # (b) E70 and E100, main case with band over other settings and surrogate B
    ax = axs[1]
    for key, lim, c, mk, lab in (("E70", 20, OI["blue"], "o", "E70"), ("E100", 46, OI["verm"], "s", "E100")):
        x, y = series(d86["A"], "case", "nrtl_R_blend", key)
        band = np.array([series(d86[s], "case", cs, key)[1] for s in d86 for cs in sorted({r["case"] for r in d86[s]})])
        ax.fill_between(x, band.min(0), band.max(0), color=c, alpha=0.15, lw=0)
        ax.plot(x, y, color=c, marker=mk, markevery=4, label=lab, clip_on=False)
        ax.axhline(lim, color=c, ls=":", lw=0.8)
    ax.set_xlabel("MTBE (vol%)"); ax.set_ylabel("Evaporated (vol%)"); ax.set_xlim(0, 20)
    ax.legend(loc="upper left"); label(ax, "(b)")
    # (c) DI
    ax = axs[2]
    x, y = series(d86["A"], "case", "nrtl_R_blend", "DI (F)")
    band = np.array([series(d86[s], "case", cs, "DI (F)")[1] for s in d86 for cs in sorted({r["case"] for r in d86[s]})])
    ax.fill_between(x, band.min(0), band.max(0), color=OI["black"], alpha=0.12, lw=0)
    ax.plot(x, y, color=OI["black"], marker="o", markevery=4, clip_on=False)
    ax.axhline(1250, color=OI["verm"], ls="--", lw=0.8)
    ax.text(20, 1252, "ASTM D4814 max.", ha="right", va="bottom", fontsize=7, color=OI["verm"])
    ax.set_xlabel("MTBE (vol%)"); ax.set_ylabel("DI (°F)"); ax.set_xlim(0, 20); label(ax, "(c)")
    fig.tight_layout(w_pad=2.0)
    save(fig, "fig2_volatility")


# ------------------------------------------------------------------ figure 3
def fig3():
    w = json.load(open(os.path.join(ROOT, "results/windows/windows_summary.json")))
    oc = json.load(open(os.path.join(ROOT, "results/octane/octane_summary.json")))["bounds"]
    o2 = w["checks"]["oxygen_limit_vol_pct"]
    rows = w["windows"][1:]
    fig, ax = plt.subplots(figsize=(COL2 * MM, 55 * MM))
    labels, ys = [], []
    for i, r in enumerate(rows):
        y = len(rows) - 1 - i; lo = float(r[1]); ys.append(y); labels.append(r[0])
        perf = r[0].startswith("Performance")
        ax.barh(y, o2 - lo, left=lo, height=0.5, color=OI["orange"] if perf else OI["sky"],
                edgecolor="black", lw=0.5, hatch=None if perf else "")
        if not perf:
            g = r[0].split(", ")[1]; p5, p95 = oc[g]["p5"], oc[g]["p95"]
            med = oc[g]["median"]
            ax.errorbar(med, y, xerr=[[max(0.0, med - p5)], [max(0.0, p95 - med)]], fmt="o", color="black",
                        ms=3, capsize=2, lw=0.8, clip_on=False)
        else:
            a, b = [float(t) for t in r[3].split("E100 ")[1].split(",")[0].split("-")]
            ax.plot([a, b], [y - 0.38, y - 0.38], color="black", lw=1.5)
            ax.text(b + 0.3, y - 0.38, "E100 range over all cases", va="center", fontsize=7)
        ax.text(lo, y + 0.3, f"{lo:g}", ha="center", va="bottom", fontsize=7)
    ax.axvline(o2, color="black", ls="--", lw=0.8)
    ax.text(o2, len(rows) - 0.4, f"O$_2$ 3.7 wt% ({o2:.1f} vol%)", ha="center", va="bottom", fontsize=7)
    ax.set_yticks(ys); ax.set_yticklabels(labels); ax.set_xlim(0, 25)
    ax.set_xlabel("MTBE (vol%)")
    save(fig, "fig4_windows")


# ------------------------------------------------------------------ figure 4
def fig4():
    cases = {s: json.load(open(os.path.join(ROOT, "results/distillation", f"d86_{s}_summary.json")))["cases"] for s in ("A", "B")}
    keys = (("dT10_10vol_C", "dT10_10vol_C", "ΔT10 at 10 vol% (°C)"),
            ("dT50_10vol_C", "dT50_10vol_C", "ΔT50 at 10 vol% (°C)"),
            ("dE70_per_vol_pct_0to15", "dE70_per_vol_pct", "ΔE70 per vol% MTBE"))
    show = (("nrtl_R_blend", "Calibrated reflux", OI["verm"], "o"), ("nrtl_R_calibrated", "Base-fit reflux", OI["blue"], "s"),
            ("nrtl_R_0", "No reflux", OI["green"], "^"))
    fig, axs = plt.subplots(1, 3, figsize=(COL2 * MM, 55 * MM))
    for ax, (mk_key, ms_key, lab), pl in zip(axs, keys, "abc"):
        meas = MEASURED[ms_key]
        ax.scatter(np.zeros(len(meas)), meas, color=OI["black"], marker="D", s=14, zorder=3, label="Measured")
        for j, (case, name, c, mk) in enumerate(show, 1):
            vals = [cases[s][case][mk_key] for s in ("A", "B")]
            ax.scatter([j] * 2, vals, color=c, marker=mk, s=14, label=name)
        ax.set_xticks(range(4)); ax.set_xticklabels(["Measured", "Cal.", "Base", "R=0"])
        ax.set_ylabel(lab); ax.set_xlim(-0.5, 3.5); label(ax, f"({pl})")
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", ncol=4, bbox_to_anchor=(0.5, 1.06))
    fig.tight_layout(w_pad=2.0)
    save(fig, "fig3_validation")


if __name__ == "__main__":
    fig1(); fig2(); fig3(); fig4()
