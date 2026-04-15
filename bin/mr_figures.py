#!/usr/bin/env python3
"""
mr_figures.py — MR Pipeline CDC 1.0.0  Stage B5

"""

import argparse
import glob
import os
import sys
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from matplotlib import rcParams

# ── Global style ────────────────────────────────────────────────
rcParams.update({
    "font.family":         "sans-serif",
    "font.sans-serif":     ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size":           11,
    "axes.linewidth":      1.0,
    "axes.spines.top":     False,
    "axes.spines.right":   False,
    "axes.grid":           False,
    "axes.labelsize":      11,
    "axes.titlesize":      12,
    "axes.titleweight":    "bold",
    "axes.titlepad":       14,
    "xtick.labelsize":     10,
    "ytick.labelsize":     10,
    "xtick.major.width":   1.0,
    "ytick.major.width":   1.0,
    "xtick.major.size":    4.5,
    "ytick.major.size":    4.5,
    "xtick.major.pad":     5,
    "ytick.major.pad":     5,
    "xtick.minor.visible": False,
    "ytick.minor.visible": False,
    "legend.frameon":      True,
    "legend.framealpha":   0.93,
    "legend.edgecolor":    "#c5d8ee",
    "legend.fontsize":     9.5,
    "legend.title_fontsize": 9.5,
    "legend.borderpad":    0.6,
    "legend.labelspacing": 0.5,
    "figure.dpi":          300,
    "savefig.dpi":         300,
    "savefig.bbox":        "tight",
    "savefig.pad_inches":  0.22,
    "pdf.fonttype":        42,
    "ps.fonttype":         42,
})

# ── Palette ─────────────────────────────────────────────────────
COL = {
    "ivw":    "#1f4e8c",
    "egger":  "#c0392b",
    "wm":     "#27ae60",
    "mvmr":   "#7d3c98",
    "het":    "#e67e22",
    "sig":    "#c0392b",
    "ns":     "#7f8c8d",
    "zero":   "#95a5a6",
}

METHOD_COL = {
    "inverse variance weighted": COL["ivw"],
    "mr egger":                  COL["egger"],
    "weighted median":           COL["wm"],
    "mvmr-ivw":                  COL["ivw"],
    "mvmr-egger":                COL["egger"],
}

def method_color(m):
    return METHOD_COL.get(str(m).strip().lower(), "#555555")


AUTHOR_STR = ""   # set from args in main()
INSTITUTE_STR = ""

def save_fig(fig, path_stem, caption=None):
    """Save figure with optional bottom caption line showing author/institute."""
    if caption:
        fig.text(
            0.5, 0.005, caption,
            ha="center", va="bottom",
            fontsize=8, color="#5a7299",
            fontstyle="italic",
            transform=fig.transFigure,
        )
    fig.savefig(f"{path_stem}.png")
    fig.savefig(f"{path_stem}.pdf")
    print(f"  Saved: {path_stem}.png / .pdf")
    plt.close(fig)


def make_caption(fig_num, description):
    parts = [f"Figure {fig_num}. {description}"]
    if AUTHOR_STR:
        parts.append(AUTHOR_STR)
    if INSTITUTE_STR:
        parts.append(INSTITUTE_STR)
    return "  |  ".join(parts)


def pfmt(p):
    try:
        p = float(p)
        if p < 1e-4:  return f"{p:.1e}"
        return f"{p:.4f}"
    except: return str(p)


def load_tsvs(pattern):
    rows = []
    for f in sorted(glob.glob(pattern)):
        try:
            df = pd.read_csv(f, sep="\t", low_memory=False)
            rows.append(df)
        except Exception as e:
            print(f"  [WARN] {f}: {e}", file=sys.stderr)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


# ────────────────────────────────────────────────────────────────
# F1  Forest plot — uvMR
# ────────────────────────────────────────────────────────────────
def figure_forest_uvmr(mr_df, out_dir):
    if mr_df.empty:
        print("  [F1] No data — skipping"); return

    # Normalise column names
    mr = mr_df.copy()
    mr.columns = [c.lower() for c in mr.columns]
    for col in ["b","beta"]:
        if col in mr.columns and "b" not in mr.columns:
            mr["b"] = mr[col]
    if "b" not in mr.columns: print("  [F1] No beta column — skipping"); return

    mr["b"]    = pd.to_numeric(mr.get("b", mr.get("beta")), errors="coerce")
    mr["se"]   = pd.to_numeric(mr["se"], errors="coerce")
    mr["pval"] = pd.to_numeric(mr.get("pval", mr.get("pval.exposure",
                  mr.get("p", np.nan))), errors="coerce")
    mr = mr.dropna(subset=["b","se"])

    # Keep only main methods
    keep = ["inverse variance weighted","mr egger","weighted median"]
    mr["method_lc"] = mr["method"].str.strip().str.lower()
    mr = mr[mr["method_lc"].isin(keep)].copy()
    if mr.empty: print("  [F1] No data after method filter"); return

    mr["lo"] = mr["b"] - 1.96 * mr["se"]
    mr["hi"] = mr["b"] + 1.96 * mr["se"]

    # Build row labels
    def row_label(r):
        pair = r.get("causal_pair","")
        dirn = str(r.get("direction","")).replace("->","→")
        return f"{pair}  ·  {dirn}"

    mr["label"] = mr.apply(row_label, axis=1)

    n_rows = len(mr)
    row_h  = 0.55     # cm per row
    fig_h  = max(8, n_rows * row_h / 2.54 + 2.5)
    fig_w  = 22 / 2.54   # wider for legend outside

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    fig.subplots_adjust(right=0.72, bottom=0.08)  # space for legend + caption

    x_all = np.concatenate([mr["lo"].values, mr["hi"].values, [0]])
    x_min = np.nanmin(x_all)
    x_max = np.nanmax(x_all)
    pad   = (x_max - x_min) * 0.12 or 0.3
    ax.set_xlim(x_min - pad, x_max + pad)

    y_pos  = np.arange(n_rows)[::-1]
    labels = []
    legend_handles = {}

    for i, (_, row) in enumerate(mr.iterrows()):
        y = y_pos[i]
        b, lo, hi = row["b"], row["lo"], row["hi"]
        p   = row.get("pval", np.nan)
        col = method_color(row["method"])
        sig = pd.notna(p) and float(p) < 0.05

        # CI line
        ax.plot([lo, hi], [y, y], color=col, lw=1.6, solid_capstyle="round",
                zorder=3)
        # Cap ticks
        ax.plot([lo, lo], [y - 0.12, y + 0.12], color=col, lw=1.2, zorder=3)
        ax.plot([hi, hi], [y - 0.12, y + 0.12], color=col, lw=1.2, zorder=3)
        # Diamond point
        d = 0.18
        diamond = plt.Polygon([[b, y-d],[b+d, y],[b, y+d],[b-d, y]],
                               closed=True, fc=col if sig else "white",
                               ec=col, lw=1.4, zorder=4)
        ax.add_patch(diamond)

        # p-value label on right
        ax.text(x_max + pad * 0.55, y, pfmt(p),
                va="center", ha="left", fontsize=8,
                color=COL["sig"] if sig else COL["ns"],
                fontweight="bold" if sig else "normal")

        labels.append(row["label"])

        # Legend handles
        mname = row["method"].strip()
        if mname not in legend_handles:
            legend_handles[mname] = mlines.Line2D(
                [], [], color=col, lw=2,
                marker="D", markersize=5,
                markerfacecolor=col, markeredgecolor=col,
                label=mname)

    # Null line
    ax.axvline(0, color=COL["zero"], lw=0.9, ls="--", zorder=1)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=8.5)
    ax.set_xlabel("Causal effect estimate (β)", fontsize=10, labelpad=8)
    ax.spines["left"].set_visible(False)
    ax.tick_params(left=False)

    ax.text(x_max + pad * 0.55, n_rows - 0.2, "P-value",
            va="center", ha="left", fontsize=8.5, fontweight="bold")

    ax.legend(handles=list(legend_handles.values()),
              loc="upper left", bbox_to_anchor=(1.01, 1.0),
              borderaxespad=0, fontsize=9.5,
              handlelength=2.2, handletextpad=0.7,
              title="Method", title_fontsize=9.5,
              frameon=True, framealpha=0.95)

    ax.set_title("Univariable MR — Causal Effect Estimates",
                 fontsize=11, fontweight="bold", pad=12)

    cap = make_caption("F1", "Univariable MR — Causal Effect Estimates (IVW · MR-Egger · Weighted Median). β ± 95% CI; filled ◆ = p<0.05.")
    save_fig(fig, os.path.join(out_dir, "F1_forest_uvmr"), caption=cap)


# ────────────────────────────────────────────────────────────────
# F2  Heterogeneity panel
# ────────────────────────────────────────────────────────────────
def figure_heterogeneity(het_df, out_dir):
    if het_df.empty:
        print("  [F2] No data — skipping"); return

    het = het_df.copy()
    het.columns = [c.lower() for c in het.columns]
    het["pvalue"]    = pd.to_numeric(het.get("pvalue"), errors="coerce")
    het["statistic"] = pd.to_numeric(het.get("statistic"), errors="coerce")
    het["i2"]        = pd.to_numeric(het.get("i2"), errors="coerce")

    # Build label
    het["label"] = (het.get("causal_pair","").astype(str) + "\n" +
                    het.get("direction","").astype(str).str.replace("->","→"))

    tests = {
        "Cochran_Q":       ("Cochran Q  (Q statistic)",  "statistic"),
        "I2":              ("I²  (%)",                   "i2"),
        "Egger_intercept": ("Egger intercept",           "statistic"),
        "Steiger":         ("Steiger pass / fail",       "statistic"),
    }

    available = [k for k in tests if k in het["test"].values]
    if not available:
        print("  [F2] No recognised tests — skipping"); return

    ncols = 2
    nrows = (len(available) + 1) // 2
    fig_w = 18 / 2.54
    fig_h = nrows * 6.5 / 2.54 + 1.0

    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(fig_w, fig_h),
                             constrained_layout=True)
    fig.get_layout_engine().set(h_pad=0.12, w_pad=0.08, hspace=0.08, wspace=0.06)
    axes = np.array(axes).flatten()

    for ax_i, test_key in enumerate(available):
        ax    = axes[ax_i]
        title, val_col = tests[test_key]
        sub   = het[het["test"] == test_key].copy()
        sub   = sub.dropna(subset=[val_col])
        if sub.empty: ax.axis("off"); continue

        vals  = sub[val_col].values
        lbls  = sub["label"].values
        pvs   = sub["pvalue"].values
        n     = len(vals)
        y     = np.arange(n)

        colors = [COL["sig"] if pd.notna(p) and float(p) < 0.05
                  else COL["ivw"] for p in pvs]

        bars = ax.barh(y, vals, color=colors, height=0.55,
                       edgecolor="white", linewidth=0.5)

        # p-value annotations
        for j, (v, p) in enumerate(zip(vals, pvs)):
            if pd.notna(p):
                ax.text(v + np.nanmax(np.abs(vals)) * 0.03, j,
                        pfmt(p), va="center", fontsize=7.5,
                        color=COL["sig"] if float(p) < 0.05 else COL["ns"])

        # Reference lines
        if test_key == "I2":
            for thresh, ls in [(25, ":"), (50, "--"), (75, "-.")]:
                ax.axvline(thresh, color="#bdc3c7", lw=0.8, ls=ls, zorder=0)
                ax.text(thresh, n - 0.3, f"{thresh}%",
                        ha="center", va="bottom", fontsize=7, color="#95a5a6")

        ax.set_yticks(y)
        ax.set_yticklabels(lbls, fontsize=7.5)
        ax.set_xlabel(val_col.replace("_"," ").title(), fontsize=9, labelpad=6)
        ax.set_title(title, fontsize=9.5, fontweight="bold", pad=8)
        ax.spines["left"].set_visible(False)
        ax.tick_params(left=False)

    # Turn off unused panels
    for ax_i in range(len(available), len(axes)):
        axes[ax_i].axis("off")

    fig.suptitle("Heterogeneity & Pleiotropy Diagnostics",
                 fontsize=12, fontweight="bold")

    cap = make_caption("F2", "Heterogeneity & Pleiotropy Diagnostics. Red bars/values = p<0.05. I² thresholds: 25% low, 50% moderate, 75% high.")
    save_fig(fig, os.path.join(out_dir, "F2_heterogeneity"), caption=cap)


# ────────────────────────────────────────────────────────────────
# F3  Scatter + funnel (one panel per analysis)
# ────────────────────────────────────────────────────────────────
def figure_scatter_funnel(wald_df, out_dir):
    if wald_df.empty:
        print("  [F3] No data — skipping"); return

    w = wald_df.copy()
    w.columns = [c.lower() for c in w.columns]
    for c in ["beta_exp","beta_out","se_exp","se_out","f_stat","wald_ratio"]:
        if c in w.columns:
            w[c] = pd.to_numeric(w[c], errors="coerce")

    # Direction label
    w["analysis"] = (w.get("exposure","").astype(str) + "->" +
                     w.get("outcome","").astype(str) + "  [" +
                     w.get("direction","").astype(str) + "]")

    analyses = w["analysis"].dropna().unique()[:10]  # cap at 10
    n        = len(analyses)
    if n == 0: print("  [F3] No analyses — skipping"); return

    cols   = 2
    rows_f = (n + cols - 1) // cols
    fig_w  = 18 / 2.54
    fig_h  = rows_f * 7.5 / 2.54 + 1.0

    # Scatter plots
    fig_s, axes_s = plt.subplots(rows_f, cols,
                                  figsize=(fig_w, fig_h),
                                  constrained_layout=True)
    axes_s = np.array(axes_s).flatten()

    fig_f, axes_f = plt.subplots(rows_f, cols,
                                  figsize=(fig_w, fig_h),
                                  constrained_layout=True)
    axes_f = np.array(axes_f).flatten()

    for i, analysis in enumerate(analyses):
        sub = w[w["analysis"] == analysis].dropna(subset=["beta_exp","beta_out"])
        if sub.empty: continue
        ax_s = axes_s[i]
        ax_f = axes_f[i]

        bx = sub["beta_exp"].values
        by = sub["beta_out"].values

        # ── Scatter ──────────────────────────────────────────
        ax_s.scatter(bx, by, s=20, color=COL["ivw"], alpha=0.75,
                     linewidths=0.4, edgecolors="white", zorder=4)

        if len(bx) >= 2:
            # IVW slope (origin-constrained)
            w_ivw = 1 / (sub["se_out"].values ** 2 + 1e-12)
            slope = np.sum(w_ivw * bx * by) / np.sum(w_ivw * bx**2)
            x_rng = np.linspace(bx.min(), bx.max(), 100)
            ax_s.plot(x_rng, slope * x_rng, color=COL["ivw"],
                      lw=1.6, label=f"IVW β={slope:.3f}", zorder=3)

        ax_s.axhline(0, color=COL["zero"], lw=0.8, ls="--")
        ax_s.axvline(0, color=COL["zero"], lw=0.8, ls="--")
        ax_s.set_xlabel("SNP–exposure β", fontsize=8.5, labelpad=5)
        ax_s.set_ylabel("SNP–outcome β",  fontsize=8.5, labelpad=5)
        ax_s.set_title(analysis, fontsize=8, fontweight="bold", pad=6)
        ax_s.tick_params(labelsize=7.5)
        if len(bx) >= 2:
            ax_s.legend(fontsize=7.5, handlelength=1.5)

        # ── Funnel ───────────────────────────────────────────
        wr = sub["wald_ratio"].values if "wald_ratio" in sub.columns else by / bx
        fs = sub["f_stat"].values if "f_stat" in sub.columns else np.ones(len(bx))
        precision = 1.0 / np.sqrt(np.clip(fs, 1e-6, None))

        ax_f.scatter(wr, precision, s=20, color=COL["het"],
                     alpha=0.75, linewidths=0.4, edgecolors="white", zorder=4)

        ivw_est = np.nanmedian(wr)
        ax_f.axvline(ivw_est, color=COL["ivw"], lw=1.2, ls="--",
                     label=f"Median = {ivw_est:.3f}")
        ax_f.axvline(0, color=COL["zero"], lw=0.8, ls=":")

        ax_f.set_xlabel("Wald ratio", fontsize=8.5, labelpad=5)
        ax_f.set_ylabel("1 / √F-stat",  fontsize=8.5, labelpad=5)
        ax_f.set_title(analysis, fontsize=8, fontweight="bold", pad=6)
        ax_f.tick_params(labelsize=7.5)
        ax_f.legend(fontsize=7.5, handlelength=1.5)

        # Invert y so small precision at bottom
        ax_f.invert_yaxis()

    for i in range(len(analyses), len(axes_s)):
        axes_s[i].axis("off")
        axes_f[i].axis("off")

    fig_s.suptitle("MR Scatter Plots — SNP–Exposure vs SNP–Outcome",
                   fontsize=10, fontweight="bold", y=1.01)
    fig_f.suptitle("Funnel Plots — Wald Ratio vs Instrument Precision",
                   fontsize=10, fontweight="bold", y=1.01)

    # Save as combined figure (scatter top, funnel bottom)
    n_panels = len(analyses)
    cols_c = 2
    rows_c = (n_panels + cols_c - 1) // cols_c
    fig_c, axes_c = plt.subplots(rows_c * 2, cols_c,
                                  figsize=(20/2.54, rows_c * 14/2.54 + 1.0),
                                  constrained_layout=True)

    for i, analysis in enumerate(analyses):
        sub = w[w["analysis"] == analysis].dropna(subset=["beta_exp","beta_out"])
        if sub.empty: continue
        row_i = (i // cols_c) * 2
        col_i = i % cols_c

        bx = sub["beta_exp"].values
        by = sub["beta_out"].values
        ax_s = axes_c[row_i, col_i]
        ax_f = axes_c[row_i + 1, col_i]

        ax_s.scatter(bx, by, s=20, color=COL["ivw"], alpha=0.75,
                     linewidths=0.4, edgecolors="white", zorder=4)
        if len(bx) >= 2:
            w_ivw = 1 / (sub["se_out"].values ** 2 + 1e-12)
            slope = np.sum(w_ivw * bx * by) / np.sum(w_ivw * bx**2)
            x_rng = np.linspace(bx.min(), bx.max(), 100)
            ax_s.plot(x_rng, slope * x_rng, color=COL["ivw"], lw=1.5,
                      label=f"IVW β={slope:.3f}")
        ax_s.axhline(0, color=COL["zero"], lw=0.7, ls="--")
        ax_s.axvline(0, color=COL["zero"], lw=0.7, ls="--")
        ax_s.set_xlabel("SNP–exposure β", fontsize=8, labelpad=4)
        ax_s.set_ylabel("SNP–outcome β", fontsize=8, labelpad=4)
        ax_s.set_title(f"Scatter — {analysis}", fontsize=8, fontweight="bold", pad=5)
        ax_s.tick_params(labelsize=7)
        if len(bx) >= 2: ax_s.legend(fontsize=7)

        wr = sub["wald_ratio"].values if "wald_ratio" in sub.columns else by / bx
        fs = sub["f_stat"].values if "f_stat" in sub.columns else np.ones(len(bx))
        prec = 1.0 / np.sqrt(np.clip(fs, 1e-6, None))
        ax_f.scatter(wr, prec, s=20, color=COL["het"], alpha=0.75,
                     linewidths=0.4, edgecolors="white", zorder=4)
        ax_f.axvline(np.nanmedian(wr), color=COL["ivw"], lw=1.2, ls="--")
        ax_f.axvline(0, color=COL["zero"], lw=0.7, ls=":")
        ax_f.invert_yaxis()
        ax_f.set_xlabel("Wald ratio", fontsize=8, labelpad=4)
        ax_f.set_ylabel("1 / √F-stat", fontsize=8, labelpad=4)
        ax_f.set_title(f"Funnel — {analysis}", fontsize=8, fontweight="bold", pad=5)
        ax_f.tick_params(labelsize=7)

    for i in range(n_panels, rows_c * cols_c):
        row_i = (i // cols_c) * 2
        col_i = i % cols_c
        axes_c[row_i, col_i].axis("off")
        axes_c[row_i + 1, col_i].axis("off")

    fig_c.suptitle("Scatter & Funnel Plots  (top row: scatter · bottom row: funnel)",
                   fontsize=12, fontweight="bold")
    cap = make_caption("F3", "Scatter plots (SNP–exposure vs SNP–outcome) and funnel plots (Wald ratio vs 1/√F-stat). Dashed line = IVW slope.")
    save_fig(fig_c, os.path.join(out_dir, "F3_scatter_funnel"), caption=cap)
    plt.close(fig_s); plt.close(fig_f)


# ────────────────────────────────────────────────────────────────
# F4  Leave-one-out forest
# ────────────────────────────────────────────────────────────────
def figure_loo(loo_df, out_dir):
    if loo_df.empty:
        print("  [F4] No data — skipping"); return

    loo = loo_df.copy()
    loo.columns = [c.lower() for c in loo.columns]
    loo["b"]  = pd.to_numeric(loo.get("b",  loo.get("beta", np.nan)), errors="coerce")
    loo["se"] = pd.to_numeric(loo.get("se", np.nan),                  errors="coerce")
    loo["p"]  = pd.to_numeric(loo.get("p",  loo.get("pval", np.nan)), errors="coerce")
    loo = loo.dropna(subset=["b","se"])

    analyses = loo.get("analysis_id", loo.get("causal_pair","")).unique()
    n_total  = len(loo)
    if n_total == 0: print("  [F4] Empty LOO — skipping"); return

    row_h  = 0.45
    fig_h  = max(7, n_total * row_h / 2.54 + 2.5)
    fig_w  = 21 / 2.54   # wider for legend outside

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    fig.subplots_adjust(right=0.75, bottom=0.07)

    loo["lo"] = loo["b"] - 1.96 * loo["se"]
    loo["hi"] = loo["b"] + 1.96 * loo["se"]

    x_all = np.concatenate([loo["lo"].values, loo["hi"].values, [0]])
    x_min = np.nanmin(x_all); x_max = np.nanmax(x_all)
    pad   = (x_max - x_min) * 0.14 or 0.3
    ax.set_xlim(x_min - pad, x_max + pad)

    y_pos  = np.arange(len(loo))[::-1]
    labels = []
    group_colors = ["#1f4e8c","#2980b9","#27ae60","#8e44ad","#c0392b",
                    "#e67e22","#16a085","#2c3e50","#7f8c8d","#6c5ce7"]

    anal_list = list(dict.fromkeys(loo.get("analysis_id",
                     loo.get("causal_pair","")).values))
    anal_col  = {a: group_colors[i % len(group_colors)]
                 for i, a in enumerate(anal_list)}

    for i, (_, row) in enumerate(loo.iterrows()):
        y   = y_pos[i]
        b, lo, hi = row["b"], row["lo"], row["hi"]
        anal = row.get("analysis_id", row.get("causal_pair",""))
        col  = anal_col.get(anal, COL["ivw"])

        ax.plot([lo, hi], [y, y], color=col, lw=1.3,
                solid_capstyle="round", zorder=3)
        ax.plot([lo, lo], [y-0.1, y+0.1], color=col, lw=1.0, zorder=3)
        ax.plot([hi, hi], [y-0.1, y+0.1], color=col, lw=1.0, zorder=3)

        # Circle point
        ax.plot(b, y, "o", color=col, ms=4.5,
                markeredgecolor="white", markeredgewidth=0.8, zorder=4)

        snp_name = str(row.get("snp", "")).replace("All - excluding","excl.")
        labels.append(snp_name[:45])

    ax.axvline(0, color=COL["zero"], lw=0.9, ls="--", zorder=1)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlabel("IVW estimate (β)", fontsize=10, labelpad=8)
    ax.spines["left"].set_visible(False)
    ax.tick_params(left=False)
    ax.set_title("Leave-One-Out Analysis — IVW Estimate Stability",
                 fontsize=11, fontweight="bold", pad=12)

    # Legend
    handles = [mpatches.Patch(color=c, label=a)
               for a, c in list(anal_col.items())[:8]]
    ax.legend(handles=handles, title="Analysis",
              loc="upper left", bbox_to_anchor=(1.01, 1.0),
              borderaxespad=0, fontsize=9, title_fontsize=9,
              handlelength=1.4, frameon=True, framealpha=0.95)

    cap = make_caption("F4", "Leave-One-Out Analysis — IVW estimate when each instrument is excluded in turn. Stability indicates absence of influential variants.")
    save_fig(fig, os.path.join(out_dir, "F4_leave_one_out"), caption=cap)


# ────────────────────────────────────────────────────────────────
# F5  Single-SNP MR forest
# ────────────────────────────────────────────────────────────────
def figure_single_snp(snp_df, out_dir):
    if snp_df.empty:
        print("  [F5] No data — skipping"); return

    snp = snp_df.copy()
    snp.columns = [c.lower() for c in snp.columns]
    snp["b"]  = pd.to_numeric(snp.get("b",  snp.get("beta", np.nan)), errors="coerce")
    snp["se"] = pd.to_numeric(snp.get("se", np.nan),                  errors="coerce")
    snp["p"]  = pd.to_numeric(snp.get("p",  snp.get("pval", np.nan)), errors="coerce")
    snp = snp.dropna(subset=["b","se"])

    # Drop IVW/Egger summary rows (they have SNP == "All - IVW" etc.)
    mask = snp.get("snp","").astype(str).str.startswith("All")
    snp  = snp[~mask]
    if snp.empty: print("  [F5] No individual SNP rows"); return

    snp["lo"] = snp["b"] - 1.96 * snp["se"]
    snp["hi"] = snp["b"] + 1.96 * snp["se"]

    n_rows = len(snp)
    row_h  = 0.45
    fig_h  = max(7, n_rows * row_h / 2.54 + 2.5)
    fig_w  = 21 / 2.54   # wider for legend outside

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    fig.subplots_adjust(right=0.76, bottom=0.07)

    x_all = np.concatenate([snp["lo"].values, snp["hi"].values, [0]])
    x_min = np.nanmin(x_all); x_max = np.nanmax(x_all)
    pad   = (x_max - x_min) * 0.14 or 0.3
    ax.set_xlim(x_min - pad, x_max + pad)

    y_pos  = np.arange(n_rows)[::-1]
    labels = []

    for i, (_, row) in enumerate(snp.iterrows()):
        y   = y_pos[i]
        b, lo, hi = row["b"], row["lo"], row["hi"]
        p   = row.get("p", np.nan)
        col = COL["sig"] if pd.notna(p) and float(p) < 0.05 else COL["ivw"]

        ax.plot([lo, hi], [y, y], color=col, lw=1.3,
                solid_capstyle="round", zorder=3)
        ax.plot([lo, lo], [y-0.1, y+0.1], color=col, lw=1.0, zorder=3)
        ax.plot([hi, hi], [y-0.1, y+0.1], color=col, lw=1.0, zorder=3)
        ax.plot(b, y, "s", color=col, ms=4,
                markeredgecolor="white", markeredgewidth=0.7, zorder=4)

        f_str = f"  F={float(row['f_stat']):.0f}" if "f_stat" in row and pd.notna(row.get("f_stat")) else ""
        lbl   = f"{row.get('snp','')}"[:40] + f_str
        labels.append(lbl)

    ax.axvline(0, color=COL["zero"], lw=0.9, ls="--", zorder=1)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=7.5)
    ax.set_xlabel("Wald ratio (β)", fontsize=10, labelpad=8)
    ax.spines["left"].set_visible(False)
    ax.tick_params(left=False)
    ax.set_title("Single-SNP MR — Per-Instrument Wald Ratios",
                 fontsize=11, fontweight="bold", pad=12)

    sig_patch = mpatches.Patch(color=COL["sig"], label="p < 0.05")
    ns_patch  = mpatches.Patch(color=COL["ivw"], label="p ≥ 0.05")
    ax.legend(handles=[sig_patch, ns_patch], title="Significance",
              loc="upper left", bbox_to_anchor=(1.01, 1.0),
              borderaxespad=0, fontsize=9.5, title_fontsize=9.5,
              handlelength=1.4, frameon=True, framealpha=0.95)

    cap = make_caption("F5", "Single-SNP MR — Per-Instrument Wald Ratios. Red = p<0.05; squares = point estimates ± 95% CI. F-statistic annotated.")
    save_fig(fig, os.path.join(out_dir, "F5_single_snp"), caption=cap)


# ────────────────────────────────────────────────────────────────
# F6  MVMR forest
# ────────────────────────────────────────────────────────────────
def figure_mvmr_forest(mvmr_df, out_dir):
    if mvmr_df.empty:
        print("  [F6] No data — skipping"); return

    mv = mvmr_df.copy()
    mv.columns = [c.lower() for c in mv.columns]
    mv["beta"]  = pd.to_numeric(mv.get("beta", mv.get("b", np.nan)), errors="coerce")
    mv["se"]    = pd.to_numeric(mv["se"], errors="coerce")
    mv["pval"]  = pd.to_numeric(mv.get("pval", np.nan), errors="coerce")
    mv["cond_f"]= pd.to_numeric(mv.get("cond_f", np.nan), errors="coerce")
    mv = mv.dropna(subset=["beta","se"])

    # Only IVW for the forest
    ivw = mv[mv["method"].str.lower().str.contains("ivw", na=False)].copy()
    if ivw.empty: print("  [F6] No MVMR-IVW rows"); return

    ivw["lo"] = ivw["beta"] - 1.96 * ivw["se"]
    ivw["hi"] = ivw["beta"] + 1.96 * ivw["se"]

    def row_label(r):
        exp = str(r.get("exposure",""))
        lbl = str(r.get("label",""))
        dirn = str(r.get("direction","")).replace("->","→")
        return f"{exp}  [{lbl}  {dirn}]"

    ivw["row_lbl"] = ivw.apply(row_label, axis=1)

    n_rows = len(ivw)
    row_h  = 0.55
    fig_h  = max(7, n_rows * row_h / 2.54 + 2.5)
    fig_w  = 22 / 2.54  # wider for legend outside

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    fig.subplots_adjust(right=0.70, bottom=0.08)

    x_all = np.concatenate([ivw["lo"].values, ivw["hi"].values, [0]])
    x_min = np.nanmin(x_all); x_max = np.nanmax(x_all)
    pad   = (x_max - x_min) * 0.14 or 0.3
    ax.set_xlim(x_min - pad, x_max + pad)

    y_pos  = np.arange(n_rows)[::-1]
    labels = []

    # Colour by outcome
    outcomes    = ivw["outcome"].unique()
    out_colors  = {o: COL[k] for o, k in
                   zip(outcomes, ["ivw","egger","wm","mvmr","het"])}

    for i, (_, row) in enumerate(ivw.iterrows()):
        y   = y_pos[i]
        b, lo, hi = row["beta"], row["lo"], row["hi"]
        p   = row.get("pval", np.nan)
        col = out_colors.get(row.get("outcome",""), COL["mvmr"])
        sig = pd.notna(p) and float(p) < 0.05

        ax.plot([lo, hi], [y, y], color=col, lw=1.8,
                solid_capstyle="round", zorder=3)
        ax.plot([lo, lo], [y-0.13, y+0.13], color=col, lw=1.2, zorder=3)
        ax.plot([hi, hi], [y-0.13, y+0.13], color=col, lw=1.2, zorder=3)

        d = 0.18
        diamond = plt.Polygon([[b, y-d],[b+d, y],[b, y+d],[b-d, y]],
                               closed=True, fc=col if sig else "white",
                               ec=col, lw=1.4, zorder=4)
        ax.add_patch(diamond)

        cf = row.get("cond_f", np.nan)
        cf_str = f"  condF={float(cf):.1f}" if pd.notna(cf) else ""
        ax.text(x_max + pad * 0.55, y,
                pfmt(p) + cf_str, va="center", ha="left", fontsize=7.5,
                color=COL["sig"] if sig else COL["ns"],
                fontweight="bold" if sig else "normal")

        labels.append(row["row_lbl"])

    ax.axvline(0, color=COL["zero"], lw=0.9, ls="--", zorder=1)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("Direct causal effect (β, MVMR-IVW)", fontsize=10, labelpad=8)
    ax.spines["left"].set_visible(False)
    ax.tick_params(left=False)

    ax.text(x_max + pad * 0.55, n_rows - 0.3, "P  condF",
            va="center", ha="left", fontsize=8.5, fontweight="bold")

    handles = [mpatches.Patch(color=c, label=o)
               for o, c in out_colors.items()]
    ax.legend(handles=handles, title="Outcome", fontsize=9.5,
              title_fontsize=9.5,
              loc="upper left", bbox_to_anchor=(1.01, 1.0),
              borderaxespad=0, handlelength=1.4, frameon=True, framealpha=0.95)

    ax.set_title("Multivariable MR — Direct Causal Effects (MVMR-IVW)",
                 fontsize=11, fontweight="bold", pad=12)

    cap = make_caption("F6", "Multivariable MR — Direct Causal Effects (MVMR-IVW). β ± 95% CI; filled ◆ = p<0.05. Conditional F-statistic and p-value annotated.")
    save_fig(fig, os.path.join(out_dir, "F6_mvmr_forest"), caption=cap)


# ────────────────────────────────────────────────────────────────
# Main
# ────────────────────────────────────────────────────────────────
def main():
    p = argparse.ArgumentParser(description="MR publication figures")
    p.add_argument("--mr_dir",    required=True)
    p.add_argument("--out_dir",   required=True)
    p.add_argument("--author",    default="")
    p.add_argument("--institute", default="")
    p.add_argument("--affiliation", default="")
    args = p.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    global AUTHOR_STR, INSTITUTE_STR
    AUTHOR_STR    = args.author
    INSTITUTE_STR = args.institute

    print("\n" + "="*60)
    print("  MR Pipeline — Generating Publication Figures")
    print(f"  Author    : {args.author}")
    print(f"  Institute : {args.institute}")
    print("="*60 + "\n")

    mr_df   = load_tsvs(os.path.join(args.mr_dir, "*_mr.tsv"))
    wald_df = load_tsvs(os.path.join(args.mr_dir, "*_wald.tsv"))
    het_df  = load_tsvs(os.path.join(args.mr_dir, "*_het.tsv"))
    loo_df  = load_tsvs(os.path.join(args.mr_dir, "*_loo.tsv"))
    snp_df  = load_tsvs(os.path.join(args.mr_dir, "*_single.tsv"))
    mvmr_df = load_tsvs(os.path.join(args.mr_dir, "*_mvmr.tsv"))

    print(f"[F1] Forest plot (uvMR)...")
    figure_forest_uvmr(mr_df, args.out_dir)

    print(f"[F2] Heterogeneity panel...")
    figure_heterogeneity(het_df, args.out_dir)

    print(f"[F3] Scatter + funnel...")
    figure_scatter_funnel(wald_df, args.out_dir)

    print(f"[F4] Leave-one-out...")
    figure_loo(loo_df, args.out_dir)

    print(f"[F5] Single-SNP MR...")
    figure_single_snp(snp_df, args.out_dir)

    print(f"[F6] MVMR forest...")
    figure_mvmr_forest(mvmr_df, args.out_dir)

    print("\n[FIGURES] Done.")


if __name__ == "__main__":
    main()
