#!/bin/bash
#
# analyse_md.sh — post-production analysis for all systems.
#
#   bash analyse_md.sh              # everything in systems.conf
#   bash analyse_md.sh fold_013     # one run
#
# Environment:
#   SKIP_NS   ns of equilibration to discard before averaging (default 20)
#   REGIONS   path to a file of "NAME START END" lines defining regions.
#             Defaults to $PIPELINE_DIR/regions.conf so region boundaries are
#             defined ONCE and cannot drift between scripts.
#
# -----------------------------------------------------------------------------
# CRITICAL: PBC correction.
# The first version of this script ran gmx rms directly on prod.xtc. Over 100 ns
# a trimer diffuses across the periodic boundary and protomers get wrapped to
# opposite sides of the box, producing large fake jumps in RMSD and Rg. Every
# trajectory must be made whole and centred BEFORE any metric is computed.
#
# CRITICAL: part files.
# Runs spanning several walltime windows write prod.xtc + prod.partNNNN.xtc.
# This script prefers prod_full.xtc (see tools/concat_trajectories.sh). Without
# it you would analyse 60 ns of one system against 100 ns of another.
# -----------------------------------------------------------------------------
set -uo pipefail

BASE="${BASE_DIR:-$HOME/md_runs}"
PIPELINE_DIR="${PIPELINE_DIR:-$HOME/md_pipeline}"
OUT="$BASE/analysis"
CONF="${CONF:-$PIPELINE_DIR/systems.conf}"
REGIONS="${REGIONS:-$PIPELINE_DIR/regions.conf}"
SKIP_NS="${SKIP_NS:-20}"

mkdir -p "$OUT"

if [ $# -ge 1 ]; then
    SYSTEMS=("$@")
else
    SYSTEMS=($(awk '!/^#/ && NF {print $1}' "$CONF"))
fi

echo "=================================================================="
echo "  Analysing ${#SYSTEMS[@]} run(s); discarding first ${SKIP_NS} ns"
echo "=================================================================="

for S in "${SYSTEMS[@]}"; do
    d="$BASE/$S"
    [ -d "$d" ] || { echo "[SKIP] $S — no run directory"; continue; }
    cd "$d" || continue

    # prefer the joined trajectory
    if   [ -s prod_full.xtc ]; then TRJ=prod_full.xtc
    elif [ -s prod.xtc ];      then TRJ=prod.xtc
        if ls prod.part*.xtc >/dev/null 2>&1; then
            echo "[WARN] $S has part files but no prod_full.xtc."
            echo "       Run tools/concat_trajectories.sh first, or you are"
            echo "       analysing only the FIRST segment of this run."
        fi
    else
        echo "[SKIP] $S — no trajectory"; continue
    fi

    echo ""
    echo "[$S]  using $TRJ"

    # ---------- PBC correction ----------
    if [ ! -s fit.xtc ] || [ "$TRJ" -nt fit.xtc ]; then
        echo "  PBC correction (make whole, centre, remove tumbling)..."
        echo "Protein System" | gmx trjconv -s prod.tpr -f "$TRJ" \
            -o whole.xtc -pbc whole -center >/dev/null 2>&1
        echo "Backbone System" | gmx trjconv -s prod.tpr -f whole.xtc \
            -o fit.xtc -fit rot+trans >/dev/null 2>&1
        rm -f whole.xtc
    fi
    [ -s fit.xtc ] || { echo "  [ERROR] PBC correction failed"; continue; }

    B=$((SKIP_NS * 1000))     # ps

    # RMSD — did the fold hold? Should plateau.
    echo "Backbone Backbone" | gmx rms -s prod.tpr -f fit.xtc \
        -o "$OUT/rmsd_${S}.xvg" -tu ns >/dev/null 2>&1

    # RMSF — per-residue flexibility. THE key metric. Computed after
    # discarding equilibration; RMSF over a drifting start is meaningless.
    echo "C-alpha" | gmx rmsf -s prod.tpr -f fit.xtc \
        -o "$OUT/rmsf_${S}.xvg" -res -b $B >/dev/null 2>&1

    # Radius of gyration — is the trimer swelling apart?
    echo "Protein" | gmx gyrate -s prod.tpr -f fit.xtc \
        -o "$OUT/rg_${S}.xvg" >/dev/null 2>&1

    # SASA — is the graft still presented on the surface?
    echo "Protein" | gmx sasa -s prod.tpr -f fit.xtc \
        -o "$OUT/sasa_${S}.xvg" -b $B >/dev/null 2>&1 || true

    # Inter-protomer H-bonds — direct probe of interface integrity.
    echo "Protein Protein" | gmx hbond -s prod.tpr -f fit.xtc \
        -num "$OUT/hbond_${S}.xvg" -b $B >/dev/null 2>&1 || true

    echo "  done"
done

echo ""
echo "Generating figures..."

SKIP_NS="$SKIP_NS" OUT="$OUT" CONF="$CONF" REGIONS="$REGIONS" python3 << 'PYEOF'
import os
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = os.environ["OUT"]
CONF = os.environ["CONF"]
REGIONS_FILE = os.environ["REGIONS"]

SYSTEMS = [l.split()[0] for l in open(CONF)
           if l.strip() and not l.startswith('#')]
# wild type first: it is the reference, drawn black and thick
SYSTEMS = ([s for s in SYSTEMS if 'wild' in s] +
           [s for s in SYSTEMS if 'wild' not in s])

# Region boundaries: read from one file so they cannot drift between scripts.
REGIONS = {}
if os.path.exists(REGIONS_FILE):
    for line in open(REGIONS_FILE):
        line = line.strip()
        if line and not line.startswith('#'):
            p = line.split()
            if len(p) >= 3:
                REGIONS[p[0]] = (int(p[1]), int(p[2]))
if not REGIONS:
    print("  [note] no regions.conf found — figures will have no region shading")

PALETTE = ['#1b7a5a', '#c1852b', '#4aa38c', '#8a5e17', '#6fbba6', '#2e8b74']
COLORS = {s: ('#2b2b2b' if 'wild' in s else PALETTE[i % len(PALETTE)])
          for i, s in enumerate(SYSTEMS)}

def load_xvg(fn):
    x, y = [], []
    for line in open(fn):
        if line.startswith(('#', '@')):
            continue
        p = line.split()
        if len(p) >= 2:
            x.append(float(p[0])); y.append(float(p[1]))
    return np.array(x), np.array(y)

plt.rcParams.update({'font.size': 10, 'axes.spines.top': False,
                     'axes.spines.right': False, 'figure.dpi': 150,
                     'savefig.dpi': 150, 'savefig.bbox': 'tight'})

def style(s):
    return dict(color=COLORS[s], lw=2.5 if 'wild' in s else 1.4,
                alpha=1.0 if 'wild' in s else 0.85, label=s)

summary = {}

# ---- FIG 1: RMSD ----
fig, ax = plt.subplots(figsize=(11, 4.5))
for s in SYSTEMS:
    fn = f"{OUT}/rmsd_{s}.xvg"
    if not os.path.exists(fn):
        continue
    t, r = load_xvg(fn)
    ax.plot(t, r * 10, **style(s))
    half = r[len(r)//2:] * 10
    summary.setdefault(s, {})['RMSD_mean'] = half.mean()
    if len(t) > 4:
        summary[s]['RMSD_drift'] = np.polyfit(t[len(t)//2:], half, 1)[0]
ax.set_xlabel("Time (ns)"); ax.set_ylabel("Backbone RMSD (Å)")
ax.set_title("Backbone RMSD — did each system settle?", fontweight='bold')
ax.legend(fontsize=9, frameon=False, ncol=3)
plt.savefig(f"{OUT}/fig1_rmsd.png"); plt.close()

# ---- FIG 2: RMSF with regions (key plot) ----
fig, ax = plt.subplots(figsize=(14, 5))
for s in SYSTEMS:
    fn = f"{OUT}/rmsf_{s}.xvg"
    if not os.path.exists(fn):
        continue
    res, f = load_xvg(fn)
    ax.plot(res, f * 10, **style(s))
ymax = ax.get_ylim()[1]
for name, (a, b) in REGIONS.items():
    ax.axvspan(a, b, alpha=0.12, color='#c1852b')
    ax.text((a + b) / 2, ymax * 0.95, name, ha='center', fontsize=7.5,
            color='#8a5e17', fontweight='bold')
ax.set_xlabel("Residue number"); ax.set_ylabel("Cα RMSF (Å)")
ax.set_title("Per-residue flexibility — engineered regions shaded (KEY PLOT)",
             fontweight='bold')
ax.legend(fontsize=8.5, frameon=False, ncol=3, loc='upper right')
plt.savefig(f"{OUT}/fig2_rmsf.png"); plt.close()

# ---- FIG 3: Rg ----
fig, ax = plt.subplots(figsize=(11, 4.5))
for s in SYSTEMS:
    fn = f"{OUT}/rg_{s}.xvg"
    if not os.path.exists(fn):
        continue
    t, g = load_xvg(fn)
    ax.plot(t / 1000, g * 10, **style(s))
    summary.setdefault(s, {})['Rg_mean'] = (g[len(g)//2:] * 10).mean()
ax.set_xlabel("Time (ns)"); ax.set_ylabel("Radius of gyration (Å)")
ax.set_title("Trimer compactness — flat is healthy", fontweight='bold')
ax.legend(fontsize=9, frameon=False, ncol=3)
plt.savefig(f"{OUT}/fig3_rg.png"); plt.close()

# ---- FIG 4: per-region mean RMSF (decision figure) ----
wt = next((s for s in SYSTEMS if 'wild' in s), None)
if wt and REGIONS and os.path.exists(f"{OUT}/rmsf_{wt}.xvg"):
    names = list(REGIONS)
    variants = [s for s in SYSTEMS
                if s != wt and os.path.exists(f"{OUT}/rmsf_{s}.xvg")]
    width = 0.8 / (len(variants) + 1)
    xs = np.arange(len(names))

    def region_mean(res, f, a, b):
        m = (res >= a) & (res <= b)
        return (f[m] * 10).mean() if m.any() else np.nan

    fig, ax = plt.subplots(figsize=(12, 4.5))
    r0, f0 = load_xvg(f"{OUT}/rmsf_{wt}.xvg")
    ax.bar(xs, [region_mean(r0, f0, *REGIONS[n]) for n in names],
           width, color='#2b2b2b', label=wt)
    for i, s in enumerate(variants):
        r, f = load_xvg(f"{OUT}/rmsf_{s}.xvg")
        ax.bar(xs + width * (i + 1),
               [region_mean(r, f, *REGIONS[n]) for n in names],
               width, color=COLORS[s], label=s, alpha=0.9)
    ax.set_xticks(xs + width * len(variants) / 2)
    ax.set_xticklabels(names)
    ax.set_ylabel("Mean Cα RMSF (Å)")
    ax.set_title("Per-region flexibility vs reference — the decision figure",
                 fontweight='bold')
    ax.legend(fontsize=8.5, frameon=False, ncol=3)
    plt.savefig(f"{OUT}/fig4_region_bars.png"); plt.close()

# ---- text summary ----
with open(f"{OUT}/summary.txt", "w") as fh:
    fh.write(f"{'run':<16}{'RMSD mean (A)':>15}{'RMSD drift (A/ns)':>20}"
             f"{'Rg mean (A)':>14}\n")
    fh.write("-" * 65 + "\n")
    for s in SYSTEMS:
        d = summary.get(s, {})
        fh.write(f"{s:<16}{d.get('RMSD_mean', float('nan')):>15.2f}"
                 f"{d.get('RMSD_drift', float('nan')):>20.4f}"
                 f"{d.get('Rg_mean', float('nan')):>14.2f}\n")
    fh.write("\nHow to read this:\n")
    fh.write("  RMSD drift near zero = equilibrated. Persistently positive =\n")
    fh.write("  still moving, so the run is not converged.\n")
    fh.write("  Rg within ~1-2% of the reference = trimer intact.\n")
    fh.write("  The pass/fail call comes from per-region RMSF (fig4).\n")

print(open(f"{OUT}/summary.txt").read())
print("Figures written to", OUT)
for f in ["fig1_rmsd.png", "fig2_rmsf.png", "fig3_rg.png", "fig4_region_bars.png"]:
    print("  ", f)
PYEOF

echo "=================================================================="
echo "  Analysis complete: $OUT"
echo "=================================================================="
