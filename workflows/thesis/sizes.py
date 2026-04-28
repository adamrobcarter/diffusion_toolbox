import matplotlib.pyplot as plt
import numpy as np
import workflows.thesis.common

fig, ax = plt.subplots(figsize=(workflows.thesis.common.textwidth, 4))

# ── Data ──────────────────────────────────────────────────────────────────────
# Each entry: (label, x_start, x_end)
particles = [
    # ("Paint pigment",         0.01,        1.0,         False),
    ("Colloids",              0.001,       10,          True),
    (r"H\textsubscript{2}O",  0.28e-3*0.9, 0.28e-3*1.1, False),
    ("DNA",                   2e-3*0.9,       2e-3*1.1,     False),
    ("Viruses",               0.02,       0.3,         True),
    ("Bacteria",              0.5,         5,          True),
    ("Human hair",            20,          180,         True),
    # ("Mist",                  5,           100,         False),
    ("Pollen",                  5,           100,         True),
    ("Red blood cells",       2,           6,          True),
    # ("Atomic radii",          0.0001,      0.0006,      False),
    ("Wheat flour",          20,           100,         True), # https://www.researchgate.net/figure/Particle-size-distribution-profile-of-flour-with-different-particle-sizes-Results-are-a_fig2_335277332
    ("Milk fats",          0.2,           15,         True),
    # ("Sand",                  40,          200,         True),
]

# ── Layout: assign each particle a y-row ──────────────────────────────────────
# Pack rows tightly by assigning to the first row where there's no overlap.
rows = []  # list of (x_end_last) per row
row_assignments = []

for label, x0, x1, fuzzy in particles:
    placed = False
    for i, row_end in enumerate(rows):
        if label != particles[0][0] and i == 0:
            continue # nothing else on colloids row
        if x0 >= row_end:          # small gap between bars
            rows[i] = x1
            row_assignments.append(i)
            placed = True
            break
    if not placed: # new row
        rows.append(x1)
        row_assignments.append(len(rows) - 1)

n_rows = len(rows)

# ── Draw bars ─────────────────────────────────────────────────────────────────
bar_height = 0.8
colors = plt.cm.tab20.colors          # 20 distinct colours

def draw_fuzzy_bar(ax, x0, x1, y, height, color, n_steps=30, fuzz_frac=0.08):
    """Draw a horizontal bar with faded left and right edges on a log axis."""
    # Solid centre
    log0, log1 = np.log10(x0), np.log10(x1)
    fuzz = (log1 - log0) * fuzz_frac          # fade zone width in log-space

    solid_x0 = 10 ** (log0 + fuzz)
    solid_x1 = 10 ** (log1 - fuzz)

    ax.barh(y, solid_x1 - solid_x0, left=solid_x0,
            height=height, color=color, linewidth=0, alpha=0.85)

    # Left and right fade strips
    for side in ("left", "right"):
        if side == "left":
            edge_log0, edge_log1 = log0, log0 + fuzz
        else:
            edge_log0, edge_log1 = log1 - fuzz, log1

        steps = np.linspace(edge_log0, edge_log1, n_steps + 1)

        for i in range(n_steps):
            sx0 = 10 ** steps[i]
            sx1 = 10 ** steps[i + 1]
            # Alpha ramps from 0 at the outer edge to 0.85 at the inner edge
            if side == "left":
                alpha = 0.85 * (i / n_steps)
            else:
                alpha = 0.85 * (1 - i / n_steps)

            ax.barh(y, sx1 - sx0, left=sx0,
                    height=height, color=color, linewidth=0, alpha=alpha)

for idx, ((label, x0, x1, fuzzy), row) in enumerate(zip(particles, row_assignments)):
    y = n_rows - row                   # flip so row-0 is at top
    color = colors[idx % len(colors)]

    # ax.barh(
    #     y, width=x1 - x0, left=x0,
    #     height=bar_height,
    #     color=color, edgecolor="white", linewidth=0.5, alpha=0.7
    # )
    draw_fuzzy_bar(ax, x0, x1, y, bar_height, color, fuzz_frac=0.08 if fuzzy else 0)

    # Label position: centre of bar if it fits, otherwise to the right
    x_centre = np.sqrt(x0 * x1)       # geometric centre (log scale)
    ax.text(
        x_centre, y, label,
        ha="center", va="center",
        fontsize=8, fontweight="bold", color="black",
        clip_on=True
    )

# # ── Range bands at the bottom ─────────────────────────────────────────────────
# range_bands = [
#     ("IONIC\nRANGE",           0.0001, 0.001,  "#b0c4de"),
#     ("MOLECULAR\nRANGE",       0.001,  0.01,   "#c8dfc8"),
#     ("MACRO MOLECULAR\nRANGE", 0.01,   0.1,    "#f0e68c"),
#     ("MICRO PARTICLE\nRANGE",  0.1,    1,      "#f4a460"),
#     ("MACRO PARTICLE\nRANGE",  1,      100,    "#cd853f"),
# ]

# band_y = 0.2          # just below the bars
# band_h = 0.5

# for label, x0, x1, color in range_bands:
#     ax.barh(band_y, width=x1 - x0, left=x0,
#             height=band_h, color=color, edgecolor="gray", linewidth=0.5)
#     x_centre = np.sqrt(x0 * x1)
#     ax.text(x_centre, band_y, label,
#             ha="center", va="center", fontsize=6, fontweight="bold")

# ── Axes formatting ───────────────────────────────────────────────────────────
ax.set_xscale("log")
# ax.set_xlim(0.0001, 100)
ax.set_ylim(0, n_rows + 1)

ax.set_xlabel("diameter")
ax.set_yticks([])                      # hide y-axis ticks / labels
ax.xaxis.set_ticks_position("top")
ax.xaxis.set_label_position("top")

# Nicer tick labels on the log axis
ax.set_xticks([1e-4, 0.001, 0.01, 0.1, 1, 10, 100])
ax.set_xticklabels([r"100\si{\pico\meter}", r"1\si{\nano\meter}", r"10\si{\nano\meter}", r"100\si{\nano\meter}", r"1\si{\micro\meter}", r"10\si{\micro\meter}", r"100\si{\micro\meter}"])


ax.spines["left"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["bottom"].set_visible(False)

ax.grid(axis="x", which="major", linestyle="--", linewidth=0.4,
        color="gray", alpha=0.4)

plt.tight_layout()
plt.savefig("workflows/thesis/figures/sizes.pdf",
            dpi=150, bbox_inches="tight")
print("Saved → workflows/thesis/figures/sizes.pdf")
