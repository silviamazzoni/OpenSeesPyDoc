"""
OpenSeesPy Eigenvalue Analysis: Mass vs Weight Units Pitfall
============================================================
A common mistake when using US customary units (kip-inch-second)
in OpenSees eigenvalue analysis: using weight (kips) instead of
mass (kip·s²/in) at nodes.

OpenSees is unit-agnostic — it won't warn you. The resulting
natural periods will be off by sqrt(g) ≈ 19.7x.

This example demonstrates the error on a realistic 3-span
continuous steel plate girder bridge (100'-150'-100' = 350' total).

Author: Michael Martello, PE
"""

import openseespy.opensees as ops
import math

# ==============================================================
# UNIT DEFINITIONS (US Customary: kip, inch, second)
# ==============================================================
kip = 1.0
inch = 1.0
sec = 1.0
ft = 12.0 * inch
ksi = kip / inch**2
g = 386.4 * inch / sec**2  # gravitational acceleration

# ==============================================================
# BRIDGE GEOMETRY
# ==============================================================
# 3-span continuous: 100' - 150' - 100' (1200" - 1800" - 1200")
spans = [100.0 * ft, 150.0 * ft, 100.0 * ft]
total_length = sum(spans)

# Node spacing for each span (10 segments per span)
n_segments_per_span = 10

# Support locations (abutments + 2 interior piers)
support_locations = [0.0, spans[0], spans[0] + spans[1], total_length]

# ==============================================================
# SECTION PROPERTIES (Typical Steel Plate Girder)
# ==============================================================
# 72" deep built-up steel plate girder (typical highway bridge)
E_steel = 29000.0 * ksi
d_girder = 72.0 * inch          # girder depth
tw = 0.5625 * inch              # web thickness
bf = 18.0 * inch                # flange width
tf = 1.25 * inch                # flange thickness

# Section properties (I-shape approximation)
A_girder = 2 * bf * tf + (d_girder - 2 * tf) * tw  # ~84 in²
Iz_girder = (bf * d_girder**3 / 12.0
             - (bf - tw) * (d_girder - 2 * tf)**3 / 12.0)  # ~97,000 in⁴

# ==============================================================
# DECK AND TRIBUTARY MASS
# ==============================================================
# 8.5" concrete deck, 40' roadway width, 5 girders @ 8' spacing
deck_thickness = 8.5 * inch
deck_width = 40.0 * ft          # tributary to one girder: ~8'
girder_spacing = 8.0 * ft
gamma_concrete = 150.0 / 1000.0 / 1728.0  # pcf -> kip/in³ (8.68e-5)
gamma_steel = 490.0 / 1000.0 / 1728.0     # pcf -> kip/in³ (2.84e-4)

# Weight per unit length (tributary to one girder)
w_deck = gamma_concrete * deck_thickness * girder_spacing   # kip/in
w_girder = gamma_steel * A_girder                           # kip/in
w_barriers = 0.45 / 1000.0 * kip / inch  # 0.45 plf -> kip/in per girder (2 barriers / 5 girders)
w_total = w_deck + w_girder + w_barriers  # kip/in (~0.085 kip/in total)


def build_and_analyze(use_correct_mass=True, label=""):
    """
    Build a 3-span bridge model and run eigenvalue analysis.

    Parameters
    ----------
    use_correct_mass : bool
        If True, divides weight by g to get mass (correct).
        If False, uses weight directly as mass (WRONG).
    label : str
        Label for output.
    """
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 3)

    # ----------------------------------------------------------
    # CREATE NODES
    # ----------------------------------------------------------
    node_id = 1
    node_x = []
    span_start = 0.0
    support_nodes = []

    for span_idx, span_len in enumerate(spans):
        dx = span_len / n_segments_per_span
        for i in range(n_segments_per_span + (1 if span_idx == len(spans) - 1 else 0)):
            x = span_start + i * dx
            ops.node(node_id, x, 0.0)
            node_x.append(x)

            # Apply mass at each node
            trib_length = dx  # tributary length
            if i == 0 or (span_idx == len(spans) - 1 and i == n_segments_per_span):
                trib_length = dx / 2.0  # half tributary at ends

            weight = w_total * trib_length  # weight in kips

            if use_correct_mass:
                # ✅ CORRECT: mass = weight / g
                mass = weight / g
            else:
                # ❌ WRONG: using weight directly as mass
                mass = weight

            ops.mass(node_id, mass, mass, 0.0)

            # Check if this is a support node
            for sup_loc in support_locations:
                if abs(x - sup_loc) < 1.0:  # within 1 inch
                    support_nodes.append(node_id)
                    break

            node_id += 1
        span_start += span_len

    total_nodes = node_id - 1

    # ----------------------------------------------------------
    # BOUNDARY CONDITIONS
    # ----------------------------------------------------------
    # Pin at first abutment, rollers at piers and far abutment
    for i, sn in enumerate(support_nodes):
        if i == 0:
            ops.fix(sn, 1, 1, 0)  # pin (Dx, Dy fixed, Rz free)
        else:
            ops.fix(sn, 0, 1, 0)  # roller (Dy fixed only)

    # ----------------------------------------------------------
    # ELEMENTS (Elastic Beam-Column)
    # ----------------------------------------------------------
    transf_tag = 1
    ops.geomTransf('Linear', transf_tag)

    for i in range(1, total_nodes):
        ops.element('elasticBeamColumn', i, i, i + 1,
                     A_girder, E_steel, Iz_girder, transf_tag)

    # ----------------------------------------------------------
    # EIGENVALUE ANALYSIS
    # ----------------------------------------------------------
    num_modes = 5
    eigenvalues = ops.eigen(num_modes)

    print(f"\n{'='*60}")
    print(f"  {label}")
    print(f"  Mass = {'W/g (correct)' if use_correct_mass else 'W (WRONG — weight used as mass)'}")
    print(f"{'='*60}")
    print(f"  {'Mode':<6} {'ω² (rad²/s²)':<18} {'T (sec)':<12} {'f (Hz)':<10}")
    print(f"  {'-'*46}")

    periods = []
    for i, ev in enumerate(eigenvalues):
        if ev > 0:
            omega = math.sqrt(ev)
            T = 2.0 * math.pi / omega
            f = 1.0 / T
            periods.append(T)
            print(f"  {i+1:<6} {ev:<18.4f} {T:<12.4f} {f:<10.4f}")
        else:
            periods.append(None)
            print(f"  {i+1:<6} {ev:<18.4f} {'N/A':<12} {'N/A':<10}")

    return periods


# ==============================================================
# RUN BOTH ANALYSES
# ==============================================================
print("\n" + "=" * 60)
print("  3-SPAN CONTINUOUS BRIDGE EIGENVALUE ANALYSIS")
print(f"  Spans: {spans[0]/ft:.0f}' - {spans[1]/ft:.0f}' - {spans[2]/ft:.0f}'")
print(f"  Girder: {d_girder:.0f}\" deep plate girder")
print(f"  Deck: {deck_thickness:.1f}\" concrete, {girder_spacing/ft:.0f}' girder spacing")
print("=" * 60)

# ❌ WRONG: weight used as mass
wrong_periods = build_and_analyze(use_correct_mass=False,
                                   label="WRONG — Weight as Mass")

# ✅ CORRECT: mass = weight / g
correct_periods = build_and_analyze(use_correct_mass=True,
                                     label="CORRECT — Mass = Weight / g")

# ==============================================================
# COMPARISON
# ==============================================================
print(f"\n{'='*60}")
print(f"  COMPARISON: Error from using weight as mass")
print(f"{'='*60}")
print(f"  {'Mode':<6} {'T_wrong (s)':<14} {'T_correct (s)':<16} {'Ratio':<10} {'Error':<10}")
print(f"  {'-'*56}")

for i in range(len(correct_periods)):
    if correct_periods[i] and wrong_periods[i]:
        ratio = wrong_periods[i] / correct_periods[i]
        error_pct = (ratio - 1.0) * 100.0
        print(f"  {i+1:<6} {wrong_periods[i]:<14.4f} {correct_periods[i]:<16.4f} {ratio:<10.2f} {error_pct:+.1f}%")

print(f"\n  Expected ratio: sqrt(g) = sqrt(386.4) = {math.sqrt(g):.2f}x")
print(f"  Using weight as mass inflates periods by ~{math.sqrt(g):.1f}x")
print(f"\n  Rule of thumb: For a {spans[1]/ft:.0f}' main span steel bridge,")
print(f"  T₁ should be roughly 0.5–2.0 seconds (vertical).")
print(f"  If you're getting T₁ > 10 seconds, check your mass units!\n")
