import numpy as np
import math

from geometry import naca_airfoil, surface_points
from grading import expansion_ratio, growth_ratio

'''levels = 0, 1, 2'''

XS = 0.3
delta = 0.1
phi = math.radians(25)
base_cells = dict(front = 40, aft = 60, wake = 50, inner = 70, outer = 50)
ds_le = 0.002
ds_te = 0.002

def distance(a, b):
    return math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2)



def normal(x, sign):
    if x >= XS - 0.0000001:
        return np.array([0.0,sign])
    if x <= 0: 
        return np.array([-1.0,0.0])
    h = 0.0000001
    slope = (naca_airfoil(x + h) - naca_airfoil(x - h)) / (2 * h)
    n = np.array([-slope, 1.0]) / math.sqrt((slope ** 2) + 1)
    return np.array([n[0], sign * n[1]])

def offset(x0, x1, sign , n=60):
    upper, lower = surface_points(n, x0, x1)
    if sign > 0:
        pts = upper
    else:
        pts = lower

    return np.array([p + delta * normal(p[0], sign) for p in pts])

def curve_length(points):
    return float(np.sum(np.hypot(*np.diff(points, axis = 0).T)))

def make_mesh(level, dy1_coarse, radius, out_file):
    r = 1.5 ** level

    n_front = round(base_cells["front"] * r)
    n_aft   = round(base_cells["aft"] * r)
    n_wake  = round(base_cells["wake"] * r)
    n_inner = round(base_cells["inner"] * r)
    n_outer = round(base_cells["outer"] * r)

    dy1 = dy1_coarse / r
    le_size = ds_le / r
    te_size = ds_te / r

    ys = naca_airfoil(XS)
    y_top = radius * math.cos(phi)
    x_a = XS - radius * math.sin(phi)
    x_t = 1 + y_top * math.tan(phi)

    pts = {}
    pts["LE"] = (0, 0) #leading edge
    pts["Su"] = (XS, ys) #upper surface at the thickest point of the airfoil
    pts["TE"] = (1, 0) #trailing edge of the airfoil
    pts["Sl"] = (XS, -ys) #lower surface at the thickest point of the airfoil
    pts["LEi"] = (-delta, 0) #delta in front of the leading edge
    pts["Sui"] = (XS, ys + delta) #delta straight above Su
    pts["Sli"] = (XS, -ys - delta) # delta straight below Sl
    pts["TEui"] = (1, delta) #delta straight above the trailing edge
    pts["TEli"] = (1, -delta) #delta straight below the trailing edge
    pts["Wm"] = (1 + radius, 0) #end of the domain at the wake line
    pts["Wui"] = (1 + radius, delta) #delta above Wm
    pts["Wli"] = (1 + radius, -delta) #delta below Wm
    pts["Af"] = (XS - radius, 0) #radius in front of the thickest airfoil point
    pts["At"] = (x_a, y_top) #top end of arc in front of airfoil
    pts["Ab"] = (x_a, -y_top) #bottom end of arc in front of airfoil
    pts["Tt"] = (x_t, y_top) #top far corner of line rising from Te tilted at phi
    pts["Tb"] = (x_t, -y_top) #same as Tt but bottom
    pts["Wt"] = (1 + radius, y_top) #top right corner of domain
    pts["Wb"] = (1 + radius, -y_top)  #bottom right corner of domain

    idx = {}
    number = 0
    for name in pts:
        idx[name] = number
        number = number + 1
    nv = len(pts)

    upper_front, lower_front = surface_points(2000, 0, XS)
    upper_aft, lower_aft = surface_points(2000, XS, 1)
    L_front = curve_length(upper_front)
    L_aft = curve_length(upper_aft)
    g_front = expansion_ratio(le_size, L_front, n_front)
    g_aft = 1 / expansion_ratio(te_size, L_aft, n_aft)
    g_wake = expansion_ratio(te_size, radius, n_wake)
    g_arc = g_front ** 0.3
    top_cell = te_size * distance(pts["At"], pts["Tt"]) / L_aft
    g_farwake = expansion_ratio(top_cell, distance(pts["Tt"], pts["Wt"]), n_wake)
    q_in = growth_ratio(dy1, delta, n_inner)
    g_in = q_in ** (n_inner - 1)
    first_out = dy1 * q_in ** n_inner
    go_le = expansion_ratio(first_out, distance(pts["LEi"], pts["Af"]), n_outer)
    go_s  = expansion_ratio(first_out, distance(pts["Sui"], pts["At"]), n_outer)
    go_te = expansion_ratio(first_out, distance(pts["TEui"], pts["Tt"]), n_outer)
    go_w  = expansion_ratio(first_out, distance(pts["Wui"], pts["Wt"]), n_outer)

    

if __name__ == "__main__":
    make_mesh(level=0, dy1_coarse=3e-6, radius=100, out_file="blockMeshDict")


