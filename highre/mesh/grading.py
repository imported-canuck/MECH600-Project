import math

def first_cell_height(y_plus, Re, chord = 1.0):
    cf = 0.026 / Re**(1/7) 
    u_tau = math.sqrt(cf / 2)
    nu = 1.0 / Re
    y = y_plus * nu / u_tau
    return (2 * y * chord)

def growth_ratio(first, length, n):
    lo, hi = 1.000000000000001, 2.0

    for _ in range(200):              
        q = 0.5 * (lo + hi)
        total = first * (q**n - 1) / (q - 1)
        if total < length:
            lo = q
        else:
            hi = q
    return q

def expansion_ratio(first, length, n):
    q = growth_ratio(first, length, n)
    return q**(n - 1)

if __name__ == "__main__":
    Re = 6e6
    d1 = first_cell_height(y_plus=1.0, Re=Re)
    print(f"first cell for y+ = 1:  {d1:.2e} chords")
    d1_wf = first_cell_height(y_plus=50, Re=Re)
    print(f"first cell for y+ = 50: {d1_wf:.2e} chords")
    q = growth_ratio(d1, 100.0, 150)
    print(f"150 cells from wall to 100 chords: growth {q:.4f}, "
          f"expansion ratio {expansion_ratio(d1, 100.0, 150):.3e}")