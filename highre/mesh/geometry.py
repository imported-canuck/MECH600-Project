import numpy as np

def naca_airfoil(x):
    return 0.594689181 * (0.298222773 * np.sqrt(x) - 0.127125232 * x - 0.357907906 * x**2 +
                        0.291984971 * x**3 - 0.105174606 * x**4)

def surface_points(n=100, x0 = 0, x1 = 1):
    theta = np.linspace(np.arccos(1 - 2 * x0), np.arccos(1 - 2 * x1), n)
    x = 0.5 * (1 - np.cos(theta))
    y = naca_airfoil(x)
    upper = np.column_stack([x, y])
    lower = np.column_stack([x, -y])
    return upper, lower

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    up, lo = surface_points()
    print("half-thickness at x=0.3:", naca_airfoil(0.3))
    print("half-thickness at x=1.0:", naca_airfoil(1.0))
    plt.plot(up[:, 0], up[:, 1], ".-", ms=2)
    plt.plot(lo[:, 0], lo[:, 1], ".-", ms=2)
    plt.axis("equal")
    plt.title("TMR NACA 0012")
    plt.savefig("airfoil.png", dpi=150)
    print("saved airfoil.png")
