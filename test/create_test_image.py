r"""
Create images to visualize the semi integration 
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.integrate import cumulative_trapezoid
from ec_tools.semi_integration import semi_integration


# define values as gaussian curve
x = np.linspace(0, 8, 1000)
y = norm.pdf(x, 4, 1)

# semi integration
d_semi = semi_integration(y, x, alg="r1", transonic_backend="numba")

dd_semi = semi_integration(
    semi_integration(y, x, alg="r1", transonic_backend="numba"),

    x,
    alg="r1",
    transonic_backend="numba",
)

# full (numerical) integration
d_full = cumulative_trapezoid(y, x, initial=0)

# Peak plot
# fig = plt.figure()
fig, ax = plt.subplots(1)
plt.plot(x, y, linewidth=5.0)
plt.xlabel("x", fontsize=18)
plt.ylabel("y", fontsize=18)
ax.set_xticks([])
ax.set_yticks([])
plt.savefig("data/images/show_peak.png", dpi=300)

# Semidiff plot
# fig = plt.figure()
fig, ax = plt.subplots(1)
plt.plot(x, d_semi, "r", linewidth=5.0)
plt.xlabel("x", fontsize=18)
plt.ylabel("y", fontsize=18)
ax.set_xticks([])
ax.set_yticks([])
plt.savefig("data/images/show_semidiff.png", dpi=300)

# Step plot
# fig = plt.figure()
fig, ax = plt.subplots(1)
plt.plot(x, d_full, linewidth=5.0)
plt.xlabel("x", fontsize=18)
plt.ylabel("y", fontsize=18)
ax.set_xticks([])
ax.set_yticks([])
plt.savefig("data/images/show_step.png", dpi=300)
