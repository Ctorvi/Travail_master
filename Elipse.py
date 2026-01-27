# python
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch

# ellipse parameters
cx, cy = 0.0, 0.0   # centre
a, b = 0.2, 0.1     # semi-major (a) and semi-minor (b)

# parametric ellipse
t = np.linspace(0, 2*np.pi, 400)
x = cx + a * np.cos(t)
y = cy + b * np.sin(t)

# choose point on ellipse where arrow will be tangent (angle in radians)
theta = 1.2  # change this to move the tangent point
px = cx + a * np.cos(theta)
py = cy + b * np.sin(theta)

# tangent vector (derivative of parametric form)
dx_dt = -a * np.sin(theta)
dy_dt =  b * np.cos(theta)
v = np.array([dx_dt, dy_dt], dtype=float)
v_norm = np.linalg.norm(v)
if v_norm == 0:
    v_unit = np.array([1.0, 0.0])
else:
    v_unit = v / v_norm

# arrow length in data units

fig, ax = plt.subplots(figsize=(6,6))
ax.plot(x, y, color='blue', lw=2)                       # ellipse
ax.scatter([cx], [cy],
           s=220,            # marker size
           facecolor='cyan',
           edgecolor='black',
           linewidth=1.5,
           zorder=30)   # center point



# small length inside the ellipse for the arrow shaft; adjust factor for desired head size
small_len = 0.03 * max(a, b)
# place start slightly inside the ellipse along the inward tangent direction
posA = (px - v_unit[0] * small_len, py - v_unit[1] * small_len)
# place end at the ellipse point (tip)
posB = (px, py)

arrow = FancyArrowPatch(
    posA=posA, posB=posB,
    arrowstyle='-|>',               # simple head
    mutation_scale=30,              # head size (adjust)
    linewidth=0,                    # hide shaft line
    color='blue',
    zorder=50,
    shrinkA=0, shrinkB=0
)
ax.add_patch(arrow)

# cosmetics
ax.set_aspect('equal', 'box')
pad = 0.2
ax.set_xlim(cx - a - pad, cx + a + pad)
ax.set_ylim(cy - b - pad, cy + b + pad)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_title('Ellipse with center and tangent arrow')

plt.show()