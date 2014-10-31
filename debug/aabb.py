from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure()
ax = fig.add_subplot(111, projection = "3d")

def draw_aabb(xmin, xmax, ymin, ymax, zmin, zmax):
    ax.plot([xmin, xmax], [ymin, ymin], [zmin, zmin], c="b")
    ax.plot([xmin, xmax], [ymin, ymin], [zmax, zmax], c="b")
    ax.plot([xmin, xmax], [ymax, ymax], [zmin, zmin], c="b")
    ax.plot([xmin, xmax], [ymax, ymax], [zmax, zmax], c="b")

    ax.plot([xmin, xmin], [ymin, ymax], [zmin, zmin], c="b")
    ax.plot([xmin, xmin], [ymin, ymax], [zmax, zmax], c="b")
    ax.plot([xmax, xmax], [ymin, ymax], [zmin, zmin], c="b")
    ax.plot([xmax, xmax], [ymin, ymax], [zmax, zmax], c="b")

    ax.plot([xmin, xmin], [ymin, ymin], [zmin, zmax], c="b")
    ax.plot([xmin, xmin], [ymax, ymax], [zmin, zmax], c="b")
    ax.plot([xmax, xmax], [ymin, ymin], [zmin, zmax], c="b")
    ax.plot([xmax, xmax], [ymax, ymax], [zmin, zmax], c="b")

def draw_ray(start, direction, m=1, c="r"):
    p1 = start
    p2 = [start[0] + m * direction[0], start[1] + m * direction[1], start[2] + m * direction[2]]
    ax.plot([p1[0], p2[0]], [p1[1], p2[1]], [p1[2], p2[2]], c=c)
    ax.scatter([start[0]], [start[1]], [start[2]], c=c)

def Ray(start=None, direction=None):
    draw_ray(start, direction, m=5)

def Sphere(center=[0, 0, 0], radius=1):
    u, v = np.mgrid[0 : 2 * np.pi : 20j, 0 : np.pi : 10j]
    x = center[0] + radius * np.cos(u) * np.sin(v)
    y = center[1] + radius * np.sin(u) * np.sin(v)
    z = -center[2] + radius * np.cos(v)
    ax.plot_wireframe(x, y, z, color="b")

# draw_aabb(-1, 2, 3, 4, -5, -2)
# draw_ray((0, -1, -1), (0, 1.5, -1), m=5, c="r")
# draw_ray((0, -1, -1), (0, 0.5, -1), m=5, c="y")
# draw_ray((0, -1, -1), (0, 2.5, -1), m=5, c="g")
# draw_ray((0, -1, -1), (0, 4, -1), m=5, c="k")
# draw_ray((0, -1, -1), (0, 2, -0.2), m=5, c="r")
# draw_ray((5, 4, 0), (1, 0.2, 1), m=3, c="r")

Sphere(center=[100, 100, 100], radius=10)
Sphere(center=[50, 80, 90], radius=10)
Sphere(center=[10, 80, 75], radius=10)
Sphere(center=[100, 100, -100], radius=10)
Sphere(center=[25, 60, -99], radius=10)
Sphere(center=[87, 23, -50], radius=10)
Sphere(center=[100, -100, 100], radius=10)
Sphere(center=[67, -43, 13], radius=10)
Sphere(center=[87, -12, 90], radius=10)
Sphere(center=[100, -100, -100], radius=10)
Sphere(center=[13, -41, -23], radius=10)
Sphere(center=[12, -12, -34], radius=10)
Sphere(center=[-100, 100, 100], radius=10)
Sphere(center=[-49, 79, 49], radius=10)
Sphere(center=[-53, 89, 28], radius=10)
Sphere(center=[-100, 100, -100], radius=10)
Sphere(center=[-34, 23, -89], radius=10)
Sphere(center=[-34, 95, -34], radius=10)
Sphere(center=[-100, -100, 100], radius=10)
Sphere(center=[-35, -9, 34], radius=10)
Sphere(center=[-34, -56, 98], radius=10)
Sphere(center=[-100, -100, -100], radius=10)
Sphere(center=[-23, -23, -10], radius=10)
Sphere(center=[-68, -10, -12], radius=10)

plt.show()
