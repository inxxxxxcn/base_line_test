import matplotlib.pyplot as plt
import numpy as np
import sympy as sym


def b(x, y, z):
    u0 = 4 * np.pi * (10 ** -7) * (10 ** 9)
    p = u0 / (4 * np.pi)
    # M
    j, k, l = 152 * np.cos(1.13) * np.cos(-0.18), 152 * np.cos(1.13) * np.sin(-0.18), 152 * np.sin(1.13)
    # U
    u, v, w = np.cos(1.1) * np.cos(-0.18), np.cos(1.1) * np.sin(-0.18), np.sin(1.1)
    # # 传感器位置
    # a, b, c = 0, 0, 0
    ba = p * ((2 * j * u - k * v - l * w) * x ** 2 + (2 * k * v - j * u - l * w) * y ** 2 + (
            2 * l * w - j * u - k * v) * z ** 2 + 3 * (k * u + j * v) * x * y + 3 * (
                      l * u + j * w) * x * z + 3 * (l * v + k * w) * y * z) / (
                 x ** 2 + y ** 2 + z ** 2) ** (5 / 2)
    return ba


def gbx0(x0, x, y, z):
    u0 = 4 * np.pi * (10 ** -7) * (10 ** 9)
    p = u0 / (4 * np.pi)
    # M
    j, k, l = 152 * np.cos(1.13) * np.cos(-0.18), 152 * np.cos(1.13) * np.sin(-0.18), 152 * np.sin(1.13)
    # U
    u, v, w = np.cos(1.1) * np.cos(-0.18), np.cos(1.1) * np.sin(-0.18), np.sin(1.1)
    a = 2 * j * u - k * v - l * w
    b = 2 * k * v - j * u - l * w
    c = 2 * l * w - j * u - k * v
    d = 3 * (k * u + j * v)
    e = 3 * (l * u + j * w)
    f = 3 * (l * v + k * w)
    # ba = p * (a * x ** 2 + b * y ** 2 + c * z ** 2 + d * x * y + e * x * z + f * y * z) * (x ** 2 + y ** 2 + z ** 2) ** (
    #         -5 / 2)
    gb0 = p * ((2 * a * (x - x0) + d * y + e * z) * ((x - x0) ** 2 + y ** 2 + z ** 2) ** (-5 / 2) - 5 * (
            a * (x - x0) ** 3 + b * (x - x0) * y ** 2 + c * (x - x0) * z ** 2 + d * (x - x0) ** 2 * y + e * (
            x - x0) ** 2 * z + f * (x - x0) * y * z) * ((x - x0) ** 2 + y ** 2 + z ** 2) ** (-7 / 2))
    return gb0


if __name__ == "__main__":
    x, y, z = sym.symbols("x, y, z")
    b = b(x, y, z)
    gb = sym.diff(b, x)
    ggb = sym.diff(gb, x)
    size = 500
    x0, y0, z0 = np.linspace(-100, 100, size), np.linspace(0, 0, size), np.linspace(0, 0, size)
    x_step = np.linspace(-100, 100, size)
    # Gp = np.zeros(shape=size)
    # for i in range(size):
    #     Gp[i] = ggb.subs(x, x0[i]).subs(y, y0[i]).subs(z, z0[i])
    # plt.plot(x0, Gp)
    # plt.show()
    print(b.subs(x, 10).subs(y, 10).subs(z, 10))

    # x, y, z, x0 = sym.symbols("x, y, z, x0")
    # gb = sym.diff(b(x, y, z), x)
    # gb0 = gbx0(x0, x, y, z)
    # ggb0 = sym.diff(gb0, x0)
    # gggb0 = sym.diff(ggb0, x0)
    # # print(gb0.subs(x0,0).subs(x, 34).subs(y, 42).subs(z, 0))
    # xt, yt, zt = np.linspace(0, 0, 101), np.linspace(42.1, 42.1, 101), np.linspace(0, 0, 101)
    # x0t = np.linspace(-100, 100, 101)
    # Gp = np.zeros(shape=101)
    # for i in range(101):
    #     Gp[i] = gb0.subs(x0, x0t[i]).subs(x, xt[i]).subs(y, yt[i]).subs(z, zt[i])
    # plt.plot(x0t, Gp)
    # plt.show()
    # print(gb.subs(x, 6.4).subs(y, -19.7).subs(z, 4.9))
    # print(gbx0(0, x, y, z).subs(x, -13).subs(y, -4).subs(z, -10))
