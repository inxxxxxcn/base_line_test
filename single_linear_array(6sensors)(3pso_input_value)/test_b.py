import numpy as np
import sympy as sym


def ba(x0, y0, z0, x1, y1, z1, m, I1, D1, I0, D0, be, bn):
    # 目标到传感器的距离
    x, y, z = x1 - x0, y1 - y0, z1 - z0
    # 系数
    u0 = 4 * np.pi * (10 ** -7) * (10 ** 9)
    p = u0 / (4 * np.pi)
    # 磁矩大小和方向
    mx, my, mz = m * sym.cos(I1) * sym.cos(D1), m * sym.cos(I1) * sym.sin(D1), m * sym.sin(I1)
    # 地磁方向
    # I0, D0 = 1.1, -0.18
    u, v, w = sym.cos(I0) * sym.cos(D0), sym.cos(I0) * sym.sin(D0), sym.sin(I0)
    # 参数变形
    a = 2 * mx * u - my * v - mz * w
    b = 2 * my * v - mx * u - mz * w
    c = 2 * mz * w - mx * u - my * v
    d = 3 * (my * u + mx * v)
    e = 3 * (mz * u + mx * w)
    f = 3 * (mz * v + my * w)
    # g = j * u + my * v + l * w
    r = (x ** 2 + y ** 2 + z ** 2) ** (1 / 2)
    # ba = p * (a * x ** 2 + b * y ** 2 + c * z ** 2 + d * x * y + e * x * z + f * y * z) * (
    #         x ** 2 + y ** 2 + z ** 2) ** (-5 / 2)
    # ba = p* (a * x ** 2 + b * y ** 2 + c * z ** 2 + d * x * y + e * x * z + f * y * z) * r ** (-5)
    # ba = p * r ** (-3) * (a * x ** 2 + b * y ** 2 + c * z ** 2 + d * x * y + e * x * z + f * y * z) * r ** (-2)
    btr = (a * x ** 2 + b * y ** 2 + c * z ** 2 + d * x * y + e * x * z + f * y * z)
    # 磁异常
    ba = p * r ** (-5) * btr
    return ba+be+bn*np.random.random(1)


if __name__ == "__main__":
    print(ba(0, 0, 0, 200, 200, 0, 20 * (10 ** 4), 1.1, -0.18, 1.1, -0.18, 000, 10 * (10 ** -3)))
