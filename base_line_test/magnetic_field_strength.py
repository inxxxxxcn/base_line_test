# -*- coding: utf-8 -*-
import numpy as np
import sympy as sym


class Target:
    def __init__(self, xt0, yt0, zt0, xt, yt, zt, t, m, I1, D1):
        self.xt0, self.yt0, self.zt0 = xt0, yt0, zt0
        self.xt, self.yt, self.zt = xt, yt, zt
        self.t = t
        self.m, self.I1, self.D1 = m, I1, D1


class Measured_Point(Target):
    def __init__(self, x0, y0, z0, x1, y1, z1, m, I1, D1, I0, D0, be, bn):
        # 测量点位置
        self.x0, self.y0, self.z0 = x0, y0, z0
        # 目标位置
        self.x1, self.y1, self.z1 = x1, y1, z1
        # 目标磁矩，磁偏角，磁倾角
        self.m, self.I1, self.D1 = m, I1, D1
        # 测量点地磁偏角，地磁倾角，地磁噪声及误差
        self.I0, self.D0, self.be, self.bn = I0, D0, be, bn

    def ba(self):
        x0, y0, z0 = self.x0, self.y0, self.z0
        x1, y1, z1 = self.x1, self.y1, self.z1
        m, I1, D1 = self.m, self.I1, self.D1
        I0, D0, be, bn = self.I0, self.D0, self.be, self.bn
        # 目标到传感器的距离
        x, y, z = x1 - x0, y1 - y0, z1 - z0
        # 系数
        u0 = 4 * np.pi * (10 ** -7) * (10 ** 9)
        p = u0 / (4 * np.pi)
        # 磁矩分量
        mx, my, mz = m * sym.cos(I1) * sym.cos(D1), m * sym.cos(I1) * sym.sin(D1), m * sym.sin(I1)
        # 地磁强度单位分量
        u, v, w = sym.cos(I0) * sym.cos(D0), sym.cos(I0) * sym.sin(D0), sym.sin(I0)
        # 参数变形
        a = 2 * mx * u - my * v - mz * w
        b = 2 * my * v - mx * u - mz * w
        c = 2 * mz * w - mx * u - my * v
        d = 3 * (my * u + mx * v)
        e = 3 * (mz * u + mx * w)
        f = 3 * (mz * v + my * w)
        r = (x ** 2 + y ** 2 + z ** 2) ** (1 / 2)
        # ba = p * (a * x ** 2 + b * y ** 2 + c * z ** 2 + d * x * y + e * x * z + f * y * z) * (
        #         x ** 2 + y ** 2 + z ** 2) ** (-5 / 2)
        # ba = p* (a * x ** 2 + b * y ** 2 + c * z ** 2 + d * x * y + e * x * z + f * y * z) * r ** (-5)
        # ba = p * r ** (-3) * (a * x ** 2 + b * y ** 2 + c * z ** 2 + d * x * y + e * x * z + f * y * z) * r ** (-2)
        btr = (a * x ** 2 + b * y ** 2 + c * z ** 2 + d * x * y + e * x * z + f * y * z)
        # 磁异常
        ba = p * r ** (-5) * btr
        return ba + be + bn * np.random.random(1)


if __name__ == "__main__":
    # print(Measured_Point(0, 0, 0, 200, 200, 0, 20 * (10 ** 4), 1.1, -0.18, 1.1, -0.18, 000, 10 * (10 ** -3)))
    t = sym.Symbol("t")
    x0 = sym.Symbol("x0")
    target = Target(-300, 200, 0, 10, 0, 0)
    point = Measured_Point(x0, 0, 0, target.position(t)[0], target.position(t)[1], target.position(t)[2],
                           20 * (10 ** 4), 1.1, -0.18, 1.1, -0.18, 000, 0 * (10 ** -3))
    print(point.ba().subs(t, 50).subs(x0, 0))
