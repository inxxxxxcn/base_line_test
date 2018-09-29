import numpy as np


class Simulation(object):
    def __init__(self, x):
        # 目标位置
        self.x = x
        # 磁矩参数
        self.m = 20 * (10 ** 4)
        self.M0 = np.array([np.cos(1.1) * np.cos(-0.18), np.cos(1.1) * np.sin(-0.18), np.sin(1.1)])
        self.M = self.m * self.M0
        # 常数系数
        self.u0 = 4 * np.pi * (10 ** -7) * (10 ** 9)
        self.p = self.u0 / (4 * np.pi)
        # 传感器位置参数
        self.l = np.array([300, 300, 300, 300, 300, 300])
        # self.s = np.array([self.l[0] / 2, self.l[0] + self.l[1] / 2, self.l[0] + self.l[1] + self.l[2] / 2,
        #                    self.l[0] + self.l[1] + self.l[2] + self.l[3] / 2,
        #                    self.l[0] + self.l[1] + self.l[2] + self.l[3] + self.l[4] / 2,
        #                    self.l[0] + self.l[1] + self.l[2] + self.l[3] + self.l[4] + self.l[5] / 2])
        self.L = np.array(
            [[0, 0, 0], [self.l[0], 0, 0], [self.l[0] + self.l[1], 0, 0], [self.l[0] + self.l[1] + self.l[2], 0, 0],
             [self.l[0] + self.l[1] + self.l[2] + self.l[3], 0, 0],
             [self.l[0] + self.l[1] + self.l[2] + self.l[3] + self.l[4], 0, 0],
             [self.l[0] + self.l[1] + self.l[2] + self.l[3] + self.l[4] + self.l[5], 0, 0]])
        # 地磁场参数
        self.U = np.array([np.cos(1.1) * np.cos(-0.18), np.cos(1.1) * np.sin(-0.18), np.sin(1.1)])
        self.be = 0
        # 噪声
        self.bn = 0 * (10 ** -3) * np.random.random(7)
        # 磁力仪测量值
        self.t = self.t()
        # 磁力仪测量值差值
        self.db = self.db()
        # 磁力仪i（i>0）与磁力仪0的差值
        self.db0 = self.db0()

    def t(self):
        R = np.zeros(shape=(7, 3))
        B = np.zeros(shape=(7, 3))
        b = np.zeros(shape=(7, 1))
        r = np.zeros(shape=(7, 1))
        x = self.x.copy()
        M = self.M
        for i in range(0, 7):
            R[i] = x - self.L[i]
            r[i] = (R[i][0] ** 2 + R[i][1] ** 2 + R[i][2] ** 2) ** (1 / 2)
            B[i] = self.p * ((3 * np.dot(R[i], M) / np.power(r[i], 5) * R[i]) - 1 / np.power(r[i], 3) * M)
            b[i] = np.dot(B[i], self.U)
        b = b.reshape(-1)
        t = b + self.be + self.bn
        return t

        # db10,20,30,40,50,60

    def db0(self):
        t = self.t
        db0 = np.zeros(shape=6)
        for i in range(0, 6):
            db0[i] = t[i + 1] - t[0]
        return db0

        # db01,db12,db23,db34,db45,db56

    def db(self):
        t = self.t
        db = np.zeros(shape=6)
        for i in range(0, 6):
            db[i] = t[i] - t[i + 1]
        return db


if __name__ == "__main__":
    a = np.array([200, 200, 0])
    s = Simulation(a)
    print(s.t[0])
