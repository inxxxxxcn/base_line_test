import numpy as np
from magnetic_field_simulation import Simulation


class Theory(object):
    def __init__(self, x):
        # 目标位置
        self.x = x
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
        # self.be = 54000
        # # 噪声
        # self.bn = 1 * (10 ** -3) * np.random.random(7)
        # 磁矩
        self.M = self.M()
        # 磁偶极子磁感应强度
        self.B = self.B()
        # 磁异常
        self.b = self.b()
        # 磁异常差值
        self.db = self.db()

    def M(self):
        # 磁偶极子磁场B=KM
        x = self.x.copy()
        C = np.zeros(shape=(3, 3))
        K = np.zeros(shape=(7, 3, 3))
        dbt = np.zeros(shape=3)
        R = np.zeros(shape=(7, 3))
        r = np.zeros(shape=(7, 1))
        E = np.eye(3)
        sdb0 = Simulation(x).db0
        for i in range(0, 7):
            R[i] = x - self.L[i]
            r[i] = (R[i][0] ** 2 + R[i][1] ** 2 + R[i][2] ** 2) ** (1 / 2)
            K[i] = 3 * np.dot(R[i].reshape(3, 1), np.array([R[i]])) / np.power(r[i], 5) - E / np.power(r[i], 3)
        for i in range(0, 3):
            C[i] = self.p * np.dot(self.U, K[i + 1] - K[0])
            dbt[i] = sdb0[i]
        C = np.mat(C).I.getA()
        M = np.dot(C, dbt)
        return M

    def B(self):
        R = np.zeros(shape=(7, 3))
        B = np.zeros(shape=(7, 3))
        r = np.zeros(shape=(7, 1))
        K = np.zeros(shape=(7, 3, 3))
        E = np.eye(3)
        x = self.x.copy()
        # 确定位置参数K
        for i in range(0, 7):
            R[i] = x - self.L[i]
            r[i] = (R[i][0] ** 2 + R[i][1] ** 2 + R[i][2] ** 2) ** (1 / 2)
            K[i] = 3 * np.dot(R[i].reshape(3, 1), np.array([R[i]])) / np.power(r[i], 5) - E / np.power(r[i], 3)
            B[i] = self.p * np.dot(self.M, K[i])
        return B

    def b(self):
        b = np.zeros(shape=(7, 1))
        # 磁异常b=UB
        for i in range(0, 7):
            b[i] = np.dot(self.U, self.B[i])
        b = b.reshape(-1)
        return b

    def db(self):
        db = np.zeros(shape=6)
        for i in range(0, 6):
            db[i] = self.b[i] - self.b[i + 1]
        return db


if __name__ == "__main__":
    a = np.array([200, 200, 0])
    t = Theory(a)
    print(t.db)
