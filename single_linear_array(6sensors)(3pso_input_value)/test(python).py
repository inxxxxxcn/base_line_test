import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from sympy import *
from magnetic_field_simulation import *
from magnetic_field_theory import *


class test(object):
    def __init__(self, p):
        self.l1 = 1

    def first(self):
        print(self.p)

    def second(self,x):
        self.p=2
        self.first()
        return self.p


if __name__ == "__main__":
    # test(2).second(2)
    # a = np.dot(np.array([[3],[2],[1]]),np.array([[1,2,3]]))
    # print(a)
    # def f(x):
    #     f = np.array([2*x[1],3*x[2],4*x[2]])
    #     return f
    #
    # def gb0(x):
    #     dx = 0.1
    #     gb0 = np.eye(3)
    #     for i in range(0, 3):
    #         xx = x
    #         xx[i] = xx[i] + dx
    #         print(xx,x)
    #         gb0[i] = (f(xx)[i]-f(x)[i])/dx
    #     return gb0
    #
    # print(gb0(np.array([1,2,3])))
    # x = range(-2, 3)
    # y = x
    # plt.plot(x,y)
    # print(np.arange(1,2,0.1))
    # def B(x):
    #     l = 0.5
    #     S = np.array([[0, 0, 0], [l, 0, 0], [-l, 0, 0], [0, l, 0], [0, -l, 0], [0, 0, l], [0, 0, -l]])
    #     R = np.zeros(shape=(7, 3))
    #     C = np.zeros(shape=(7, 3))
    #     T = x
    #     # for i in range(0, 7):
    #     #     R[i] = T - S[i]
    #     C=T-S[1]
    #     R[1]=C
    # #     print(R)
    #
    # #
    # x1,x2,x3=symbols("x1 x2 x3")
    # x=np.array([x1-1,x2,x3])
    #
    # z=np.zeros(shape=(3,3),dtype=Symbol)
    # z[1]=x-
    # print(z)

    # a = np.arange(-2, 3.5, 0.1)
    # # (-2, 3.4+0.1, 0.1)
    # b = np.zeros(shape=(6, 54))
    # i = 0
    # for x in a:
    #     t = np.array([x, 4.21, 0])
    #     s = Simulation(t)
    #     b[:][i] = s.db0()
    #     i += 1
    # plt.plot(a, b)
    # plt.show()

    # def y_calculate(x):
    #     a = [x, 4.21, 0]
    #     s = Simulation(a)
    #     y0 = s.db0()[1]/(s.l*2)
    #     return y0
    #
    # y = np.zeros(shape=54)
    # x = np.arange(-2, 3.4, 0.1)
    # for i in range(0, 54):
    #     y[i] = y_calculate(x[i])
    # print(y)
    # plt.plot(x, y)
    # # plt.show()
    # a=np.array([3.4,4.21,0])
    # print(Theory(a).gb0())
    # print(Simulation(a).gb0())
    # # print(Theory(a).gb0()-Simulation(a).gb0())
    # fig = plt.figure()  # an empty figure with no axes
    # ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    # # fig, ax_lst = plt.subplots(2, 2)  # a figure with a 2x2 grid of Axes
    # plt.show()
    # 绘制普通图像
    # from matplotlib.font_manager import FontProperties
    #
    # roman = FontProperties(fname=r'C:\Windows\Fonts\times.ttf', size=10)  # Times new roman
    # x = np.linspace(-1, 1, 50)
    # y1 = 2 * x + 1
    # y2 = x ** 2
    #
    # plt.figure()
    # # 在绘制时设置lable, 逗号是必须的
    # l1, = plt.plot(x, y1, label='line')
    # l2, = plt.plot(x, y2, label=r'$parabola$', color='red', linestyle="--", marker="o")
    # print(type(l1,))
    #
    # # 设置坐标轴的取值范围
    # plt.xlim((-1, 1))
    # plt.ylim((0, 2))
    #
    # # 设置坐标轴的lable
    # plt.xlabel("X axis", fontproperties=roman)
    # plt.ylabel('Y axis')
    #
    # # 设置x坐标轴刻度, 原来为0.25, 修改后为0.5
    # plt.xticks(np.linspace(-1, 1, 5))
    # # 设置y坐标轴刻度及标签, $$是设置字体
    # plt.yticks([0, 0.5], ['$minimum$', 'normal'])
    #
    # # 设置legend
    # plt.legend(handles=[l1, l2, ], loc='best', fancybox=False)
    # plt.show()

    # x1 = np.linspace(-50, 50, 50, endpoint=True)
    # x2 = np.linspace(-50, 50, 50, endpoint=True)
    # x3 = np.linspace(0, 0, 50, endpoint=True)
    # a = np.array([[1, 1, 1]])
    # for i in x2:
    #     x02 = np.linspace(i, i, 101, endpoint=True)
    #     a0 = np.array([x1, x02, x3]).swapaxes(1, 0)
    #     a = np.concatenate((a, a0), axis=0)
    # s = np.zeros(shape=(10202, 6))
    # for i in range(10202):
    #     s[i] = Simulation(a[i]).gb()
    # print(s)
    a=np.array([[1,2,3],[1,2,3]])
    print(a.reshape(2,3))