import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

from magnetic_field_PSO import *


class Plotting():
    def __init__(self):
        self.simsun = FontProperties(fname=r'C:\Windows\Fonts\simsun.ttc', size=10)  # 宋体
        self.roman = FontProperties(fname=r'C:\Windows\Fonts\times.ttf', size=10)  # Times new roman

    def fig1(self):
        x1 = np.linspace(-2, 3.4, 256, endpoint=True)
        x2 = np.linspace(-2, 3.4, 12, endpoint=True)
        y1 = np.zeros(shape=x1.size)
        y2 = np.zeros(shape=x2.size)
        for i in range(0, x1.size):
            a1 = np.array([x1[i], 4.21, 0])
            y1[i] = Theory(a1).db0()[1]
        for i in range(0, x2.size):
            a2 = np.array([x2[i], 4.21, 0])
            y2[i] = Simulation(a2).db0()[1]
        plt.plot(x1, y1, linestyle="-.", color="black", linewidth=0.5, label="磁法理论值")
        plt.plot(x2, y2, linestyle="none", marker="o", markerfacecolor="none", label="仿真测量值")
        plt.xlim(-2, 3.4)
        plt.xticks(np.arange(-2, 4, 1), rotation=0)
        plt.ylim(-40, 40)
        plt.yticks(np.arange(-40, 60, 20))
        plt.xlabel(r"$x_{A}$/ m")
        plt.ylabel(r"$\Delta$$T$"r"$_{20}$/ nT")
        plt.legend(loc="upper left", frameon=True, prop=self.simsun)
        ax = plt.gca()
        ax.tick_params(which='both', direction="in", top=True, right=True)

    def fig2(self):
        x1 = np.linspace(-2, 3.4, 256, endpoint=True)
        x2 = np.linspace(-2, 3.4, 12, endpoint=True)
        y1 = np.zeros(shape=x1.size)
        y2 = np.zeros(shape=x2.size)
        for i in range(0, x1.size):
            a1 = np.array([x1[i], 4.21, 0])
            y1[i] = Theory(a1).db0()[3]
        for i in range(0, x2.size):
            a2 = np.array([x2[i], 4.21, 0])
            y2[i] = Simulation(a2).db0()[3]
        plt.plot(x1, y1, linestyle="-.", color="black", linewidth=0.5, label="磁法理论值")
        plt.plot(x2, y2, linestyle="none", marker="o", markerfacecolor="none", label="仿真测量值")
        plt.xlim(-2, 3.4)
        plt.xticks(np.arange(-2, 4, 1), rotation=0)
        plt.ylim(10, 70)
        plt.yticks(np.arange(10, 90, 20))
        plt.xlabel(r"$x_{A}$/ m")
        plt.ylabel(r"$\Delta$$T$"r"$_{40}$/ nT")
        plt.legend(loc="upper left", frameon=True, prop=self.simsun)
        ax = plt.gca()
        ax.tick_params(which='both', direction="in", top=True, right=True)

    def fig3(self):
        x1 = np.linspace(-2, 3.4, 256, endpoint=True)
        x2 = np.linspace(-2, 3.4, 12, endpoint=True)
        y1 = np.zeros(shape=x1.size)
        y2 = np.zeros(shape=x2.size)
        for i in range(0, x1.size):
            a1 = np.array([x1[i], 4.21, 0])
            y1[i] = Theory(a1).db0()[5]
        for i in range(0, x2.size):
            a2 = np.array([x2[i], 4.21, 0])
            y2[i] = Simulation(a2).db0()[5]
        plt.plot(x1, y1, linestyle="-.", color="black", linewidth=0.5, label="磁法理论值")
        plt.plot(x2, y2, linestyle="none", marker="o", markerfacecolor="none", label="仿真测量值")
        plt.xlim(-2, 3.4)
        plt.xticks(np.arange(-2, 4, 1), rotation=0)
        plt.ylim(-20, 20)
        plt.yticks(np.arange(-20, 30, 10))
        plt.xlabel(r"$x_{A}$/ m")
        plt.ylabel(r"$\Delta$$T$"r"$_{60}$/ nT")
        plt.legend(loc="upper left", frameon=True, prop=self.simsun, fancybox=False)
        ax = plt.gca()
        ax.tick_params(which='both', direction="in", top=True, right=True)

    def fig4(self):
        size = 10
        x1 = np.linspace(-300, 300, size, endpoint=True)
        y1 = np.linspace(200, 200, size, endpoint=True)
        z1 = np.linspace(0, 0, size, endpoint=True)
        x2 = np.zeros(shape=size)
        y2 = np.zeros(shape=size)
        for i in range(0, size):
            a2 = np.array([x1[i], y1[i], z1[i]])
            print(a2)
            xyz = PSO(a2).solution()
            x2[i] = xyz[0]
            y2[i] = xyz[1]
        plt.plot(x1, y1, linestyle="--", color="black", marker="o", markerfacecolor="none", linewidth=0.5, label="实际位置")
        plt.plot(x2, y2, linestyle="-", color="red", marker="o", markerfacecolor="none", linewidth=0.5, label="磁法定位")
        plt.xlim(-350, 350)
        plt.xticks(np.arange(-350, 400, 100), rotation=0)
        plt.ylim(0, 350)
        plt.yticks(np.arange(0, 400, 50))
        plt.xlabel(r"$x_{A}$/ m")
        plt.ylabel(r"$y_{A}$/ m")
        plt.legend(loc="lower right", frameon=True, prop=self.simsun, fancybox=False)
        ax = plt.gca()
        ax.tick_params(which='both', direction="in", top=True, right=True)


if __name__ == "__main__":
    Plotting().fig4()
    plt.show()
