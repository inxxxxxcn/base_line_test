from magnetic_field_strength import *
import matplotlib.pyplot as plt


class base_line_x:
    def __init__(self):
        t = sym.Symbol("t")
        x0 = sym.Symbol("x0")
        # 目标运动参数
        self.target = Target(-300, 200, 0, 10, 0, 0)
        target = self.target
        self.point = Measured_Point(x0, 0, 0, target.position(t)[0], target.position(t)[1], target.position(t)[2],
                                    20 * (10 ** 4), 1.1, -0.18, 1.1, -0.18, 000, 0 * (10 ** -3))

    def ba(self, t_input, x0_input):
        t = sym.Symbol("t")
        x0 = sym.Symbol("x0")
        return self.point.ba().subs(t, t_input).subs(x0, x0_input)

    def dbat(self, t1_input, t2_input):
        t = sym.Symbol("t")
        x0 = sym.Symbol("x0")
        return self.point.ba().subs(t, t1_input) - self.point.ba().subs(t, t2_input)

    def gbax(self, x1_input, x2_input, t1_input, t2_input):
        t = sym.Symbol("t")
        x0 = sym.Symbol("x0")
        x2 = sym.Symbol("x2")
        t2 = sym.Symbol("t2")
        dbat = self.dbat(t1_input, t2_input)
        return sym.diff((dbat.subs(x0, x1_input) - dbat.subs(x0, x2)), x2).subs(x2, x2_input)


if __name__ == "__main__":
    print(base_line_x().gbax(0, 1, 50, 0))
    size = 200
    x = np.linspace(0, 1000, size)
    y = np.zeros(size)
    for i in range(size):
        y[i] = base_line_x().gbax(0, x[i], 50, 0)
    plt.plot(x,y)
    plt.show()

