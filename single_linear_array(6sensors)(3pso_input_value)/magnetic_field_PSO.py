from magnetic_field_theory import *


class PSO(object):
    def __init__(self, x0):
        self.x0 = x0
        # 空间维度
        self.dim = 3
        # 粒子群大小
        self.size = 300
        # 最大迭代次数
        self.steps = 300
        # 加速因子
        self.c1 = self.c2 = 2
        # 初始值
        self.x = self.produce_x()
        self.v = np.random.rand(self.size, self.dim)
        fitness = self.fitness()
        self.p = self.x
        self.pg = self.x[np.argmin(fitness)]
        self.individual_best_fitness = fitness
        self.global_best_fitness = np.min(fitness)

    def produce_x(self):
        dim = self.dim
        size = self.size
        x_bound = np.zeros(shape=(dim, 2))
        # 位置范围
        x_bound[0] = [-300, 300]
        x_bound[1] = [0, 200]
        x_bound[2] = [-10, 10]
        # # 磁矩范围
        # x_bound[3] = []
        # x_bound[4] = []
        # x_bound[5] = []
        x = np.zeros(shape=(size, dim))
        for i in range(0, dim):
            yy = np.random.uniform(x_bound[i][0], x_bound[i][1], (size, 1))
            x[:, i] = yy.reshape(-1)
        return x

    def fitness(self):
        x = self.x
        x0 = self.x0
        sub = np.zeros(shape=(self.size, 6))
        for i in range(0, self.size):
            sub[i] = Theory(x[i]).db - Simulation(x0).db
        fitness = np.sum(np.square(sub), axis=1)
        return fitness

    def solution(self):
        for step in range(self.steps):
            r1 = np.random.rand(self.size, self.dim)
            r2 = np.random.rand(self.size, self.dim)
            # 惯性因子
            w = 0.9 - 0.3 / (self.steps - 1) * step
            self.v = w * self.v + self.c1 * r1 * (self.p - self.x) + self.c2 * r2 * (self.pg - self.x)
            self.x = self.v + self.x
            # plt.clf()
            # plt.scatter(self.x[:, 0], self.x[:, 1], s=30, color='k')
            # plt.pause(10)
            fitness = self.fitness()
            update_id = np.greater(self.individual_best_fitness, fitness)
            self.p[update_id] = self.x[update_id]
            self.individual_best_fitness[update_id] = fitness[update_id]
            if np.min(fitness) < self.global_best_fitness:
                self.pg = self.x[np.argmin(fitness)]
                self.global_best_fitness = np.min(fitness)
            print(self.pg)
            print(self.global_best_fitness)
        print("\n")
        print(self.global_best_fitness)
        print("\n")
        return self.pg


if __name__ == "__main__":
    a = np.array([-300, 200, 0])
    print(PSO(a).solution())
    # plt.show()

