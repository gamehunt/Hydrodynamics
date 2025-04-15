import math
import tabulate
import numpy as np
import matplotlib.pyplot as plt

# Метод Сопряженных Градиентов

def mu(i, j, grid):
    x = i * grid.h
    y = j * grid.k

    if i == 0 and j >= grid.m // 4:
        return y ** 2

    if i == grid.n:
        return y ** 2 + 4

    if j == grid.m // 4 and i <= grid.n // 2:
        return x ** 2 + 0.0625

    if j == grid.m:
        return x ** 2 + 1

    if i == grid.n // 2 and j < grid.m // 4:
        return y ** 2 + 1

    if j == 0 and i >= grid.n // 2:
        return x ** 2

    if j < grid.m // 4 and i < grid.n // 2:
        return np.nan

    return 0

def f(i, j, grid):
    return -4

def u(x, y):
    return x**2 + y**2

class Grid:
    def __init__(self, w, h, n, m, f, mu, u = None):
        self.width = w
        self.height = h
        self.n = n
        self.m = m
        self.h = w / n
        self.k = h / m
        self.k_sq = self.k ** 2
        self.h_star = 1 / (self.h ** 2)
        self.k_star = 1 / self.k_sq
        self.a_star = -2 * (self.h_star + self.k_star)
        self.grid = []
        self.f = f
        self.u = u

        self.accuracy = []

        self.x = np.linspace(0, self.width,  self.n + 1)
        self.y = np.linspace(0, self.height, self.m + 1)

        self.X, self.Y  = np.meshgrid(self.x, self.y)

        if u:
            self.uv = np.array(u(np.ravel(self.X), np.ravel(self.Y)))
        
        for j in range(m + 1):
            self.grid.append([0] * (n + 1))
            for i in range(n + 1):
                self.set(i, j, mu(i, j, self))

    def __str__(self):
        return self.grid.__str__()

    def get(self, i, j):
        return self.grid[j][i]

    def set(self, i, j, v):
        self.grid[j][i] = v

    def __update(self, i, j):
        prev = self.get(i, j)
        left  = self.h_star * self.get(i - 1, j)
        right = self.h_star * self.get(i + 1, j)
        up    = self.k_star * self.get(i, j + 1)
        down  = self.k_star * self.get(i, j - 1)
        f     = self.f(i, j, self)
        new   = (f + left + right + up + down) / -self.a_star
        self.set(i, j, new)
        return math.fabs(prev - new)

    def __opt_update(self, i, j):
        prev  = self.get(i, j)
        left  = self.get(i - 1, j)
        right = self.get(i + 1, j)
        up    = self.get(i, j + 1)
        down  = self.get(i, j - 1)
        f     = self.f(i, j, self)
        new   = 0.4 * (f * self.k_sq + 0.25 * (left + right) + (up + down))
        self.set(i, j, new)
        return math.fabs(prev - new)

    def solve(self, precision, max_iterations, optimized = False):
        self.accuracy = []

        max_delta = precision + 1
        iterations = 0

        update_function = self.__update
        if optimized:
            update_function = self.__opt_update

        while iterations < max_iterations and max_delta > precision:
            max_delta = -1
            for j in range(1, self.m // 4 + 1):
                for i in range(self.n // 2 + 1, self.n):
                    d = update_function(i, j)
                    if d > max_delta:
                        max_delta = d
            for j in range(self.m // 4 + 1, self.m):
                for i in range(1, self.n):
                    d = update_function(i, j)
                    if d > max_delta:
                        max_delta = d
            iterations += 1
            self.accuracy.append(max_delta)
            print(f'i = {iterations}, z = {max_delta}')

        error = None
        if self.u:
            z     = np.ravel(self.grid)
            diffs = np.absolute(np.subtract(z, self.uv))
            error     = max(diffs[~np.isnan(diffs)])

        return iterations, max_delta, error

    def plot(self):
        fig  = plt.figure(1)

        ax   = fig.add_subplot(1, 1, 1)

        z    = np.ravel(self.grid)
        Z    = z.reshape(self.X.shape)

        cs = ax.contourf(self.X, self.Y, Z)     
        fig.colorbar(cs)

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.xaxis.set_major_locator(plt.MultipleLocator(self.h))
        ax.yaxis.set_major_locator(plt.MultipleLocator(self.k))

        plt.legend()
        plt.grid()


        rect = plt.Rectangle((0, 0), 1, 0.25, linewidth=1, edgecolor='none', facecolor='black')
        ax.add_patch(rect)

        plt.show()

if __name__ == '__main__':
    grid = Grid(2, 1, 20, 20, f, mu, u)
    max_iterations = 1000
    iterations, acc, error = grid.solve(1e-14, max_iterations, True)
    print(tabulate.tabulate([['Количество итераций', f'{iterations}/{max_iterations}'], 
                             ['Точность', acc], 
                             ['Погрешность', error]], 
                            tablefmt='simple_grid', colalign=('left','right')))
    # print(grid)
    grid.plot()
