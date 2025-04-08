import math
import tabulate
import numpy as np
import matplotlib.pyplot as plt
import pyopencl as cl

# Метод сопряженных градиентов

def mu(i, j, h, k):
    x = i * h
    y = j * k
    if i == 0:
        return y
    if j == 4:
        return 2
    return 0

if __name__ == '__main__':
    width = 4
    height = 2
    m = 4
    n = 8
    h = width / n
    k = height / m
    precision = 1e-14
    max_iterations = 1000
    k_sq = k ** 2
    h_star = 1 / (h ** 2)
    k_star = 1 / k_sq
    a_star = -2 * (h_star + k_star)
    grid = [] 
    bad_nodes = [(6, 1), (7, 2)]
    invalid_nodes = [(7, 1)]

    for i in range(m + 1):
        grid.append([0] * (n + 1))
        for j in range(n + 1):
            grid[i][j] = mu(j, i, h, k)

    max_delta = precision + 1
    iterations = 0
    while iterations < max_iterations and max_delta > precision:
        max_delta = -1
        for j in range(1, m):
            for i in range(1, n):
                if (i, j) in invalid_nodes:
                    continue
                left  = h_star * grid[j][i - 1]
                right = h_star * grid[j][i + 1]
                up    = k_star * grid[j + 1][i]
                down  = k_star * grid[j - 1][i]
                f     = 0
                old   = grid[j][i]
                if (i, j) in bad_nodes:
                    if i == 6:
                        xw = math.sqrt(1 - (j * k) ** 2) + 4
                        x  = i * h
                        alpha = (xw - x) / h
                        a_star_corrected = -2 * (k_star + h_star / alpha )
                        grid[j][i] = (up + down + 2 * (left / (1 + alpha) + right / (alpha * (1 + alpha)))) / -a_star_corrected
                    else:
                        yw = math.sqrt(1 - (i * h - 4) ** 2)
                        y  = j * k
                        beta = (y - yw) / k
                        a_star_corrected = -2 * (h_star + k_star / beta)
                        grid[j][i] = (left + right + 2 * (up / (1 + beta) + down / (beta * (1 + beta)))) / -a_star_corrected
                else:
                    grid[j][i] = (f + left + right + up + down) / -a_star
                d = math.fabs(grid[j][i] - old)
                if d > max_delta:
                    max_delta = d
        iterations += 1

    x = np.linspace(0, width,  n + 1)
    y = np.linspace(0, height, m + 1)

    X, Y  = np.meshgrid(x, y)

    fig  = plt.figure(1)
    ax   = fig.add_subplot(1, 1, 1, projection='3d')

    z    = np.ravel(grid)
    Z    = z.reshape(X.shape)

    ax.plot_wireframe(X, Y, Z, color='b', rstride = 1, cstride = 1, label='~u(x, y)')     
    ax.scatter(X, Y, Z, color = 'red')

    plt.legend()
    plt.show()
