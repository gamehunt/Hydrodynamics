from math import sin, pi
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation

# Начальные условия и генератор cетки 
def create_grid(length, time_slice, n, m, u0, mu):
    grid = []

    h = length / n
    t = time_slice / m

    for j in range(m):
        grid.append([0] * n)

    for i in range(n):
        grid[0][i] = u0(i * h)

    for j in range(m):
        grid[j][0] = mu(j * t)

    return grid, h, t

def u0_0(x):
    return 1.0 if 1 <= x <= 2 else 0.0

def u0_1(x):
    if x < 1 or x > 2:
        return 0
    elif x >= 1 and x < 1.5:
        return 2 * (x - 1)
    else:
        return 1 - 2 * (x - 1.5)

def u0_2(x):
    if x < 1 or x > 2:
        return 0
    else:
        return 0.5 * (1 + sin(2 * pi * (x - 1) - pi / 2))

def u0_3(x):
    if x < -1 or x > 1:
        return 0
    else:
        return 1 - x**2  # функция для половины параболы
    
def mu3(t):
    global a
    return u0_3(a * t)

def mu(t):
    return 0

# Решение разностной схемой
#     ()
#     |
# * - *
def solve1(grid, a, h, tau):
    for t in range(1, len(grid)):
        for x in range(1, len(grid[t])):
            grid[t][x] = grid[t - 1][x] - a * tau * ((grid[t - 1][x] - grid[t - 1][x - 1]) / h)
    return grid


# Решение разностной схемой
#     ()
#     |
# * - * - *
def solve2(grid, a, h, tau):
    for t in range(1, len(grid)):
        for x in range(1, len(grid[t]) - 1):
            grid[t][x] = grid[t - 1][x] - a * tau * ((grid[t - 1][x + 1] - grid[t - 1][x - 1]) / (2 * h))
    return grid

# Итерация прогонки
def solve_single(t, grid, a, h, tau):
    prev = grid[t - 1]
    curr = grid[t]

    ti   = 1 / tau

    mu1  = curr[0]
    mu2  = 0

    kappa1 = 0
    kappa2 = 0

    ai = [kappa1]
    bi = [mu1]

    C = -ti
    A = -a / (4 * h)
    B =  a / (4 * h)

    for x in range(1, len(curr) - 1):
        phi  = prev[x] * ti - a / (4 * h) * (prev[x + 1] - prev[x - 1])
        alpha = B / (C - A * ai[x - 1])
        beta  = (-phi + A * bi[x - 1]) / (C - A * ai[x - 1])
        ai.append(alpha)
        bi.append(beta)

    curr[-1] = (kappa2 * bi[-1] + mu2) / (1 - kappa2 * ai[-1])
    for x in range(len(curr) - 2, -1, -1):
        curr[x] = ai[x] * curr[x + 1] + bi[x]


# Решение разностной схемой
# * - () - *
#     |
# * - * - *
def solve3(grid, a, h, tau):
    for t in range(1, len(grid)):
        solve_single(t, grid, a, h, tau)
    return grid

# Генератор для точного решения
def generate_precise(u0, a, i, h, n, tau):
    grid = []
    for x in range(n):
        grid.append(u0(x * h - i * tau * a))
    return grid

a = 2
w = 15 # x
h = 10 # t
n = w * 10
m = h * 20

targets = [u0_0, u0_1, u0_2, u0_3]
solvers = [solve1, solve2, solve3]

# tau / h < 0.5

if __name__ == '__main__':
    s = int(input('u0: '))

    # Выбор н.у
    if s < 1 or s > len(targets):
        print('u0 should be 1, 2, 3 or 4')
        exit(1)

    target = None
    _mu = mu

    # Коррекция mu для последнего случая
    target = targets[s - 1]
    if s == 4:
        _mu = mu3

    grid, hx, tau = create_grid(w, h, n, m, target, _mu)

    # Выбор метода
    method = int(input('Method: '))
    if method < 1 or method > 3:
        print('Method should be 1, 2 or 3')
        exit(1)

    # Проверка условия Куррента (3 метод сходится абсолютно - его исключаем)
    if method != 3:
        c = abs(a) * tau / hx
        if c > 1:
            print('Doesnt converge: c = ', c)
            print(tau, hx)
            exit(1)

    grid_solved = solvers[method - 1](grid, a, hx, tau)

    fig, ax = plt.subplots()
    x = np.linspace(0, w, n)
    ax.xaxis.set_major_locator(plt.MultipleLocator(1))

    # Отрисовка конкретного фрейма
    def update(frame):
        plt.cla()
        plt.title(f"t = {tau * frame:.2f}")
        numeric = grid_solved[frame]
        precise = generate_precise(target, a, frame, hx, n, tau)
        ax.plot(x, numeric, '-', label = 'num. method')
        ax.plot(x, precise, '-', label = 'precise')
        plt.legend()

        eps = np.absolute(np.subtract(numeric, precise)).max()
        print(f"t = {frame * tau}, eps (Погрешность) = {eps}")
        print(numeric)
        print(precise)

    # Ничего интересного - если ввели пустое время, то рисуем мультик, иначе 
    # рисуем фрейм заданного времени
    fr = input('Time: ')
    anim = len(fr) == 0
    if not anim:
        fr = int(float(fr) / tau)

    if anim:
        ani = animation.FuncAnimation(fig, update, frames=m, interval=30, repeat = False)
    else:
        update(fr)

    plt.show()
