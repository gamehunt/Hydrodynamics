import numpy as np
import math
import matplotlib.pyplot as plt

def u(x):
    return 10 + 90 * (x ** 2)

def phi(x):
    return -2110 + 450 * (x ** 2)

def solve(n):
    mu1 = 10
    mu2 = 100

    kappa1 = kappa2 = 0

    h = 1 / n

    A = 12 / (h ** 2)
    B = 12 / (h ** 2)
    C = 24 / (h ** 2) + 5

    ai = [kappa1]
    bi = [mu1]

    xi = [0]

    for i in range(1, n):
        _xi = i * h
        xi.append(_xi)
        a = B / (C - A * ai[i - 1])
        b = (phi(_xi) + A * bi[i - 1]) / (C - A * ai[i-1])

        ai.append(a)
        bi.append(b)

    xi.append(n * h)

    yn = (kappa2 * bi[n - 1] + mu2) / (1 - kappa2 * ai[n - 1])
    yi = [yn]

    for i in range(n - 1, -1, -1):
        yi.append(ai[i] * yi[n - i - 1] + bi[i])

    yi = yi[::-1]

    return xi, yi

def solve_opt1(n):
    mu1 = 10
    mu2 = 100

    kappa1 = kappa2 = 0

    h = 1 / n

    hs = h ** 2 / 12

    A = 12 / (h ** 2)
    B = 12 / (h ** 2)
    C = 24 / (h ** 2) + 5

    ai = [kappa1]
    bi = [mu1]

    xi = [0]

    for i in range(1, n):
        _xi = i * h
        xi.append(_xi)
        a = 1 / ((hs * 5 - ai[i - 1]) + 2)
        # a = B / (C - A * ai[i - 1])
        b = (hs * phi(_xi) + bi[i - 1]) * a
        ai.append(a)
        bi.append(b)

    xi.append(n * h)

    yn = mu2
    yi = [yn]

    for i in range(n - 1, -1, -1):
        yi.append(ai[i] * yi[n - i - 1] + bi[i])

    yi = yi[::-1]

    return xi, yi

def make_plot(i, xs, ys, title):
    plt.subplot(3, 1, i)
    plt.title(title)

    # plt.xticks(np.arange(0, 1, 1 / n))
    plt.yticks(np.arange(0, 110, 10))

    plt.grid()
    plt.plot(np.array(xs), np.array(ys))

def z(xs, ys):
    d = []
    for i in range(len(xs)):
        d.append(math.fabs(u(xs[i]) - ys[i]))
    return max(d)

if __name__ == '__main__':
    s = input('Enter grid size: ')

    if len(s) > 0:
        n = int(s)
        
        uxs = np.linspace(0, 1, n)
        uys = u(uxs)

        vxs, vys = solve(n)
        print(f'z = {z(vxs, vys)}')

        vxs_opt, vys_opt = solve_opt1(n)
        print(f'z_opt = {z(vxs_opt, vys_opt)}')

        make_plot(1, uxs, uys, 'U(x)')
        make_plot(2, vxs, vys, 'V(x)')
        make_plot(3, vxs_opt, vys_opt, 'V(x) (Opt 1)')

        plt.show()
    else:
        n = 10
        ns = []
        zs = []
        zs_opt = []
        while n <= 10000000:
            vxs, vys = solve(n)
            vxs_opt, vys_opt = solve_opt1(n)
            zv = z(vxs, vys)
            zv_opt = z(vxs_opt, vys_opt)
            print(f'n = {n}, z = {zv}, z_opt = {zv_opt}')
            ns.append(n)
            zs.append(zv)
            zs_opt.append(zv_opt)
            n = n * 10

        plt.yscale('log')
        plt.plot(np.array(ns), np.array(zs), label = 'original')
        plt.plot(np.array(ns), np.array(zs_opt), label = 'optimized')
        plt.legend(loc="upper left")
        plt.grid()
        plt.show()




