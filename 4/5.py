import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def gauss_seidel(u, f, h, bc_func, max_iter=10000, tol=1e-14):
    u = u.copy()
    for iteration in range(max_iter):
        u_old = u.copy()
        for i in range(1, u.shape[0]-1):
            for j in range(1, u.shape[1]-1):
                u[i,j] = 0.25 * (u[i+1,j] + u[i-1,j] + u[i,j+1] + u[i,j-1] - h**2 * f[i,j])
        u = bc_func(u, h)
        e = np.abs(u - u_old)
        if np.max(e) < tol:
            print(f"Количество итераций: {iteration+1}")
            print(f"Точность: {np.max(e)}")
            break
    return u

def problem5_neumann_bc():
    Nx, Ny = 20, 20
    x = np.linspace(0, 1, Nx)
    y = np.linspace(0, 1, Ny)
    h = x[1] - x[0]
    X, Y = np.meshgrid(x, y, indexing='ij')
    f = np.ones_like(X)
    u0 = np.zeros_like(f)

    def bc_func(u, h):
        u[0, :] = y**2  # mu1(y)
        u[:, 0] = x**2  # mu3(x)
        u[:, -1] = x**2 + 1  # mu4(x)
        u[-1, :] = u[-2, :] + 2 * h  # du/dx = 2x (Neumann) Т.к x = 1, то сводится к этому. Я хз правильно или нет))))
        return u

    u0 = bc_func(u0, h)
    return x, y, u0, f, h, bc_func

def plot_solutions(x, y, solutions, titles):
    fig = plt.figure()

    for idx, (u, title) in enumerate(zip(solutions, titles)):
        X, Y = np.meshgrid(x, y, indexing='ij')


        ax = fig.add_subplot(1, 1, 1)
        contour = ax.contourf(X, Y, u, levels=20, cmap='viridis')
        fig.colorbar(contour, ax=ax, shrink=0.5, aspect=10)
        #ax.set_title(f'Top View: {title}', fontsize=12)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.grid(True, linestyle='--', alpha=0.5)

    plt.tight_layout()
    plt.show()

def solve_problem5_all_methods():

    x, y, u0, f, h, bc_func = problem5_neumann_bc()


    methods = {
        'Метод Зейделя': gauss_seidel,
    }

    solutions = []
    titles = []

    for name, method in methods.items():
        u = u0.copy()
        u = method(u, f, h, bc_func)
        solutions.append(u)
        titles.append(f"Метод: {name}")

    plot_solutions(x, y, solutions, titles)

if __name__ == '__main__':
    solve_problem5_all_methods()
