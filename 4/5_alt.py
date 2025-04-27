import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def solve_poisson_gauss_seidel():
    # Параметры сетки
    nx = 20  # количество узлов по x
    ny = 20  # количество узлов по y
    dx = 1.0 / (nx - 1)
    dy = 1.0 / (ny - 1)

    h_sq = 1 / dx**2
    k_sq = 1 / dy**2
    a_star  = 2 * (h_sq + k_sq)
    
    # Инициализация сетки
    u = np.zeros((ny, nx))
    x = np.linspace(0, 1, nx)
    y = np.linspace(0, 1, ny)
    
    # Установка граничных условий
    # Левая граница: u(0, y) = y^2 (Дирихле)
    u[:, 0] = y**2
    
    # Правая граница: du/dn = 2x (Нейман)
    # Используем одностороннюю разность: (u[:, -1] - u[:, -2])/dx = 2x | (x = 1) = 2
    u[:, -1] = u[:, -2] + 2 * dx
    
    # Верхняя граница: u(x, 1) = x^2 + 1 (Дирихле)
    u[-1, :] = x**2 + 1
    
    # Нижняя граница: u(x, 0) = x^2 (Дирихле)
    u[0, :] = x**2

    # Правая часть уравнения Пуассона (f(x,y) = 4)
    f = 4 * np.ones((ny, nx))
    
    # Параметры итерационного процесса
    max_iter = 10000
    tolerance = 1e-14
    
    # Итерационный процесс (метод Гаусса-Зейделя без релаксации)
    for iteration in range(max_iter):
        u_old = u.copy()
        
        for j in range(1, ny-1):
            for i in range(1, nx-1):
                u[j, i] = ((u[j, i+1] + u[j, i-1]) * h_sq + (u[j+1, i] + u[j-1, i]) * k_sq - f[j, i]) / a_star
        
        # Обновление граничного условия Неймана на правой границе
        u[:, -1] = u[:, -2] + 2 * dx 
        
        # Проверка сходимости
        residual = np.linalg.norm(u - u_old)
        if residual < tolerance:
            print(f"Сходимость достигнута на итерации {iteration}")
            print(f"Точность: {residual}")
            break
    
    # Точное решение
    exact = np.zeros((ny, nx))
    for j in range(ny):
        for i in range(nx):
            exact[j, i] = x[i]**2 + y[j]**2
    
    error = np.abs(u - exact).max()
    print(f"Погрешность: {error}")
    
    # Визуализация
    X, Y = np.meshgrid(x, y)
    
    fig = plt.figure()#figsize=(12, 5))
    
    ax1 = fig.add_subplot(111)
    ax1.set_box_aspect(1.0)
    cm = ax1.contourf(X, Y, u, levels = 20)
    ax1.set_title('Численное решение (Гаусс-Зейдель)')
    
    # ax2 = fig.add_subplot(122)
    # cm = ax2.contourf(X, Y, exact)
    # ax2.set_title('Точное решение')

    fig.colorbar(cm)
    
    plt.grid()
    plt.tight_layout()
    plt.show()
    
    return u, exact

# Запуск решения
numerical_solution, exact_solution = solve_poisson_gauss_seidel() 
