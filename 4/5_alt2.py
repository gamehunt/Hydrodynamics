import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def solve_poisson_equation():
    # Параметры сетки
    N = 20  # Количество узлов по x и y
    h = 1.0 / (N - 1)  # Шаг сетки
    x = np.linspace(0, 1, N)
    y = np.linspace(0, 1, N)
    X, Y = np.meshgrid(x, y)
    
    # Инициализация решения
    u = np.zeros((N, N))
    
    # Задание граничных условий
    # Левая граница: u(0, y) = y^2
    u[0, :] = y**2
    
    # Нижняя граница: u(x, 0) = x^2
    u[:, 0] = x**2
    
    # Верхняя граница: u(x, 1) = x^2 + 1
    u[:, -1] = x**2 + 1
    
    # Правая граница: du/dn = 2x (Неймана)
    # Используем одностороннюю разность для производной
    
    # Итерационный метод (Якоби или Гаусса-Зейделя)
    max_iter = 10000
    tolerance = 1e-6
    f = 4.0  # Правая часть уравнения Пуассона
    
    for _ in range(max_iter):
        u_old = u.copy()
        for i in range(1, N-1):
            for j in range(1, N-1):
                u[i, j] = 0.25 * (u_old[i+1, j] + u_old[i-1, j] + 
                                  u_old[i, j+1] + u_old[i, j-1] - h**2 * f)
        
        # Обновление правой границы (Неймана)
        for j in range(1, N-1):
            u[-1, j] = u[-2, j] + 2 * x[-1] * h  # du/dx = 2x => u_{N} = u_{N-1} + 2x h
        
        # Проверка сходимости
        if np.max(np.abs(u - u_old)) < tolerance:
            break
    
    # Точное решение
    exact_solution = X**2 + Y**2
    
    # Вычисление максимальной погрешности
    error = np.max(np.abs(u - exact_solution))
    print(f"Максимальная погрешность: {error:.6f}")
    
    # Визуализация
    fig = plt.figure()
    
    # Численное решение
    ax1 = fig.add_subplot(111)
    ax1.contourf(X, Y, u, cmap='viridis')
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    # ax1.set_zlabel('u(x, y)')
    ax1.set_title('Численное решение')
    ax1.set_box_aspect(1.0)
    
    # Точное решение
    # ax2 = fig.add_subplot(122, projection='3d')
    # ax2.plot_surface(X, Y, exact_solution, cmap='plasma')
    # ax2.set_xlabel('x')
    # ax2.set_ylabel('y')
    # ax2.set_zlabel('u(x, y)')
    # ax2.set_title('Точное решение')
    
    plt.tight_layout()
    plt.show()
    
    return u, exact_solution, error

# Запуск решения
numerical_solution, exact_solution, max_error = solve_poisson_equation()
