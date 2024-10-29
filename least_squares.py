from math import log, sin
from Gauss import Gauss_method

a = 0
b = 3
h = (b - a)/10

X = []  # Точки X_i
Y = []  # Значения функции в точке X_i
n = 11

for i in range(0, n):
    X.append(a + i * h)
    Y.append(log(2 * X[i] + 1) + 2 * sin(3 * X[i]))  # Функция


def least_square_method(degree):  # МНК

    matrixBC = []  # Матрица СЛАУ

    for i in range(0, degree + 1):
        row = []
        for j in range(0, degree + 2):
            row.append(0)
        matrixBC.append(row)

    for p in range(0, degree + 1):  # Заполняем матрицу по формулам B_pq и C_p

        for q in range(p, degree + 1):
            for i in range(0, n):
                matrixBC[p][q] = matrixBC[p][q] + X[i] ** (p + q)   # B_pq
            matrixBC[q][p] = matrixBC[p][q]

        for i in range(0, n):
            matrixBC[p][degree + 1] = matrixBC[p][degree + 1] + Y[i] * X[i] ** p   # C_p

    return matrixBC


def find_phi_values(coefficients):  # Значения искомых функций Ф

    phi1 = []

    for i in range(0, n):
        phi = 0
        for j in range(0, len(coefficients)):
            phi = phi + coefficients[j] * X[i] ** j
        phi1.append(float(phi))

    return phi1


def calculate_epsilon(phi_values):  # Вычисление невязок
    epsilon = []
    for i in range(0, n):
        epsilon.append(phi_values[i] - Y[i])

    return epsilon


# Линейная функция
SLAU = least_square_method(1)  # Получаем СЛАУ, решение которого — коэффициенты искомой функции
coefficients_of_polinoms = Gauss_method(SLAU)  # Решаем СЛАУ методом Гаусса — получаем коэффициенты искомой функции
phi1 = find_phi_values(coefficients_of_polinoms)  # Поиск значений найденный функции в точках X_i
eps1 = calculate_epsilon(phi1)  # Вычисление невязок
print(coefficients_of_polinoms)


# Квадратная функция
SLAU = least_square_method(2)
coefficients_of_polinoms = Gauss_method(SLAU)
phi2 = find_phi_values(coefficients_of_polinoms)
eps2 = calculate_epsilon(phi2)
print(coefficients_of_polinoms)

# Кубическая функция
SLAU = least_square_method(3)
coefficients_of_polinoms = Gauss_method(SLAU)
phi3 = find_phi_values(coefficients_of_polinoms)
eps3 = calculate_epsilon(phi3)
print(coefficients_of_polinoms)

# Вывод таблицы
print('  X       Y        Ф1       Ф2       Ф3      Eps1     Eps2     Eps3')
for i in range(0, n):
    print(f'{X[i]: .2f} {Y[i]: .5f} {phi1[i]: .5f} {phi2[i]: .5f} {phi3[i]: .5f} {eps1[i]: .5f} {eps2[i]: .5f} {eps3[i]: .5f}')
print('                  Сумма квадратов невязок:')

