import matplotlib.pyplot as plt
import tkinter as tk
import numpy as np
import math

def function(x):
    return x * (1 + x) ** (1 / 3)

def tabulate(a = 1, b = 9, n = 8):
    step = (b - a) / n
    results = []
    for i in range(n + 1):            # <-- исправлено: цикл по счётчику
        x = a + i * step
        results.append((x, function(x)))
    return results

graphic = tabulate()
x_values = [point[0] for point in graphic]
y_values = [point[1] for point in graphic]

def mnk(x_values, y_values, m):
    n = len(x_values)
    A = np.zeros((n, m + 1))
    for i in range(n):
        for j in range(m + 1):
            A[i][j] = x_values[i] ** j
    y = np.array(y_values)
    Y = y[:, np.newaxis]

    AT = A.T
    ATA = AT @ A
    ATY = AT @ Y
    coeffs = np.linalg.solve(ATA, ATY)

    func = f'{coeffs[0][0]:.4f}'
    for i in range(1, len(coeffs)):
        func += f' + {coeffs[i][0]:.4f} * x ^ {i}'
    print(func)

    global y_apr_value
    y_apr_value = []
    apr_y = 0
    for x in x_values:
        for j in range(len(coeffs)):
            apr_y += coeffs[j][0] * x ** j
        y_apr_value.append(apr_y)
        apr_y = 0

    print(y_apr_value)

def lagrange_basis(x_values, i, x):
    p = 1.0
    for j in range(len(x_values)):
        if i != j:
            p *= (x - x_values[j]) / (x_values[i] - x_values[j])
    return p

def lagrange(x_values, y_values, x):
    result = 0.0
    for i in range(len(x_values)):
        result += y_values[i] * lagrange_basis(x_values, i, x)
    return result


def exp_aprox(x_values, y_values, N = 11):
    sumX = 0
    sumY = 0
    sumX2 = 0
    sumxY = 0
    for i in range(len(x_values)):
        Y = math.log(y_values[i])
        sumX += x_values[i]
        sumY += Y
        sumX2 += x_values[i] * x_values[i]
        sumxY += x_values[i] * Y

    b_exp = (N * sumxY - sumX * sumY) / (N * sumX2 - sumX**2)

    A = (sumY - b_exp * sumX) / N
    a_exp = math.exp(A)

    return a_exp, b_exp

print(x_values, y_values)
lagr_val = [lagrange(x_values, y_values, x) for x in x_values]
mnk(x_values, y_values, m = 3)
aex, bex = exp_aprox(x_values, y_values)
yexp = []
for x in x_values:
    yexp.append(aex * math.exp(bex * x))
plt.plot(x_values, yexp, color='green', linestyle=':', linewidth=2, marker='^', label='Функция экспоненты')
plt.plot(x_values, y_values, color='red', linestyle='-', linewidth=2, marker='o', label='f(x)')
# plt.plot(x_values, y_apr_value, color='blue', linestyle='--', linewidth=2, marker='s', label='МНК')
# plt.plot(x_values, lagr_val, color='green', linestyle=':', linewidth=2, marker='^', label='Лагранж')
plt.legend()
plt.grid(True)
plt.show()