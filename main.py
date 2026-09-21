import matplotlib.pyplot as plt
import tkinter as tk
import numpy as np

def function(x):
    return x * (1 + x) ** (1 / 3)



def tabulate(a = 1, b = 9, n = 8):
    step = (b - a) / n
    results = []
    while a != b:
        x, y = a, function(a)
        results.append((a, function(a)))
        a += step
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


def lagrange(x_valuse, y_values):
    znam = x_values[1] - x_values[0]
    shab = []
    for i in range(len(x_values)):
        p = 1
        for j in range(len(x_values)):
            if i != j:
                p*= (x_values[i] - x_values[j])
        shab.append(p)
        p = 1
    print(shab)

print(x_values, y_values)
lagrange(x_values, y_values)
# mnk(x_values, y_values, m = 3)
# plt.plot(x_values, y_values, color = 'red')
# print(x_values, y_apr_value)
# plt.plot(x_values, y_apr_value, 'o--', color = 'blue' , markerfacecolor = 'green')
# plt.grid(True)
# plt.show()