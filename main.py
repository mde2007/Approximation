import math
import tkinter as tk

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg


def function(x):
    return x * (1 + x) ** (1 / 3)


def tabulate(a=1, b=9, n=8):
    step = (b - a) / n
    results = []
    for i in range(n + 1):
        x = a + i * step
        results.append((x, function(x)))
    return results


graphic = tabulate()
x_values = [point[0] for point in graphic]
y_values = [point[1] for point in graphic]
a_tab = 1
b_tab = 9
n_tab = 8
h_tab = (b_tab - a_tab) / n_tab
x_mid = [a_tab + h_tab * (j + 0.5) for j in range(n_tab)]
y_mid_exact = [function(x) for x in x_mid]


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

    y_apr_value = []
    for x in x_values:
        apr_y = 0
        for j in range(len(coeffs)):
            apr_y += coeffs[j][0] * x ** j
        y_apr_value.append(apr_y)

    coeffs_flat = [c[0] for c in coeffs]
    return y_apr_value, func, coeffs_flat


def eval_poly(coeffs, x):
    return sum(c * x ** j for j, c in enumerate(coeffs))


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


def exp_aprox(x_values, y_values, N=9):
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

    b_exp = (N * sumxY - sumX * sumY) / (N * sumX2 - sumX ** 2)
    A = (sumY - b_exp * sumX) / N
    a_exp = math.exp(A)

    return a_exp, b_exp



lagr_val = [lagrange(x_values, y_values, x) for x in x_values]
aex, bex = exp_aprox(x_values, y_values)
yexp = [aex * math.exp(bex * x) for x in x_values]
lagr_mid_val = [lagrange(x_values, y_values, x) for x in x_mid]
lagr_mid_dev = [abs(lagr_mid_val[j] - y_mid_exact[j]) for j in range(n_tab)]
lagr_max_dev = max(lagr_mid_dev)

canvas1 = None
canvas2 = None
mnk_dev_labels = []
lagr_dev_labels = []


def draw_first_graph(m):
    global canvas1
    if canvas1 is not None:
        canvas1.get_tk_widget().destroy()
    y_apr_value, func, coeffs = mnk(x_values, y_values, m)

    mnk_mid_val = [eval_poly(coeffs, x) for x in x_mid]
    mnk_mid_dev = [abs(mnk_mid_val[j] - y_mid_exact[j]) for j in range(n_tab)]
    mnk_max_dev = max(mnk_mid_dev)

    for lbl in mnk_dev_labels:
        lbl.destroy()
    mnk_dev_labels.clear()

    header = tk.Label(
        win,
        text=f'Погрешность МНК (m={m}), макс. отклонение = {mnk_max_dev:.6f}',
        font=('Arial', 10, 'bold')
    )
    header.place(x=0, y=580)
    mnk_dev_labels.append(header)

    for j in range(n_tab):
        lbl = tk.Label(
            win,
            text=(
                f'x={x_mid[j]:.4f}  f(x)={y_mid_exact[j]:.4f}  '
                f'P(x)={mnk_mid_val[j]:.4f}  |откл.|={mnk_mid_dev[j]:.6f}'
            ),
            font=('Arial', 9)
        )
        lbl.place(x=0, y=605 + 18 * j)
        mnk_dev_labels.append(lbl)

    fig1, ax1 = plt.subplots(figsize=(6, 4))
    ax1.plot(x_values, y_values,'r-o', linewidth=2, markersize=6, label='f(x)')
    ax1.plot(x_values, y_apr_value,'b--s', linewidth=2, markersize=6, label=f'МНК (m={m})')
    ax1.plot(x_values, lagr_val,'g:^', linewidth=2, markersize=6, label='Лагранж')
    ax1.set_xlabel('x', fontsize=11)
    ax1.set_ylabel('y', fontsize=11)
    ax1.set_title('Исходная функция, МНК и Лагранж', fontsize=12)
    ax1.grid(True, alpha=0.3)
    ax1.legend(loc='best')
    canvas1 = FigureCanvasTkAgg(fig1, master=win)
    canvas1.draw()
    canvas1.get_tk_widget().place(x=0, y=150, width=580, height=420)
    func_label.config(text=f'Полином МНК степени {m}:\n{func}')


def draw_lagrange_deviation():
    for lbl in lagr_dev_labels:
        lbl.destroy()
    lagr_dev_labels.clear()

    header = tk.Label(
        win,
        text=f'Погрешность Лагранжа, макс. отклонение = {lagr_max_dev:.6f}',
        font=('Arial', 10, 'bold')
    )
    header.place(x=600, y=580)
    lagr_dev_labels.append(header)

    for j in range(n_tab):
        lbl = tk.Label(
            win,
            text=(
                f'x={x_mid[j]:.4f}  f(x)={y_mid_exact[j]:.4f}  '
                f'L(x)={lagr_mid_val[j]:.4f}  |откл.|={lagr_mid_dev[j]:.6f}'
            ),
            font=('Arial', 9)
        )
        lbl.place(x=600, y=605 + 18 * j)
        lagr_dev_labels.append(lbl)


def draw_second_graph():
    global canvas2
    fig2, ax2 = plt.subplots(figsize=(6, 4))
    ax2.plot(x_values, y_values,'r-o', linewidth=2, markersize=6, label='f(x)')
    ax2.plot(x_values, yexp, 'g:^', linewidth=2, markersize=6, label='Функция экспоненты')
    ax2.set_xlabel('x', fontsize=11)
    ax2.set_ylabel('y', fontsize=11)
    ax2.set_title('Исходная функция и экспоненциальная аппроксимация', fontsize=12)
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc='best')
    canvas2 = FigureCanvasTkAgg(fig2, master=win)
    canvas2.draw()
    canvas2.get_tk_widget().place(x=600, y=150, width=580, height=420)


def button_function():
    m = int(ment.get())
    draw_first_graph(m)
    draw_second_graph()
    draw_lagrange_deviation()


win = tk.Tk()
win.geometry('1200x950')
win.title('Лаба')
title = tk.Label(win, text='Аппроксимация функции f(x) = x*(1+x)^(1/3)', font=('Arial', 14))
title.place(x=0, y=0)
text = tk.Label(win, text='Введите степень аппроксимации m для МНК', font=('Arial', 12))
text.place(x=0, y=40)
ment = tk.Entry(win)
ment.insert(0, '3')
ment.place(x=0, y=70)
but = tk.Button(win, text='Построить МНК', command=button_function)
but.place(x=150, y=68)
func_label = tk.Label(win, text='', font=('Arial', 10), justify='left')
func_label.place(x=0, y=100)

win.mainloop()