import math

import numpy as np
from tabulate import tabulate
import matplotlib.pyplot as plt


def func(x, y, z):
    return -z / x


# def func(x, y, z):
#    return 2 * x * z / ((x ** 2) + 1)


# def test_func(x):
#     return (-x ** 3) + 3 * x + 1

def test_func(x):
    return 1 + np.log1p(np.abs(x) - 1)


def L_def(x0, y0, z0, h):
    return h * func(x0, y0, z0)


def K_def(z, h):
    return h * z


def L_calc(x0, y0, z0, h, p):
    # h - step,
    # x0 - first x, y0=y(x0) z0=y`(x0)
    # p - max value of order

    if p == 1:
        return L_def(x0, y0, z0, h)
    if p == 4:
        return h * func(x0 + h,
                        y0 + K_calc(x0, y0, z0, h, p - 1),
                        z0 + L_calc(x0, y0, z0, h, p - 1))
    return h * func(x0 + h / 2,
                    y0 + K_calc(x0, y0, z0, h, p - 1) / 2,
                    z0 + L_calc(x0, y0, z0, h, p - 1) / 2)


def K_calc(x0, y0, z0, h, p):
    # h - step,
    # x0 - first x, y0=y(x0) z0=y`(x0)
    # p - max value of order

    if p == 1:
        return K_def(z0, h)
    if p == 4:
        return h * (z0 + L_calc(x0, y0, z0, h, p - 1))
    return h * (z0 + L_calc(x0, y0, z0, h, p - 1) / 2)


def y_delta(x0, y0, z0, h, p):
    return 1 / 6 * (
            K_calc(x0, y0, z0, h, 1) +
            2 * K_calc(x0, y0, z0, h, 2) +
            2 * K_calc(x0, y0, z0, h, 3) +
            K_calc(x0, y0, z0, h, 4))


def z_delta(x0, y0, z0, h, p):
    return 1 / 6 * (
            L_calc(x0, y0, z0, h, 1) +
            2 * L_calc(x0, y0, z0, h, 2) +
            2 * L_calc(x0, y0, z0, h, 3) +
            L_calc(x0, y0, z0, h, 4))


def y_point_def(x, y, z, h, p):
    return y + y_delta(x, y, z, h, p)


def z_point_def(x, y, z, h, p):
    return z + z_delta(x, y, z, h, p)


def k_list_gen(x0, xk, h):
    list = []
    a = 0
    for i in np.arange(x0, xk + h, h):
        list.append(a)
        a += 1
    return list


def x_list_gen(x0, xk, h):
    list = []
    for i in np.arange(x0, xk + h, h):
        list.append(i)
    return list


def main_list_gen(x0, y0, z0, h, p, x_list, k_list):
    list = [[], [], [], [], []]
    list[0].append(y0)
    list[1].append(z0)
    for i in k_list:
        list[0].append(y_point_def(x_list[i], list[0][i], list[1][i], h, p))
        list[1].append(z_point_def(x_list[i], list[0][i], list[1][i], h, p))
        list[2].append(y_delta(x_list[i], list[0][i], list[1][i], h, p))
        list[3].append(z_delta(x_list[i], list[0][i], list[1][i], h, p))
        list[4].append(test_func(x_list[i]))
    return list


def check(gen_list):
    list = []
    for i in range(0, len(gen_list[0]) - 2, 1):
        list.append((gen_list[0][i] - gen_list[4][i]) / gen_list[0][i])
    return list


def euler(x0, y0, z0, xk, h):
    xList = []
    yList = [y0]
    zList = [z0]
    i = 0
    while x0 <= xk:
        xList.append(x0)
        x0 += h

        y = yList[i] + h * zList[i]
        yList.append(y)

        z = zList[i] + h * (func(xList[i],yList[i], zList[i])+func(xList[i],yList[i+1], zList[i]))/2
        zList.append(z)
        i += 1
    xList.append(x0)
    return [xList, yList]


def graphic(x0, y0, z0, xk, h, mList, genlist):
    plt.title("y=1+ln|x|")
    plt.xlabel("X")
    plt.ylabel("f(x)")
    plt.plot(mList[0], mList[2], color='green', linestyle='solid', label='f(x)=1+nl|x|')
    plt.plot(mList[0], mList[1], color="blue", label='Метод Ейлера')
    plt.plot(mList[0], genlist[0], color="orange", label='Метод Рунге Кутте')

    plt.legend(loc='upper left')
    plt.grid(True)
    plt.show()


xk = 2
x0 = 1
y0 = 1
z0 = 1
h = 0.1
p=4



k_list = k_list_gen(x0, xk, h)
x_list = x_list_gen(x0, xk, h)

gen_list = main_list_gen(x0, y0, z0, h, p, x_list, k_list)
del gen_list[0][-1]
check_list = check(gen_list)

index = ['x', 'y', 'z', 'D_y', 'D_z', 'test', 'wrong']
index_res = ['x', 'y_accurate','y_runge', 'y_eulera']
table = [k_list, x_list, gen_list[0], gen_list[1], gen_list[2], gen_list[3], gen_list[4], check_list]
print(tabulate(table, headers='firstrow', tablefmt='grid', showindex=index))




mList = euler(x0, y0, z0, xk, h)
mList.append(gen_list[4])
graphic(x0, y0, z0, xk, h, mList, gen_list)


print("\n\n\n")
table_res = [k_list, x_list,gen_list[4], gen_list[0], mList[1]]
print(tabulate(table_res, headers='firstrow', tablefmt='grid', showindex=index_res))