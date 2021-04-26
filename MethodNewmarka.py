import math
import matplotlib.pyplot as plt
import numpy as np


def func(u, ut, t):
    return math.sin(omega * h * t)


def newmark(dt, m, alpha, c):
    CONST = 1 / (m / (a1 * dt * dt) + (alpha * a2) / (a1 * dt) + c)

    u = 0
    u_t = 0
    u_2t = -(alpha * u_t) / m - (c * u) / m + k * func(u, u_t, 0) / m

    result_u.append([u])
    result_u_t.append([u_t])
    result_u_2t.append([u_2t])

    for i in range(int(n) - 1):
        u_n = CONST * m * (u / (a1 * dt * dt) + u_t / (a1 * dt) - u_2t * (1 - 1 / (2 * a1))) \
              + CONST * alpha * ((a2 * u) / (a1 * dt) + u_t * (1 - a2 / a1) - dt / 2 * (2 - a2 / a1) * u_2t) \
              + func(u, u_t, i+1) * CONST
        u_t_n = a2 * (u_n - u) / (a1 * dt) + (1 - a2 / a1) * u_t + u_2t * dt * (2 - a2 / a1) / 2
        u_2t_n = (u_n - u) / (a1 * dt * dt) - u_t / (a1 * dt) + (1 - 1 / (2 * a1)) * u_2t

        u = u_n
        u_t = u_t_n
        u_2t = u_2t_n

        result_u.append([u])
        result_u_t.append([u_t])
        result_u_2t.append([u_2t])


result_u = []
result_u_t = []
result_u_2t = []

a1 = 0.25
a2 = 0.50
k = 1
m = 5  # 5
a = k / m
c = 300  # 20
alpha = 0
beta = alpha / (2 * m)
w2 = c / m
w = math.sqrt(w2)
omega = 0.01  # w-0.01

t_start = 0
t_end = 10
n = 1000
h = (t_end - t_start) / n

newmark(h, m, alpha, c)
t_array = []
for i in range(int(n)):
    t_array.append(h * i)

fig, ax = plt.subplots()
ax.plot(t_array, result_u)

plt.show()
