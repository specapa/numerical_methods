import math

from lab1.lu import lu_decompose, lu_solve
from lab1.matops import vec_norm_inf, vec_sub


def make_trig_system(a):
    a = float(a)

    def F(x):
        return [
            x[0] - math.cos(x[1]) - a,
            x[1] - math.sin(x[0]) - a,
        ]

    def J(x):
        return [
            [1.0, math.sin(x[1])],
            [-math.cos(x[0]), 1.0],
        ]

    def phi(x):
        return [
            math.cos(x[1]) + a,
            math.sin(x[0]) + a,
        ]

    return F, J, phi


def simple_iteration(phi, x0, eps, max_iter=10000):
    x = [float(v) for v in x0]
    hist = []
    for k in range(1, max_iter + 1):
        nxt = phi(x)
        err = vec_norm_inf(vec_sub(nxt, x))
        hist.append(err)
        x = nxt
        if err < eps:
            return x, k, hist
    return x, max_iter, hist


def newton(F, J, x0, eps, max_iter=10000):
    x = [float(v) for v in x0]
    hist = []
    for k in range(1, max_iter + 1):
        fx = F(x)
        jx = J(x)
        L, U, P, _ = lu_decompose(jx)
        delta = lu_solve(L, U, P, fx)
        nxt = [x[i] - delta[i] for i in range(len(x))]
        err = vec_norm_inf(vec_sub(nxt, x))
        hist.append(err)
        x = nxt
        if err < eps:
            return x, k, hist
    return x, max_iter, hist
