import math


def copy_vec(x):
    return list(x)


def add_scaled(a, b, s):
    return [a[i] + s * b[i] for i in range(len(a))]


def axpy(a, sa, b, sb):
    return [sa * a[i] + sb * b[i] for i in range(len(a))]


def vec_add(*vecs):
    n = len(vecs[0])
    out = [0.0] * n
    for v in vecs:
        for i in range(n):
            out[i] += v[i]
    return out


def scale(v, s):
    return [s * x for x in v]


def max_abs_diff(a, b):
    return max(abs(a[i] - b[i]) for i in range(len(a)))


def max_abs(v):
    return max(abs(x) for x in v)


def linspace(a, b, h):
    """Сетка от a до b; h > 0. Если b < a, узлы идут влево с шагом -h."""
    if h <= 0:
        raise ValueError("шаг h должен быть положительным")
    direction = 1.0 if b >= a else -1.0
    step = direction * h
    n = int(round(abs(b - a) / h))
    if abs(a + n * step - b) > 1e-9 * max(1.0, abs(a), abs(b)):
        raise ValueError(f"отрезок [{a}, {b}] не делится на шаг h={h}")
    return [a + i * step for i in range(n + 1)], n


def runge_romberg(y_fine, y_coarse, p):
    """Refine fine-grid solution using coarse (step 2h) values at shared nodes.

    y_fine: values with step h (length 2m+1)
    y_coarse: values with step 2h (length m+1)
    Returns refined values on the coarse nodes and error estimates R.
    """
    if p < 1:
        raise ValueError("порядок p должен быть >= 1")
    denom = 2 ** p - 1
    refined = []
    errs = []
    for i, yc in enumerate(y_coarse):
        yf = y_fine[2 * i]
        r = (yf - yc) / denom
        refined.append(yf + r)
        errs.append(r)
    return refined, errs


def runge_romberg_component(ys_fine, ys_coarse, p, component=0):
    fine = [row[component] for row in ys_fine]
    coarse = [row[component] for row in ys_coarse]
    return runge_romberg(fine, coarse, p)


def gauss_solve(A, b):
    n = len(b)
    M = [A[i][:] + [b[i]] for i in range(n)]
    for k in range(n):
        piv = k
        mx = abs(M[k][k])
        for i in range(k + 1, n):
            v = abs(M[i][k])
            if v > mx:
                mx = v
                piv = i
        if mx == 0:
            raise ValueError("вырожденная система")
        if piv != k:
            M[k], M[piv] = M[piv], M[k]
        akk = M[k][k]
        for j in range(k, n + 1):
            M[k][j] /= akk
        for i in range(n):
            if i == k:
                continue
            factor = M[i][k]
            if factor == 0:
                continue
            for j in range(k, n + 1):
                M[i][j] -= factor * M[k][j]
    return [M[i][n] for i in range(n)]


def fmt_row(xs, width=14, prec=8):
    parts = []
    for x in xs:
        if isinstance(x, float):
            parts.append(f"{x:{width}.{prec}g}")
        else:
            parts.append(f"{x!s:>{width}}")
    return " ".join(parts)
