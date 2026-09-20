from lab3.linsolve import gauss_solve


def monomial_powers(x, degree):
    row = [1.0]
    p = 1.0
    for _ in range(degree):
        p *= x
        row.append(p)
    return row


def least_squares_poly(xs, ys, degree):
    """
    Приближающий многочлен степени degree через нормальную систему МНК.
    Возвращает коэффициенты a0..adegree (по возрастанию степени).
    """
    n_coef = degree + 1
    N = len(xs)
    if N < n_coef:
        raise ValueError("точек меньше, чем коэффициентов")

    A = [[0.0] * n_coef for _ in range(n_coef)]
    b = [0.0] * n_coef
    for j in range(N):
        phi = monomial_powers(xs[j], degree)
        for k in range(n_coef):
            b[k] += ys[j] * phi[k]
            for i in range(n_coef):
                A[k][i] += phi[k] * phi[i]
    coeffs = gauss_solve(A, b)
    return coeffs, A, b


def poly_eval(coeffs, x):
    s = 0.0
    p = 1.0
    for a in coeffs:
        s += a * p
        p *= x
    return s


def sse(xs, ys, coeffs):
    s = 0.0
    for i in range(len(xs)):
        r = ys[i] - poly_eval(coeffs, xs[i])
        s += r * r
    return s


def poly_string(coeffs):
    parts = []
    for k, a in enumerate(coeffs):
        if k == 0:
            parts.append(f"{a:.10g}")
            continue
        sign = "+" if a >= 0 else "-"
        if k == 1:
            parts.append(f" {sign} {abs(a):.10g}*x")
        else:
            parts.append(f" {sign} {abs(a):.10g}*x^{k}")
    return "".join(parts)
