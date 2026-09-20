import math


def arccot(x):
    """arcctg(x) ∈ (0, π)."""
    return math.atan2(1.0, x)


def sample_function(xs, f=None):
    if f is None:
        f = arccot
    return [f(x) for x in xs]


def lagrange_basis(xs, i, x):
    n = len(xs)
    li = 1.0
    xi = xs[i]
    for j in range(n):
        if j == i:
            continue
        li *= (x - xs[j]) / (xi - xs[j])
    return li


def lagrange_eval(xs, ys, x):
    s = 0.0
    for i in range(len(xs)):
        s += ys[i] * lagrange_basis(xs, i, x)
    return s


def divided_differences(xs, ys):
    n = len(xs)
    table = [[0.0] * n for _ in range(n)]
    for i in range(n):
        table[i][0] = ys[i]
    for j in range(1, n):
        for i in range(n - j):
            table[i][j] = (table[i + 1][j - 1] - table[i][j - 1]) / (xs[i + j] - xs[i])
    coeffs = [table[0][j] for j in range(n)]
    return coeffs, table


def newton_eval(xs, coeffs, x):
    n = len(coeffs)
    s = coeffs[0]
    prod = 1.0
    for k in range(1, n):
        prod *= x - xs[k - 1]
        s += coeffs[k] * prod
    return s


def newton_build(xs, ys):
    coeffs, table = divided_differences(xs, ys)
    return coeffs, table


def interpolation_error(true_y, approx_y):
    return abs(true_y - approx_y)


def poly_string_newton(xs, coeffs):
    parts = [f"{coeffs[0]:.10g}"]
    for k in range(1, len(coeffs)):
        factors = "".join(f"(x{xs[j]:+.10g})" for j in range(k))
        c = coeffs[k]
        sign = "+" if c >= 0 else "-"
        parts.append(f" {sign} {abs(c):.10g}*{factors}")
    return "".join(parts)


def poly_string_lagrange(xs, ys):
    terms = []
    n = len(xs)
    for i in range(n):
        den = 1.0
        factors = []
        for j in range(n):
            if j == i:
                continue
            den *= xs[i] - xs[j]
            factors.append(f"(x{xs[j]:+.10g})")
        coef = ys[i] / den
        sign = "+" if coef >= 0 else "-"
        terms.append(f" {sign} {abs(coef):.10g}*{'*'.join(factors)}")
    s = "".join(terms).lstrip()
    if s.startswith("+ "):
        s = s[2:]
    return s
