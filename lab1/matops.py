import math


def zeros(n, m=None):
    if m is None:
        m = n
    return [[0.0] * m for _ in range(n)]


def identity(n):
    E = zeros(n)
    for i in range(n):
        E[i][i] = 1.0
    return E


def copy_mat(A):
    return [row[:] for row in A]


def copy_vec(x):
    return x[:]


def mat_mul(A, B):
    n = len(A)
    m = len(B[0])
    k = len(B)
    C = zeros(n, m)
    for i in range(n):
        for t in range(k):
            a = A[i][t]
            if a == 0:
                continue
            for j in range(m):
                C[i][j] += a * B[t][j]
    return C


def mat_vec(A, x):
    n = len(A)
    y = [0.0] * n
    for i in range(n):
        s = 0.0
        row = A[i]
        for j in range(len(row)):
            s += row[j] * x[j]
        y[i] = s
    return y


def transpose(A):
    n = len(A)
    m = len(A[0])
    return [[A[i][j] for i in range(n)] for j in range(m)]


def mat_sub(A, B):
    n = len(A)
    m = len(A[0])
    return [[A[i][j] - B[i][j] for j in range(m)] for i in range(n)]


def vec_sub(a, b):
    return [a[i] - b[i] for i in range(len(a))]


def vec_norm_inf(x):
    m = 0.0
    for v in x:
        av = abs(v)
        if av > m:
            m = av
    return m


def vec_norm2(x):
    s = 0.0
    for v in x:
        s += v * v
    return math.sqrt(s)


def mat_norm_inf(A):
    m = 0.0
    for row in A:
        s = 0.0
        for v in row:
            s += abs(v)
        if s > m:
            m = s
    return m


def max_abs_mat(A):
    m = 0.0
    for row in A:
        for v in row:
            av = abs(v)
            if av > m:
                m = av
    return m


def diag_from(vals):
    n = len(vals)
    D = zeros(n)
    for i in range(n):
        D[i][i] = vals[i]
    return D


def fmt_num(x, nd=8):
    if isinstance(x, complex):
        if abs(x.imag) < 1e-12:
            return f"{x.real:.{nd}g}"
        sign = "+" if x.imag >= 0 else "-"
        return f"{x.real:.{nd}g}{sign}{abs(x.imag):.{nd}g}j"
    return f"{x:.{nd}g}"


def fmt_vec(x, nd=8):
    return "[" + ", ".join(fmt_num(v, nd) for v in x) + "]"


def fmt_mat(A, nd=8):
    w = 14
    lines = []
    for row in A:
        parts = []
        for v in row:
            parts.append(f"{fmt_num(v, nd):>{w}}")
        lines.append(" ".join(parts))
    return "\n".join(lines)


def almost_equal_mat(A, B, tol=1e-8):
    n = len(A)
    m = len(A[0])
    for i in range(n):
        for j in range(m):
            if abs(A[i][j] - B[i][j]) > tol:
                return False
    return True
